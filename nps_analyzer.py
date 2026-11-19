#!/usr/bin/env python3
"""NPS & CSAT Survey CSV/JSON Analyzer CLI (Lite).

Zero-dependency Python 3.9+ standard-library CLI to compute Net Promoter Score
(NPS = %Promoters - %Detractors), CSAT (1-5), Wilson confidence intervals,
per-branch / per-staff breakdowns and rankings, period-over-period deltas, and
urgent detractor alert flags from CSV or JSON survey exports.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import math
import re
import statistics
import sys
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

__version__ = "1.0.0"

DEFAULT_CONFIDENCE = 0.95
DEFAULT_MIN_SAMPLE = 5

DEFAULT_URGENT_KEYWORDS: Tuple[str, ...] = (
    "cancel",
    "unsubscribe",
    "refund",
    "chargeback",
    "overcharged",
    "double charged",
    "fraud",
    "scam",
    "unauthorized",
    "billing error",
    "lawyer",
    "attorney",
    "lawsuit",
    "legal",
    "sue",
    "gdpr",
    "breach",
    "unsafe",
    "hazard",
    "broken",
    "outage",
    "crash",
    "unacceptable",
    "rude",
    "harass",
    "manager",
    "supervisor",
    "escalate",
    "furious",
    "horrible",
    "terrible",
    "worst",
    "never again",
    "last time",
    "switching",
    "churn",
)

ID_ALIASES: Tuple[str, ...] = ("id", "response_id", "survey_id", "ticket_id", "ref")
DATE_ALIASES: Tuple[str, ...] = (
    "date",
    "submitted_at",
    "created_at",
    "responded_at",
    "timestamp",
    "survey_date",
    "day",
)
PERIOD_ALIASES: Tuple[str, ...] = ("period", "cohort", "month", "quarter", "week")
NPS_ALIASES: Tuple[str, ...] = (
    "nps",
    "nps_score",
    "nps_rating",
    "recommend_score",
    "promoter_score",
)
CSAT_ALIASES: Tuple[str, ...] = (
    "csat",
    "csat_score",
    "csat_rating",
    "satisfaction",
    "satisfaction_score",
    "stars",
)
GENERIC_SCORE_ALIASES: Tuple[str, ...] = ("score", "rating", "value")
BRANCH_ALIASES: Tuple[str, ...] = (
    "branch",
    "location",
    "store",
    "clinic",
    "office",
    "site",
    "region",
    "department",
)
STAFF_ALIASES: Tuple[str, ...] = (
    "staff",
    "agent",
    "employee",
    "rep",
    "specialist",
    "advisor",
    "operator",
    "owner",
)
CHANNEL_ALIASES: Tuple[str, ...] = ("channel", "source", "platform", "touchpoint")
COMMENT_ALIASES: Tuple[str, ...] = (
    "comment",
    "feedback",
    "comments",
    "review",
    "text",
    "note",
    "notes",
    "verbatim",
    "reason",
)
CUSTOMER_ALIASES: Tuple[str, ...] = (
    "customer",
    "client",
    "customer_name",
    "client_name",
    "respondent",
    "email",
    "user",
)


@dataclass(frozen=True)
class SurveyResponse:
    """Normalized survey response record."""

    response_id: str
    response_date: Optional[date]
    period_label: str
    nps_score: Optional[float]
    csat_score: Optional[float]
    branch: str
    staff: str
    channel: str
    comment: str
    customer: str
    raw_fields: Dict[str, str] = field(default_factory=dict)
    source_index: int = 0

    def primary_score(self, mode: str) -> Optional[float]:
        return self.nps_score if mode == "nps10" else self.csat_score


@dataclass(frozen=True)
class WilsonInterval:
    """Confidence interval bounds [low, high] expressed as percentages or NPS points."""

    low: Optional[float]
    high: Optional[float]

    def to_dict(self, precision: int = 2) -> Dict[str, Optional[float]]:
        return {
            "low": round(self.low, precision) if self.low is not None else None,
            "high": round(self.high, precision) if self.high is not None else None,
        }


@dataclass(frozen=True)
class NpsMetrics:
    """Computed Net Promoter Score and trinomial breakdown metrics."""

    count: int
    promoters: int
    passives: int
    detractors: int
    promoter_pct: Optional[float]
    passive_pct: Optional[float]
    detractor_pct: Optional[float]
    nps: Optional[float]
    nps_ci: WilsonInterval
    promoter_ci: WilsonInterval
    passive_ci: WilsonInterval
    detractor_ci: WilsonInterval
    mean_score: Optional[float]
    median_score: Optional[float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "count": self.count,
            "promoters": self.promoters,
            "passives": self.passives,
            "detractors": self.detractors,
            "promoter_pct": _r(self.promoter_pct),
            "passive_pct": _r(self.passive_pct),
            "detractor_pct": _r(self.detractor_pct),
            "nps": _r(self.nps),
            "nps_wilson_ci": self.nps_ci.to_dict(2),
            "promoter_wilson_ci": self.promoter_ci.to_dict(2),
            "passive_wilson_ci": self.passive_ci.to_dict(2),
            "detractor_wilson_ci": self.detractor_ci.to_dict(2),
            "mean_score": _r(self.mean_score),
            "median_score": _r(self.median_score),
        }


@dataclass(frozen=True)
class CsatMetrics:
    """Computed 1-5 CSAT summary metrics."""

    count: int
    mean: Optional[float]
    median: Optional[float]
    csat_pct: Optional[float]
    csat_pct_ci: WilsonInterval
    satisfied_count: int
    neutral_count: int
    dissatisfied_count: int
    distribution: Dict[str, int]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "count": self.count,
            "mean": _r(self.mean),
            "median": _r(self.median),
            "csat_pct": _r(self.csat_pct),
            "csat_pct_wilson_ci": self.csat_pct_ci.to_dict(2),
            "satisfied_count": self.satisfied_count,
            "neutral_count": self.neutral_count,
            "dissatisfied_count": self.dissatisfied_count,
            "distribution": dict(self.distribution),
        }


@dataclass(frozen=True)
class DetractorAlert:
    """Flagged urgent detractor/dissatisfied response."""

    response_id: str
    response_date: Optional[date]
    branch: str
    staff: str
    channel: str
    customer: str
    score: float
    severity: str  # CRITICAL | HIGH
    reasons: Tuple[str, ...]
    matched_keywords: Tuple[str, ...]
    comment: str
    source_index: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "response_id": self.response_id,
            "date": self.response_date.isoformat() if self.response_date else None,
            "branch": self.branch,
            "staff": self.staff,
            "channel": self.channel,
            "customer": self.customer,
            "score": _r(self.score),
            "severity": self.severity,
            "reasons": list(self.reasons),
            "matched_keywords": list(self.matched_keywords),
            "comment": self.comment,
        }


@dataclass(frozen=True)
class PeriodSnapshot:
    """Metrics snapshot for a single period."""

    label: str
    count: int
    nps: Optional[float]
    promoter_pct: Optional[float]
    detractor_pct: Optional[float]
    csat_mean: Optional[float]
    csat_pct: Optional[float]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "label": self.label,
            "count": self.count,
            "nps": _r(self.nps),
            "promoter_pct": _r(self.promoter_pct),
            "detractor_pct": _r(self.detractor_pct),
            "csat_mean": _r(self.csat_mean),
            "csat_pct": _r(self.csat_pct),
        }


@dataclass(frozen=True)
class PeriodDelta:
    """Period-over-period comparison between previous and latest periods."""

    available: bool
    granularity: str
    previous: Optional[PeriodSnapshot]
    current: Optional[PeriodSnapshot]
    delta_nps: Optional[float]
    delta_promoter_pct: Optional[float]
    delta_detractor_pct: Optional[float]
    delta_csat_mean: Optional[float]
    delta_csat_pct: Optional[float]
    trend: str

    def to_dict(self) -> Dict[str, Any]:
        return {
            "available": self.available,
            "granularity": self.granularity,
            "previous_period": self.previous.to_dict() if self.previous else None,
            "current_period": self.current.to_dict() if self.current else None,
            "delta_nps": _r(self.delta_nps),
            "delta_promoter_pct": _r(self.delta_promoter_pct),
            "delta_detractor_pct": _r(self.delta_detractor_pct),
            "delta_csat_mean": _r(self.delta_csat_mean),
            "delta_csat_pct": _r(self.delta_csat_pct),
            "trend": self.trend,
        }


@dataclass(frozen=True)
class GroupBreakdown:
    """Per-group metrics, ranking, Wilson CI, alerts, and period delta."""

    rank: int
    group: str
    count: int
    low_sample: bool
    nps_metrics: Optional[NpsMetrics]
    csat_metrics: Optional[CsatMetrics]
    delta_nps: Optional[float]
    delta_csat_mean: Optional[float]
    delta_csat_pct: Optional[float]
    urgent_alert_count: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rank": self.rank,
            "group": self.group,
            "count": self.count,
            "low_sample": self.low_sample,
            "nps": self.nps_metrics.to_dict() if self.nps_metrics else None,
            "csat": self.csat_metrics.to_dict() if self.csat_metrics else None,
            "delta_nps": _r(self.delta_nps),
            "delta_csat_mean": _r(self.delta_csat_mean),
            "delta_csat_pct": _r(self.delta_csat_pct),
            "urgent_alert_count": self.urgent_alert_count,
        }


@dataclass(frozen=True)
class ThresholdCheck:
    """CI/ops threshold verification outcome."""

    enabled: bool
    fail_below: Optional[float]
    metric: Optional[str]
    actual: Optional[float]
    passed: bool

    def to_dict(self) -> Dict[str, Any]:
        return {
            "enabled": self.enabled,
            "fail_below": _r(self.fail_below),
            "metric": self.metric,
            "actual": _r(self.actual),
            "passed": self.passed,
        }


@dataclass(frozen=True)
class AnalysisReport:
    """Complete survey analysis output."""

    mode: str
    group_by: str
    confidence: float
    z_score: float
    min_sample: int
    total_rows: int
    analyzed_count: int
    skipped_missing_score: int
    nps_metrics: Optional[NpsMetrics]
    csat_metrics: Optional[CsatMetrics]
    period_delta: PeriodDelta
    groups: Tuple[GroupBreakdown, ...]
    total_detractors: int
    alerts: Tuple[DetractorAlert, ...]
    threshold_check: ThresholdCheck

    def to_dict(self) -> Dict[str, Any]:
        critical_count = sum(1 for a in self.alerts if a.severity == "CRITICAL")
        high_count = sum(1 for a in self.alerts if a.severity == "HIGH")
        return {
            "mode": self.mode,
            "group_by": self.group_by,
            "confidence": self.confidence,
            "z_score": round(self.z_score, 6),
            "min_sample": self.min_sample,
            "summary": {
                "total_rows": self.total_rows,
                "analyzed_count": self.analyzed_count,
                "skipped_missing_score": self.skipped_missing_score,
                "nps": self.nps_metrics.to_dict() if self.nps_metrics else None,
                "csat": self.csat_metrics.to_dict() if self.csat_metrics else None,
            },
            "period_delta": self.period_delta.to_dict(),
            "groups": [g.to_dict() for g in self.groups],
            "detractor_alerts": {
                "total_detractors": self.total_detractors,
                "urgent_count": len(self.alerts),
                "critical_count": critical_count,
                "high_count": high_count,
                "alerts": [a.to_dict() for a in self.alerts],
            },
            "threshold_check": self.threshold_check.to_dict(),
        }


def _r(value: Optional[float], digits: int = 2) -> Optional[float]:
    if value is None:
        return None
    return round(float(value), digits)


def confidence_to_z(confidence: float) -> float:
    """Convert confidence level (0 < confidence < 1) into standard normal two-sided z-score."""
    if not (0.0 < confidence < 1.0):
        raise ValueError(f"Confidence must be strictly between 0 and 1, got {confidence}.")
    return statistics.NormalDist().inv_cdf((1.0 + confidence) / 2.0)


def wilson_proportion_interval(successes: int, total: int, z: float) -> WilsonInterval:
    """Compute Wilson score confidence interval (0..100%) for a binomial proportion."""
    if total <= 0:
        return WilsonInterval(None, None)
    if successes < 0 or successes > total:
        raise ValueError(f"Invalid proportion counts: successes={successes}, total={total}.")

    p = successes / total
    z2 = z * z
    denom = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / denom
    margin = (z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * total)) / total)) / denom
    low = 0.0 if successes == 0 else max(0.0, (center - margin) * 100.0)
    high = 100.0 if successes == total else min(100.0, (center + margin) * 100.0)
    return WilsonInterval(low=low, high=high)


def wilson_nps_interval(promoters: int, passives: int, detractors: int, z: float) -> WilsonInterval:
    """Compute Newcombe-Wilson score CI (-100..+100) for trinomial NPS = p_P - p_D.

    Uses Wilson score bounds for Promoter and Detractor proportions paired with
    Newcombe's Method 10 correlation adjustment for mutually exclusive outcomes.
    """
    total = promoters + passives + detractors
    if total <= 0:
        return WilsonInterval(None, None)

    p_prom = promoters / total
    p_det = detractors / total

    ci_prom = wilson_proportion_interval(promoters, total, z)
    ci_det = wilson_proportion_interval(detractors, total, z)
    assert ci_prom.low is not None and ci_prom.high is not None
    assert ci_det.low is not None and ci_det.high is not None

    l_p = ci_prom.low / 100.0
    u_p = ci_prom.high / 100.0
    l_d = ci_det.low / 100.0
    u_d = ci_det.high / 100.0

    # Paired 2x2 contingency table where (promoter=1, detractor=1) cell is 0.
    denom_sq = promoters * (passives + detractors) * detractors * (passives + promoters)
    phi = (-promoters * detractors / math.sqrt(denom_sq)) if denom_sq > 0 else 0.0

    delta_low_sq = max(
        0.0,
        (p_prom - l_p) ** 2 - 2.0 * phi * (p_prom - l_p) * (u_d - p_det) + (u_d - p_det) ** 2,
    )
    delta_high_sq = max(
        0.0,
        (u_p - p_prom) ** 2 - 2.0 * phi * (u_p - p_prom) * (p_det - l_d) + (p_det - l_d) ** 2,
    )

    nps_prop = p_prom - p_det
    low = max(-100.0, (nps_prop - math.sqrt(delta_low_sq)) * 100.0)
    high = min(100.0, (nps_prop + math.sqrt(delta_high_sq)) * 100.0)
    return WilsonInterval(low=low, high=high)


def parse_survey_date(raw_value: str) -> Optional[date]:
    """Parse common ISO-8601 and regional date formats cleanly for Python 3.9+."""
    cleaned = raw_value.strip()
    if not cleaned:
        return None

    # Handle ISO timestamps such as 2026-04-15T14:22:00Z or 2026-04-15 14:22:00+02:00
    date_candidate = cleaned
    if "T" in cleaned:
        date_candidate = cleaned.split("T", 1)[0].strip()
    elif " " in cleaned and (
        re.match(r"^\d{4}[-/]\d{2}[-/]\d{2}\s+", cleaned)
        or re.match(r"^\d{2}[./]\d{2}[./]\d{4}\s+", cleaned)
    ):
        date_candidate = cleaned.split(" ", 1)[0].strip()

    try:
        return date.fromisoformat(date_candidate)
    except ValueError:
        pass

    for fmt in ("%Y/%m/%d", "%d.%m.%Y", "%m/%d/%Y", "%d/%m/%Y"):
        try:
            return datetime.strptime(date_candidate, fmt).date()
        except ValueError:
            continue

    raise ValueError(f"Unrecognized date format '{raw_value}'. Expected YYYY-MM-DD or ISO-8601 timestamp.")


def _detect_csv_delimiter(clean_text: str) -> str:
    """Detect CSV delimiter from the header row alone to avoid punctuation in comments."""
    first_line = ""
    for line in clean_text.splitlines():
        if line.strip():
            first_line = line
            break
    if ";" in first_line and "," not in first_line:
        return ";"
    if "\t" in first_line and "," not in first_line and ";" not in first_line:
        return "\t"
    return ","


def _first_present(norm_row: Dict[str, Any], aliases: Sequence[str], default: str = "") -> str:
    for key in aliases:
        if key in norm_row and norm_row[key] is not None:
            val = str(norm_row[key]).strip()
            if val != "":
                return val
    return default


def _parse_numeric_score(raw: str, row_ref: str, min_val: float, max_val: float, label: str) -> float:
    cleaned = raw.strip()
    try:
        val = float(cleaned)
    except ValueError as exc:
        raise ValueError(f"Survey '{row_ref}' has non-numeric {label} score '{raw}'.") from exc
    if math.isnan(val) or math.isinf(val) or not (min_val <= val <= max_val):
        raise ValueError(
            f"Survey '{row_ref}' has {label} score {raw!r} outside valid range [{int(min_val)}, {int(max_val)}]."
        )
    return val


def _parse_survey_record(
    raw: Dict[str, Any],
    index: int,
    mode: str,
    score_column: Optional[str] = None,
    csat_column: Optional[str] = None,
    date_column: Optional[str] = None,
    comment_column: Optional[str] = None,
) -> SurveyResponse:
    norm: Dict[str, Any] = {str(k).strip().lower(): v for k, v in raw.items() if k is not None}
    raw_fields: Dict[str, str] = {
        str(k).strip().lower(): ("" if v is None else str(v).strip())
        for k, v in raw.items()
        if k is not None
    }

    response_id = _first_present(norm, ID_ALIASES, f"ROW-{index + 1}")

    date_keys = (date_column.strip().lower(),) if date_column else DATE_ALIASES
    date_raw = _first_present(norm, date_keys)
    parsed_date = parse_survey_date(date_raw) if date_raw else None
    period_label = _first_present(norm, PERIOD_ALIASES)

    comment_keys = (comment_column.strip().lower(),) if comment_column else COMMENT_ALIASES
    comment = _first_present(norm, comment_keys)

    branch = _first_present(norm, BRANCH_ALIASES, "Unassigned")
    staff = _first_present(norm, STAFF_ALIASES, "Unassigned")
    channel = _first_present(norm, CHANNEL_ALIASES, "Unknown")
    customer = _first_present(norm, CUSTOMER_ALIASES)

    # Resolve NPS and CSAT columns cleanly according to mode and explicit overrides
    nps_score: Optional[float] = None
    csat_score: Optional[float] = None

    explicit_nps_col = score_column.strip().lower() if (score_column and mode == "nps10") else None
    explicit_csat_col = csat_column.strip().lower() if csat_column else (
        score_column.strip().lower() if (score_column and mode == "csat5") else None
    )

    nps_aliases = (explicit_nps_col,) if explicit_nps_col else (
        NPS_ALIASES + GENERIC_SCORE_ALIASES if mode == "nps10" else NPS_ALIASES
    )
    csat_aliases = (explicit_csat_col,) if explicit_csat_col else (
        CSAT_ALIASES + GENERIC_SCORE_ALIASES if mode == "csat5" else CSAT_ALIASES
    )

    nps_raw = _first_present(norm, nps_aliases)
    if nps_raw != "":
        nps_score = _parse_numeric_score(nps_raw, response_id, 0.0, 10.0, "NPS")

    csat_raw = _first_present(norm, csat_aliases)
    if csat_raw != "":
        csat_score = _parse_numeric_score(csat_raw, response_id, 1.0, 5.0, "CSAT")

    return SurveyResponse(
        response_id=response_id,
        response_date=parsed_date,
        period_label=period_label,
        nps_score=nps_score,
        csat_score=csat_score,
        branch=branch,
        staff=staff,
        channel=channel,
        comment=comment,
        customer=customer,
        raw_fields=raw_fields,
        source_index=index,
    )


def load_surveys_text(
    text: str,
    mode: str = "nps10",
    format_hint: str = "auto",
    score_column: Optional[str] = None,
    csat_column: Optional[str] = None,
    date_column: Optional[str] = None,
    comment_column: Optional[str] = None,
) -> List[SurveyResponse]:
    """Load survey responses from CSV or JSON text."""
    stripped = text.lstrip("\ufeff \t\r\n")
    if not stripped:
        raise ValueError("Survey input is empty.")

    fmt = format_hint.lower()
    if fmt == "auto":
        fmt = "json" if stripped[0] in ("{", "[") else "csv"

    responses: List[SurveyResponse] = []

    if fmt == "json":
        try:
            payload = json.loads(stripped)
        except json.JSONDecodeError as exc:
            raise ValueError(f"Invalid JSON survey input: {exc}") from exc

        if isinstance(payload, list):
            raw_rows = payload
        elif isinstance(payload, dict):
            raw_rows = None
            for key in ("surveys", "responses", "rows", "items", "data", "records"):
                if key in payload:
                    raw_rows = payload[key]
                    break
            if raw_rows is None:
                raise ValueError(
                    "JSON survey object must contain a list under 'surveys', 'responses', 'rows', 'items', or 'data'."
                )
            if not isinstance(raw_rows, list):
                raise ValueError("Survey JSON collection field must be a list.")
        else:
            raise ValueError("Survey JSON root must be an array or an object containing survey records.")

        for idx, item in enumerate(raw_rows):
            if not isinstance(item, dict):
                raise ValueError(f"Survey record at index {idx} must be a JSON object.")
            responses.append(
                _parse_survey_record(
                    item,
                    idx,
                    mode=mode,
                    score_column=score_column,
                    csat_column=csat_column,
                    date_column=date_column,
                    comment_column=comment_column,
                )
            )
        return responses

    if fmt == "csv":
        clean_text = text.lstrip("\ufeff")
        delimiter = _detect_csv_delimiter(clean_text)
        reader = csv.DictReader(io.StringIO(clean_text), delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("CSV survey file has no header row.")
        for idx, row in enumerate(reader):
            if row is None:
                continue
            if all((v is None or str(v).strip() == "") for v in row.values()):
                continue
            responses.append(
                _parse_survey_record(
                    row,
                    idx,
                    mode=mode,
                    score_column=score_column,
                    csat_column=csat_column,
                    date_column=date_column,
                    comment_column=comment_column,
                )
            )
        return responses

    raise ValueError(f"Unsupported input format hint '{format_hint}'.")


def load_surveys_file(
    path: Path,
    mode: str = "nps10",
    score_column: Optional[str] = None,
    csat_column: Optional[str] = None,
    date_column: Optional[str] = None,
    comment_column: Optional[str] = None,
) -> List[SurveyResponse]:
    """Load survey records from file path or '-' for stdin."""
    if str(path) == "-":
        return load_surveys_text(
            sys.stdin.read(),
            mode=mode,
            format_hint="auto",
            score_column=score_column,
            csat_column=csat_column,
            date_column=date_column,
            comment_column=comment_column,
        )
    if not path.is_file():
        raise ValueError(f"Survey file not found: {path}")
    text = path.read_text(encoding="utf-8-sig")
    suffix = path.suffix.lower()
    hint = "json" if suffix == ".json" else "csv" if suffix == ".csv" else "auto"
    return load_surveys_text(
        text,
        mode=mode,
        format_hint=hint,
        score_column=score_column,
        csat_column=csat_column,
        date_column=date_column,
        comment_column=comment_column,
    )


def compute_nps_metrics(scores: Sequence[float], z: float) -> Optional[NpsMetrics]:
    """Compute NPS (-100..+100), splits, and Wilson confidence intervals."""
    n = len(scores)
    if n == 0:
        return None

    promoters = sum(1 for s in scores if s >= 9.0)
    passives = sum(1 for s in scores if 7.0 <= s < 9.0)
    detractors = sum(1 for s in scores if s < 7.0)

    prom_pct = (promoters / n) * 100.0
    pass_pct = (passives / n) * 100.0
    det_pct = (detractors / n) * 100.0
    nps_val = prom_pct - det_pct

    return NpsMetrics(
        count=n,
        promoters=promoters,
        passives=passives,
        detractors=detractors,
        promoter_pct=prom_pct,
        passive_pct=pass_pct,
        detractor_pct=det_pct,
        nps=nps_val,
        nps_ci=wilson_nps_interval(promoters, passives, detractors, z),
        promoter_ci=wilson_proportion_interval(promoters, n, z),
        passive_ci=wilson_proportion_interval(passives, n, z),
        detractor_ci=wilson_proportion_interval(detractors, n, z),
        mean_score=statistics.fmean(scores),
        median_score=statistics.median(scores),
    )


def compute_csat_metrics(scores: Sequence[float], z: float) -> Optional[CsatMetrics]:
    """Compute 1-5 CSAT mean, median, Top-2-Box satisfaction %, and Wilson CI."""
    n = len(scores)
    if n == 0:
        return None

    satisfied = sum(1 for s in scores if s >= 4.0)
    neutral = sum(1 for s in scores if 3.0 <= s < 4.0)
    dissatisfied = sum(1 for s in scores if s < 3.0)
    csat_pct = (satisfied / n) * 100.0

    dist: Dict[str, int] = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0}
    for s in scores:
        bucket = min(5, max(1, int(round(s))))
        dist[str(bucket)] += 1

    return CsatMetrics(
        count=n,
        mean=statistics.fmean(scores),
        median=statistics.median(scores),
        csat_pct=csat_pct,
        csat_pct_ci=wilson_proportion_interval(satisfied, n, z),
        satisfied_count=satisfied,
        neutral_count=neutral,
        dissatisfied_count=dissatisfied,
        distribution=dist,
    )


def _match_keywords(comment: str, keywords: Sequence[str]) -> Tuple[str, ...]:
    if not comment.strip():
        return ()
    lower_comment = comment.lower()
    matched: List[str] = []
    for kw in keywords:
        token = kw.strip().lower()
        if not token:
            continue
        if " " in token:
            if token in lower_comment:
                matched.append(token)
        else:
            if re.search(rf"\b{re.escape(token)}\w*\b", lower_comment):
                matched.append(token)
    # Deduplicate preserving order
    seen = set()
    unique: List[str] = []
    for item in matched:
        if item not in seen:
            seen.add(item)
            unique.append(item)
    return tuple(unique)


def detect_detractor_alerts(
    responses: Sequence[SurveyResponse],
    mode: str,
    urgent_keywords: Sequence[str] = DEFAULT_URGENT_KEYWORDS,
    alert_max_score: Optional[float] = None,
) -> Tuple[int, Tuple[DetractorAlert, ...]]:
    """Identify urgent detractor comments with transparent severity & reason codes."""
    if mode == "nps10":
        detractor_cutoff = 6.0 if alert_max_score is None else alert_max_score
        severe_cutoff = min(3.0, detractor_cutoff)
    else:
        detractor_cutoff = 2.0 if alert_max_score is None else alert_max_score
        severe_cutoff = min(1.0, detractor_cutoff)

    total_detractors = 0
    alerts: List[DetractorAlert] = []

    for r in responses:
        score = r.primary_score(mode)
        if score is None or score > detractor_cutoff:
            continue
        total_detractors += 1

        is_severe_score = score <= severe_cutoff
        matched_kws = _match_keywords(r.comment, urgent_keywords)

        if not is_severe_score and not matched_kws:
            continue

        reasons: List[str] = []
        if is_severe_score:
            reasons.append(f"severe_score<={int(severe_cutoff)}")
        if matched_kws:
            reasons.append("urgent_keyword")

        if is_severe_score and matched_kws:
            severity = "CRITICAL"
        elif is_severe_score and score <= 1.0:
            severity = "CRITICAL"
        else:
            severity = "HIGH"

        alerts.append(
            DetractorAlert(
                response_id=r.response_id,
                response_date=r.response_date,
                branch=r.branch,
                staff=r.staff,
                channel=r.channel,
                customer=r.customer,
                score=score,
                severity=severity,
                reasons=tuple(reasons),
                matched_keywords=matched_kws,
                comment=r.comment,
                source_index=r.source_index,
            )
        )

    severity_rank = {"CRITICAL": 0, "HIGH": 1}
    alerts.sort(
        key=lambda a: (
            severity_rank.get(a.severity, 9),
            a.score,
            -(a.response_date.toordinal() if a.response_date else 0),
            a.response_id,
        )
    )
    return total_detractors, tuple(alerts)


def _bucket_period_key(resp: SurveyResponse, granularity: str) -> Optional[str]:
    if granularity == "explicit":
        return resp.period_label or None
    if resp.response_date is None:
        return None
    if granularity == "month":
        return f"{resp.response_date.year:04d}-{resp.response_date.month:02d}"
    if granularity == "week":
        iso_year, iso_week, _ = resp.response_date.isocalendar()
        return f"{iso_year:04d}-W{iso_week:02d}"
    if granularity == "day":
        return resp.response_date.isoformat()
    return None


def _resolve_period_buckets(
    responses: Sequence[SurveyResponse], period_mode: str
) -> Tuple[str, Dict[str, List[SurveyResponse]], List[str]]:
    """Bucket responses by period granularity and return ordered period keys."""
    if period_mode == "none" or not responses:
        return "none", {}, []

    if period_mode == "split":
        dated = sorted(
            [r for r in responses if r.response_date is not None],
            key=lambda r: (r.response_date, r.source_index),  # type: ignore[arg-type]
        )
        if len(dated) < 2:
            return "split", {}, []
        mid = len(dated) // 2
        buckets = {"First Half": dated[:mid], "Second Half": dated[mid:]}
        return "split", buckets, ["First Half", "Second Half"]

    if period_mode == "auto":
        dated_count = sum(1 for r in responses if r.response_date is not None)
        explicit_count = sum(1 for r in responses if r.period_label)
        if dated_count >= 2:
            for gran in ("month", "week", "day"):
                buckets: Dict[str, List[SurveyResponse]] = {}
                for r in responses:
                    key = _bucket_period_key(r, gran)
                    if key is not None:
                        buckets.setdefault(key, []).append(r)
                if len(buckets) >= 2:
                    return gran, buckets, sorted(buckets.keys())
        if explicit_count >= 2:
            buckets = {}
            for r in responses:
                if r.period_label:
                    buckets.setdefault(r.period_label, []).append(r)
            if len(buckets) >= 2:
                return "explicit", buckets, sorted(buckets.keys())
        return "auto", {}, []

    # Explicit granularity: month, week, or day
    buckets = {}
    for r in responses:
        key = _bucket_period_key(r, period_mode)
        if key is not None:
            buckets.setdefault(key, []).append(r)
    if len(buckets) < 2:
        return period_mode, {}, []
    return period_mode, buckets, sorted(buckets.keys())


def _snapshot_for_responses(label: str, items: Sequence[SurveyResponse], z: float) -> PeriodSnapshot:
    nps_scores = [r.nps_score for r in items if r.nps_score is not None]
    csat_scores = [r.csat_score for r in items if r.csat_score is not None]
    nps_m = compute_nps_metrics(nps_scores, z) if nps_scores else None
    csat_m = compute_csat_metrics(csat_scores, z) if csat_scores else None
    return PeriodSnapshot(
        label=label,
        count=len(items),
        nps=nps_m.nps if nps_m else None,
        promoter_pct=nps_m.promoter_pct if nps_m else None,
        detractor_pct=nps_m.detractor_pct if nps_m else None,
        csat_mean=csat_m.mean if csat_m else None,
        csat_pct=csat_m.csat_pct if csat_m else None,
    )


def compute_period_delta(
    responses: Sequence[SurveyResponse],
    mode: str,
    period_mode: str,
    z: float,
    previous_key: Optional[str] = None,
    current_key: Optional[str] = None,
) -> Tuple[PeriodDelta, str, Optional[str], Optional[str]]:
    """Compute overall period-over-period delta and return resolved keys for group deltas."""
    granularity, buckets, ordered_keys = _resolve_period_buckets(responses, period_mode)
    if not ordered_keys or len(ordered_keys) < 2:
        return (
            PeriodDelta(
                available=False,
                granularity=granularity,
                previous=None,
                current=None,
                delta_nps=None,
                delta_promoter_pct=None,
                delta_detractor_pct=None,
                delta_csat_mean=None,
                delta_csat_pct=None,
                trend="unavailable",
            ),
            granularity,
            None,
            None,
        )

    prev_label = previous_key if previous_key else ordered_keys[-2]
    curr_label = current_key if current_key else ordered_keys[-1]
    if prev_label not in buckets:
        raise ValueError(f"Previous period '{prev_label}' not found in dataset ({', '.join(ordered_keys)}).")
    if curr_label not in buckets:
        raise ValueError(f"Current period '{curr_label}' not found in dataset ({', '.join(ordered_keys)}).")

    prev_snap = _snapshot_for_responses(prev_label, buckets[prev_label], z)
    curr_snap = _snapshot_for_responses(curr_label, buckets[curr_label], z)

    delta_nps = (curr_snap.nps - prev_snap.nps) if (curr_snap.nps is not None and prev_snap.nps is not None) else None
    delta_prom = (
        (curr_snap.promoter_pct - prev_snap.promoter_pct)
        if (curr_snap.promoter_pct is not None and prev_snap.promoter_pct is not None)
        else None
    )
    delta_det = (
        (curr_snap.detractor_pct - prev_snap.detractor_pct)
        if (curr_snap.detractor_pct is not None and prev_snap.detractor_pct is not None)
        else None
    )
    delta_csat_mean = (
        (curr_snap.csat_mean - prev_snap.csat_mean)
        if (curr_snap.csat_mean is not None and prev_snap.csat_mean is not None)
        else None
    )
    delta_csat_pct = (
        (curr_snap.csat_pct - prev_snap.csat_pct)
        if (curr_snap.csat_pct is not None and prev_snap.csat_pct is not None)
        else None
    )

    primary_delta = delta_nps if mode == "nps10" else delta_csat_mean
    if primary_delta is None:
        trend = "unavailable"
    elif primary_delta > 0.05:
        trend = "improving"
    elif primary_delta < -0.05:
        trend = "declining"
    else:
        trend = "flat"

    return (
        PeriodDelta(
            available=True,
            granularity=granularity,
            previous=prev_snap,
            current=curr_snap,
            delta_nps=delta_nps,
            delta_promoter_pct=delta_prom,
            delta_detractor_pct=delta_det,
            delta_csat_mean=delta_csat_mean,
            delta_csat_pct=delta_csat_pct,
            trend=trend,
        ),
        granularity,
        prev_label,
        curr_label,
    )


def _group_value(resp: SurveyResponse, group_by: str) -> str:
    dims = [d.strip().lower() for d in group_by.split(",") if d.strip()]
    parts: List[str] = []
    for dim in dims:
        if dim == "branch":
            parts.append(resp.branch)
        elif dim == "staff":
            parts.append(resp.staff)
        elif dim == "channel":
            parts.append(resp.channel)
        else:
            val = resp.raw_fields.get(dim, "").strip()
            parts.append(val if val else "Unassigned")
    return " / ".join(parts)


def analyze_surveys(
    responses: Iterable[SurveyResponse],
    mode: str = "nps10",
    group_by: str = "branch",
    confidence: float = DEFAULT_CONFIDENCE,
    period_mode: str = "auto",
    previous_period: Optional[str] = None,
    current_period: Optional[str] = None,
    min_sample: int = DEFAULT_MIN_SAMPLE,
    branch_filter: Optional[str] = None,
    staff_filter: Optional[str] = None,
    urgent_keywords: Sequence[str] = DEFAULT_URGENT_KEYWORDS,
    alert_max_score: Optional[float] = None,
    fail_below: Optional[float] = None,
    csat_threshold_metric: str = "auto",
) -> AnalysisReport:
    """Analyze survey responses and produce full report."""
    if mode not in ("nps10", "csat5"):
        raise ValueError(f"Invalid mode '{mode}', expected 'nps10' or 'csat5'.")
    if min_sample < 1:
        raise ValueError("--min-sample must be >= 1.")

    z = confidence_to_z(confidence)

    filtered = list(responses)
    if branch_filter:
        b_norm = branch_filter.strip().lower()
        filtered = [r for r in filtered if r.branch.strip().lower() == b_norm]
    if staff_filter:
        s_norm = staff_filter.strip().lower()
        filtered = [r for r in filtered if r.staff.strip().lower() == s_norm]

    total_rows = len(filtered)
    usable = [r for r in filtered if r.primary_score(mode) is not None]
    skipped_missing = total_rows - len(usable)

    if not usable:
        expected_label = "NPS (0-10)" if mode == "nps10" else "CSAT (1-5)"
        raise ValueError(f"No valid {expected_label} survey scores found in the selected input.")

    all_nps_scores = [r.nps_score for r in usable if r.nps_score is not None]
    all_csat_scores = [r.csat_score for r in usable if r.csat_score is not None]

    overall_nps = compute_nps_metrics(all_nps_scores, z) if all_nps_scores else None
    overall_csat = compute_csat_metrics(all_csat_scores, z) if all_csat_scores else None

    period_delta, resolved_granularity, prev_key, curr_key = compute_period_delta(
        usable,
        mode=mode,
        period_mode=period_mode,
        z=z,
        previous_key=previous_period,
        current_key=current_period,
    )

    total_detractors, alerts = detect_detractor_alerts(
        usable,
        mode=mode,
        urgent_keywords=urgent_keywords,
        alert_max_score=alert_max_score,
    )

    group_by_norm = group_by.strip().lower()
    alert_source_indices = {a.source_index for a in alerts}
    alert_counts_by_group: Dict[str, int] = {}
    if group_by_norm != "none":
        for r in usable:
            if r.source_index in alert_source_indices:
                g_key = _group_value(r, group_by)
                alert_counts_by_group[g_key] = alert_counts_by_group.get(g_key, 0) + 1

    groups: List[GroupBreakdown] = []
    if group_by_norm != "none":
        grouped_map: Dict[str, List[SurveyResponse]] = {}
        for r in usable:
            g_name = _group_value(r, group_by)
            grouped_map.setdefault(g_name, []).append(r)

        unsorted_groups: List[GroupBreakdown] = []
        for g_name, g_rows in grouped_map.items():
            g_nps_scores = [r.nps_score for r in g_rows if r.nps_score is not None]
            g_csat_scores = [r.csat_score for r in g_rows if r.csat_score is not None]
            g_nps = compute_nps_metrics(g_nps_scores, z) if g_nps_scores else None
            g_csat = compute_csat_metrics(g_csat_scores, z) if g_csat_scores else None

            g_delta_nps: Optional[float] = None
            g_delta_csat_mean: Optional[float] = None
            g_delta_csat_pct: Optional[float] = None
            if period_delta.available and prev_key and curr_key:
                if resolved_granularity == "split":
                    g_dated = sorted(
                        [r for r in g_rows if r.response_date is not None],
                        key=lambda r: (r.response_date, r.source_index),  # type: ignore[arg-type]
                    )
                    if len(g_dated) >= 2:
                        mid = len(g_dated) // 2
                        p_snap = _snapshot_for_responses(prev_key, g_dated[:mid], z)
                        c_snap = _snapshot_for_responses(curr_key, g_dated[mid:], z)
                    else:
                        p_snap = c_snap = None
                else:
                    p_rows = [r for r in g_rows if _bucket_period_key(r, resolved_granularity) == prev_key]
                    c_rows = [r for r in g_rows if _bucket_period_key(r, resolved_granularity) == curr_key]
                    p_snap = _snapshot_for_responses(prev_key, p_rows, z) if p_rows else None
                    c_snap = _snapshot_for_responses(curr_key, c_rows, z) if c_rows else None

                if p_snap and c_snap:
                    if p_snap.nps is not None and c_snap.nps is not None:
                        g_delta_nps = c_snap.nps - p_snap.nps
                    if p_snap.csat_mean is not None and c_snap.csat_mean is not None:
                        g_delta_csat_mean = c_snap.csat_mean - p_snap.csat_mean
                    if p_snap.csat_pct is not None and c_snap.csat_pct is not None:
                        g_delta_csat_pct = c_snap.csat_pct - p_snap.csat_pct

            unsorted_groups.append(
                GroupBreakdown(
                    rank=0,
                    group=g_name,
                    count=len(g_rows),
                    low_sample=len(g_rows) < min_sample,
                    nps_metrics=g_nps,
                    csat_metrics=g_csat,
                    delta_nps=g_delta_nps,
                    delta_csat_mean=g_delta_csat_mean,
                    delta_csat_pct=g_delta_csat_pct,
                    urgent_alert_count=alert_counts_by_group.get(g_name, 0),
                )
            )

        def _group_sort_key(item: GroupBreakdown) -> Tuple[float, float, int, str]:
            if mode == "nps10":
                primary = item.nps_metrics.nps if (item.nps_metrics and item.nps_metrics.nps is not None) else -999.0
                secondary = (
                    item.nps_metrics.promoter_pct
                    if (item.nps_metrics and item.nps_metrics.promoter_pct is not None)
                    else -999.0
                )
            else:
                primary = item.csat_metrics.mean if (item.csat_metrics and item.csat_metrics.mean is not None) else -999.0
                secondary = (
                    item.csat_metrics.csat_pct
                    if (item.csat_metrics and item.csat_metrics.csat_pct is not None)
                    else -999.0
                )
            return (-round(primary, 6), -round(secondary, 6), -item.count, item.group.lower())

        unsorted_groups.sort(key=_group_sort_key)
        for idx, g in enumerate(unsorted_groups, start=1):
            groups.append(
                GroupBreakdown(
                    rank=idx,
                    group=g.group,
                    count=g.count,
                    low_sample=g.low_sample,
                    nps_metrics=g.nps_metrics,
                    csat_metrics=g.csat_metrics,
                    delta_nps=g.delta_nps,
                    delta_csat_mean=g.delta_csat_mean,
                    delta_csat_pct=g.delta_csat_pct,
                    urgent_alert_count=g.urgent_alert_count,
                )
            )

    # Threshold check
    threshold_check = ThresholdCheck(enabled=False, fail_below=None, metric=None, actual=None, passed=True)
    if fail_below is not None:
        if mode == "nps10":
            if not (-100.0 <= fail_below <= 100.0):
                raise ValueError("--fail-below in nps10 mode must be between -100 and 100.")
            assert overall_nps is not None and overall_nps.nps is not None
            metric_name = "nps"
            actual_val = overall_nps.nps
        else:
            if not (0.0 <= fail_below <= 100.0):
                raise ValueError("--fail-below in csat5 mode must be between 0 and 100.")
            assert overall_csat is not None
            use_pct = (
                csat_threshold_metric == "pct"
                or (csat_threshold_metric == "auto" and fail_below > 5.0)
            )
            if use_pct:
                assert overall_csat.csat_pct is not None
                metric_name = "csat_pct"
                actual_val = overall_csat.csat_pct
            else:
                assert overall_csat.mean is not None
                metric_name = "csat_mean"
                actual_val = overall_csat.mean
        threshold_check = ThresholdCheck(
            enabled=True,
            fail_below=fail_below,
            metric=metric_name,
            actual=actual_val,
            passed=actual_val >= fail_below,
        )

    return AnalysisReport(
        mode=mode,
        group_by=group_by,
        confidence=confidence,
        z_score=z,
        min_sample=min_sample,
        total_rows=total_rows,
        analyzed_count=len(usable),
        skipped_missing_score=skipped_missing,
        nps_metrics=overall_nps,
        csat_metrics=overall_csat,
        period_delta=period_delta,
        groups=tuple(groups),
        total_detractors=total_detractors,
        alerts=alerts,
        threshold_check=threshold_check,
    )


def _fmt_signed(val: Optional[float], suffix: str = "") -> str:
    if val is None:
        return "n/a"
    sign = "+" if val > 0 else ""
    return f"{sign}{val:.1f}{suffix}"


def _fmt_ci(ci: WilsonInterval, signed: bool = False, suffix: str = "") -> str:
    if ci.low is None or ci.high is None:
        return "n/a"
    if signed:
        l_s = f"+{ci.low:.1f}" if ci.low > 0 else f"{ci.low:.1f}"
        h_s = f"+{ci.high:.1f}" if ci.high > 0 else f"{ci.high:.1f}"
        return f"[{l_s}, {h_s}]{suffix}"
    return f"[{ci.low:.1f}%, {ci.high:.1f}%]"


def _format_table(headers: Sequence[str], rows: Sequence[Sequence[str]], right_align: Sequence[int] = ()) -> str:
    widths = [len(h) for h in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(cell))

    def _align(text: str, idx: int) -> str:
        return text.rjust(widths[idx]) if idx in right_align else text.ljust(widths[idx])

    header_line = " | ".join(_align(h, i) for i, h in enumerate(headers)).rstrip()
    sep_line = "-+-".join("-" * widths[i] for i in range(len(headers)))
    body_lines = [
        " | ".join(_align(cell, i) for i, cell in enumerate(row)).rstrip() for row in rows
    ]
    return "\n".join([header_line, sep_line, *body_lines])


def format_report_text(report: AnalysisReport, max_alerts: int = 20) -> str:
    """Format AnalysisReport as crisp aligned terminal tables and summaries."""
    lines: List[str] = []
    conf_pct = int(round(report.confidence * 100))

    status_tag = ""
    if report.threshold_check.enabled:
        status_tag = " | Status: PASS" if report.threshold_check.passed else " | Status: FAIL"

    skipped_note = f" (skipped missing: {report.skipped_missing_score})" if report.skipped_missing_score else ""
    lines.append(
        f"Survey Analysis ({report.mode.upper()}) | Responses: {report.analyzed_count}{skipped_note} | "
        f"CI: {conf_pct}% Wilson{status_tag}"
    )

    if report.mode == "nps10" and report.nps_metrics is not None:
        m = report.nps_metrics
        lines.append(
            f"Overall NPS: {_fmt_signed(m.nps)}  "
            f"(Wilson {conf_pct}% CI: {_fmt_ci(m.nps_ci, signed=True)}) | "
            f"Mean: {m.mean_score:.2f}/10 | Median: {m.median_score:.1f}/10"
        )
        lines.append(
            f"Splits     : "
            f"Promoters (9-10): {m.promoters}/{m.count} ({m.promoter_pct:.1f}% {_fmt_ci(m.promoter_ci)}) | "
            f"Passives (7-8): {m.passives}/{m.count} ({m.passive_pct:.1f}% {_fmt_ci(m.passive_ci)}) | "
            f"Detractors (0-6): {m.detractors}/{m.count} ({m.detractor_pct:.1f}% {_fmt_ci(m.detractor_ci)})"
        )

    if report.csat_metrics is not None:
        c = report.csat_metrics
        lines.append(
            f"CSAT (1-5) : Mean {c.mean:.2f}/5.00 | Median: {c.median:.1f} | "
            f"Top-2-Box (>=4): {c.satisfied_count}/{c.count} ({c.csat_pct:.1f}% "
            f"Wilson {conf_pct}% CI: {_fmt_ci(c.csat_pct_ci)})"
        )
        if report.mode == "csat5":
            d = c.distribution
            lines.append(
                f"Ratings    : 5*: {d['5']} | 4*: {d['4']} | 3*: {d['3']} | 2*: {d['2']} | 1*: {d['1']}"
            )

    if report.threshold_check.enabled:
        tc = report.threshold_check
        state = "PASS" if tc.passed else "FAIL"
        lines.append(
            f"CI Check   : {state} ({tc.metric} = {tc.actual:.2f}, required >= {tc.fail_below:.2f})"
        )

    pd = report.period_delta
    if pd.available and pd.previous and pd.current:
        if report.mode == "nps10":
            lines.append(
                f"Period Delta ({pd.granularity}): "
                f"{pd.previous.label} ({_fmt_signed(pd.previous.nps)}, n={pd.previous.count}) -> "
                f"{pd.current.label} ({_fmt_signed(pd.current.nps)}, n={pd.current.count}) | "
                f"Delta: {_fmt_signed(pd.delta_nps, ' pts')} [{pd.trend.upper()}]"
            )
        else:
            lines.append(
                f"Period Delta ({pd.granularity}): "
                f"{pd.previous.label} ({pd.previous.csat_mean:.2f}/5, n={pd.previous.count}) -> "
                f"{pd.current.label} ({pd.current.csat_mean:.2f}/5, n={pd.current.count}) | "
                f"Delta: {_fmt_signed(pd.delta_csat_mean, ' pts')} "
                f"(CSAT%: {_fmt_signed(pd.delta_csat_pct, ' pp')}) [{pd.trend.upper()}]"
            )
    else:
        lines.append("Period Delta: unavailable (need >= 2 distinct periods)")

    if report.groups:
        lines.append("")
        lines.append(f"Rankings by {report.group_by}:")
        if report.mode == "nps10":
            headers = (
                "Rank",
                report.group_by.title(),
                "N",
                "NPS",
                f"Wilson {conf_pct}% CI",
                "Prom%",
                "Pass%",
                "Det%",
                "CSAT",
                "Delta",
                "Alerts",
            )
            rows: List[List[str]] = []
            for g in report.groups:
                gm = g.nps_metrics
                assert gm is not None
                label = f"{g.group}*" if g.low_sample else g.group
                csat_cell = f"{g.csat_metrics.mean:.2f}" if g.csat_metrics and g.csat_metrics.mean is not None else "-"
                rows.append(
                    [
                        str(g.rank),
                        label,
                        str(g.count),
                        _fmt_signed(gm.nps),
                        _fmt_ci(gm.nps_ci, signed=True),
                        f"{gm.promoter_pct:.1f}%",
                        f"{gm.passive_pct:.1f}%",
                        f"{gm.detractor_pct:.1f}%",
                        csat_cell,
                        _fmt_signed(g.delta_nps),
                        str(g.urgent_alert_count),
                    ]
                )
            lines.append(_format_table(headers, rows, right_align=(0, 2, 3, 5, 6, 7, 8, 9, 10)))
        else:
            headers = (
                "Rank",
                report.group_by.title(),
                "N",
                "CSAT Avg",
                "Top2Box%",
                f"Wilson {conf_pct}% CI",
                "Sat",
                "Neut",
                "Dissat",
                "Delta",
                "Alerts",
            )
            rows = []
            for g in report.groups:
                cm = g.csat_metrics
                assert cm is not None
                label = f"{g.group}*" if g.low_sample else g.group
                rows.append(
                    [
                        str(g.rank),
                        label,
                        str(g.count),
                        f"{cm.mean:.2f}",
                        f"{cm.csat_pct:.1f}%",
                        _fmt_ci(cm.csat_pct_ci),
                        str(cm.satisfied_count),
                        str(cm.neutral_count),
                        str(cm.dissatisfied_count),
                        _fmt_signed(g.delta_csat_mean),
                        str(g.urgent_alert_count),
                    ]
                )
            lines.append(_format_table(headers, rows, right_align=(0, 2, 3, 4, 6, 7, 8, 9, 10)))

        if any(g.low_sample for g in report.groups):
            lines.append(f"* Low sample size (< {report.min_sample} responses); interpret wide Wilson CI cautiously.")

    lines.append("")
    if not report.alerts:
        lines.append(f"Urgent Detractor Alerts: 0 (out of {report.total_detractors} detractors)")
    else:
        shown = report.alerts[:max_alerts]
        lines.append(
            f"Urgent Detractor Alerts: {len(report.alerts)} flagged (out of {report.total_detractors} detractors)"
        )
        alert_headers = ("Severity", "ID", "Date", "Branch", "Staff", "Score", "Signals", "Comment")
        alert_rows: List[List[str]] = []
        for a in shown:
            signals = ", ".join(a.matched_keywords) if a.matched_keywords else "severe_score"
            comment_preview = a.comment if len(a.comment) <= 72 else a.comment[:69] + "..."
            alert_rows.append(
                [
                    a.severity,
                    a.response_id,
                    a.response_date.isoformat() if a.response_date else "-",
                    a.branch,
                    a.staff,
                    f"{a.score:g}",
                    signals,
                    comment_preview or "(no comment)",
                ]
            )
        lines.append(_format_table(alert_headers, alert_rows, right_align=(5,)))
        if len(report.alerts) > len(shown):
            lines.append(f"... and {len(report.alerts) - len(shown)} more alerts (use --format json for full list).")

    return "\n".join(lines)


def format_report_csv(report: AnalysisReport) -> str:
    """Format per-group breakdown (or OVERALL if ungrouped) as RFC-4180 CSV."""
    fieldnames = (
        "rank",
        "group",
        "responses",
        "low_sample",
        "nps",
        "nps_ci_low",
        "nps_ci_high",
        "promoters",
        "passives",
        "detractors",
        "promoter_pct",
        "passive_pct",
        "detractor_pct",
        "csat_mean",
        "csat_pct",
        "csat_ci_low",
        "csat_ci_high",
        "delta_nps",
        "delta_csat_mean",
        "urgent_alerts",
    )
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()

    def _row_from(
        rank: int,
        group: str,
        count: int,
        low_sample: bool,
        nps_m: Optional[NpsMetrics],
        csat_m: Optional[CsatMetrics],
        delta_nps: Optional[float],
        delta_csat_mean: Optional[float],
        urgent_alerts: int,
    ) -> Dict[str, Any]:
        return {
            "rank": rank,
            "group": group,
            "responses": count,
            "low_sample": str(low_sample).lower(),
            "nps": _r(nps_m.nps) if nps_m else "",
            "nps_ci_low": _r(nps_m.nps_ci.low) if (nps_m and nps_m.nps_ci.low is not None) else "",
            "nps_ci_high": _r(nps_m.nps_ci.high) if (nps_m and nps_m.nps_ci.high is not None) else "",
            "promoters": nps_m.promoters if nps_m else "",
            "passives": nps_m.passives if nps_m else "",
            "detractors": nps_m.detractors if nps_m else "",
            "promoter_pct": _r(nps_m.promoter_pct) if nps_m else "",
            "passive_pct": _r(nps_m.passive_pct) if nps_m else "",
            "detractor_pct": _r(nps_m.detractor_pct) if nps_m else "",
            "csat_mean": _r(csat_m.mean) if csat_m else "",
            "csat_pct": _r(csat_m.csat_pct) if csat_m else "",
            "csat_ci_low": _r(csat_m.csat_pct_ci.low) if (csat_m and csat_m.csat_pct_ci.low is not None) else "",
            "csat_ci_high": _r(csat_m.csat_pct_ci.high) if (csat_m and csat_m.csat_pct_ci.high is not None) else "",
            "delta_nps": _r(delta_nps) if delta_nps is not None else "",
            "delta_csat_mean": _r(delta_csat_mean) if delta_csat_mean is not None else "",
            "urgent_alerts": urgent_alerts,
        }

    if report.groups:
        for g in report.groups:
            writer.writerow(
                _row_from(
                    g.rank,
                    g.group,
                    g.count,
                    g.low_sample,
                    g.nps_metrics,
                    g.csat_metrics,
                    g.delta_nps,
                    g.delta_csat_mean,
                    g.urgent_alert_count,
                )
            )
    else:
        writer.writerow(
            _row_from(
                1,
                "OVERALL",
                report.analyzed_count,
                report.analyzed_count < report.min_sample,
                report.nps_metrics,
                report.csat_metrics,
                report.period_delta.delta_nps,
                report.period_delta.delta_csat_mean,
                len(report.alerts),
            )
        )

    return buf.getvalue().rstrip("\n")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="nps_analyzer.py",
        description=(
            "Zero-dependency Python CLI to compute NPS (%Promoters - %Detractors), CSAT (1-5), "
            "Wilson confidence intervals, per-branch/per-staff rankings, period-over-period deltas, "
            "and urgent detractor comment alerts from CSV/JSON survey exports."
        ),
    )
    parser.add_argument(
        "input_path",
        nargs="?",
        help="Path to survey CSV or JSON file (use '-' for stdin).",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_opt",
        metavar="PATH",
        help="Path to survey CSV or JSON file (alternative to positional argument).",
    )
    parser.add_argument(
        "--mode",
        choices=("nps10", "csat5"),
        default="nps10",
        help="Primary scoring mode: 'nps10' (0-10 NPS, default) or 'csat5' (1-5 CSAT).",
    )
    parser.add_argument(
        "--group-by",
        default="branch",
        metavar="DIM",
        help="Breakdown dimension: 'branch' (default), 'staff', 'channel', 'branch,staff', 'none', or custom column.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format: text (default), json, or csv.",
    )
    parser.add_argument(
        "--confidence",
        type=float,
        default=DEFAULT_CONFIDENCE,
        metavar="FLOAT",
        help="Confidence level for Wilson intervals (e.g. 0.90, 0.95, 0.99; default: 0.95).",
    )
    parser.add_argument(
        "--period",
        choices=("auto", "month", "week", "day", "split", "none"),
        default="auto",
        help="Period-over-period delta granularity (default: auto).",
    )
    parser.add_argument(
        "--previous-period",
        metavar="LABEL",
        help="Explicit baseline period label (e.g. 2026-03 or 2026-W14).",
    )
    parser.add_argument(
        "--current-period",
        metavar="LABEL",
        help="Explicit current period label (e.g. 2026-04 or 2026-W15).",
    )
    parser.add_argument(
        "--min-sample",
        type=int,
        default=DEFAULT_MIN_SAMPLE,
        metavar="N",
        help="Minimum sample size before flagging group as low-sample (default: 5).",
    )
    parser.add_argument(
        "--branch",
        metavar="NAME",
        help="Filter responses to a specific branch name (case-insensitive).",
    )
    parser.add_argument(
        "--staff",
        metavar="NAME",
        help="Filter responses to a specific staff/agent name (case-insensitive).",
    )
    parser.add_argument(
        "--urgent-keywords",
        metavar="CSV",
        help="Additional comma-separated keywords to trigger urgent detractor alerts.",
    )
    parser.add_argument(
        "--alert-max-score",
        type=float,
        default=None,
        metavar="SCORE",
        help="Maximum score considered for detractor alert scanning (default: 6 for nps10, 2 for csat5).",
    )
    parser.add_argument(
        "--max-alerts",
        type=int,
        default=20,
        metavar="N",
        help="Maximum detractor alerts printed in text output (default: 20).",
    )
    parser.add_argument(
        "--score-column",
        metavar="COL",
        help="Explicit primary score column name in input CSV/JSON.",
    )
    parser.add_argument(
        "--csat-column",
        metavar="COL",
        help="Explicit CSAT score column name in input CSV/JSON.",
    )
    parser.add_argument(
        "--date-column",
        metavar="COL",
        help="Explicit date/timestamp column name in input CSV/JSON.",
    )
    parser.add_argument(
        "--comment-column",
        metavar="COL",
        help="Explicit comment/feedback column name in input CSV/JSON.",
    )
    parser.add_argument(
        "--fail-below",
        type=float,
        default=None,
        metavar="SCORE",
        help="Exit with code 1 if overall NPS (nps10) or CSAT (csat5: <=5 checks mean, >5 checks Top2Box %%) falls below SCORE.",
    )
    parser.add_argument(
        "--csat-threshold-metric",
        choices=("auto", "mean", "pct"),
        default="auto",
        help="In csat5 mode, force --fail-below to evaluate 1-5 mean ('mean') or Top-2-Box %% ('pct').",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.input_path and args.input_opt and args.input_path != args.input_opt:
        print("Error: Provide input via positional path OR --input, not two different paths.", file=sys.stderr)
        return 2
    target_path = args.input_opt or args.input_path
    if not target_path:
        parser.print_usage(sys.stderr)
        print("Error: Survey file path is required (use '-' for stdin).", file=sys.stderr)
        return 2

    extra_keywords: List[str] = list(DEFAULT_URGENT_KEYWORDS)
    if args.urgent_keywords:
        extra_keywords.extend(kw.strip() for kw in args.urgent_keywords.split(",") if kw.strip())

    try:
        responses = load_surveys_file(
            Path(target_path) if target_path != "-" else Path("-"),
            mode=args.mode,
            score_column=args.score_column,
            csat_column=args.csat_column,
            date_column=args.date_column,
            comment_column=args.comment_column,
        )
        report = analyze_surveys(
            responses=responses,
            mode=args.mode,
            group_by=args.group_by,
            confidence=args.confidence,
            period_mode=args.period,
            previous_period=args.previous_period,
            current_period=args.current_period,
            min_sample=args.min_sample,
            branch_filter=args.branch,
            staff_filter=args.staff,
            urgent_keywords=extra_keywords,
            alert_max_score=args.alert_max_score,
            fail_below=args.fail_below,
            csat_threshold_metric=args.csat_threshold_metric,
        )
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2

    if args.format == "json":
        print(json.dumps(report.to_dict(), indent=2))
    elif args.format == "csv":
        print(format_report_csv(report))
    else:
        print(format_report_text(report, max_alerts=max(1, args.max_alerts)))

    if report.threshold_check.enabled and not report.threshold_check.passed:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
