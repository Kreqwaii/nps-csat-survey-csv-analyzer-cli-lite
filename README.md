# nps-csat-survey-csv-analyzer-cli-lite

Zero-dependency Python 3.9+ CLI to compute Net Promoter Score (`NPS = %Promoters - %Detractors`), CSAT (1–5), Wilson score confidence intervals, per-branch / per-staff breakdowns and rankings, period-over-period deltas, and urgent detractor alert flags from CSV or JSON survey exports.

---

## Features

- **Zero external dependencies**: Pure Python 3.9+ standard library (`argparse`, `csv`, `datetime`, `json`, `math`, `re`, `statistics`). No `pip install`, pandas, or virtual environment required.
- **Dual Scoring Modes (`--mode nps10|csat5`)**:
  - `nps10`: Computes NPS (`-100` to `+100`), Promoter (`9–10`), Passive (`7–8`), and Detractor (`0–6`) splits, plus side-by-side CSAT summary if a CSAT column is also present.
  - `csat5`: Computes mean & median 1–5 CSAT, star distribution (`5*` down to `1*`), and Top-2-Box satisfaction percentage (`>= 4`).
- **Wilson Score Confidence Intervals**:
  - Computes exact binomial Wilson score intervals for `%Promoters`, `%Passives`, `%Detractors`, and CSAT Top-2-Box `%` (default `95%`, configurable via `--confidence`).
  - Computes Newcombe-Wilson paired-trinomial confidence intervals (`[-100, +100]`) for NPS so small-sample and 100%/0% branches never collapse to zero-width bounds.
- **Per-Branch / Per-Staff Breakdowns (`--group-by`)**: Ranks groups by primary score (`branch`, `staff`, `channel`, composite `branch,staff`, `none`, or any custom CSV/JSON column) and flags small samples (`--min-sample`, default `5`).
- **Period-over-Period Delta (`--period`)**: Automatically compares latest vs previous month, ISO week, calendar day, explicit period label, or chronological half-split (`--period auto|month|week|day|split|none`), both overall and per group.
- **Urgent Detractor Alert Flags**: Scans detractors (`<= 6` for NPS, `<= 2` for CSAT) for severe ratings and escalation keywords (`refund`, `chargeback`, `overcharged`, `cancel`, `broken`, `crash`, `legal`, `rude`, `escalate`, etc.), tagging `CRITICAL` vs `HIGH` severity with matched signals while ignoring promoter false positives.
- **CI / Ops Automation (`--format`, `--fail-below`)**: Outputs aligned terminal tables (`text`), machine-readable `json`, or RFC-4180 `csv`, and returns exit code `1` when scores drop below `--fail-below`.

---

## Quick Start

### 1. Run NPS + Branch Breakdown + Detractor Alerts on `sample_surveys.csv`

```bash
py nps_analyzer.py sample_surveys.csv
```

Output:

```text
Survey Analysis (NPS10) | Responses: 38 | CI: 95% Wilson
Overall NPS: +23.7  (Wilson 95% CI: [-3.9, +47.4]) | Mean: 7.58/10 | Median: 8.5/10
Splits     : Promoters (9-10): 19/38 (50.0% [34.8%, 65.2%]) | Passives (7-8): 9/38 (23.7% [13.0%, 39.2%]) | Detractors (0-6): 10/38 (26.3% [15.0%, 42.0%])
CSAT (1-5) : Mean 3.97/5.00 | Median: 4.0 | Top-2-Box (>=4): 28/38 (73.7% Wilson 95% CI: [58.0%, 85.0%])
Period Delta (month): 2026-03 (+11.1, n=18) -> 2026-04 (+35.0, n=20) | Delta: +23.9 pts [IMPROVING]

Rankings by branch:
Rank | Branch   |  N |   NPS | Wilson 95% CI  | Prom% | Pass% |  Det% | CSAT | Delta | Alerts
-----+----------+----+-------+----------------+-------+-------+-------+------+-------+-------
   1 | Downtown | 11 | +63.6 | [+13.2, +86.0] | 72.7% | 18.2% |  9.1% | 4.55 | +43.3 |      0
   2 | Online   |  7 | +57.1 | [-8.7, +86.5]  | 71.4% | 14.3% | 14.3% | 4.57 | -16.7 |      0
   3 | Westside | 10 | +30.0 | [-19.7, +65.7] | 50.0% | 30.0% | 20.0% | 4.00 | +60.0 |      1
   4 | Airport  | 10 | -50.0 | [-77.6, -0.4]  | 10.0% | 30.0% | 60.0% | 2.90 | -20.0 |      4

Urgent Detractor Alerts: 5 flagged (out of 10 detractors)
Severity | ID       | Date       | Branch   | Staff      | Score | Signals                                     | Comment
---------+----------+------------+----------+------------+-------+---------------------------------------------+-------------------------------------------------------------------------
CRITICAL | SRV-1015 | 2026-03-25 | Airport  | Noah Vance |     1 | unacceptable, rude, supervisor, never again | Rude supervisor at counter, unacceptable attitude. Never again!
CRITICAL | SRV-1032 | 2026-04-11 | Airport  | Noah Vance |     2 | chargeback, overcharged, manager, escalate  | Overcharged $45 on invoice, escalate to billing manager or I will fil...
CRITICAL | SRV-1009 | 2026-03-20 | Westside | Maya Lin   |     2 | refund, double charged                      | Double charged on my card and nobody answered support. Need a refund ...
CRITICAL | SRV-1030 | 2026-04-05 | Airport  | Noah Vance |     3 | cancel, broken                              | Portal was broken all morning and customer service could not cancel m...
HIGH     | SRV-1013 | 2026-03-12 | Airport  | Chloe Shaw |     4 | crash                                       | App crashed twice during checkout and lost my queue spot.
```

### 2. Analyze CSAT (1–5), Staff Rankings, JSON/CSV Exports & CI Gates

```bash
# 1-5 CSAT mode with star distribution and Top-2-Box Wilson CI
py nps_analyzer.py sample_surveys.csv --mode csat5

# Rank individual staff members instead of branches
py nps_analyzer.py sample_surveys.csv --group-by staff

# Export structured JSON report (for dashboards or alerts pipelines)
py nps_analyzer.py sample_surveys.csv --format json

# Export branch rankings as RFC-4180 CSV
py nps_analyzer.py sample_surveys.csv --group-by branch --format csv

# CI / scheduled ops gate: exit 1 if overall NPS drops below +30
py nps_analyzer.py sample_surveys.csv --fail-below 30

# Read JSON survey array from stdin and extend urgent keywords
cat surveys.json | py nps_analyzer.py - --urgent-keywords "latency,outage,unresponsive"
```

---

## Input Formats

### CSV Exports
Supports UTF-8 (with or without BOM), comma, semicolon, or tab delimiters (detected from header row so commas inside customer comments never break parsing).

Flexible column aliases are detected automatically (or overridden via `--score-column`, `--csat-column`, `--date-column`, `--comment-column`):
- **ID**: `id`, `response_id`, `survey_id`, `ticket_id`, `ref`
- **Date**: `date`, `submitted_at`, `created_at`, `responded_at`, `timestamp`, `survey_date`, `day` (ISO-8601 dates/timestamps, `YYYY/MM/DD`, `DD.MM.YYYY`)
- **NPS (0–10)**: `nps`, `nps_score`, `nps_rating`, `recommend_score`, `promoter_score` (or generic `score`/`rating` in `nps10` mode)
- **CSAT (1–5)**: `csat`, `csat_score`, `csat_rating`, `satisfaction`, `satisfaction_score`, `stars` (or generic `score`/`rating` in `csat5` mode)
- **Dimensions**: `branch` (`location`, `store`, `clinic`, `office`, `site`, `region`, `department`), `staff` (`agent`, `employee`, `rep`, `specialist`), `channel` (`source`, `platform`, `touchpoint`)
- **Comment**: `comment`, `feedback`, `comments`, `review`, `text`, `note`, `verbatim`, `reason`

### JSON Exports
Accepts either a top-level array of response objects or an object wrapping records under `surveys`, `responses`, `rows`, `items`, `data`, or `records`:

```json
{
  "responses": [
    {
      "id": "SRV-201",
      "submitted_at": "2026-04-12T10:30:00Z",
      "branch": "Downtown",
      "staff": "Elena Ross",
      "nps_score": 10,
      "csat_score": 5,
      "comment": "Fast service and helpful staff."
    }
  ]
}
```

---

## CLI Options & Exit Codes

| Flag | Default | Description |
| :--- | :--- | :--- |
| `input_path` / `-i`, `--input` | — | Path to survey `.csv` or `.json` file (`-` for stdin) |
| `--mode` | `nps10` | Scoring mode: `nps10` (0–10 NPS) or `csat5` (1–5 CSAT) |
| `--group-by` | `branch` | Breakdown dimension: `branch`, `staff`, `channel`, `branch,staff`, `none`, or custom column |
| `--format` | `text` | Output format: `text`, `json`, or `csv` |
| `--confidence` | `0.95` | Wilson confidence interval level (`0.90`, `0.95`, `0.99`, etc.) |
| `--period` | `auto` | Period comparison granularity: `auto`, `month`, `week`, `day`, `split`, or `none` |
| `--previous-period` / `--current-period` | `auto` | Explicit period labels to compare (e.g. `2026-03` vs `2026-04`) |
| `--min-sample` | `5` | Minimum responses before flagging group as low-sample (`*`) |
| `--branch` / `--staff` | `—` | Case-insensitive filter by branch or staff member |
| `--urgent-keywords` | `—` | Extra comma-separated keywords for urgent detractor flagging |
| `--alert-max-score` | `6` / `2` | Max score scanned for detractor alerts (`nps10` / `csat5`) |
| `--fail-below` | `—` | Exit `1` if overall NPS (`nps10`) or CSAT (`<=5`: mean, `>5`: Top-2-Box `%`) `< SCORE` |
| `--csat-threshold-metric` | `auto` | Force `--fail-below` in `csat5` mode to check `mean` (`1–5`) or `pct` (`0–100%`) |

**Exit Codes**: `0` = success (within threshold policy), `1` = threshold failed under `--fail-below`, `2` = invalid arguments, unreadable input, or out-of-range scores.

---

## Lite CLI vs Full Telegram Client Feedback, NPS Survey & Review Routing Bot + Light Web Dashboard Pro

This open-source Lite CLI analyzes historical CSV/JSON survey exports offline from your terminal or CI pipeline. If you want an automated **live Telegram feedback collection, positive review routing, instant detractor escalation, and light web dashboard** ready to deploy for branches and service teams, compare both editions below:

| Capability | Lite CLI (Free MIT) | Full Telegram Client Feedback, NPS Survey & Review Routing Bot + Light Web Dashboard Pro |
| :--- | :--- | :--- |
| Offline CSV/JSON NPS & CSAT calculation + Wilson CIs | Yes | Yes (live analytics + exportable CSV/JSON reports) |
| Per-branch / per-staff rankings & period-over-period delta | Yes | Yes (real-time rolling KPIs & branch comparisons) |
| Urgent detractor comment detection | Offline batch scan | Instant real-time Telegram manager/admin escalation alerts |
| Interactive Telegram survey collection flow (NPS/CSAT + comment) | — | Yes (inline button rating flow + QR/deep-link per branch/staff) |
| Smart positive review routing (Promoters -> Google/Maps/2GIS/Trust links) | — | Yes (automated review funnel for 9–10 / 4–5 star clients) |
| Detractor recovery workflow & ticket status tracking | — | Yes (private apology/follow-up capture before public review) |
| Light Web Dashboard + SQLite persistence | — | Yes (browser KPI overview, branch leaderboard, alert inbox) |
| Production deployment templates & configuration guide | — | Yes |

- **Get the Pro Edition on FunPay:** https://funpay.com/lots/offer?id=82604181

---

## Running Tests

```bash
py -m unittest discover -s tests -v
```

---

## License

MIT License — see [LICENSE](LICENSE).
