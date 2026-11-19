"""Offline unit and CLI integration tests for nps_analyzer.py."""

from __future__ import annotations

import contextlib
import csv
import io
import json
import sys
import tempfile
import unittest
from datetime import date
from pathlib import Path

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.insert(0, str(ROOT_DIR))

import nps_analyzer as na  # noqa: E402

SAMPLE_CSV = ROOT_DIR / "sample_surveys.csv"


class TestStatisticalMath(unittest.TestCase):
    def test_wilson_proportion_interval_standard_and_extremes(self) -> None:
        z95 = na.confidence_to_z(0.95)
        ci_mid = na.wilson_proportion_interval(50, 100, z95)
        self.assertIsNotNone(ci_mid.low)
        self.assertIsNotNone(ci_mid.high)
        assert ci_mid.low is not None and ci_mid.high is not None
        self.assertAlmostEqual(ci_mid.low, 40.38, places=2)
        self.assertAlmostEqual(ci_mid.high, 59.62, places=2)

        # Boundary cases must remain non-degenerate and bounded within [0, 100]
        ci_zero = na.wilson_proportion_interval(0, 10, z95)
        assert ci_zero.low is not None and ci_zero.high is not None
        self.assertEqual(ci_zero.low, 0.0)
        self.assertGreater(ci_zero.high, 20.0)
        self.assertLess(ci_zero.high, 35.0)

        ci_full = na.wilson_proportion_interval(10, 10, z95)
        assert ci_full.low is not None and ci_full.high is not None
        self.assertGreater(ci_full.low, 70.0)
        self.assertEqual(ci_full.high, 100.0)

        ci_empty = na.wilson_proportion_interval(0, 0, z95)
        self.assertIsNone(ci_empty.low)
        self.assertIsNone(ci_empty.high)

    def test_wilson_nps_interval_bounds_and_confidence_scaling(self) -> None:
        z90 = na.confidence_to_z(0.90)
        z95 = na.confidence_to_z(0.95)
        z99 = na.confidence_to_z(0.99)

        ci_all_prom = na.wilson_nps_interval(10, 0, 0, z95)
        assert ci_all_prom.low is not None and ci_all_prom.high is not None
        self.assertLess(ci_all_prom.low, 100.0)
        self.assertGreater(ci_all_prom.low, 50.0)
        self.assertEqual(ci_all_prom.high, 100.0)

        ci_all_det = na.wilson_nps_interval(0, 0, 10, z95)
        assert ci_all_det.low is not None and ci_all_det.high is not None
        self.assertEqual(ci_all_det.low, -100.0)
        self.assertGreater(ci_all_det.high, -100.0)

        ci90 = na.wilson_nps_interval(50, 30, 20, z90)
        ci95 = na.wilson_nps_interval(50, 30, 20, z95)
        ci99 = na.wilson_nps_interval(50, 30, 20, z99)
        assert ci90.low is not None and ci95.low is not None and ci99.low is not None
        assert ci90.high is not None and ci95.high is not None and ci99.high is not None
        self.assertLess(ci99.low, ci95.low)
        self.assertLess(ci95.low, ci90.low)
        self.assertGreater(ci99.high, ci95.high)
        self.assertGreater(ci95.high, ci90.high)

    def test_confidence_to_z_values_and_validation(self) -> None:
        self.assertAlmostEqual(na.confidence_to_z(0.95), 1.959964, places=5)
        self.assertAlmostEqual(na.confidence_to_z(0.90), 1.644854, places=5)
        self.assertAlmostEqual(na.confidence_to_z(0.99), 2.575829, places=5)
        for bad in (0.0, 1.0, -0.5, 1.2):
            with self.assertRaises(ValueError):
                na.confidence_to_z(bad)

    def test_compute_nps_metrics_splits_and_boundaries(self) -> None:
        scores = [10, 9, 8, 7, 6, 5, 0, 9, 10, 6]
        z95 = na.confidence_to_z(0.95)
        m = na.compute_nps_metrics(scores, z95)
        assert m is not None
        self.assertEqual(m.count, 10)
        self.assertEqual(m.promoters, 4)
        self.assertEqual(m.passives, 2)
        self.assertEqual(m.detractors, 4)
        self.assertAlmostEqual(m.promoter_pct or 0.0, 40.0)
        self.assertAlmostEqual(m.passive_pct or 0.0, 20.0)
        self.assertAlmostEqual(m.detractor_pct or 0.0, 40.0)
        self.assertAlmostEqual(m.nps or 0.0, 0.0)
        self.assertAlmostEqual(m.mean_score or 0.0, 7.0)
        self.assertAlmostEqual(m.median_score or 0.0, 7.5)

    def test_compute_csat_metrics_top2box_and_distribution(self) -> None:
        scores = [5, 5, 4, 4, 3, 2, 1, 5]
        z95 = na.confidence_to_z(0.95)
        c = na.compute_csat_metrics(scores, z95)
        assert c is not None
        self.assertEqual(c.count, 8)
        self.assertEqual(c.satisfied_count, 5)
        self.assertEqual(c.neutral_count, 1)
        self.assertEqual(c.dissatisfied_count, 2)
        self.assertAlmostEqual(c.csat_pct or 0.0, 62.5)
        self.assertAlmostEqual(c.mean or 0.0, 3.625)
        self.assertAlmostEqual(c.median or 0.0, 4.0)
        self.assertEqual(c.distribution, {"1": 1, "2": 1, "3": 1, "4": 2, "5": 3})


class TestParsingAndValidation(unittest.TestCase):
    def test_load_csv_utf8_bom_semicolon_and_quoted_commas(self) -> None:
        raw = (
            "\ufeffid;date;branch;staff;nps_score;comment\n"
            'R1;2026-04-01;North;Alice;10;"Great service, fast delivery"\n'
            ";;;;;\n"
            'R2;2026-04-02;North;Bob;5;"Slow checkout, needed help"\n'
        )
        rows = na.load_surveys_text(raw, mode="nps10", format_hint="csv")
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0].response_id, "R1")
        self.assertEqual(rows[0].comment, "Great service, fast delivery")
        self.assertEqual(rows[1].nps_score, 5.0)

    def test_load_json_array_and_wrapped_object_with_iso_z(self) -> None:
        arr_json = json.dumps(
            [
                {
                    "response_id": "J1",
                    "submitted_at": "2026-04-10T09:15:00Z",
                    "location": "East",
                    "agent": "Sam",
                    "nps": 9,
                    "feedback": "Awesome!",
                }
            ]
        )
        rows1 = na.load_surveys_text(arr_json, mode="nps10")
        self.assertEqual(len(rows1), 1)
        self.assertEqual(rows1[0].response_date, date(2026, 4, 10))
        self.assertEqual(rows1[0].branch, "East")
        self.assertEqual(rows1[0].staff, "Sam")

        wrapped_json = json.dumps(
            {
                "responses": [
                    {"id": "J2", "date": "2026-04-11 16:40:00+02:00", "csat": 4, "branch": "West"}
                ]
            }
        )
        rows2 = na.load_surveys_text(wrapped_json, mode="csat5")
        self.assertEqual(len(rows2), 1)
        self.assertEqual(rows2[0].response_date, date(2026, 4, 11))
        self.assertEqual(rows2[0].csat_score, 4.0)

    def test_generic_score_alias_and_missing_score_skipping(self) -> None:
        csv_text = (
            "id,date,branch,score\n"
            "S1,2026-04-01,Main,5\n"
            "S2,2026-04-02,Main,\n"
            "S3,2026-04-03,Main,4\n"
        )
        rows_csat = na.load_surveys_text(csv_text, mode="csat5")
        report_csat = na.analyze_surveys(rows_csat, mode="csat5")
        self.assertEqual(report_csat.total_rows, 3)
        self.assertEqual(report_csat.analyzed_count, 2)
        self.assertEqual(report_csat.skipped_missing_score, 1)
        assert report_csat.csat_metrics is not None
        self.assertAlmostEqual(report_csat.csat_metrics.mean or 0.0, 4.5)

        rows_nps = na.load_surveys_text(csv_text, mode="nps10")
        report_nps = na.analyze_surveys(rows_nps, mode="nps10")
        assert report_nps.nps_metrics is not None
        self.assertAlmostEqual(report_nps.nps_metrics.nps or 0.0, -100.0)

    def test_custom_column_overrides(self) -> None:
        csv_text = (
            "ticket,when,team,loyalty,happiness,verbatim_note\n"
            "T1,2026-04-01,Ops,10,5,Zero issues\n"
            "T2,2026-04-02,Ops,2,1,Billing error refund requested\n"
        )
        rows = na.load_surveys_text(
            csv_text,
            mode="nps10",
            score_column="loyalty",
            csat_column="happiness",
            date_column="when",
            comment_column="verbatim_note",
        )
        report = na.analyze_surveys(rows, mode="nps10", group_by="team")
        self.assertEqual(report.analyzed_count, 2)
        assert report.nps_metrics is not None and report.csat_metrics is not None
        self.assertAlmostEqual(report.nps_metrics.nps or 0.0, 0.0)
        self.assertAlmostEqual(report.csat_metrics.mean or 0.0, 3.0)
        self.assertEqual(len(report.alerts), 1)
        self.assertEqual(report.groups[0].group, "Ops")

    def test_date_format_variants_and_errors(self) -> None:
        self.assertEqual(na.parse_survey_date("2026-04-15"), date(2026, 4, 15))
        self.assertEqual(na.parse_survey_date("2026/04/15"), date(2026, 4, 15))
        self.assertEqual(na.parse_survey_date("15.04.2026"), date(2026, 4, 15))
        self.assertEqual(na.parse_survey_date("2026-04-15T12:34:56Z"), date(2026, 4, 15))
        self.assertIsNone(na.parse_survey_date("   "))
        with self.assertRaises(ValueError):
            na.parse_survey_date("not-a-date")

    def test_invalid_inputs_and_out_of_range_scores(self) -> None:
        with self.assertRaises(ValueError):
            na.load_surveys_text("   \n  ")
        with self.assertRaises(ValueError):
            na.load_surveys_text("{bad json")
        with self.assertRaises(ValueError):
            na.load_surveys_text('{"foo": []}')
        with self.assertRaises(ValueError):
            na.load_surveys_text("id,nps\n1,11\n", mode="nps10")
        with self.assertRaises(ValueError):
            na.load_surveys_text("id,nps\n1,-1\n", mode="nps10")
        with self.assertRaises(ValueError):
            na.load_surveys_text("id,csat\n1,6\n", mode="csat5")
        with self.assertRaises(ValueError):
            na.load_surveys_text("id,nps\n1,abc\n", mode="nps10")


class TestGroupingPeriodsAndAlerts(unittest.TestCase):
    def test_group_by_dimensions_ranking_and_low_sample(self) -> None:
        csv_text = (
            "id,date,branch,staff,channel,nps\n"
            "1,2026-04-01,Alpha,Ann,Web,10\n"
            "2,2026-04-02,Alpha,Ann,Web,9\n"
            "3,2026-04-03,Alpha,Ben,Telegram,7\n"
            "4,2026-04-04,Beta,Ann,Web,10\n"
            "5,2026-04-05,Beta,Ben,Email,2\n"
        )
        rows = na.load_surveys_text(csv_text, mode="nps10")
        rep_branch = na.analyze_surveys(rows, mode="nps10", group_by="branch", min_sample=3)
        self.assertEqual([g.group for g in rep_branch.groups], ["Alpha", "Beta"])
        self.assertEqual([g.rank for g in rep_branch.groups], [1, 2])
        self.assertFalse(rep_branch.groups[0].low_sample)
        self.assertTrue(rep_branch.groups[1].low_sample)

        rep_comp = na.analyze_surveys(rows, mode="nps10", group_by="branch,staff")
        self.assertEqual(len(rep_comp.groups), 4)
        self.assertEqual(rep_comp.groups[0].group, "Alpha / Ann")

        rep_none = na.analyze_surveys(rows, mode="nps10", group_by="none")
        self.assertEqual(len(rep_none.groups), 0)

    def test_branch_and_staff_filters(self) -> None:
        rows = na.load_surveys_file(SAMPLE_CSV, mode="nps10")
        rep = na.analyze_surveys(rows, mode="nps10", branch_filter="downtown", staff_filter="elena ross")
        self.assertEqual(rep.analyzed_count, 6)
        assert rep.nps_metrics is not None
        self.assertAlmostEqual(rep.nps_metrics.nps or 0.0, 83.33, places=2)
        self.assertEqual(len(rep.alerts), 0)

    def test_period_delta_month_week_split_and_explicit(self) -> None:
        csv_text = (
            "id,date,branch,nps,csat\n"
            "1,2026-03-02,HQ,6,3\n"
            "2,2026-03-09,HQ,7,4\n"
            "3,2026-04-02,HQ,10,5\n"
            "4,2026-04-16,HQ,10,5\n"
        )
        rows = na.load_surveys_text(csv_text, mode="nps10")

        rep_month = na.analyze_surveys(rows, mode="nps10", period_mode="month")
        self.assertTrue(rep_month.period_delta.available)
        self.assertEqual(rep_month.period_delta.granularity, "month")
        assert rep_month.period_delta.previous is not None and rep_month.period_delta.current is not None
        self.assertEqual(rep_month.period_delta.previous.label, "2026-03")
        self.assertEqual(rep_month.period_delta.current.label, "2026-04")
        self.assertAlmostEqual(rep_month.period_delta.delta_nps or 0.0, 150.0)
        self.assertEqual(rep_month.period_delta.trend, "improving")

        rep_week = na.analyze_surveys(rows, mode="nps10", period_mode="week")
        self.assertTrue(rep_week.period_delta.available)
        self.assertEqual(rep_week.period_delta.granularity, "week")

        rep_split = na.analyze_surveys(rows, mode="nps10", period_mode="split")
        self.assertTrue(rep_split.period_delta.available)
        self.assertEqual(rep_split.period_delta.granularity, "split")

    def test_period_delta_unavailable_for_single_period(self) -> None:
        csv_text = "id,date,branch,nps\n1,2026-04-01,HQ,9\n2,2026-04-01,HQ,10\n"
        rows = na.load_surveys_text(csv_text, mode="nps10")
        rep = na.analyze_surveys(rows, mode="nps10", period_mode="auto")
        self.assertFalse(rep.period_delta.available)
        self.assertEqual(rep.period_delta.trend, "unavailable")

    def test_urgent_detractor_alert_logic_and_no_promoter_false_positives(self) -> None:
        csv_text = (
            "id,date,branch,staff,nps,comment\n"
            "A1,2026-04-01,HQ,Ann,10,Never going to cancel this plan!\n"
            "A2,2026-04-01,HQ,Ann,6,Wait time was slightly long.\n"
            "A3,2026-04-01,HQ,Ann,5,Want a refund after billing error.\n"
            "A4,2026-04-01,HQ,Ann,3,Portal was slow.\n"
            "A5,2026-04-01,HQ,Ann,2,System crashed and I will cancel.\n"
            "A6,2026-04-01,HQ,Ann,0,Bad.\n"
        )
        rows = na.load_surveys_text(csv_text, mode="nps10")
        total_det, alerts = na.detect_detractor_alerts(rows, mode="nps10")
        self.assertEqual(total_det, 5)  # A1 (score 10) is Promoter, ignored!
        by_id = {a.response_id: a for a in alerts}
        self.assertNotIn("A1", by_id)
        self.assertNotIn("A2", by_id)  # Mild detractor without urgent keyword
        self.assertEqual(by_id["A3"].severity, "HIGH")
        self.assertEqual(by_id["A3"].matched_keywords, ("refund", "billing error"))
        self.assertEqual(by_id["A4"].severity, "HIGH")
        self.assertEqual(by_id["A5"].severity, "CRITICAL")
        self.assertEqual(by_id["A6"].severity, "CRITICAL")

    def test_custom_urgent_keywords_and_alert_max_score(self) -> None:
        csv_text = (
            "id,date,branch,nps,comment\n"
            "K1,2026-04-01,HQ,7,Latency issue during peak hours\n"
            "K2,2026-04-01,HQ,6,Minor typo\n"
        )
        rows = na.load_surveys_text(csv_text, mode="nps10")
        _, alerts = na.detect_detractor_alerts(
            rows,
            mode="nps10",
            urgent_keywords=[*na.DEFAULT_URGENT_KEYWORDS, "latency"],
            alert_max_score=7.0,
        )
        self.assertEqual(len(alerts), 1)
        self.assertEqual(alerts[0].response_id, "K1")
        self.assertEqual(alerts[0].matched_keywords, ("latency",))


class TestCLIIntegration(unittest.TestCase):
    def _run_cli(self, argv: list[str], stdin_text: str = "") -> tuple[int, str, str]:
        out = io.StringIO()
        err = io.StringIO()
        old_stdin = sys.stdin
        try:
            sys.stdin = io.StringIO(stdin_text)
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                code = na.main(argv)
        finally:
            sys.stdin = old_stdin
        return code, out.getvalue(), err.getvalue()

    def test_fail_below_nps10_pass_and_fail(self) -> None:
        code_pass, out_pass, _ = self._run_cli([str(SAMPLE_CSV), "--fail-below", "20"])
        self.assertEqual(code_pass, 0)
        self.assertIn("Status: PASS", out_pass)

        code_fail, out_fail, _ = self._run_cli([str(SAMPLE_CSV), "--fail-below", "30"])
        self.assertEqual(code_fail, 1)
        self.assertIn("Status: FAIL", out_fail)
        self.assertIn("Urgent Detractor Alerts", out_fail)

    def test_fail_below_csat5_mean_and_pct(self) -> None:
        # Sample CSAT mean is 3.97, Top-2-Box % is 73.68%
        c1, out1, _ = self._run_cli([str(SAMPLE_CSV), "--mode", "csat5", "--fail-below", "3.9"])
        self.assertEqual(c1, 0)
        self.assertIn("csat_mean", out1)

        c2, out2, _ = self._run_cli([str(SAMPLE_CSV), "--mode", "csat5", "--fail-below", "4.2"])
        self.assertEqual(c2, 1)
        self.assertIn("Status: FAIL", out2)

        c3, out3, _ = self._run_cli([str(SAMPLE_CSV), "--mode", "csat5", "--fail-below", "80"])
        self.assertEqual(c3, 1)
        self.assertIn("csat_pct", out3)

        c4, out4, _ = self._run_cli(
            [str(SAMPLE_CSV), "--mode", "csat5", "--fail-below", "4.0", "--csat-threshold-metric", "pct"]
        )
        self.assertEqual(c4, 0)
        self.assertIn("csat_pct", out4)

    def test_cli_json_format_on_sample_csv(self) -> None:
        code, out, err = self._run_cli([str(SAMPLE_CSV), "--format", "json"])
        self.assertEqual(code, 0, msg=err)
        payload = json.loads(out)
        self.assertEqual(payload["mode"], "nps10")
        self.assertEqual(payload["summary"]["analyzed_count"], 38)
        self.assertAlmostEqual(payload["summary"]["nps"]["nps"], 23.68, places=2)
        self.assertAlmostEqual(payload["summary"]["csat"]["mean"], 3.97, places=2)
        self.assertTrue(payload["period_delta"]["available"])
        self.assertEqual(payload["period_delta"]["previous_period"]["label"], "2026-03")
        self.assertEqual(payload["period_delta"]["current_period"]["label"], "2026-04")
        self.assertEqual([g["group"] for g in payload["groups"]], ["Downtown", "Online", "Westside", "Airport"])
        self.assertEqual(payload["detractor_alerts"]["urgent_count"], 5)

    def test_cli_csv_format_grouped_and_ungrouped(self) -> None:
        code_grp, out_grp, _ = self._run_cli([str(SAMPLE_CSV), "--group-by", "branch", "--format", "csv"])
        self.assertEqual(code_grp, 0)
        rows_grp = list(csv.DictReader(io.StringIO(out_grp)))
        self.assertEqual(len(rows_grp), 4)
        self.assertEqual(rows_grp[0]["group"], "Downtown")
        self.assertEqual(rows_grp[0]["rank"], "1")
        self.assertEqual(rows_grp[3]["group"], "Airport")
        self.assertEqual(rows_grp[3]["urgent_alerts"], "4")

        code_none, out_none, _ = self._run_cli([str(SAMPLE_CSV), "--group-by", "none", "--format", "csv"])
        self.assertEqual(code_none, 0)
        rows_none = list(csv.DictReader(io.StringIO(out_none)))
        self.assertEqual(len(rows_none), 1)
        self.assertEqual(rows_none[0]["group"], "OVERALL")
        self.assertEqual(rows_none[0]["responses"], "38")

    def test_cli_stdin_json_and_error_exit_codes(self) -> None:
        stdin_payload = json.dumps(
            [
                {"id": "X1", "date": "2026-04-01", "branch": "Central", "nps": 10, "comment": "Great"},
                {"id": "X2", "date": "2026-04-02", "branch": "Central", "nps": 1, "comment": "Need refund"},
            ]
        )
        code, out, err = self._run_cli(["-", "--format", "json"], stdin_text=stdin_payload)
        self.assertEqual(code, 0, msg=err)
        parsed = json.loads(out)
        self.assertEqual(parsed["summary"]["nps"]["nps"], 0.0)
        self.assertEqual(parsed["detractor_alerts"]["urgent_count"], 1)

        with tempfile.TemporaryDirectory() as tmpdir:
            missing = Path(tmpdir) / "does_not_exist.csv"
            c_missing, _, e_missing = self._run_cli([str(missing)])
            self.assertEqual(c_missing, 2)
            self.assertIn("not found", e_missing.lower())

        c_no_arg, _, e_no_arg = self._run_cli([])
        self.assertEqual(c_no_arg, 2)
        self.assertIn("required", e_no_arg.lower())

        c_conflict, _, _ = self._run_cli([str(SAMPLE_CSV), "--input", "other.csv"])
        self.assertEqual(c_conflict, 2)


if __name__ == "__main__":
    unittest.main()
