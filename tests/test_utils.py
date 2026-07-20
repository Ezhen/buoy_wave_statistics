"""
Regression tests for utils.get_dt_hours - the era-aware sampling-
interval accessor built to fix a live correctness gap in Stages 11b,
13, and 24: each read `sampling_interval_hours` directly from Stage
01's load_summary.json, which only reports the DOMINANT era's value -
silently wrong for any segment actually analyzed from a non-dominant
era on a multi-era buoy.

Run with: pytest tests/test_utils.py -v
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import utils  # noqa: E402


def test_single_era_summary_falls_back_to_top_level_field():
    """The common case (11 of 19 buoys) - no 'eras' list, or a
    single-entry one - must behave exactly as the pre-fix code did."""
    summary_no_eras = {"sampling_interval_hours": 0.5}
    assert utils.get_dt_hours(summary_no_eras, "2020-01-01") == 0.5

    summary_one_era = {
        "sampling_interval_hours": 1.0,
        "eras": [{"start": "2000-01-01", "end": "2020-01-01", "inferred_freq_hours": 1.0}],
    }
    assert utils.get_dt_hours(summary_one_era, "2010-06-15") == 1.0


def test_multi_era_returns_correct_era_specific_value():
    """The core fix: a segment sitting in the non-dominant era must get
    THAT era's dt, not the dominant era's."""
    summary = {
        "sampling_interval_hours": 0.25,  # dominant era's value
        "eras": [
            {"start": "2010-01-01", "end": "2019-01-01", "inferred_freq_hours": 0.25},
            {"start": "2019-01-01", "end": "2024-12-31", "inferred_freq_hours": 0.5},
        ],
    }
    # segment starting inside the dominant era
    assert utils.get_dt_hours(summary, "2012-05-01") == 0.25
    # segment starting inside the NON-dominant era - this is the case
    # the old dominant-era-only lookup got silently wrong
    assert utils.get_dt_hours(summary, "2020-03-01") == 0.5


def test_era_aware_gap_detection_catches_within_era_gap():
    """A genuine gap WITHIN a single era must still be caught using
    that era's own dt - not missed by an inappropriately coarse shared
    threshold."""
    e1a = pd.date_range("2010-01-01", "2010-06-01 00:00", freq="15min")
    e1b = pd.date_range("2010-06-01 00:40", "2011-01-01", freq="15min")  # 40min gap
    idx = e1a.append(e1b)
    s = pd.Series(np.arange(len(idx)), index=idx)
    summary = {"eras": [{"start": str(e1a[0]), "end": str(e1b[-1]), "inferred_freq_hours": 0.25}]}

    segs = utils.segments_by_time_gap_era_aware(s, summary)
    assert len(segs) == 2
    assert segs[0].index[-1] == e1a[-1]
    assert segs[1].index[0] == e1b[0]


def test_era_aware_gap_detection_always_breaks_at_era_boundary():
    """The core risk this function exists to avoid: a small transition
    seam that would NOT clear a naive shared-threshold check (using the
    coarser era's dt) must still be caught, since gap-detection never
    even evaluates data across two eras in the same call."""
    e1 = pd.date_range("2015-01-01", "2015-06-01 00:00", freq="15min")
    # 20-minute seam - would NOT exceed 1.5 * 30min = 45min if era2's
    # dt were (wrongly) used as a single shared threshold
    e2 = pd.date_range("2015-06-01 00:20", "2016-01-01", freq="30min")
    idx = e1.append(e2)
    s = pd.Series(np.arange(len(idx)), index=idx)
    summary = {
        "eras": [
            {"start": str(e1[0]), "end": str(e1[-1]), "inferred_freq_hours": 0.25},
            {"start": str(e2[0]), "end": str(e2[-1]), "inferred_freq_hours": 0.5},
        ]
    }

    segs = utils.segments_by_time_gap_era_aware(s, summary)
    assert len(segs) == 2
    assert segs[0].index[-1] == e1[-1]
    assert segs[1].index[0] == e2[0]


def test_era_aware_gap_detection_single_era_fallback():
    """Single-era summaries must behave exactly like the plain
    (non-era-aware) function - the common case, unchanged."""
    idx = pd.date_range("2010-01-01", periods=1000, freq="h")
    s = pd.Series(np.arange(1000), index=idx)
    summary = {"sampling_interval_hours": 1.0}  # no 'eras' key at all

    segs_era_aware = utils.segments_by_time_gap_era_aware(s, summary)
    segs_plain = utils.segments_by_time_gap(s, 1.0)
    assert len(segs_era_aware) == len(segs_plain) == 1


def test_emit_warning_structures_and_prints(capsys):
    """The structured-warnings convention (2026-07-20 session): one call
    site produces both the console message and a machine-readable entry,
    instead of stages maintaining two separate code paths."""
    warnings = []
    entry = utils.emit_warning(warnings, "warning", "test_code", "a test message", value=42)

    captured = capsys.readouterr()
    assert "WARNING: a test message" in captured.out
    assert len(warnings) == 1
    assert warnings[0] == entry
    assert entry["severity"] == "warning"
    assert entry["code"] == "test_code"
    assert entry["context"] == {"value": 42}


def test_emit_warning_omits_context_key_when_empty():
    """No spurious empty 'context' key when no extra values are passed."""
    warnings = []
    entry = utils.emit_warning(warnings, "note", "code", "message")
    assert "context" not in entry


def test_raises_on_timestamp_outside_every_era():
    """A silently wrong dt_hours would corrupt every downstream lag-
    based calculation - fail loudly instead of guessing. Needs 2+ eras
    in the fixture: a single-era summary correctly takes the fallback
    branch (tested above) rather than ever reaching this check."""
    summary = {
        "sampling_interval_hours": 0.25,
        "eras": [
            {"start": "2010-01-01", "end": "2015-01-01", "inferred_freq_hours": 0.25},
            {"start": "2016-01-01", "end": "2019-01-01", "inferred_freq_hours": 0.5},
        ],
    }
    with pytest.raises(ValueError):
        utils.get_dt_hours(summary, "2025-01-01")


def test_end_to_end_divergence_from_dominant_era_only_lookup():
    """Full pipeline: build a real multi-era series where the dominant
    era (by sample count) is heavily fragmented and the non-dominant
    era is the one longest_contiguous_segment actually selects -
    proving get_dt_hours changes the real answer, not just a value in
    isolation."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("stage01", REPO_ROOT / "01_load_clean.py")
    stage01 = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stage01)

    rng = np.random.default_rng(52)

    era1_idx = pd.date_range("2010-01-01", "2018-12-31 23:45", freq="15min")
    n1 = len(era1_idx)
    vhm0_1 = np.clip(1.0 + rng.normal(0, 0.2, n1), 0.05, None)
    block = n1 // 61
    for i in range(60):
        start = block * (i + 1) - 500
        vhm0_1[start : start + 1000] = np.nan

    era2_idx = pd.date_range("2019-01-01", "2024-12-31 23:30", freq="30min")
    n2 = len(era2_idx)
    vhm0_2 = np.clip(1.0 + rng.normal(0, 0.2, n2), 0.05, None)

    idx = era1_idx.append(era2_idx)
    vhm0 = np.concatenate([vhm0_1, vhm0_2])
    s = pd.Series(vhm0, index=idx, name="VHM0")

    s_clean, report = stage01.regularize_and_clean(s, "VHM0")
    segment, meta = utils.longest_contiguous_segment(s_clean)

    # Sanity: this test is only meaningful if the segment actually
    # landed in the non-dominant era - assert that setup assumption
    # explicitly so a future change to the fixture can't silently
    # make this test pass for the wrong reason.
    assert pd.Timestamp(meta["segment_start"]) >= pd.Timestamp("2019-01-01")

    correct_dt = utils.get_dt_hours(report, meta["segment_start"])
    dominant_era_dt = report["sampling_interval_hours"]
    assert correct_dt == 0.5
    assert dominant_era_dt == 0.25
    assert correct_dt != dominant_era_dt
