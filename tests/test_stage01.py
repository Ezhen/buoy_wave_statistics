"""
Regression tests for 01_load_clean.py's era-aware regularization.

These transcribe the synthetic known-answer cases built and manually
validated during the session that found and fixed two real bugs in
Stage 01:
  1. Era-mismatch fabrication - a single assumed native frequency for
     the whole record silently fabricated interpolated data for any
     buoy whose true rate changed mid-record (confirmed on real
     Zeebrugge data: 23.8% of its record was fabricated).
  2. Timestamp-jitter false missingness - exact-timestamp reindex
     silently mismatched real samples carrying ordinary sub-minute
     telemetry jitter (confirmed synthetically: +/-3s of jitter alone
     produced 85.9% false missingness with zero real gap).
  3. A THIRD bug found while fixing #1 - a genuine long gap between two
     stretches of the SAME native rate was mis-split into two separate
     "eras," silently erasing the gap period (fewer output rows,
     shrunk record_years) instead of representing it as NaN.

Run with: pytest tests/test_stage01.py -v

Without this file, catching a regression in any of the three fixes
above requires manually rebuilding these synthetic cases by hand again -
exactly what happened three times in the session that produced this
file. The whole point of transcribing them here is that the next
regression gets caught automatically instead.
"""

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent


def _load_stage01_module():
    """01_load_clean.py starts with a digit, so it isn't a valid module
    name for a plain `import` - load it directly by file path instead.

    spec_from_file_location does NOT add the loaded file's own
    directory to sys.path the way running a script directly does - and
    01_load_clean.py itself does `from utils import ...`, so without
    this, that import fails with ModuleNotFoundError as soon as
    exec_module runs, regardless of where pytest is invoked from
    (pytest's own rootdir insertion puts tests/ on sys.path when
    there's no __init__.py here, not the repo root)."""
    if str(REPO_ROOT) not in sys.path:
        sys.path.insert(0, str(REPO_ROOT))
    spec = importlib.util.spec_from_file_location(
        "stage01", REPO_ROOT / "01_load_clean.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["stage01"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def stage01():
    return _load_stage01_module()


# ---------------------------------------------------------------------
# Case 1: single-era buoy with ordinary scattered short gaps (VALUES
# missing at otherwise-regular timestamps, not missing timestamps) -
# the baseline case that must keep working exactly as before.
# ---------------------------------------------------------------------
def _single_era_series(seed=11):
    np.random.seed(seed)
    idx = pd.date_range("2005-01-01", "2020-12-31 23:00", freq="h")
    n = len(idx)
    vhm0 = np.clip(1.0 + np.random.normal(0, 0.2, n), 0.05, None)
    gap_starts = np.random.choice(n - 50, 20, replace=False)
    for gs in gap_starts:
        length = np.random.randint(1, 10)
        vhm0[gs : gs + length] = np.nan
    return pd.Series(vhm0, index=idx, name="VHM0")


def test_single_era_no_false_era_split(stage01):
    """Ordinary scattered short gaps must never be mistaken for a
    sampling-rate change - this is the most common real-world case and
    must stay a single era."""
    s = _single_era_series()
    s_clean, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["multi_era_detected"] is False
    assert len(report["eras"]) == 1
    assert report["eras"][0]["inferred_freq_hours"] == 1.0


def test_single_era_gap_classification_unchanged(stage01):
    """Exact counts from this session's manual validation - short gaps
    (<=3 samples) interpolated, longer gaps left as NaN, nothing
    fabricated or dropped."""
    s = _single_era_series()
    s_clean, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["n_samples_regularized"] == len(s)
    assert report["n_short_gap_interpolated"] == 8
    assert report["n_long_gap_left_as_nan"] == 101
    assert report["longest_gap_samples"] == 9


# ---------------------------------------------------------------------
# Case 2: pure rate change, no real gap - must split into exactly 2
# eras, each regularized at its own true rate, with ZERO fabricated
# interpolation (the original bug this fix targets).
# ---------------------------------------------------------------------
def _rate_change_series(seed=12):
    rng = np.random.default_rng(seed)
    fine = pd.date_range("2009-01-01", "2017-12-31 23:45", freq="15min")
    coarse = pd.date_range("2018-01-01", "2023-12-31 23:30", freq="30min")
    idx = fine.append(coarse)
    n = len(idx)
    vhm0 = np.clip(1.0 + rng.normal(0, 0.2, n), 0.05, None)
    return pd.Series(vhm0, index=idx, name="VHM0")


def test_rate_change_detected_as_two_eras(stage01):
    s = _rate_change_series()
    _, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["multi_era_detected"] is True
    assert len(report["eras"]) == 2
    freqs = sorted(e["inferred_freq_hours"] for e in report["eras"])
    assert freqs == pytest.approx([0.25, 0.5])


def test_rate_change_no_fabrication(stage01):
    """The core bug: era-mismatched data must NOT be silently
    interpolated. Zero fabrication expected since both eras are
    internally clean (no real gaps injected in this fixture)."""
    s = _rate_change_series()
    _, report = stage01.regularize_and_clean(s, "VHM0")
    for era in report["eras"]:
        assert era["n_short_gap_interpolated"] == 0


# ---------------------------------------------------------------------
# Case 3: a genuine long gap BETWEEN two stretches of the SAME native
# rate, plus a genuine rate change elsewhere - must merge the two
# same-rate stretches into ONE era (gap represented as NaN inside it),
# not silently erase the gap period between two artificial eras.
# ---------------------------------------------------------------------
def _gap_plus_rate_change_series(seed=20):
    rng = np.random.default_rng(seed)
    era1a = pd.date_range("2009-01-01", "2010-06-30 23:45", freq="15min")
    era1b = pd.date_range("2012-01-01", "2017-12-31 23:45", freq="15min")
    era2 = pd.date_range("2018-01-01", "2020-12-31 23:30", freq="30min")
    idx = era1a.append(era1b).append(era2)
    n = len(idx)
    vhm0 = np.clip(1.0 + rng.normal(0, 0.2, n), 0.05, None)
    true_gap_hours = (era1b[0] - era1a[-1]).total_seconds() / 3600.0
    return pd.Series(vhm0, index=idx, name="VHM0"), true_gap_hours


def test_gap_between_same_rate_stretches_merges_to_one_era(stage01):
    s, _ = _gap_plus_rate_change_series()
    _, report = stage01.regularize_and_clean(s, "VHM0")
    # 2 eras: the merged 15-min stretch (with the gap inside it) + the
    # separate 30-min era - NOT 3 (which would mean the gap incorrectly
    # split same-rate data into two separate eras).
    assert len(report["eras"]) == 2


def test_gap_represented_as_nan_not_erased(stage01):
    """The bug this specifically fixes: the gap period must survive as
    an explicit long-NaN run, not vanish (fewer rows, shrunk span)
    between two artificially-separated eras."""
    s, true_gap_hours = _gap_plus_rate_change_series()
    s_clean, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["longest_gap_hours"] == pytest.approx(true_gap_hours, rel=0.01)
    true_span_years = (s.index.max() - s.index.min()).days / 365.25
    observed_span_years = (s_clean.index.max() - s_clean.index.min()).days / 365.25
    assert observed_span_years == pytest.approx(true_span_years, rel=0.01)


# ---------------------------------------------------------------------
# Case 4: ordinary telemetry jitter must NOT cause false missingness -
# the second bug found this session, in the base reindex logic itself
# (not era-detection).
# ---------------------------------------------------------------------
def _jittered_series(jitter_seconds, seed=31, n=50000):
    rng = np.random.default_rng(seed)
    base = pd.date_range("2018-01-01", periods=n, freq="30min")
    jitter = rng.integers(-jitter_seconds, jitter_seconds + 1, n)
    idx = base + pd.to_timedelta(jitter, unit="s")
    vhm0 = np.clip(1.0 + rng.normal(0, 0.2, n), 0.05, None)
    return pd.Series(vhm0, index=idx, name="VHM0")


def test_small_jitter_does_not_cause_false_missingness(stage01):
    """+/-3s jitter reproduced 85.9% false missingness before the
    snap-to-grid fix. After the fix, real samples must land on their
    intended grid slot regardless of a few seconds of jitter."""
    s = _jittered_series(jitter_seconds=3)
    _, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["pct_missing_after_clean"] < 1.0


def test_snap_collisions_are_reported_not_silent(stage01):
    """Large jitter (comparable to half the sampling interval) WILL
    legitimately collide two samples onto the same grid slot - that's
    expected, but must be counted, not silently dropped without a
    trace."""
    s = _jittered_series(jitter_seconds=600)  # +/-10min, on a 30min grid
    _, report = stage01.regularize_and_clean(s, "VHM0")
    total_collisions = sum(e["n_snap_collisions"] for e in report["eras"])
    assert total_collisions >= 0  # field must exist and be well-formed
    assert "n_snap_collisions" in report["eras"][0]


# ---------------------------------------------------------------------
# Case 5: a real rate-change seam with only a SHORT natural gap (no
# deliberately injected long gap) must still always break contiguity -
# found via direct empirical test, not theory: a rate change with a
# short seam concatenated with ZERO NaN marker by default, meaning
# longest_contiguous_segment (used by Stages 11b, 13, 24, 26) would
# silently splice two different sampling rates into one "contiguous"
# segment. Reproduces the exact real Zeebrugge seam timing (era1 ends
# 09:30, era2 begins 10:15 - a 45-minute natural gap).
# ---------------------------------------------------------------------
def _short_seam_rate_change_series(seed=40):
    rng = np.random.default_rng(seed)
    era1 = pd.date_range("2017-01-01 00:00", "2017-10-16 09:30:00", freq="15min")
    era2 = pd.date_range("2017-10-16 10:15:00", "2018-06-01 00:00:00", freq="30min")
    idx = era1.append(era2)
    vhm0 = np.clip(1.0 + rng.normal(0, 0.2, len(idx)), 0.05, None)
    return pd.Series(vhm0, index=idx, name="VHM0")


def test_era_boundary_always_breaks_contiguity(stage01):
    """The core bug: a short natural seam between two eras must not be
    silently bridged. Must insert an explicit break even when no long
    gap already separates the two rates."""
    s = _short_seam_rate_change_series()
    s_clean, report = stage01.regularize_and_clean(s, "VHM0")
    assert report["n_era_boundary_markers"] >= 1

    seam_window = s_clean["2017-10-16 09:00:00":"2017-10-16 11:00:00"]
    assert seam_window.isna().any(), (
        "no NaN break found spanning the era transition - a contiguous-"
        "segment check would silently splice two different sampling "
        "rates together"
    )


def test_era_boundary_respected_by_longest_contiguous_segment(stage01):
    """The actual consumer-facing guarantee: longest_contiguous_segment
    must stop at the era boundary, not just 'a NaN exists somewhere
    nearby' - this is what 11b/13/24/26 actually rely on."""
    sys.path.insert(0, str(REPO_ROOT))
    import utils  # noqa: E402 - path bootstrap must run first

    s = _short_seam_rate_change_series()
    s_clean, _ = stage01.regularize_and_clean(s, "VHM0")
    segment, meta = utils.longest_contiguous_segment(s_clean)

    era_boundary_start = pd.Timestamp("2017-10-16 09:30:00")
    era_boundary_end = pd.Timestamp("2017-10-16 10:15:00")
    spans_boundary = (
        segment.index.min() <= era_boundary_start
        and segment.index.max() >= era_boundary_end
    )
    assert not spans_boundary, (
        "longest_contiguous_segment returned a segment spanning the era "
        "boundary - it would mix two different sampling rates under one "
        "dt_hours assumption"
    )
