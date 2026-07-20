# Session Handoff — Belgian Coastal Zone Wave Buoy Pipeline

Paste this at the start of a new chat to continue exactly where this session left off.
Supersedes the previous handoff (2026-07-12) - a lot changed since then, including
two real correctness bugs found in the foundational cleaning stage.

## Project identity

Statistical characterization + forecasting pipeline for 19 CMEMS in-situ wave buoys
in the Belgian Coastal Zone (BCZ), built on ULiège/CECI HPC infrastructure. User
addresses Claude as "Twin." Repo lives at `~/wave_statistics/` on the HPC
(Python 3.10.20). Sandbox dev/test environment used during Claude sessions runs
Python 3.12.3 - this version gap is a real, proven source of bugs (never write a
multi-line expression inside an f-string's `{}` - fine on 3.12, `SyntaxError` on 3.10).

Key reference files, all current as of this session:
- `PLAN_next_session.md` — the living priority/task tracker; **read this first**,
  it now has 1883 lines covering everything from Priority 1 through this session's
  Stage 01 fixes and repository restructuring. The final section, "Status as of
  2026-07-20, end of session," lists every open thread in priority order.
- `README.md` — project overview, stage table (execution order), and a NEW
  "What question each stage answers" table (grouped by scientific question)
- `CHANGELOG.md` — full bug-fix and build history, dated entries
- `METHODS.md` — formal methods writeup, now includes Section 7 (Stages 24-26)
  and a Section 8 limitations note on the Stage 01 data-integrity history
- `tests/README.md` — how to run the pytest suite (19 tests, `pytest tests/ -v`)

## What changed this session (2026-07-19/20) - the short version

Started as "run SSA on Zeebrugge" (Stage 26). Surfaced three real, previously-
invisible bugs in `01_load_clean.py` along the way - found through internal
consistency checks, not a priori review, same pattern as prior sessions' bugs:

1. **Era-mismatch fabrication** (8 of 19 buoys) - a single assumed native sampling
   frequency per buoy silently fabricated interpolated data for any buoy whose true
   rate changed mid-record. Up to 23.8% of Zeebrugge's record was fabricated data,
   concentrated in its best-looking recent decade.
2. **Timestamp-jitter false missingness** (network-wide) - exact-timestamp reindex
   silently mismatched real samples carrying ordinary telemetry jitter.
3. **Era-boundary silent splicing** - a rate change with only a short natural seam
   concatenated with zero gap marker, letting downstream gap-aware tools silently
   mix two different sampling rates into one "contiguous" segment.

All three fixed and validated against synthetic ground truth. Follow-on fix: Stages
11b/13/24 were reading a single (dominant-era-only) `sampling_interval_hours` value
as a lag/dt parameter - wrong for any non-dominant-era segment. Fixed via
`utils.get_dt_hours()` and `utils.segments_by_time_gap_era_aware()`.

Then did a 4-phase repository restructuring (buoy/stage registries + a generic
`tools/run_stage.py`, replacing 3 independent eligibility-discovery implementations
with one shared source of truth), plus started a structured-warnings convention
(`utils.emit_warning()`) on 2 of 26 stages as proof-of-concept.

**Full detail on all of the above: see `PLAN_next_session.md`'s final sections.**

## Working discipline (carried forward + reinforced this session)

All prior discipline (never trust a new library API without an isolated test;
validate every new method on synthetic data with a KNOWN answer BEFORE trusting
real results; verify schemas against actual source, never guess; gap-fragmentation
requires `longest_contiguous_segment`/`all_contiguous_segments`/`segments_by_time_gap`,
never naive `dropna()`) still applies, reinforced by this session specifically:

- **A synthetic test that "passes" can still miss the real bug** if it doesn't
  exercise the actual failure mode - the first Stage 01 regression test used
  scattered NaN VALUES on an already-regular raw timestamp index, which never
  exercised a genuine raw-timestamp GAP, missing the exact bug it was meant to
  catch. Rebuilt with a real gap+rate-change combination once this was caught.
- **A fix can introduce a new bug while fixing the old one** - the era-merge fix
  (for era-mismatch fabrication) initially caused a DIFFERENT bug (silently
  erasing a real multi-year gap between two same-rate stretches) - caught via a
  `record_years` sanity check that didn't match a known-correct value, not by
  the original validation suite.
- **Building on a system correctly means checking, not assuming, that the system
  actually guarantees what you're about to rely on** - the `get_dt_hours` accessor
  was only safe to build because `longest_contiguous_segment` was confirmed (by
  direct empirical test, not assumption) to never span an era boundary - which
  required fixing bug #3 above first, discovered specifically while checking this.

## Current open threads, in priority order (full detail in PLAN_next_session.md)

1. **AkkaertSouthwestBuoy's snap-collision anomaly** (25,389 collisions, ~2.9% of
   its record, 10-100x every other buoy) - flagged repeatedly, never directly
   investigated. Natural next step: `tools/check_sampling_interval_history.py`
   against it specifically.
2. **Westhinder's 2001-2006 change-point** - every artifact explanation tested has
   failed (coverage, QC, network cross-reference, file history). Next step
   identified, not executed: direct inquiry to MDK about deployment/mooring/sensor
   history, citing the sand-dynamics report's "problems with the wave buoy"
   precedent.
3. **Stage 24 non-dominant-era reruns** for the 8 multi-era buoys - deferred by
   explicit request. Real unused data (up to 42% of valid samples for Zeebrugge).
4. **`run_changepoint_batch.py` migration** to wrap `tools/run_stage.py` - not
   urgent, current tool works fine standalone.
5. **`requirements.txt` placeholders** for copernicusmarine/cdsapi - still open
   from before this session.
6. **Broader literature grounding** - only Caires 2011 checked so far.
7. **Structured-warnings convention** - only 2 of 26 stages retrofitted
   (Stage 15, Stage 25); apply opportunistically going forward, not as a
   standalone pass.

## Repo structure quick-reference

```
01-26_*.py          numbered pipeline stages (see README's two stage tables)
utils.py            shared library - schema helpers, gap-aware segmentation,
                     get_dt_hours, stage_eligible, emit_warning
tools/               reactive diagnostics + registry-driven runner
  build_buoy_registry.py / build_stage_registry.py / run_stage.py
  check_sampling_interval_history.py / inspect_netcdf_metadata.py
  run_changepoint_batch.py
tests/               pytest suite - 19 tests, `pytest tests/ -v`
run_all_buoys.py     default Core batch (Stages 02-13 only, deliberately)
buoy_registry.json / stage_registry.json   generated - regenerate via tools/build_*.py
pipeline_out/        all stage outputs, one subfolder per stage
```
