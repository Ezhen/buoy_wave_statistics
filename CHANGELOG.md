# Changelog

All notable changes to this repository, in the order they actually
happened. Dates reflect real batch-run timestamps where available.

## 2026-07-19 — Stage 01 data-integrity fixes (era-mismatch fabrication, jitter-driven false missingness); Stage 26 (SSA)

### Fixed — `01_load_clean.py`
- **Era-mismatch fabrication (8 of 19 buoys affected).** `regularize_and_clean`
  previously assumed ONE native sampling frequency for a buoy's entire
  record (the global mode of raw time diffs). Any buoy whose native rate
  changed mid-record (found via a network-wide audit: A2Buoy,
  AkkaertSouthwestBuoy, NieuwpoortBuoy, OstendEasternPalisadeBuoy,
  ScheurWielingenBuoy, TrapegeerBuoy, WandelaarBuoy,
  ZeebruggeZandopvangkadeBuoy — mostly a shared rate decrease around
  2016-2021, consistent with a single network-wide operator telemetry
  change rather than 8 independent events) had every real reading in
  its non-dominant era silently fabricated via linear interpolation
  across the mismatched grid, reported as 0% missing. Confirmed on real
  Zeebrugge data: 23.8% of its full 605,091-sample record was fabricated
  interpolation, concentrated almost entirely in 2018-2026 - the era that
  looked best-covered. Fixed via `detect_sampling_eras()` +
  per-era `regularize_and_clean_one_era()`, with adjacent same-rate runs
  explicitly merged so a genuine long GAP between two same-rate stretches
  (e.g. Zeebrugge's real ~2.93yr 2010-2012 outage) isn't misread as a
  rate change and silently erased instead of represented as NaN - this
  exact bug was caught and fixed within today's session, via a
  `record_years` sanity check that didn't match the known-correct value.
- **Timestamp-jitter false missingness (network-wide, not era-specific).**
  The regularize step's `reindex` requires exact-second timestamp
  matching; ordinary real-world telemetry jitter of even a few seconds
  caused real samples to silently miss their grid slot and register as
  missing/interpolated. Reproduced synthetically: ±3s of random jitter
  alone produced 85.9% false missingness with zero real gap. Fixed by
  snapping each raw timestamp to its nearest grid point
  (`.index.round(freq)`) before reindexing; collisions are counted
  (`n_snap_collisions`) rather than silently dropped. AkkaertSouthwestBuoy
  showed an outlier collision rate (25,389, ~2.9% of its record) worth a
  follow-up look - possibly irregular native cadence rather than simple
  jitter.
- Both fixes validated against synthetic known-answer cases (single-era
  regression - byte-identical output to the pre-fix code; pure rate
  change; gap + rate change combined; small vs large jitter) before
  trusting real-data results, per this project's established discipline.
  All 19 buoys re-run; every downstream stage depending on Stage 01
  output for the 8 affected buoys (and, to a lesser extent, all 19, given
  the jitter fix's broader scope) re-run and re-validated.

### Added
- `26_ssa_decomposition.py` — Singular Spectrum Analysis (trajectory-
  matrix/Hankel SVD, truncated via `scipy.sparse.linalg.svds`), the third
  and last post-priority signal-processing extension (alongside Stage 24
  HMM and Stage 25 change-point), scoped for Zeebrugge's unresolved
  tidal-notch anomaly. Validated against synthetic ground truth (single
  sinusoid recovery, M2/S2 near-degenerate frequency separation, pure-
  noise scree control, full synthetic M2+S2+MS4 compound-tide scenario)
  before running on real data. Real result (post-fix, 2015-01 to 2016-09
  contiguous segment, verified clear of both bugs above): no clean
  M2/S2/MS4 pair recovered despite three independent runs; one robust,
  reproducible non-tidal ~49.7h (~2.07 day) oscillatory pair recovered
  consistently across three different segment selections (49.671h,
  49.671h, 49.700h) - a genuine, well-resolved finding, not a fluke of
  one window. Negative tidal-constituent result, now free of the data-
  integrity confounds above, lends more weight to the time-varying-
  tidal-parameter hypothesis over the compound-tide (MS4) hypothesis for
  Zeebrugge's notch failure.
- `check_sampling_interval_history.py` — audits a buoy's raw `.nc` TIME
  array for native-interval changes by year, independent of Stage 01's
  own regularization; this is what first surfaced the era-mismatch bug.
- `inspect_netcdf_metadata.py` — dumps global/variable NetCDF attributes,
  QC flag distributions and their correspondence to NaN in the value
  array, and per-timestep position (mooring relocation check). Used to
  rule out (for Westhinder) a QC-methodology or history-trail explanation
  for the 2001-2006 change-point investigation; came back clean on all
  five axes checked (annual + storm-season coverage, network cross-
  reference, file history, QC-to-NaN correspondence).
- `run_changepoint_batch.py` — discovers long-record buoys from Stage
  01's own output (no hardcoded buoy list existed anywhere in the repo)
  and batch-runs Stage 25 full-record + edge-trimmed (relative to each
  buoy's own record start, not a shared calendar year).
- Stage 25 (`25_changepoint_detection.py`) extended with storm-season
  (Oct-Mar) coverage checking (annual coverage alone doesn't discriminate
  a calm-season gap from a storm-season gap - confirmed on real
  Westhinder data, where two segments had near-identical annual coverage
  but only one showed an anomalous mean) and segment-interior scanning
  (not just years adjacent to a change-point boundary - the previous
  boundary-only check missed 2003, the single most extreme year in
  Westhinder's record, sitting mid-segment).

## 2026-07-09 — Initial pipeline (Westhinder pilot)

### Added
- `download_belgian_wave_buoys.py` — CMEMS in-situ wave buoy download for
  the Belgian Coastal Zone, product `INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036`
  (auto-discovers the exact dataset id via `copernicusmarine.describe()`
  rather than hardcoding a guess).
- `plot_belgian_buoys.py` — per-buoy time series + combined location map,
  no basemap dependency (works offline on the HPC).
- Core characterization pipeline, stages 01-06: `01_load_clean.py`
  (grid regularization, gap interpolation, sanity bounds),
  `02_eda_diagnostics.py` (rolling stats, ACF/PACF, periodogram with M2
  marker), `03_stationarity_tests.py` (ADF + KPSS side by side),
  `04_transform_detrend.py` (Box-Cox + differencing),
  `05_whiteness_check.py` (Ljung-Box), `06_distribution_fit.py`
  (Rayleigh/Weibull/log-normal + KS + Q-Q plots).
- `03b_tidal_notch.py` — harmonic regression at the M2 frequency,
  inserted between stationarity testing and differencing after Westhinder's
  periodogram showed a dominant tidal peak (43x baseline power). `04`
  updated to auto-consume 03b's detided output when present.
- `07_arch_lm_test.py` — Engle's ARCH-LM test for volatility clustering.
- `08_extreme_value_analysis.py` — Peaks-over-Threshold + GPD fit, with
  mean-residual-life threshold diagnostics and return-level estimation.
- `run_all_buoys.py` — batch orchestrator, continues past a failed buoy
  instead of aborting the whole run.
- `summarize_results.py` — builds one cross-buoy comparison table from
  each stage's summary JSON/CSV output.
- `rerun_eva_all_buoys.py` — re-run only Stage 08 at a different
  threshold, without redoing the full pipeline.
- `README.md` — pipeline documentation, updated iteratively through the day.

### Fixed
- Stage 6 (distribution fit) and Stage 8 (EVA) deliberately use the raw
  cleaned level series, not the Stage 3 residual - Rayleigh/Weibull/GPD
  describe Hs itself, not a mean-zero differenced residual.

### Findings
- Westhinder: M2 tidal ratio 43x (before notch) -> ~10x (after, 2
  harmonics - never fully clean). ADF/KPSS disagree. Ljung-Box fails at
  every tested lag - diagnosed as genuine storm persistence, not a
  pipeline artifact. Weibull best fit. ARCH effects confirmed.
- Full 19-buoy run: sampling-rate heterogeneity discovered (CadzandBoei
  and Deurlo at 10-min vs. 30-min for the rest). M2 contamination has a
  real spatial gradient (6x-135x, worst at Zeebrugge). Ljung-Box fails
  network-wide, confirming the physical (not artifact) explanation.
  Weibull wins 17/19 buoys. EVA unreliable for 18/19 at the 95th
  percentile threshold (too few storm peaks) - fixed by lowering to the
  85th percentile / 24h separation, after which 19/19 buoys were
  reliable and every GPD shape parameter came out negative (bounded
  upper tail - physically sensible for a fetch-limited shelf sea).

### Changed
- `01_load_clean.py`: added `sampling_interval_hours`, `longest_gap_hours`,
  `n_duplicate_timestamps_raw` to the Stage 0 report.
- `summarize_results.py`: added `ljungbox_lags_hours` (real-time
  equivalents, since sampling intervals differ) and an explicit warning
  when sampling intervals differ across buoys.

## 2026-07-12 — Pipeline expansion, bug fixes, multi-year data, meteo

### Added - new analysis stages
- `09_cross_variable_analysis.py` — Hs/Tp/Tm02/direction correlation
  (VMDR decomposed to sin/cos, since it's circular), lagged
  cross-correlation, PCA.
- `10_regime_identification.py` — GMM clustering into
  calm/moderate/energetic/storm regimes; `--include-period` gracefully
  falls back to Hs-only if VTPK is absent (fixed a crash - see below).
- `11_spatial_statistics.py` — network-wide pairwise correlation,
  correlation-vs-distance (haversine), hierarchical buoy clustering.
- `11b_dependence_structure.py` — per-buoy integral (persistence)
  timescale from the ACF, feeding block-bootstrap block length and EVA
  declustering-window justification downstream.
- `12_confidence_intervals.py` — block bootstrap on Hs mean/quantiles
  (block length from 11b), bootstrap on GPD xi using Stage 08's
  declustered peaks directly, CI band on the Weibull fit.
- `12b_correlation_confidence.py` — Fisher z CI on Stage 11's pairwise
  correlations, with a Dawdy-Matalas AR(1) effective-N correction for
  autocorrelation (using each buoy's lag-1 ACF from 11b).
- `13_stability_analysis.py` — moving-window distribution stability,
  drop-biggest-storm jackknife (cross-referenced against Stage 12's CI
  width), block-bootstrap on Stage 10's regime fractions.
- `PLAN_next_session.md` — priority-ordered roadmap, updated throughout
  the day as each priority was completed.

### Added - data quality / infrastructure
- `01_load_clean.py`: added `available_variables` and `record_years` to
  the Stage 0 report (`utils.detect_available_variables()`).
- `utils.py`: `load_buoy_dataframe()` (multi-variable loader),
  `haversine_km()`, `count_raw_duplicate_timestamps()`,
  `resolve_block_length()` (shared block-length fallback logic, used by
  both Stage 12 and Stage 13), `detect_available_variables()`.
- `run_all_buoys.py`: rewritten with an explicit Core/Advanced tiering
  gate (`stage_eligible()`) - replaces three previously separate ad hoc
  patterns (Stage 09's manual loop exclusion, Stage 10's silent
  fallback, a stage crashing on a missing variable) with one declared
  mechanism. Stage 0 now always runs first per buoy specifically to
  populate the gate's eligibility checks.

### Fixed
- **`10_regime_identification.py`**: crashed with `KeyError: 'VTPK'`
  when `--include-period` was used on a buoy without that variable
  (discovered on Westhinder, which turned out to be VHM0-only). Now
  falls back to Hs-only with an explicit warning instead of crashing.
- **`11b_dependence_structure.py` (major)**: v1 computed the ACF on the
  raw, still tide-contaminated series. A strongly periodic component
  makes the ACF cross zero at fractions of the tidal period - not
  evidence of losing memory. Symptom: Zeebrugge (most tidally
  contaminated buoy in the network) showed the *shortest* apparent
  persistence, backwards from physical expectation. Fixed: now reads
  Stage 03b's detided series (falls back to raw with a loud warning if
  03b hasn't run), and replaced the single zero-crossing criterion with
  a significance-band test (held for 5 consecutive lags) to avoid both
  spurious early crossings and silently reporting a lower-bound-at-the-
  search-ceiling as if it were a real measurement. Reordered in
  `run_all_buoys.py` to run after `03b` instead of right after `02`.
- **`run_all_buoys.py`**: caught before shipping - `13_stability_analysis.py`
  was wired with a `--data-dir` argument it doesn't accept (would have
  crashed every buoy at that stage). Fixed prior to the first real run.

### Investigated, not fully resolved
- **Zeebrugge tidal contamination**: `--harmonics 3` barely moved either
  the residual M2 ratio (35.99 -> 36.32) or the persistence estimate
  (56.95h -> 57.16h). Ruled out "just needs more harmonics." Working
  theory: the notch's fixed-M2-exact-frequency assumption likely
  doesn't hold at this specific shallow/harbor site (possible S2
  constituent or shoaling-shifted effective tidal frequency). Confirmed
  independently by Stage 13's moving-window check: Zeebrugge's
  best-fit-distribution verdict is unstable across sub-windows (75%
  agreement, not 100%) - a third, independent signal pointing at the
  same unresolved data-quality issue. Flagged as a scoped future fix
  (fit the tidal frequency rather than assume it), not chased further
  this session.

### Findings
- Full-network Stage 12/12b run: GPD xi CIs vary hugely by peak count -
  a 7-peak buoy's CI spanned [-2.9, 0.5] (crosses zero, uninformative
  about tail boundedness) vs. a 13-peak buoy's [-2.2, -0.06] (decisive).
  Stage 12b found a ~99% effective-sample-size reduction from
  autocorrelation on real data (not just the synthetic test fixture) -
  direct confirmation that an uncorrected Fisher z CI would have been
  dramatically overconfident.
- Variable availability is structurally uneven: 17/19 buoys report only
  `VHM0`; only CadzandBoei (+VTPK) and Deurlo (+VTPK, VMDR) carry the
  fuller sensor set - the same two buoys that also run at 10-min instead
  of 30-min sampling, suggesting a genuinely different instrument class,
  not a coincidental data gap.

### Multi-year data - three corrections in sequence
1. First attempt guessed a separate delayed-mode product id
   (`INSITU_NWS_PHYBGCWAV_DISCRETE_MY_013_036`) - doesn't exist. Unlike
   the global product, NWS has no separate MY product; it's a single
   combined MYNRT product (the one already in use).
2. Second attempt requested a wide `start_datetime` via `subset()` on
   the same product - still capped at 30 days. Root cause: `subset()`
   only ever serves the `latest` dataset part for in-situ products,
   regardless of requested date range. The full archive (`history` part,
   one NetCDF per platform) requires the Files service
   (`copernicusmarine.get()`) instead - a different API entirely.
3. `download_belgian_wave_buoys_history.py` — added, using `get()` with
   `dataset_part="history"`. v1 guessed `show_outputnames=True` (doesn't
   exist) and `force_download` (deprecated) - corrected against the
   actual installed signature via `inspect.signature()` and
   `copernicusmarine get --help`. Final version: `--list` (uses the
   confirmed `create_file_list` option to preview filenames with zero
   download) then `--download --file-list <path>` (exact-path download,
   no regex ambiguity). Also filters the raw download to the Belgian
   bbox by reading each file's actual LATITUDE/LONGITUDE as a safety net.

### Multi-year - result
- Full NWS `history` listing: 9778 files, whole North West Shelf region.
  Grepped by known buoy names -> all 19 matched exactly under
  `history/MO/NO_TS_MO_<name>.nc` (a handful of *unrelated* files at the
  same sites were also present - `*Wind`, `*Weather`, `*Tide`, `*MP`
  variants - correctly excluded from the exact-match file list).
- Real download: 81.39 MB for all 19 buoys (vs. 21.38 GB for the
  unfiltered whole-region pull).
- **Real per-buoy coverage is highly uneven**: 6 buoys go back to
  1990-1997 (Westhinder longest, from 1990-07-19); the remaining 13
  range from 2009 to 2021. All extend to 2026-06-30.

### Meteo (ERA5)
- `download_era5_meteo.py` — added. Variables: 10m u/v wind (primary
  wave-generation driver), mean sea-level pressure (leading indicator),
  2m air temperature + SST (their difference affects the drag
  coefficient - a second-order correction, not a standalone driver).
- v1 (yearly chunks, hourly) hit CDS's request-size limit on a single
  year. Fixed: monthly chunking + 3-hourly default sampling (~96x
  request-size reduction) - justified by Stage 11b's own finding that
  storm persistence operates on a 50-100h scale, so 3-hourly resolves
  it comfortably without requesting resolution nothing downstream uses.
- ERA5 range decided directly from the real multi-year buoy coverage
  above, not a round number: 2010-2026 (17 years) as the primary pull -
  covers 17/19 buoys near-completely; extending back to 1990 deferred
  until a specific long-record analysis (e.g. Mann-Kendall on the 6
  oldest buoys) actually needs it.

### Infrastructure
- `.gitignore` added - excludes `pipeline_out/`, all raw data
  directories, intermediate download file-lists, Python cache, and
  credential files.
- Consolidated `download_belgian_wave_buoys.py` and `plot_belgian_buoys.py`
  into `wave_pipeline/` (previously sitting outside it, predating the
  folder's creation).

## 2026-07-13 — Multi-year data live, gap-splicing fixes, performance, Mann-Kendall

### Fixed — gap-splicing bugs, found once real multi-year (gappy) data arrived
All of these share one root cause: a naive `.dropna()` silently
concatenates the sample before a gap with the sample after it, treating
unrelated calendar periods as temporally adjacent - fine for point
statistics (distribution fits, EVA), invalid for anything lag/order-based.

- **`utils.py`**: three new helpers - `longest_contiguous_segment()`,
  `all_contiguous_segments()` (length-weighted aggregation across every
  qualifying segment, not just the longest), `segments_by_time_gap()`
  (for series where gap rows are entirely absent rather than
  NaN-marked, e.g. Stage 10's regime labels).
- **`resolve_coord_name()`** added to `utils.py` - case-insensitive
  lookup for `TIME`/`LATITUDE`/`LONGITUDE`, since the NRT (`subset()`)
  and multi-year (`get()`) downloads returned different capitalization
  for the same logical fields. Applied throughout `utils.py`,
  `download_belgian_wave_buoys_history.py`, `plot_belgian_buoys.py`.
- **`03b_tidal_notch.py`**: rewritten to fit the harmonic regression on
  valid data only but evaluate/output on the FULL grid, so gap
  positions survive as NaN into every downstream stage instead of being
  silently collapsed - this was the actual source of the erasure
  everything else inherited.
- **`02_eda_diagnostics.py`, `03_stationarity_tests.py`,
  `04_transform_detrend.py`**: now restrict their lag-based computation
  (ACF/PACF/periodogram, ADF/KPSS, differencing) to the longest
  contiguous segment instead of the gap-spliced full record.
- **`11b_dependence_structure.py` (v3, major rewrite)**: v2's "use the
  longest segment" discarded ~93% of valid data on Westhinder's
  fragmented 36-year record (longest single segment = 6.9% of valid
  samples). Now aggregates the integral timescale across every
  qualifying segment (length-weighted mean) - Westhinder went from
  using 6.9% of data to 61.4%. `--max-lag` default raised 500->2000
  (applied per-segment).
- **`13_stability_analysis.py` `[A]`**: was splitting the post-dropna()
  COLLAPSED array positionally - a "window" on a fragmented record
  could be a patchwork of disjoint calendar periods, and printed window
  start/end dates were actively wrong. Now splits by real calendar
  time; window sample counts correctly come out unequal reflecting real
  gap-coverage differences across eras.
- **`13_stability_analysis.py` `[C]`**: same root issue, different
  detection method needed (Stage 10's regime labels have gap rows
  entirely absent, not NaN-marked) - added `segments_by_time_gap()` and
  bootstraps within each detected segment separately. Cross-validated:
  detected the same 1905 segments on Westhinder that the NaN-based
  methods found independently on the raw Hs series.
- **`13_stability_analysis.py` `[C]` plotting crash**: `ValueError:
  'yerr' must not contain negative values` - crashed CadzandBoei and
  Deurlo outright in the first full batch run. Root cause: the
  bootstrap CI's lower bound can end up ABOVE the point estimate when
  many segments are shorter than the block length (they resample
  near-deterministically, biasing the bootstrap distribution - a
  limitation already flagged as a known caveat, this is its first
  concrete real-data consequence). Fixed: proper asymmetric error bars
  from both `ci_low`/`ci_high` (old code only used `ci_low`, applied
  symmetrically), negative bar lengths clipped to 0 with an explicit
  warning naming which regime(s) it happened to, instead of crashing.
  Reproduced the exact crash with a crafted CI before fixing, confirmed
  fixed after.

### Fixed — batch orchestration
- **`run_all_buoys.py`**: Stage 08 was never connected to Stage 11b's
  persistence estimate in the batch path (only in manual single-buoy
  runs) - every buoy in the first full multi-year batch used the
  generic 48h declustering default regardless of what 11b had already
  computed for that specific buoy. Added `build_stage08_args()`, same
  dynamic-dispatch pattern as Stage 10's `--include-period`. Verified
  via actual logged command lines that each buoy gets its own window
  (confirmed different per buoy in a smoke test: 68.9h/66.4h/60.4h).
  Concrete symptom this fixed: A2Buoy's GPD xi CI was [-0.198, 0.065]
  (crosses zero) at the wrong window.

### Fixed — performance at multi-year scale
- **`12_confidence_intervals.py`**: Weibull bootstrap refit used full
  MLE (`scipy.stats.weibull_min.fit()`) 1000x per buoy - ~2s/fit at
  500k+ samples meant 30-40+ minutes per buoy, silently (orchestrator
  buffers stdout until a stage exits), stalling the first full batch
  run outright. Fixed: `fast_weibull_moment_fit()` - closed-form
  method-of-moments via the coefficient-of-variation relation, used
  only inside the bootstrap loop (Stage 06's point estimate still uses
  full MLE) - validated at 0.02% relative error vs. MLE, ~630x faster.
  `--n-bootstrap` also now auto-scales down for large records (200 at
  >=500k samples, 500 at >=50k, 1000 otherwise). A 600k-sample test
  that would have taken ~33 minutes now runs in ~11 seconds.

### Added
- **`14_mann_kendall_trend.py`** - Mann-Kendall trend test with Hamed-Rao
  variance correction for autocorrelation (not the textbook version -
  Ljung-Box's confirmed serial correlation at every buoy is exactly the
  condition known to inflate plain Mann-Kendall's false-positive rate).
  Runs on ANNUAL aggregates (mean and p95 Hs separately - typical
  climate vs. storm intensity), not raw high-frequency data - both a
  computational necessity (raw MK is O(n^2), infeasible at 500k+
  samples) and standard climate-trend-detection practice. Validated
  three ways before shipping: recovered a known injected +0.02 m/year
  trend exactly (+0.0200 m/year, p=0.0043); correctly found no trend on
  an autocorrelated no-trend control; correctly refuses to run below
  `--min-years` (default 10) rather than producing an unreliable result.
  Standalone script, not yet wired into `run_all_buoys.py` - only
  meaningful on the 6 buoys with 30+ year records.
- **`download_belgian_wave_buoys_multiyear.py`** superseded by
  `download_belgian_wave_buoys_history.py` after discovering
  `copernicusmarine.subset()` only ever serves the `latest` (30-day)
  dataset part for in-situ products, regardless of requested date
  range - the full archive needs the Files service
  (`copernicusmarine.get()` with `dataset_part="history"`), a different
  API entirely. Kept in the repo for the record of what was tried.

### Multi-year data - acquired
- Full 19-buoy history download via `get()` + `dataset_part="history"`:
  9778 files across the whole North West Shelf region in the raw
  listing, grepped down to the 19 Belgian buoys by name
  (`history/MO/NO_TS_MO_<name>.nc` naming, distinguished from several
  unrelated Wind/Weather/Tide/MP files at the same sites) - 81.39 MB for
  all 19 vs. 21.38 GB for the unfiltered whole-region pull.
- **Real per-buoy coverage, now confirmed**: 6 buoys back to 1990-1997
  (Westhinder longest, 1990-07-19); 13 buoys from 2009-2021. All extend
  to 2026-06-30.
- ERA5 meteo (wind/pressure/temperature, 2010-2026, 3-hourly,
  monthly-chunked) downloaded in parallel - not yet paired with the
  buoy data in any analysis.

### Findings - full 19-buoy batch run against real multi-year data
- Findings reproduce on completely different data than the original
  2-month NRT window: Zeebrugge's singleton spatial cluster, ~99%
  effective-N reduction from autocorrelation, correlation-decreasing-
  with-distance (Spearman r=-0.506) - all close to the NRT-window
  values.
- **Westhinder shows two independent signals of possible genuine
  multi-decade drift** (persistence-timescale spread; Stage 13 mixture
  artifact with a monotonic +3.2% climb in era means) - **Zeebrugge
  shows none** (100% window agreement) - the contrast suggests this is
  a real, site-specific phenomenon, not a pipeline artifact common to
  all long records. Not yet resolved.
- A2Buoy's tail is close to exponential (xi ~ 0) even after the Stage
  08/11b window correction - a genuine finding now, not a methodology
  artifact, though its updated confidence interval still needs checking.

## 2026-07-14 — Gap-splicing generalized, Priority 4 closed, wind coupling, full forecasting pipeline (Stages 14-21)

### Fixed — gap-splicing, extended beyond the lag-based stages
- **`13_stability_analysis.py` `[C]` crash**: `ValueError: 'yerr' must
  not contain negative values` - crashed CadzandBoei/Deurlo outright.
  Root cause: bootstrap CI's lower bound can land above the point
  estimate when many segments are shorter than the block length
  (near-deterministic resampling biases the distribution). Fixed with
  proper asymmetric error bars (old code only used `ci_low`, applied
  symmetrically) and negative-width clipping with an explicit warning.
  Reproduced the exact crash with a crafted CI before fixing.
- **`utils.segments_by_time_gap()`** added - for series where gap
  positions are entirely ABSENT (missing rows) rather than NaN-marked,
  e.g. Stage 10's regime labels. Cross-validated: found the same 1905
  segments on Westhinder that NaN-based detection found independently
  on the raw Hs series.
- **`utils.integral_timescale()`** and **`select_order()`** both moved
  out of their original stage scripts (`11b_dependence_structure.py`,
  `18_forecast_arma.py`) into shared modules (`utils.py`,
  `forecast_utils.py`) so later stages could reuse them without
  duplicating logic - Python module names can't start with a digit, so
  direct import from a numbered stage file isn't possible. Verified no
  regression in both source stages after each move.

### Priority 4 (Zeebrugge tidal-frequency) — CLOSED
- `03b_tidal_notch.py` extended with `--fit-frequency`: finds the
  actual dominant period near M2 instead of assuming it exactly, via a
  periodogram search on the longest contiguous segment. Validated on
  synthetic data with an injected S2 (not M2) signal - fixed-M2 made
  things worse (669162->669434), frequency-fitting recovered the true
  period exactly and cleaned the signal ~59,000x.
- **Real result: hypothesis rejected.** Fitted period on real Zeebrugge
  data: 12.4163h - essentially exact M2 (delta -0.0043h), nowhere near
  S2. Even with the correct frequency, the notch measurably worsened
  (1142->1592). Rules out "wrong fundamental frequency" as the
  explanation. Reverted to the harmonics=2 non-fitted default.
- Joint M2+diurnal notch also tried (`--diurnal-harmonics N`) after a
  separate diurnal-signal investigation (below) - real but partial,
  inconsistent improvement on 3 test buoys (5-80% reduction, never
  below the clean threshold) - not adopted, flag left in the codebase
  as validated-but-unused infrastructure.

### Priority 5 (wind-wave coupling) — CLOSED, full network run
- `utils.load_era5_for_buoy()` added - nearest-grid-cell extraction,
  concatenated across monthly ERA5 files, derived wind_speed and
  wind_dir_from_deg (meteorological convention, hand-verified against
  all 4 cardinal directions before trusting it).
- `16_wind_wave_coupling.py` added - wind/MSLP CCF with Hs, R² summary,
  directional alignment by regime (VMDR buoys only). **Caught and fixed
  a real CCF lag-sign bug before it ran on real data** - empirically
  verified (known-lag synthetic test) that positive lag means x leads
  y, opposite of the first draft's assumption.
- **Full 19-buoy run**: consistent signature at 18/19 buoys (wind leads
  0-3h, R² 0.41-0.84, MSLP correctly leads). Zeebrugge collapses
  (R²=0.093, MSLP mechanism inverts) - see Zeebrugge section above.
- **Diurnal signal investigation**: Stage 02 extended with a general
  24h-band power-ratio check (validated on synthetic injected-signal
  data) - found a real ~24h Hs signal at every buoy tested (85-289x
  baseline), not Zeebrugge-specific. Checked whether ERA5 wind itself
  shows the same cycle (land-sea-breeze hypothesis) - confirmed at all
  4 buoys tested (75-1986x). Building the wind persistence-timescale
  check to explain a related ARMAX finding surfaced two more real bugs
  in sequence: raw wind gave an implausible 343h timescale (seasonal
  contamination, same category of mistake 11b itself originally made);
  first fix (rolling-mean deseasonalizing) gave a mathematically
  impossible NEGATIVE timescale (spurious autocorrelation artifact of
  rolling-mean subtraction); second fix (harmonic regression at the
  annual period, reusing Stage 03b's already-trusted method) gave
  57.3h - validated on a synthetic AR(1)+seasonal test with a known
  answer (57.0h true) before trusting the real result.

### Priority 6 (forecasting) — Stages A-F all built, tested, closed
- **`forecast_utils.py`** (Stage B, shared harness) +
  **`17_forecast_baseline.py`** (Stage A, persistence). Caught a real
  off-by-one (origin semantics) via an analytical ground-truth test on
  a deterministic ramp - fixed before it ran on real data. Separately
  caught and fixed a real O(n^2) performance bug on the actual
  Westhinder run (5+ minutes, still not finished) - traced to passing
  the full ever-growing history to a model that only needed a small
  bounded window; fixed via a `history_window` parameter, verified
  17s at 630k samples after the fix (down from 5+ min).
- **`18_forecast_arma.py`** (Stage C) - beats persistence at every
  horizon on Westhinder (0.10-0.30 skill, growing with horizon).
  Caught a real design bug: hardcoded `d=1` (copied from Stage 04,
  tuned for a different series) gave NEGATIVE skill at every horizon
  beyond 1h on A2Buoy; fixed by searching `d` instead of assuming it -
  default changed to `d=0` only after confirming the fix reversed
  every negative result to positive.
- **`19_forecast_arma_garch.py`** (Stage D) - point forecast + 95% PI.
  Real Westhinder fit landed almost exactly at the IGARCH boundary
  (alpha+beta=1.0000000000273) - required two rounds of fixing the
  multi-step variance formula for this regime (first attempt used a
  flat forecast, which was itself WRONG and caused a real regression,
  caught by re-running on real data; corrected to the proper linear-
  growth IGARCH limit, verified via exact match against direct
  numerical iteration of the true recursion, not just re-derived by
  eye a second time).
- **`20_forecast_armax.py`** (Stage E) - ARMA + lagged wind exog.
  Explicit honest-forecasting design: exog only genuinely known within
  the wind-Hs lag window, persisted (not fabricated) beyond it -
  validated on synthetic data showing the predicted signature (skill
  peaks at the lag boundary, reverses negative beyond it). Real
  Westhinder result diverged informatively: skill stayed positive
  through 24h rather than reversing, explained (and confirmed) by
  wind's own 57.3h persistence timescale extending well beyond the
  tested horizon range.
- **`21_forecast_exceedance.py`** (Stage F) - probabilistic "will Hs
  exceed threshold X in next Nh" classification (Brier score, ROC-AUC),
  optional wind feature. Validated on two synthetic scenarios (Hs-
  driven and wind-driven exceedance) before trusting real results.
  Real result on Westhinder: skill +0.42 to +0.68, largely reflecting
  the buoy's already-known extreme persistence recast into a
  classification frame. Wind-augmentation tested on 2 real buoys
  (Westhinder, A2Buoy) - helps consistently at every threshold/horizon
  (12/12 positive deltas), benefit does not fade with horizon at either
  buoy - though the specific mechanistic explanation proposed for
  Westhinder's flat pattern didn't survive testing against A2Buoy
  (which showed a GROWING delta, opposite of the prediction) - revised
  to a more general, better-supported explanation involving the
  relative information decay of Hs's own autoregressive features vs.
  wind's broader synoptic relevance.
