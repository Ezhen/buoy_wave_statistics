# Changelog

All notable changes to this repository, in the order they actually
happened. Dates reflect real batch-run timestamps where available.

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
