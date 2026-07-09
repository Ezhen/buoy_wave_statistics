# Belgian Coastal Wave Buoy — Statistical Characterization Pipeline

Characterizes the statistical regime of Belgian coastal wave buoy data
(CMEMS in-situ, NWS regional product) — stationarity, tidal contamination,
Hs distribution, volatility clustering, storm extremes, regime structure,
and spatial coherence across the network. Validated on all **19** BCZ
buoys (pilot: Westhinder, offshore reference, cleanest signal).

This does **not** yet include forecasting (AR/ARMA/GARCH) — see
"Where this stops" below. Next steps are tracked in
`PLAN_next_session.md`.

## Data

Downloaded via `download_belgian_wave_buoys.py` from CMEMS product
`INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036` (requires a free
data.marine.copernicus.eu account — separate from any CDSE/Sentinel
credentials). Produces one `.nc` file per buoy in `data/`.

**Variable availability is NOT uniform across the network** — confirmed
via `ncdump -h` across all 19 files:
- **17/19 buoys**: `VHM0` only.
- **CadzandBoei**: `VHM0`, `VTPK`.
- **Deurlo**: `VHM0`, `VTPK`, `VMDR`.

These same two buoys also run at **10-min sampling** vs. everyone else's
30-min (see Stage 0 findings below) — three independent signals (variable
richness, sampling rate, and whatever else differs upstream) point at
Cadzand/Deurlo being a genuinely different instrument class, not just two
stations with a data gap. This is why Stage 09 (cross-variable analysis)
is excluded from the default batch loop — see Stage 09 entry below.

## Directory layout

```
wave_pipeline/
├── utils.py                      # shared buoy-loading helpers
├── 01_load_clean.py
├── 02_eda_diagnostics.py
├── 03_stationarity_tests.py
├── 03b_tidal_notch.py
├── 04_transform_detrend.py
├── 05_whiteness_check.py
├── 06_distribution_fit.py
├── 07_arch_lm_test.py
├── 08_extreme_value_analysis.py
├── 09_cross_variable_analysis.py  # NOT in the default batch - see notes
├── 10_regime_identification.py
├── 11_spatial_statistics.py       # runs once, network-wide, not per-buoy
├── run_all_buoys.py               # orchestrator: runs 01-08,10 per buoy, then 11 once
├── summarize_results.py           # builds one cross-buoy comparison table
├── rerun_eva_all_buoys.py         # re-run ONLY Stage 08 at a different threshold
├── PLAN_next_session.md           # roadmap for what's next
data/                              # input .nc files (one per buoy)
pipeline_out/                      # all outputs, one subfolder per stage
├── 01_load_clean/ ... 11_spatial_statistics/
├── batch_run.log                  # appended, timestamped, from run_all_buoys.py
└── bcz_comparison_summary.csv     # one row per buoy, from summarize_results.py
```

Each stage reads the previous stage's output from `pipeline_out/` and
writes its own. Stages 01, 02, 03b, 04, and 08 also write a small
`*_summary.json` alongside their plots — that's what `summarize_results.py`
reads to build the comparison table.

## Running it — single buoy

```bash
python 01_load_clean.py       --data-dir data --buoy WesthinderBuoy --var VHM0
python 02_eda_diagnostics.py  --buoy WesthinderBuoy --var VHM0
python 03_stationarity_tests.py --buoy WesthinderBuoy --var VHM0
python 03b_tidal_notch.py     --buoy WesthinderBuoy --var VHM0 --harmonics 2
python 04_transform_detrend.py --buoy WesthinderBuoy --var VHM0 --diff-order 1
python 05_whiteness_check.py  --buoy WesthinderBuoy --var VHM0
python 06_distribution_fit.py --buoy WesthinderBuoy --var VHM0
python 07_arch_lm_test.py     --buoy WesthinderBuoy --var VHM0
python 08_extreme_value_analysis.py --buoy WesthinderBuoy --var VHM0 --threshold-percentile 85
python 10_regime_identification.py --data-dir data --buoy WesthinderBuoy --var VHM0
```

Stage 09 only produces real output for Cadzand/Deurlo:
```bash
python 09_cross_variable_analysis.py --data-dir data --buoy CadzandBoei
python 09_cross_variable_analysis.py --data-dir data --buoy Deurlo
```

## Running it — all buoys at once

```bash
python run_all_buoys.py --data-dir data --var VHM0
python summarize_results.py --var VHM0
```

`run_all_buoys.py` runs stages **01→08 and 10** per buoy (09 intentionally
excluded — see Data section), then runs **Stage 11 once** across the whole
network after the per-buoy loop finishes (only if ≥3 buoys succeeded). If
a stage fails for a given buoy, it logs the failure and moves to the next
buoy rather than aborting the batch. Everything is appended, timestamped,
to `pipeline_out/batch_run.log`.

`summarize_results.py` reads every stage's summary file across all buoys
into `pipeline_out/bcz_comparison_summary.csv` — sampling interval, M2
ratio before/after notch, ADF/KPSS agreement, Ljung-Box whiteness
(including real-time `ljungbox_lags_hours`, since sampling intervals
differ), best-fit distribution, ARCH significance, and EVA shape
parameter with a `fit_reliable` flag. It also **warns explicitly** if
sampling intervals differ across buoys, since that breaks direct
Ljung-Box comparability between buoys.

`rerun_eva_all_buoys.py` re-runs *only* Stage 08 at a different threshold
across all buoys with existing Stage 0 output — no need to redo the full
pipeline just to retune EVA's threshold percentile.

## What each stage does

| Stage | Script | Purpose | Key output |
|---|---|---|---|
| 0 | `01_load_clean.py` | Regularize sampling grid, interpolate short gaps, flag long gaps as NaN, sanity-bound VHM0 ≥ 0, report longest gap and duplicate-timestamp count | `_clean.csv` + `_load_summary.json` |
| 1 | `02_eda_diagnostics.py` | Rolling mean/variance, ACF/PACF, periodogram with M2 tidal frequency marked | plots + `_eda_summary.json` |
| 2 | `03_stationarity_tests.py` | ADF + KPSS on the raw series, side by side | `_stationarity.json` |
| 2b | `03b_tidal_notch.py` | Box-Cox, then harmonic regression at M2 (+ harmonics) to remove the tidal signature before differencing | `_detided_boxcox.csv` + `_notch_summary.json` |
| 3 | `04_transform_detrend.py` | Order-N differencing on the detided series | `_residual.csv` + `_detrend_summary.json` |
| 4 | `05_whiteness_check.py` | Ljung-Box on the Stage 3 residual | `_ljungbox.csv` |
| 4b | `06_distribution_fit.py` | Rayleigh/Weibull/log-normal fit to the **raw level series** | `_fit_summary.csv` + plots |
| 5 | `07_arch_lm_test.py` | Engle's ARCH-LM test on the Stage 3 residual (volatility clustering) | `_arch_lm.csv` |
| 6 | `08_extreme_value_analysis.py` | POT: decluster storm peaks, fit GPD, extrapolate return levels | `_return_levels.csv` + `_eva_summary.json` |
| 7 | `09_cross_variable_analysis.py` | Hs/Tp/Tm02/direction correlation, lagged cross-correlation, PCA — **only meaningful for CadzandBoei/Deurlo** | correlation matrix, CCF, PCA loadings |
| 8 | `10_regime_identification.py` | GMM clustering into calm/moderate/energetic/storm; `--include-period` falls back to Hs-only gracefully if VTPK is absent | `_regime_labels.csv` + `_regime_summary.csv` |
| 9 | `11_spatial_statistics.py` | Network-wide: pairwise Hs correlation, correlation-vs-distance (haversine), hierarchical buoy clustering | `_pairwise_correlation.csv`, `_dendrogram.png`, `_buoy_clusters.csv` |

**Why distribution/EVA (06/08) use the raw level series, not the Stage 3
residual:** Rayleigh/Weibull/log-normal/GPD all describe Hs itself
(positive, right-skewed). A differenced residual is mean-zero and can go
negative — those distributions don't apply to it.

## Findings — network-wide (19 buoys)

- **Sampling heterogeneity**: 17/19 buoys at 30-min; **CadzandBoei and
  Deurlo at 10-min** — same two buoys with richer variable sets (see
  Data section). `summarize_results.py` flags this automatically since it
  breaks direct Ljung-Box lag comparison across buoys.
- **M2 tidal contamination has a real spatial gradient**: raw ratio from
  ~6 (OstendEasternPalisade) to **135-211× (ZeebruggeZandopvangkade,
  varies by run)** — physically consistent with harbor/shallow-water
  shoaling. The 2-harmonic notch never fully cleans any buoy, worst at
  Zeebrugge — consistent with shallow-water overtides needing more
  harmonics than 2 there specifically.
- **Stationarity**: ADF/KPSS disagree at every buoy (ADF: stationary,
  KPSS: not) — consistent with strong periodicity, not genuine drift.
  Several buoys' KPSS statistics exceed the lookup table range entirely
  (`InterpolationWarning`) — stronger version of the same finding, not a
  bug.
- **Ljung-Box fails at all 19 buoys, no exceptions** — confirms storm
  persistence (hours-to-days AR memory) is a network-wide physical
  property, not a single-station artifact. Don't chase this to whiteness
  via more differencing/harmonics.
- **Distribution fit**: Weibull wins 17/19, log-normal wins 2
  (Raversijde2, Zeebrugge).
- **ARCH-LM**: significant at all 19 buoys — volatility clustering is
  universal, not a Westhinder quirk.
- **EVA**: at the default 95th-percentile threshold, only 1/19 buoys had
  enough peaks (>=10) for a reliable fit; shape parameters swung wildly
  (-3.68 to +0.14 - small-sample noise). **Fixed** by lowering to the
  85th percentile / 24h min separation (`rerun_eva_all_buoys.py`):
  **19/19 buoys now reliable**, xi ranges narrowly from -1.31 to -0.04,
  **all negative** (bounded upper tail everywhere - physically sensible
  for a depth/fetch-limited shelf sea). Still treat any return level
  beyond ~4 months (2x record length) as illustrative regardless of
  threshold tuning.
- **Cross-variable analysis (Stage 09)**: only ever produces real output
  for CadzandBoei and Deurlo - the other 17 buoys structurally lack the
  variables. Not a bug; excluded from the default batch loop for this
  reason.
- **Spatial statistics (Stage 11)**: correlation decreases with distance
  (Spearman r ≈ -0.47) - textbook spatial coherence. 4-cluster structure
  is geographically sensible; **Zeebrugge consistently forms its own
  singleton cluster** - independently corroborated by its extreme M2
  ratio and its (coincidentally) reliable EVA fit at the old threshold.
  Three unrelated analyses agreeing Zeebrugge is dynamically distinct
  from the rest of the coast is a solid, cross-validated finding.

## Where this stops

This pipeline characterizes the **regime** - it does not forecast, and
does not yet report uncertainty on any of its point estimates (planned
next: see `PLAN_next_session.md`, Priority 1 - bootstrap/CI stages,
particularly for the EVA shape parameter, before trusting cross-buoy
comparisons of it too literally). The AR/ARMA(-GARCH) forecasting stage
discussed separately would use the Stage 3 residual differently (as
something to model, not whiten away) and is intentionally not started
until the uncertainty stages exist.
