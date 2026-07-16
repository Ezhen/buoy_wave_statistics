# Belgian Coastal Wave Buoy — Statistical Characterization Pipeline

Characterizes the statistical regime of Belgian coastal wave buoy data
(CMEMS in-situ, NWS regional product) across all **19** BCZ buoys:
stationarity, tidal contamination, Hs distribution, volatility
clustering, storm extremes, dependence structure, uncertainty
(confidence intervals + stability), regime structure, and spatial
coherence.

Does **not** forecast — see "Where this stops" below.
Full build history: `CHANGELOG.md`. Open priorities and next steps:
`PLAN_next_session.md`.

## Data

Three separate data sources, each with its own download script and its
own caveats — worth reading before assuming any of them "just work" the
way you'd expect:

**1. NRT wave data** (`download_belgian_wave_buoys.py`) — the original,
fast download. Product `INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036` via
`copernicusmarine.subset()`. **Only ever returns a rolling ~30-day
window**, regardless of date range requested — this is a hard limitation
of the `subset()` API for in-situ products (it only serves the `latest`
dataset part), not a bug in the script. Saves to `data/`.

**2. Multi-year wave history** (`download_belgian_wave_buoys_history.py`)
— the real fix for #1. Uses `copernicusmarine.get()` with
`dataset_part="history"` instead of `subset()` — a different API/service
entirely. Two-step workflow:
```bash
python download_belgian_wave_buoys_history.py --list       # preview filenames, 0 download
python download_belgian_wave_buoys_history.py --download --file-list belgian_19_file_list.txt
```
The `history` part covers the *whole* North West Shelf region (~9800
files, ~21 GB) — don't run `--download` without `--file-list` or
`--regex` scoping it down first. The script also runs its own lat/lon
safety filter after downloading. Saves to `data_multiyear/`.
**Real per-buoy coverage is highly uneven**: 6 buoys go back to
1990-1997 (Westhinder longest, from 1990-07-19); the other 13 range from
2009 to 2021. All extend to 2026-06-30.

**3. ERA5 meteo** (`download_era5_meteo.py`) — wind (u/v @ 10m, the
primary wave-generation driver), mean sea-level pressure (leading
indicator), 2m air temp + SST (their difference affects drag
coefficient — second-order). Via CDS (`cds.climate.copernicus.eu` — a
**third**, separate Copernicus portal/account, distinct from both CDSE
and the CMEMS account used above). Chunked by month, defaults to
3-hourly (justified by Stage 11b's ~50-100h persistence finding — hourly
would be unused resolution). Saves to `meteo_era5/`.

**Variable availability is NOT uniform across the network**: 17/19
buoys report only `VHM0`. Only **CadzandBoei** (+`VTPK`) and **Deurlo**
(+`VTPK`, `VMDR`) carry the fuller set — the same two buoys that also
run at 10-min sampling instead of 30-min, suggesting a genuinely
different instrument class, not a coincidental gap.

## Directory layout

```
wave_pipeline/
├── utils.py                         # shared helpers
├── download_belgian_wave_buoys.py           # NRT (~30-day window)
├── download_belgian_wave_buoys_history.py   # multi-year (real fix)
├── download_belgian_wave_buoys_multiyear.py # superseded, kept for history - see CHANGELOG
├── download_era5_meteo.py                   # meteo
├── plot_belgian_buoys.py            # per-buoy plots + location map
├── 01_load_clean.py  ... 14_mann_kendall_trend.py   # pipeline stages, see table below
├── run_all_buoys.py                 # orchestrator with tiering gate
├── summarize_results.py             # cross-buoy comparison table
├── rerun_eva_all_buoys.py           # re-run ONLY Stage 08 at a new threshold
├── CHANGELOG.md                     # full build history
├── PLAN_next_session.md             # open priorities / roadmap
data/                                # NRT download (~30-day window)
data_multiyear/                      # multi-year download (real history)
meteo_era5/                          # ERA5 download
pipeline_out/                        # all stage outputs, one subfolder per stage
├── 01_load_clean/ ... 13_stability_analysis/
├── batch_run.log                    # appended, timestamped
└── bcz_comparison_summary.csv       # one row per buoy
```

## Running it

**All buoys, full pipeline:**
```bash
python run_all_buoys.py --data-dir data --var VHM0
python summarize_results.py --var VHM0
```
Point `--data-dir` at `data_multiyear` once you're ready to run against
real historical coverage instead of the 30-day NRT window.

**One buoy, stage by stage** (useful for debugging a specific stage):
```bash
python 01_load_clean.py       --data-dir data --buoy WesthinderBuoy --var VHM0
python 02_eda_diagnostics.py  --buoy WesthinderBuoy --var VHM0
python 03_stationarity_tests.py --buoy WesthinderBuoy --var VHM0
python 03b_tidal_notch.py     --buoy WesthinderBuoy --var VHM0 --harmonics 2
python 11b_dependence_structure.py --buoy WesthinderBuoy --var VHM0
python 04_transform_detrend.py --buoy WesthinderBuoy --var VHM0 --diff-order 1
python 05_whiteness_check.py  --buoy WesthinderBuoy --var VHM0
python 06_distribution_fit.py --buoy WesthinderBuoy --var VHM0
python 07_arch_lm_test.py     --buoy WesthinderBuoy --var VHM0
python 08_extreme_value_analysis.py --buoy WesthinderBuoy --var VHM0 --threshold-percentile 85
python 12_confidence_intervals.py --buoy WesthinderBuoy --var VHM0
python 10_regime_identification.py --data-dir data --buoy WesthinderBuoy --var VHM0
python 13_stability_analysis.py --buoy WesthinderBuoy --var VHM0
```
Network-wide, after the per-buoy loop:
```bash
python 11_spatial_statistics.py --data-dir data --var VHM0
python 12b_correlation_confidence.py --var VHM0
```
Only for CadzandBoei/Deurlo (the two buoys with VTPK/VMDR):
```bash
python 09_cross_variable_analysis.py --data-dir data --buoy CadzandBoei
```

**Mann-Kendall trend test** (Stage 14, standalone — not yet wired into
`run_all_buoys.py`; only meaningful on the 6 buoys with 30+ year records):
```bash
python 14_mann_kendall_trend.py --buoy WesthinderBuoy --var VHM0
```

`run_all_buoys.py` runs Stage 0 first, unconditionally, per buoy — it's
what populates the tiering gate's eligibility checks (see below) before
deciding which stages to run for that buoy.

## Tiering gate

`run_all_buoys.py` checks each stage's declared requirement against
what Stage 0 found for that buoy, rather than a stage crashing or
silently degrading:
- **Core** (no requirement): everything except 09 and 10's
  `--include-period` flag.
- **Advanced** (`variables_any`/`variables_all`): Stage 09 needs
  `VTPK` or `VMDR`; Stage 10's `--include-period` needs `VTPK`. Both
  decided by the gate, not the script's own fallback.
- **Dynamic argument injection** (not a tier, but the same "decide
  automatically, don't hardcode" principle): Stage 08's
  `--min-separation-hours` is set per-buoy from Stage 11b's own
  persistence estimate, not a shared default — found missing after the
  first full multi-year batch run used the same 48h window for all 19
  buoys despite 11b already having computed a real per-buoy-justified
  number (e.g. 184.8h for A2Buoy). See `build_stage08_args()`.
- **Multi-year** (`min_record_years`): mechanism ready, unused so far —
  no stage currently requires a minimum record length.

A `SKIPPED` stage (requirement not met) is logged distinctly from a
`FAILED` one in both console output and `batch_run.log`.

## Gap handling

Multi-year buoy records have real gaps (sensor outages, maintenance) —
Stage 0 leaves long gaps as NaN rather than bridging them (bridging
would fabricate structure across a gap that could be months long).
Every stage that computes a **lag/order-based** statistic (ACF,
differencing, unit-root tests, block bootstrap) has to actively avoid
treating the sample right before a gap as "adjacent" to the sample
right after it — a naive `.dropna()` silently splices them together as
if temporally adjacent, which corrupts cumulative statistics far more
than it corrupts a point estimate. This bit multiple stages before
being fixed (see `CHANGELOG.md`): 11b's persistence estimate came out
backwards on the 2-month data, then wrong again in a different way on
the 36-year data; Stage 13's moving-window split produced misleading
window boundaries; Stage 13's regime-fraction bootstrap needed a
different detection method entirely, since Stage 10's label output has
gap rows entirely absent rather than NaN-marked.

Three helpers in `utils.py` handle this, used differently depending on
what a stage needs:
- **`longest_contiguous_segment()`** — NaN-gridded series (Stage 0's
  `_clean.csv` and anything derived from it while preserving the grid,
  e.g. Stage 03b's detided output). Used by 02/03/04 for their
  one-shot lag-based computations.
- **`all_contiguous_segments()`** — same NaN-gridded input, but returns
  every qualifying segment instead of just the longest. Used by 11b to
  aggregate the persistence estimate across a fragmented record
  (length-weighted mean) rather than discarding everything except one
  segment — on Westhinder, the single longest segment covered only
  6.9% of valid data; aggregating across 50 segments covered 61.4%.
- **`segments_by_time_gap()`** — for series where gap positions are
  entirely ABSENT (missing rows) rather than NaN-marked, e.g. Stage
  10's regime labels. Detects a break wherever consecutive timestamps
  are further apart than expected, rather than relying on NaN.

A real cross-check this produced: `segments_by_time_gap` on
Westhinder's regime labels found **exactly 1905 segments** — the same
count the NaN-based methods independently found on the raw Hs series.
Two different detection methods agreeing on a fragmented 36-year record
is a solid confirmation both are correctly recovering the same gap
structure, not coincidence.

## Performance at multi-year scale

Stage 12's block-bootstrap Weibull refit originally used full
`scipy.stats.weibull_min.fit()` (MLE) inside the 1000-iteration loop —
fine at ~2850 samples (2-month NRT window), but ~2s/fit at 500k+
samples meant 30-40+ minutes per buoy once real multi-year data
arrived, silently (no visible output mid-stage, since the orchestrator
buffers stdout until a stage exits) — this stalled the first full
19-buoy batch run outright. Fixed two ways: (1) `fast_weibull_moment_fit()`
— closed-form method-of-moments via the coefficient-of-variation
relation, used only inside the bootstrap loop (Stage 06's actual point
estimate still uses full MLE) — validated at 0.02% relative error vs.
MLE, ~630x faster; (2) `--n-bootstrap` auto-scales down for large
records (200 at ≥500k samples, 500 at ≥50k, 1000 otherwise), since
bootstrap CI precision barely improves past a few hundred resamples
once the underlying sample size is already huge. A 600k-sample test
that would have taken ~33 minutes now runs in ~11 seconds.

## What each stage does

| Stage | Script | Purpose |
|---|---|---|
| 0 | `01_load_clean.py` | Regularize grid, interpolate short gaps, sanity-bound, record available variables/record length |
| 1 | `02_eda_diagnostics.py` | Rolling stats, ACF/PACF, periodogram with M2 marked |
| 2 | `03_stationarity_tests.py` | ADF + KPSS side by side |
| 2b | `03b_tidal_notch.py` | Box-Cox + harmonic regression at M2 (+harmonics) |
| — | `11b_dependence_structure.py` | Integral (persistence) timescale from the detided series' ACF — feeds block-bootstrap length and EVA declustering justification |
| 3 | `04_transform_detrend.py` | Order-N differencing on the detided series |
| 4 | `05_whiteness_check.py` | Ljung-Box on the residual |
| 4b | `06_distribution_fit.py` | Rayleigh/Weibull/log-normal on the **raw level series** |
| 5 | `07_arch_lm_test.py` | ARCH-LM volatility clustering test |
| 6 | `08_extreme_value_analysis.py` | POT + GPD fit, return levels |
| — | `12_confidence_intervals.py` | Block bootstrap on Hs mean/quantiles; bootstrap on GPD xi; CI band on Weibull fit |
| 7 | `09_cross_variable_analysis.py` | Hs/Tp/direction correlation, lag-CCF, PCA (Advanced tier) |
| 8 | `10_regime_identification.py` | GMM regime clustering |
| — | `13_stability_analysis.py` | Moving-window stability, drop-biggest-storm jackknife, regime-fraction bootstrap |
| 9 | `11_spatial_statistics.py` | Network-wide correlation, correlation-vs-distance, clustering |
| — | `12b_correlation_confidence.py` | Fisher z CI on Stage 11's correlations, autocorrelation-corrected |
| — | `14_mann_kendall_trend.py` | Hamed-Rao-corrected trend test on annual mean/p95 Hs — standalone, not in `run_all_buoys.py`; only meaningful on 30+ year records (6 of 19 buoys) |
| — | `15_seasonal_decomposition.py` | STL seasonal decomposition on monthly mean/p95 Hs — same 6-buoy scope as Stage 14 |
| — | `16_wind_wave_coupling.py` | ERA5 wind/MSLP cross-correlation with Hs, diurnal + persistence-timescale checks on wind itself, directional alignment (VMDR buoys only) |
| — | `17_forecast_baseline.py` | Persistence forecast baseline + shared rolling-origin backtest harness (`forecast_utils.py`) |
| — | `18_forecast_arma.py` | ARMA point forecast, backtested against local persistence |
| — | `19_forecast_arma_garch.py` | ARMA-GARCH: point forecast + calibrated prediction interval |
| — | `20_forecast_armax.py` | ARMAX: ARMA + lagged ERA5 wind speed as exogenous regressor |
| — | `21_forecast_exceedance.py` | Probabilistic "will Hs exceed threshold X in the next Nh" classification, optional wind feature |

**Why 06/08/12 use the raw level series, not the residual**:
Rayleigh/Weibull/GPD describe Hs itself (positive, right-skewed); a
differenced residual is mean-zero and can go negative.

**Why 11b runs on the detided (03b) series, not raw**: a strongly
periodic (tidal) component makes the ACF cross zero at fractions of the
tidal period — that's the tide's own oscillation, not evidence of losing
memory. (v1 of this stage got this wrong — see `CHANGELOG.md`.)

**Why 12/13's uncertainty methods are autocorrelation-aware, not
textbook-default**: Ljung-Box fails at every buoy, confirmed
network-wide — an ordinary/IID bootstrap or a classical Fisher z CI on
this data would understate uncertainty. Block length and the
effective-N correction both come from Stage 11b.

## Findings summary

- **Tidal contamination** has a real spatial gradient (6x-135x raw M2
  ratio), worst at Zeebrugge — whose notch never fully cleans regardless
  of harmonic count (see Zeebrugge section below).
- **Ljung-Box fails at all 19 buoys** — network-wide storm persistence,
  confirmed physical, not a pipeline artifact.
- **Weibull wins 17/19** buoys for Hs distribution; log-normal the other 2.
- **ARCH effects significant at all 19 buoys** — volatility clustering
  is universal.
- **EVA declustering window is per-buoy, not a fixed default** — each
  buoy's Stage 08 run uses its own Stage 11b persistence-based window
  (typically 100-270h on the real multi-year data), not a generic
  round-number default. This was a real fix made after discovering the
  first full multi-year batch had used the same 48h window network-wide
  regardless of what 11b had already computed per buoy — see
  `CHANGELOG.md`.
- **GPD shape parameters, corrected network run: -0.544 to +0.327** —
  tightened dramatically after the declustering fix above (previously a
  much wider, partly-artifactual range under the old fixed-window
  approach). The reliable core of the network (many storm peaks, narrow
  CIs — Trapegeer, ScheurWielingen, AkkaertSouthwest, Wandelaar)
  clusters tightly around **-0.27 to -0.32** (bounded upper tail),
  roughly 2x more negative than a published shallow-water North Sea
  reference (Caires 2011, xi≈-0.12 to -0.13) — plausibly because BCZ
  buoys are shallower/more depth-limited than that 19m reference site,
  not a methodology gap.
  **Blankenberge is a real outlier** — the only positive xi in the
  network (+0.327, CI crosses zero), with both the fewest storm peaks
  (43) and shortest record (2.77 years) of any buoy. Positive xi
  implies an unbounded tail, physically implausible for depth-limited
  shallow-sea waves — near-certainly a small-sample artifact, not a
  genuine finding. Don't trust this buoy's xi at face value.
- **GPD xi confidence intervals vary hugely by peak count** — some
  cross zero (uninformative about tail boundedness, e.g. Blankenberge
  above), others are wide-but-entirely-negative (imprecise but still
  confidently bounded, e.g. Raversijde1Buoy: CI [-0.929, -0.271] from
  only 67 peaks — different from "unreliable," worth distinguishing).
  Don't compare xi point estimates across buoys without checking their
  CIs. `summarize_results.py` now includes this cross-check
  automatically (sorted by |xi|, flags wide/zero-crossing CIs).
- **Spatial correlation decreases with distance** (Spearman r ≈ -0.47);
  4-cluster structure is geographically sensible, with **Zeebrugge as a
  singleton cluster** — independently corroborated by its extreme M2
  ratio, its distinct EVA peak count, and (once Stage 13 existed) its
  unstable moving-window distribution verdict. Three unrelated analyses
  agreeing makes this a solid finding, not a coincidence.
- **Real multi-year coverage is highly uneven** (see Data section) —
  worth checking before assuming any cross-buoy multi-year comparison is
  apples-to-apples.
- **Findings reproduce on real multi-year data, not just the 2-month
  NRT window**: full 19-buoy batch run against `data_multiyear`
  reproduced Zeebrugge's singleton spatial cluster, ~99% effective-N
  reduction from autocorrelation (362,817 raw → 2,783 effective, network
  mean), and correlation-decreasing-with-distance (Spearman r=-0.506) —
  all close to the NRT-window values, on a completely different dataset.
- **Westhinder shows possible genuine multi-decade drift**: two
  independent signals, not one — (1) 11b's persistence estimate has a
  wide per-segment spread (26-253h) that doesn't cleanly resolve to
  noise vs. real change; (2) Stage 13's calendar-windowed stability
  check shows all 4 eras individually agree on Weibull, but the pooled
  full record says lognormal (a mixture artifact from combining
  slowly-drifting sub-periods) — and window means climb monotonically
  +3.2% oldest to newest era. **Zeebrugge, by contrast, shows no such
  drift** (100% window agreement on lognormal, every era) — the
  contrast between two long-record buoys is itself informative: this
  looks like a real, site-specific phenomenon at Westhinder, not a
  pipeline artifact common to all long records. Not yet resolved — see
  `PLAN_next_session.md`.

### Zeebrugge — five independent lines of evidence it's structurally distinct

Zeebrugge's tidal notch never fully cleans regardless of approach tried
(2 vs. 3 harmonics; fitting the actual dominant frequency instead of
assuming M2). **The frequency hypothesis was tested directly and
rejected**: fitted period came back 12.4163h — essentially exact M2
(delta -0.0043h), nowhere near S2. Even with the *correct* frequency,
the notch made things measurably *worse*, ruling out "wrong fundamental
frequency" as the explanation. Two untested remaining hypotheses
(compound/shallow-water tide like MS4; time-varying tidal parameters
over the 14-year record) — not pursued further, this isn't a
parameter-tuning problem.

A joint M2+diurnal harmonic notch was also tried and **not adopted** —
partial, inconsistent real-data improvement (5-80% across 3 test
buoys, never below the "clean" threshold), doesn't justify the
downstream re-validation cost.

**Five independent analyses now agree Zeebrugge is dynamically
distinct**, not just noisier: worst tidal contamination; its own
singleton spatial cluster (reproduced on both NRT and multi-year data);
unstable distribution-fit verdict across sub-windows; and — new —
**wind-wave coupling collapses** (R²=0.093 vs. network median 0.565,
and the MSLP leading-indicator mechanism inverts sign entirely, lagging
instead of leading). Plausible mechanism: as the most sheltered,
harbor-interior site in the network, local wave conditions are likely
dominated by harbor geometry/vessel wake rather than direct open-water
wind forcing — a genuinely different causal pathway. Treated as a
structurally different site (harbor case study) rather than a
characterization target to keep forcing into the same framework as the
other 18 buoys.

**A sixth check was tried and came back ambiguous, not confirmatory —
worth recording honestly rather than omitting.** Statistical-fingerprint
clustering (distribution/tail/persistence/regime similarity, a
different lens than the correlation-based spatial clustering above)
initially showed Zeebrugge as its own cluster too — but this dissolved
once record-length/missingness metadata (14.3yr, 20.8% missing — both
real outliers among the 19 buoys) was excluded from the feature set.
Zeebrugge's fingerprint-based distinctiveness looks like it was
substantially driven by those metadata differences, not by its core
behavioral properties — this specific piece of evidence is downgraded
accordingly. The other five lines above are unaffected (none depend on
this feature set), so the overall case for Zeebrugge being distinct
doesn't weaken, but it rests on five confirmed lines, not six.

### Wind-wave coupling (ERA5, Stages 16/20/21)

Strong, consistent signature across 18 of 19 buoys: wind leads Hs by
0-3h, R² up to 0.71-0.77 (Westhinder), MSLP consistently leads with the
physically correct pressure-drop-precedes-Hs-rise sign. A genuine ~24h
diurnal cycle exists in Hs at every buoy tested (85-289x baseline
power) — confirmed to be substantially wind-driven (land-sea breeze),
not primarily an unmodeled tidal constituent, by checking ERA5 wind's
own diurnal signature independently.

Wind adds real forecasting value beyond what Hs's own autoregressive
structure already captures: ARMAX beats plain ARMA at every horizon
tested on both Westhinder and A2Buoy, and the benefit does not fade
with horizon — traced to the fact that Hs's own short-term features
become progressively less informative as forecast horizon extends,
while wind (reflecting broader synoptic state) retains more relative
relevance. Same pattern holds independently for exceedance
classification (Stage 21).

## Forecasting (Stages 17-21)

Built after the characterization pipeline's uncertainty machinery
(Stage 12/13) was trusted — deliberately, not started earlier, since
validating a forecast against point estimates of unknown uncertainty
would have been premature. Six stages, each validated on synthetic data
with a known answer before trusting it on real data (see
`CHANGELOG.md` for the specific bugs each one caught along the way —
several were real, not superficial):

- **Persistence baseline + shared backtest harness** (`forecast_utils.py`,
  `17_forecast_baseline.py`) — every later model reuses this harness.
- **ARMA** (`18_forecast_arma.py`) — beats persistence at every horizon
  tested on Westhinder (skill 0.10-0.30, growing with horizon).
- **ARMA-GARCH** (`19_forecast_arma_garch.py`) — adds a calibrated
  prediction interval (95% coverage, verified on real data), not just a
  point forecast. Westhinder's volatility is near-IGARCH
  (alpha+beta≈1) — variance shocks barely decay within any
  practically forecastable horizon.
- **ARMAX** (`20_forecast_armax.py`) — adds lagged ERA5 wind as an
  exogenous regressor. Handles the honest-forecasting constraint
  explicitly: real future wind isn't known any more than real future
  Hs is, so exog is only genuinely known within the wind-Hs lag window
  and persisted (not fabricated) beyond it.
- **Exceedance forecasting** (`21_forecast_exceedance.py`) — reframes
  as "will Hs exceed threshold X in the next N hours" (probabilistic
  classification, Brier score + ROC-AUC), arguably more operationally
  relevant than point forecasts. Optional wind feature improves skill
  consistently across both buoys tested, and the benefit does not fade
  with horizon.

All of this is standalone (buoy-by-buoy, run manually), not wired into
`run_all_buoys.py` — same pattern as Stages 14-16.

## Where this stops

Point/probabilistic forecasting is now built (see above), but this is
still fundamentally a **characterization + forecasting toolkit**, not
an operational system — no automated retraining, no real-time data
feed, no alerting layer. The parked industrial-monitoring ideas
(CUSUM/EWMA) from an earlier external review would sit on top of Stage
21's exceedance probabilities if built. See `PLAN_next_session.md` for
the current priority queue and every open question that's still
genuinely unresolved (Zeebrugge's tidal-frequency mystery, the
diurnal-signal investigation, per-buoy wind-lag calibration beyond
Westhinder/A2Buoy).
