# Pipeline Expansion — Plan for Next Session

Status as of today: 19-buoy characterization pipeline (stages 01-08, 10,
11, 11b) runs clean end-to-end. Stage 09 removed from the default batch
loop (only CadzandBoei/Deurlo support it — see README/batch log for why).
Stage 11b (dependence structure / integral timescale) added ahead of the
originally-planned uncertainty stages, per the point that naive CIs on
autocorrelated data would be invalid without it — see Priority 1.

## Priority 1 — Dependence structure (BUILT — do this reading before Stage 12/13)

**Status: built as `11b_dependence_structure.py`, wired into `run_all_buoys.py`
right after Stage 0/02.**

This moved ahead of confidence intervals for a concrete reason, not just
tidiness: bootstrap CIs and Fisher z CIs on autocorrelated data (which is
every series in this pipeline — Ljung-Box fails at all 19 buoys) need an
autocorrelation-aware correction, and both corrections need the same
number — the integral (decorrelation) timescale — as an input. Building
CIs first would have meant guessing that number twice, in two different
stages, with no shared justification.

What it computes, per buoy: lag-1 ACF, the lag where ACF first crosses
zero, and the integral timescale (`tau = dt * (1 + 2*sum(rho))` up to that
crossing). From that it prints (not yet auto-applies):
- a suggested block-bootstrap block length, for Priority 2 below
- a suggested EVA declustering window (2x tau) — **compare this against
  Stage 08's 24-48h default before trusting existing xi estimates**; if a
  buoy's suggested window exceeds 48h, some of that buoy's "independent"
  storm peaks may be the same storm counted twice, which biases xi
  (usually toward appearing more bounded than reality)

**Action before building Stage 12**: run `11b` across all 19 real buoys,
check how far its suggested decluster windows land from Stage 08's
current 24-48h default. If they diverge a lot, re-run
`rerun_eva_all_buoys.py` with `--min-separation-hours` set from 11b's
per-buoy suggestion before trusting the current EVA table further.

Note: this is a simplified version of the idea behind Politis-White
optimal block length (first-zero-crossing summation), not the full
algorithm — stated in the script's own docstring so it isn't over-trusted
as more rigorous than it is.

## Priority 2 — Uncertainty & stability, now with the dependence fix baked in

**Stage 12 — Statistical confidence** (revised from the original plan —
naive/IID methods are now explicitly out)
- Hs mean/quantile CI: **block bootstrap** (moving-block or circular
  block), block length from Stage 11b's per-buoy suggestion — NOT a
  plain IID bootstrap, which would understate uncertainty given the
  confirmed autocorrelation
- GPD xi CI: bootstrap on the *declustered* exceedances — closer to
  valid than the Hs bootstrap since declustering already targets
  quasi-independence, but only as good as the declustering window itself
  (see Priority 1's action item — fix the window first)
- Fisher z CI on Stage 11 pairwise correlations: apply the Bartlett-type
  effective-N correction (N_eff ≈ N / (1 + 2·Σ ρ_x(k)·ρ_y(k))) using
  Stage 11b's per-buoy ACF instead of the raw sample size — otherwise
  the "correlation decreases with distance" CIs will look artificially
  tight
- CI bands on the fitted Weibull/lognormal PDF from Stage 06 — and add a
  one-line caveat wherever Stage 06's KS p-values are reported: they were
  never valid in the classical sense, independent of the large-n floor
  already noted, because the KS test itself assumes iid samples and the
  raw level series is autocorrelated. Point estimates (shape/scale) are
  probably fine; the p-values aren't a real p-value in the textbook sense
  either way.

**Stage 13 — Stability analysis**
- Moving-window recomputation of key stats (mean Hs, best-fit
  distribution) to check if conclusions hold across sub-windows
- Jackknife or "drop the biggest storm" sensitivity check on the Weibull
  fit and the EVA shape parameter specifically
- Bootstrap resampling of the regime fractions from Stage 10

Both stages read existing Stage 0/06/08/11/11b output — no new
downloads needed.

## Priority 3 — Pipeline architecture: tiering gate

Formalize what's currently ad hoc (Stage 10's silent Hs-only fallback,
Stage 09's manual exclusion from the batch loop) into an explicit,
declared gate in `run_all_buoys.py`:

- **Core** (always runs): 01, 02, 03, 03b, 04, 05, 06, 08
- **Advanced** (runs if the buoy's data supports it): 09 (needs
  VTPK/VMDR), 10 `--include-period` (needs VTPK), 07/10/11 already fit
  here without extra gating
- **Multi-year** (only activates once record length crosses a threshold):
  Mann-Kendall trend test, seasonal STL, robust EVA (see Priority 5)

Concretely: each stage script should declare what it needs (e.g. "requires
VTPK") and the orchestrator checks that against what Stage 0 found for
that buoy, instead of the stage crashing or silently degrading.

## Priority 4 — Decide and act on: multi-year re-download

Discussed today but not yet decided. If yes:
- Re-run `download_belgian_wave_buoys.py` with a wider time window (check
  how far back INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036 actually goes,
  or whether the MY delayed-mode product is needed instead of NRT for a
  multi-year historical pull)
- This unlocks: robust EVA (30-50+ storm peaks instead of 5-30), Mann-Kendall
  trend test, real seasonal decomposition, regime dwell-time distributions
- Decision needed: how many years, and NRT vs MY product

## Priority 5 — Forecasting pipeline (separate track, only after Priority 1-3)

Skeleton agreed on previously:
1. Persistence baseline (mandatory reference point)
2. Feature engineering — lagged Hs/Tp, cross-buoy leading indicators
   (Stage 11's spatial correlation tells you which upstream buoy to use)
3. ARMA on Stage 04's residual → ARMA-GARCH (justified by Stage 07's
   universal ARCH finding) → ARMAX once wind data exists
4. Rolling-origin backtest, skill score vs. persistence, error-vs-horizon
   curve
5. Consider exceedance forecasting (binary/probabilistic) as the
   operationally relevant framing given the DEME context

Do not start this until Priority 2 (uncertainty) exists — validating a forecast against
point estimates whose own uncertainty is unknown is premature.

## Priority 6 — Meteo data source

Check co-located sensors first (free, no new pipeline):
```bash
ncdump -h data/WesthinderBuoy.nc | grep -iE "WSPD|WDIR|ATMS|PRES"
```
If absent: ERA5 via Copernicus **Climate Data Store** (`cds.climate.copernicus.eu`,
`cdsapi` client) — a third, separate Copernicus portal from CDSE and CMEMS,
own account. Only pursue once Priority 5 actually needs an exogenous
regressor (ARMAX step).

## Lower priority / parked

- **Stage 09 fate**: keep as a standalone script for CadzandBoei/Deurlo
  only. Not worth gating into the main loop given 17/19 buoys can't use it.
- **Buoy statistical fingerprint + clustering (from the GPT review)**: PCA
  or classical MDS, not UMAP (n=19 is too small for UMAP to mean anything).
  Build after Priority 2, since a fingerprint without uncertainty bounds
  on its components is exactly the kind of point-estimate-only artifact
  Priority 1 is meant to fix.
- **Entropy/Hurst measures**: still deprioritized, low marginal value over
  existing ACF-based diagnostics.
- **Seasonality beyond M2**: blocked on Priority 4 (needs multi-year data).

## Suggested order for tomorrow specifically

1. Run `11b_dependence_structure.py` across all 19 real buoys (via
   `run_all_buoys.py` — it's already wired in). Check its suggested
   decluster windows against Stage 08's current 24-48h default before
   anything else — this gates whether the existing EVA table needs a
   third re-run.
2. Stage 12 (confidence intervals) using 11b's numbers — start with the
   EVA xi bootstrap, since that's been the least stable number across
   every run so far, and now has an actual declustering-window
   justification behind it instead of a round default
3. Stage 13 (stability/jackknife) if time allows
4. Decide Priority 4 (multi-year — yes/no, how far back) so it can run
   in the background while Stage 12/13 get built
