# Pipeline Expansion — Plan for Next Session

## Status as of today (2026-07-12, end of day)

Priorities 1-3 from the previous version of this plan are **done**:
dependence structure (11b, fixed), uncertainty/stability (12, 12b, 13),
and the tiering gate. Full detail on all of that: `CHANGELOG.md`.

**New today, not yet acted on:**
- Multi-year wave history downloaded successfully: `data_multiyear/`,
  all 19 buoys, 81.39 MB total. Real per-buoy coverage is highly uneven
  - 6 buoys back to 1990-1997 (Westhinder longest, 1990-07-19), the
  other 13 from 2009-2021. All extend to 2026-06-30.
- ERA5 meteo download launched (2010-2026, 3-hourly, monthly-chunked) -
  **status unknown as of writing this**, check it finished before
  relying on it.
- The characterization pipeline has **not yet been run against
  `data_multiyear/`** - everything built and validated so far, including
  all of Stage 12/13's uncertainty work, was tested against the ~2-month
  NRT window only.

## Priority 1 — Verify the ERA5 download completed

204 monthly requests each hit CDS's queue independently. Check for
failed months before assuming the range is complete:
```bash
ls meteo_era5/ | wc -l   # should be 17*12 = 204
```
If any are missing, `download_era5_meteo.py --start-year Y --end-year Y`
re-run for just the affected year(s) - it skips files that already
exist, so this is safe to re-run broadly if unsure.

## Priority 2 — Point the pipeline at real multi-year data, carefully

**Do not just run `run_all_buoys.py --data-dir data_multiyear` blindly.**
Several defaults were tuned against the ~2850-sample, 59-day NRT window
and may not behave sensibly on decade-plus records. Spot-check one
long-record buoy (Westhinder - 1990-2026, the longest) stage-by-stage
first, the same discipline used all day today before trusting a batch
run:

```bash
python 01_load_clean.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
python 02_eda_diagnostics.py --buoy WesthinderBuoy --var VHM0
python 03b_tidal_notch.py --buoy WesthinderBuoy --var VHM0 --harmonics 2
python 11b_dependence_structure.py --buoy WesthinderBuoy --var VHM0
python 08_extreme_value_analysis.py --buoy WesthinderBuoy --var VHM0
```

Specific things likely to need adjustment, in rough order of how likely
they are to actually break something (not just look different):

1. **Stage 08's EVA threshold/declustering defaults.** The 85th
   percentile / 24h fix solved "too few peaks in 2 months." A 36-year
   record has the opposite risk: at 85th percentile, Westhinder alone
   could have hundreds of "storm peaks," many of which may not be
   independent storms at 24h separation once the record is this long.
   May need the threshold raised back toward 95th (or higher) now that
   sample size is no longer the constraint - check n_peaks before
   trusting any multi-year xi estimate.
2. **11b's `--max-lag` default (500 samples).** Was already occasionally
   hitting its ceiling on the ~2850-sample record. On a 36-year record
   this is a tiny window - will almost certainly need scaling relative
   to record length, not a fixed sample count, or the persistence
   estimate on long-record buoys will silently be a lower bound again
   (the exact bug that was fixed in 11b's raw-vs-detided issue, just a
   different cause this time).
3. **Stage 13's fixed 4-window split.** Sensible when each window was
   ~2 weeks. On 30 years, 4 windows means ~7.5 years each - may hide
   real short-term instability. Consider `--n-windows` scaled to record
   length, or a fixed window duration instead of a fixed window count.

## Priority 3 — Multi-year-only analyses (now genuinely unblocked)

These were explicitly parked pending real multi-year data - that data
exists now. Only meaningful for the **6 buoys with 30+ year records**
(Westhinder, TrapegeerBuoy, BolVanHeistBuoy, ScheurWielingenBuoy,
WandelaarBuoy, OstendEasternPalisadeBuoy) - the other 13 don't have
enough history yet for a credible version of either:

- **Mann-Kendall trend test** - is there a real long-term drift in wave
  climate at these 6 sites, separate from storm-scale and tidal
  dynamics. Literature convention wants 10-20+ years; these buoys clear
  that now.
- **Seasonal decomposition (STL)** - separate real winter-storm-season
  vs. summer-calm seasonality from the M2 tide, which is the only
  "seasonality" the pipeline has been able to detect so far. Needs at
  least 1-2 full annual cycles; these 6 buoys have 20-30+.

New stage(s) needed - not yet built. Suggest a new script (e.g.
`14_multiyear_climatology.py`) gated via the tiering gate's
`min_record_years` requirement type (already built in Priority 3 from
today, currently unused - this is exactly what it's for).

## Priority 4 — Revisit the Zeebrugge tidal-frequency issue

With 17 years of real Zeebrugge data now available (vs. 2 months
before), a proper frequency-fitting approach - estimating the actual
dominant tidal frequency near M2 rather than assuming it's exactly
12.4206h - has a much better chance of separating true constituents
(M2 vs. a nearby S2 or shoaling-shifted frequency) than it did on a
2-month record where two similar-period constituents can't be
resolved apart. Worth attempting once Priority 2's multi-year pipeline
run is validated and stable.

## Priority 5 — Pair ERA5 with the buoys

First real cross-variable work with actual wind data, once both
downloads are confirmed complete (Priority 1) and the multi-year
pipeline run is trusted (Priority 2). Extends Stage 09's cross-variable
framework - genuinely useful across most of the network now, not just
CadzandBoei/Deurlo, since wind is external data rather than something
only 2 buoys happen to measure.

Natural first questions: does wind speed lead Hs by a predictable lag
(cross-correlation, same method as Stage 09's Hs/Tp CCF)? Does wind
direction align with wave direction (VMDR) during storms at
Cadzand/Deurlo specifically, or do they decouple (swell vs. wind-sea)?

## Priority 6 — Forecasting pipeline (still parked)

Reasoning from before still holds: don't build on top of the
characterization pipeline until the multi-year version of it (Priority
2) has been run and trusted. Skeleton unchanged from prior planning -
persistence baseline, then ARMA/ARMA-GARCH on Stage 04's residual, then
ARMAX once ERA5 wind is validated and paired (Priority 5), then
rolling-origin backtest against a skill score, not raw RMSE.

## Priority 7 — Provenance metadata (cheap, do whenever there's a gap)

Not urgent relative to Priorities 1-3, but cheap enough to slot in
anytime, including mid-way through tomorrow if there's a natural pause
(e.g. while a batch run or CDS download is in progress). From today's
external review — the single most directly actionable suggestion in it.

Attach to every stage's output JSON: git commit hash, Python version,
key package versions (statsmodels/scipy/pandas), input file hash or
mtime, random seed used, execution timestamp. Given today included
three separate wrong-API-guess corrections against `copernicusmarine`,
having exact package-version provenance on every run would have made
pinning down "which version's signature am I actually looking at"
faster. Implementation: a small `utils.get_provenance()` helper, called
once per script, merged into whatever summary JSON that stage already
writes — not a new file per stage.

## Priority 8 — Physical interpretation stage

Cheap and low-risk: templated text conditioned on parameters already
computed (not new inference), translating e.g. `xi < 0` into "bounded
tail, consistent with fetch-limited shelf sea" the way this plan and the
README already do in prose. Genuinely useful for anyone reading Stage 08
output who isn't fitting GPDs for a living. Could live as an extra
printed block in `08_extreme_value_analysis.py` and `06_distribution_fit.py`
rather than a new stage.

## Priority 9 — Stage X: four-question diagnostics report

From the external review, distinct from the fingerprint idea (Priority
10) — a narrative summary, not a numeric table. Four questions per
buoy: can I trust the data (QC/gaps), can I trust the assumptions
(stationarity/dependence/fit), can I trust the estimates (CI/stability),
what did I learn about the physics (storms/persistence/regimes/spatial).
Complements the tiering/CI machinery rather than duplicating it. Bigger
lift than 7/8 — do after Priority 2/3 (multi-year validation) are done,
since it should report on real multi-year results, not the 2-month
window.

## Priority 10 — Statistical fingerprint (now unblocked)

Previously parked pending uncertainty quantification — that prerequisite
is satisfied now that Stage 12/13 exist. One-page-per-buoy: distribution,
storm stats, persistence, CI, stability, all in one place. Still behind
7-9 in priority. Clustering across buoy fingerprints: PCA or classical
MDS, not UMAP (n=19 too small for UMAP to mean anything - unchanged
from the earlier assessment of this idea).

## Priority 11 — Assumptions-per-test summary

Lower value than it initially looks: significant overlap with
information already scattered across existing summary JSONs (Stage 03's
ADF/KPSS agreement flag, Stage 06's KS caveat, Stage 08's peak-count
reliability flag). Only worth building if it actually *consolidates*
those into one place rather than adding a fourth place the same
information lives — otherwise it's restating existing output with extra
steps.

## Explicitly rejected from the external review — do not build

- **"Uncertainty budget" decomposition** (sampling / autocorrelation /
  measurement / model-assumption / parameter-estimation as separate
  quantified slices). Standard practice in metrology *when each source
  can be independently characterized* - e.g. a sensor's calibration
  spec gives a real measurement-uncertainty number. We have no sensor
  accuracy spec or calibration history for these buoys. Decomposing
  "total uncertainty" into a measurement-uncertainty slice would mean
  inventing a number with no basis, dressed in engineering language -
  exactly the failure mode the KS caveat, the EVA threshold fix, and
  the Fisher z correction all exist to prevent. Do not build this.
- **Stage renumbering** (11 -> 11.1/11.2/11.3, etc.). Purely cosmetic,
  and actively risky on a live, actively-iterated pipeline - would
  require updating `run_all_buoys.py`, the README, the CHANGELOG, and
  breaking two days of accumulated command muscle-memory, for zero
  functional gain. Only worth revisiting if/when this becomes an actual
  versioned public release, not before.
- **"~80% of the framework would transfer to other environmental time
  series"** - not a task, just an unsupported speculative claim (no
  actual cross-domain testing behind the number). Worth remembering as
  a narrative point about the architecture's generality, not something
  to act on or repeat as if it were measured.

## Lower priority / parked

- **Stage 09 fate**: still standalone for CadzandBoei/Deurlo only, via
  the tiering gate now rather than manual exclusion.
- **Buoy statistical fingerprint + clustering**: PCA/MDS, not UMAP
  (n=19 too small). Build once Priority 2/3 give something worth
  fingerprinting beyond the current 2-month snapshot.
- **Entropy/Hurst measures**: still deprioritized.
- **`download_belgian_wave_buoys_multiyear.py`**: superseded by
  `download_belgian_wave_buoys_history.py`. Kept in the repo for the
  CHANGELOG's record of what was tried and why it didn't work
  (subset() vs. get() - see CHANGELOG), not meant to be run.

## Suggested order for tomorrow specifically

1. Priority 1 - confirm ERA5 finished (2 minutes, do this first)
2. Priority 2 - Westhinder-only spot check against `data_multiyear`,
   fix whatever breaks (EVA threshold almost certainly; 11b max-lag
   probably; Stage 13 windowing maybe) before touching the other 18 buoys
3. Once Westhinder looks right, full `data_multiyear` batch run
4. Priority 3 (Mann-Kendall/STL) on the 6 long-record buoys, if time
   allows - otherwise this is the natural start for the session after
