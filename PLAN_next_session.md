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

**Fixed today, blocking issue found mid-batch-run**: Stage 12's Weibull
bootstrap refit loop used full `scipy.stats.weibull_min.fit()` (MLE,
numerical optimizer over the whole resampled array) 1000 times per
buoy. At Westhinder's scale (~2s/fit x 1000) this was ~30-40 minutes;
A2Buoy (1.08M samples, ~2x Westhinder) would have been proportionally
worse - the batch run got stuck silently on this (no visible output
during a stage, since `run_all_buoys.py` buffers via
`capture_output=True`) and had to be killed. Fixed two ways: (1) added
`fast_weibull_moment_fit()` - closed-form method-of-moments via the
coefficient-of-variation relation, used inside the bootstrap loop only
(Stage 06's actual point estimate still uses full MLE) - validated at
0.02% relative error vs. full MLE, 630x faster; (2) `--n-bootstrap` now
auto-scales down for large records (200 at >=500k samples, 500 at
>=50k, 1000 otherwise) since bootstrap CI precision barely improves
past a few hundred resamples once the underlying n is already huge.
Combined: a 600k-sample synthetic test that would have taken ~33+
minutes now runs in ~11 seconds. Batch run needs to be restarted from
scratch with the fixed script.

**Fixed today, found from the first full multi-year batch run's own
results**: `run_all_buoys.py` never connected Stage 08 to Stage 11b's
persistence estimate the way Westhinder's manual spot-check did -
every one of the 19 buoys' EVA ran on the generic 48h default,
regardless of what 11b had already computed for that specific buoy
(e.g. A2Buoy: 184.8h suggested, 48h used). Concrete symptom: A2Buoy's
GPD xi CI at the wrong window was [-0.198, 0.065] - crosses zero,
uninformative, same failure signature as the too-few-peaks cases but
caused by too-short-window this time. Added `build_stage08_args()`,
same dynamic-dispatch pattern as Stage 10's `--include-period` -
verified via the actual logged command lines that each buoy now gets
its own persistence-justified `--min-separation-hours` (confirmed
different per buoy in a 3-buoy smoke test: 68.9h/66.4h/60.4h). **The
completed 19-buoy multi-year batch predates this fix - its Stage
08/12 xi results should be treated as provisional until re-run.**

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

1. **DONE, tested on real Westhinder data - Stage 08's EVA declustering
   window connected to Stage 11b's persistence estimate.** Ran
   `08_extreme_value_analysis.py --min-separation-hours 231.4` (11b's
   suggested window, 2x its 115.70h weighted-mean persistence estimate)
   against the default 48h baseline. Result: n_peaks 1031->483 (-53%,
   less than the ~5x naive scaling would predict - most of the original
   peaks were already reasonably independent), xi -0.2315->-0.3629 and
   sigma 0.8929->1.2729 (both shifted meaningfully, ~57% on xi), but
   **every return level moved by under 0.15m even at the 50-year
   horizon** - shape/scale changes largely cancel in the return-level
   formula, so the practically-relevant number is robust to this
   methodological choice even though the individual GPD parameters
   aren't. Independent plausibility check: 483 peaks/35.95yr = 13.4
   storms/year, physically sensible for North Sea autumn/winter storm
   seasonality; the old 48h window's 28.7/year implied an implausibly
   high independent-storm rate for this coast. KS p-value dropped
   (0.91->0.28, expected from fewer points, still non-rejecting) - not
   a concern given the standing caveat on this test's validity anyway.
   **Net: the connection was worth making** - solves the "same problem
   twice" issue the plan flagged, and the result is now both internally
   validated (return-level robustness) and physically validated
   (storm-rate plausibility) rather than just theoretically motivated.
   Do this same `--min-separation-hours` connection for the other 18
   buoys once their own multi-year 11b runs exist, rather than reusing
   Westhinder's 231.4h number - each buoy's persistence timescale is
   its own.
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
4. **Stage 11's per-pair overlap isn't reported (not a bug, a
   transparency gap).** `df.corr(method="pearson", min_periods=20)`
   already computes each buoy pair's correlation over their true
   overlapping window only (pandas' default pairwise-complete-
   observation behavior) - it does NOT compare full-record-each, so
   there's no correctness issue here. But with record lengths now
   spanning 5 years (Raversijde1/2) to 36 years (Westhinder), a
   correlation between a short- and long-record buoy gets computed
   correctly but silently over just the short buoy's window - nothing
   in the output currently shows that. Add an `n_overlap` column to
   `_correlation_vs_distance.csv` so a 5-year-overlap pair isn't read
   with the same confidence as a 30-year one.

## Priority 3 — Multi-year-only analyses (now genuinely unblocked)

These were explicitly parked pending real multi-year data - that data
exists now. Only meaningful for the **6 buoys with 30+ year records**
(Westhinder, TrapegeerBuoy, BolVanHeistBuoy, ScheurWielingenBuoy,
WandelaarBuoy, OstendEasternPalisadeBuoy) - the other 13 don't have
enough history yet for a credible version of either:

- **Mann-Kendall trend test** - is there a real long-term drift in wave
  climate at these 6 sites, separate from storm-scale and tidal
  dynamics. Literature convention wants 10-20+ years; these buoys clear
  that now. **Use a modified test, not the textbook one**: Ljung-Box
  confirms strong serial correlation out to at least 24-48h lags, which
  is exactly the condition known to inflate plain Mann-Kendall's
  false-positive rate on trend detection. Use Hamed-Rao variance
  correction (effective-sample-size adjustment, same spirit as the
  Fisher z correction already applied in Stage 12b) or Yue-Pilon
  prewhitening (remove AR(1) structure before testing). Running the
  textbook version here would be a real regression in rigor relative to
  everything else in this pipeline.
- **Seasonal decomposition (STL)** - separate real winter-storm-season
  vs. summer-calm seasonality from the M2 tide, which is the only
  "seasonality" the pipeline has been able to detect so far. Needs at
  least 1-2 full annual cycles; these 6 buoys have 20-30+.

New stage(s) needed - not yet built. Suggest a new script (e.g.
`14_multiyear_climatology.py`) gated via the tiering gate's
`min_record_years` requirement type (already built in Priority 3 from
today, currently unused - this is exactly what it's for).

## Flagged for later — Westhinder persistence drift question

11b's v3 (segment-aggregated) run on Westhinder shows a very wide
per-segment tau spread: 26.4h to 252.9h across 50 segments (weighted
mean 115.70h). Two competing explanations, not yet distinguished:
- **Noise**: shorter qualifying segments (down to 200 samples = 100h)
  produce much noisier individual ACF/tau estimates than the 36,170-
  sample longest segment.
- **Genuine drift**: a 36-year deployment could plausibly have real
  persistence changes over time - mooring relocations, sensor
  upgrades, different storm-climate eras.

Diagnostic to run when this gets picked back up: plot per-segment tau
against (a) segment length (tests the noise hypothesis) and (b) segment
midpoint date (tests the drift hypothesis). If (b) shows a real trend,
that's a genuinely interesting finding in its own right, not just
something to average away - connects directly to Priority 3's
Mann-Kendall work (persistence itself trending would be a different,
possibly more fundamental, finding than a simple Hs-level trend).

Also worth revisiting once this is checked: `--max-segments=50` (out of
345 found) currently leaves ~39% of valid Westhinder data unused. Try
raising the cap and see whether the aggregate estimate stabilizes or
keeps shifting - stabilizing would suggest 50 was already enough;
shifting would suggest it wasn't.

**Update - second independent signal now points the same direction.**
Stage 13's `[A]` moving-window check was fixed (see below) and rerun on
Westhinder: all 4 calendar-era windows individually agree on Weibull,
but the pooled full-record fit (Stage 06) says lognormal - a classic
mixture artifact of combining several Weibull sub-periods with slowly
drifting parameters. Window means also climb monotonically: 1.033 ->
1.048 -> 1.059 -> 1.066 m (+3.2% oldest to newest era). Two separate
analyses (persistence-timescale spread, distribution-shape pooling
artifact) now both suggest possible genuine long-term drift at this
buoy, not independent curiosities. This is a concrete lead for
Priority 3's Mann-Kendall work, not just a generic checklist item -
bumping that priority up given this.

**Fixed today**: Stage 13's `[A]` was splitting the post-dropna()
COLLAPSED array positionally - on Westhinder's 345+-segment fragmented
record, a "window" could be a patchwork of disjoint calendar periods,
and printed window start/end dates were actively wrong. Now splits by
real calendar time; window sample counts correctly come out unequal
(95314/130267/149287/152752) reflecting real differences in gap
coverage across eras (plausible: older decades = more sensor downtime).

**Fixed**: `10_regime_identification.py`'s regime labels have gap
positions entirely ABSENT (missing rows), not NaN-marked on a regular
grid - `utils.all_contiguous_segments()` doesn't apply directly (it
needs NaN gridding). Added `utils.segments_by_time_gap()` instead -
splits a series by actual elapsed time between consecutive index
entries, for exactly this "compacted, gap-rows-removed" case. Stage
13's `[C]` now bootstraps within each detected segment separately
(never draws a block across a real gap seam), concatenating per-segment
resamples to build each bootstrap iteration's full-length surrogate.
Tested on a synthetic fragmented regime sequence (15 injected gaps ->
9 detected segments after adjacent merges, 100% of labels used,
sensible ~25% CIs on uniform-random test labels) before running on
real data. **Result on Westhinder**: 1905 segments detected by
`segments_by_time_gap` on the regime labels - exactly matching the
1905 segments Stage 02/03/03b/11b independently found via NaN-based
detection on the same underlying Hs series (strong cross-validation
that both gap-detection methods agree). 668 segments (96.3% of labels)
used; point fractions unchanged (correct - they never depended on the
bootstrap method), CI widths broadly similar but regime 0's CI became
asymmetric ([39.2%, 41.2%] vs the old [38.2%, 40.5%]).

**Known limitation, not fixed (arguably not a bug)**: many of the 668
used segments are shorter than the 232-sample block length, which makes
`moving_block_bootstrap_resample` degenerate to a single deterministic
block (no real position to vary) - those segments contribute identical
values to every bootstrap iteration rather than genuine resampled
variability, which is what's producing the asymmetric/narrowed CI.
Defensible reading: a segment shorter than the persistence timescale
doesn't have much internal decorrelation to resample in the first
place, so the degeneracy may honestly reflect low information content
rather than being wrong. Left as a known caveat rather than fixed given
it's genuinely ambiguous whether "more resampling" would even be more
correct here.

**Confirmed consequential, not just theoretical**: this exact bias
crashed `13_stability_analysis.py` outright on CadzandBoei and Deurlo
in the first full multi-year batch run - `ValueError: 'yerr' must not
contain negative values`, because the bootstrap CI's lower bound ended
up ABOVE the point estimate for at least one regime on each buoy (both
are the finest-sampling, richest-variable-set buoys in the network -
plausibly a different segment-length distribution relative to block
length than the other 17). **Fixed** (plotting only, not the underlying
bootstrap bias itself): reproduced the crash directly with a crafted
CI where `ci_low > point_fraction`, confirmed the fix handles it -
proper asymmetric error bars from both `ci_low`/`ci_high` (the old code
only used `ci_low`, applied symmetrically, silently ignoring `ci_high`
entirely), with negative error-bar lengths clipped to 0 and an explicit
console warning naming which regime(s) it happened to, rather than
crashing silently mid-batch. The underlying bootstrap bias itself is
still the same open question as above - this only stops it from
crashing the pipeline.

## Priority 4 — Revisit the Zeebrugge tidal-frequency issue

With 17 years of real Zeebrugge data now available (vs. 2 months
before), a proper frequency-fitting approach - estimating the actual
dominant tidal frequency near M2 rather than assuming it's exactly
12.4206h - has a much better chance of separating true constituents
(M2 vs. a nearby S2 or shoaling-shifted frequency) than it did on a
2-month record where two similar-period constituents can't be
resolved apart. Worth attempting once Priority 2's multi-year pipeline
run is validated and stable.

**Confirm before spending time on the refit, don't assume**: frequency
resolution scales as 1/T, so going from 2 months to 17 years shrinks it
by roughly two orders of magnitude - M2 (12.4206h) and S2 (12.0000h)
should be comfortably resolvable at that record length. Worth a quick
check (compute the actual achievable frequency resolution at 17 years
and confirm it's well under the M2/S2 spacing) before investing time in
the harmonic refit itself - this is the kind of assumption the rest of
this pipeline has been careful about verifying rather than taking on
faith.

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
