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

## Priority 3 — Multi-year-only analyses — CLOSED

Both pieces done, run on all 6 qualifying buoys (Westhinder,
TrapegeerBuoy, BolVanHeistBuoy, ScheurWielingenBuoy, WandelaarBuoy,
OstendEasternPalisadeBuoy - the only ones with 30+ year records; the
other 13 don't have enough history for either test to be credible).

- **Mann-Kendall trend test** - DONE, run on all 6 qualifying buoys.
  Result: **5/6 show no significant trend** in mean or p95 Hs.
  OstendEasternPalisadeBuoy shows a marginal p95 (storm intensity)
  increase (p=0.046, +0.005 m/year over 30 years) - but this doesn't
  survive multiple-testing correction (6 buoys x 2 metrics = 12 tests;
  Bonferroni threshold ~0.0042, this doesn't clear it). Treat as a lead
  to revisit once ERA5 wind is paired, not an established finding.
  Hamed-Rao correction factors mostly ~1.00-1.04 (minimal correction
  needed) except Trapegeer p95 (1.59) and OstendEasternPalisade mean
  (2.28) - consistent with the synthetic no-trend control test, which
  showed real annual aggregates naturally decorrelate well.
  **Resolves the Westhinder persistence-drift question from earlier**:
  Mann-Kendall found NO significant trend at Westhinder (p=0.456 mean,
  p=0.917 p95), despite Stage 13's `[A]` showing a monotonic +3.2% climb
  across 4 eras. Read: Stage 13's signal was very likely a coverage-
  imbalance artifact of its own windowing (window sample counts were
  wildly uneven - 95k/130k/149k/153k - so gap-heavy seasons landing
  unevenly across windows could produce an apparent trend with no
  significance test behind it), not genuine multi-decade drift. Annual
  Mann-Kendall (n=37, proper variance correction) is a more rigorous
  test than a 4-bin eyeball comparison and it found nothing. Downgraded
  from "open question" to "investigated, resolved."

- **Seasonal decomposition (STL)** - DONE, `15_seasonal_decomposition.py`
  built, validated on synthetic data (known-signal recovery, no-signal
  control, min-years gate), run on all 6 buoys. **Strong, consistent,
  physically real result** - passes both discriminators set up in
  advance:
  - Peak month is winter (Nov/Dec/Jan) at all 6 buoys, both mean and
    p95 - consistent across independent buoys, matching known North
    Sea storm seasonality.
  - Amplitudes (0.49-1.23m mean, 1.3-1.8m p95) are an order of
    magnitude above the synthetic no-seasonality control (0.04-0.06m) -
    real signal, not STL fitting noise to a forced period-12 term.
  - **New finding, not visible on shorter records**: storm-intensity
    (p95) seasonal amplitude is consistently ~2-3x the mean-Hs seasonal
    amplitude at every single buoy (e.g. Trapegeer 0.49m mean vs. 1.33m
    p95; Ostend 0.50m vs. 1.50m). Winter doesn't just raise the average
    sea state, it specifically inflates the extreme end - consistent
    with the already-established universal ARCH/volatility-clustering
    finding, now with a seasonal driver identified. This is genuinely
    new: M2 tidal notching was the only "seasonality" this pipeline
    could detect before the multi-year data existed.
  - Caveat, not a problem: variance-fraction attributable to seasonality
    sits at 0.37-0.57 (vs. 0.99 in the clean synthetic test) - a real
    chunk of month-to-month variance is unexplained by seasonality
    alone, consistent with everything else already known (storm
    persistence, ARCH, general non-stationarity) rather than
    seasonality being the single dominant driver.

**Combined read across both tests**: this long-record subnetwork shows
no robust long-term trend (Mann-Kendall), but a strong, consistent,
physically sensible winter-storm seasonal cycle (STL) that specifically
concentrates in storm intensity rather than typical conditions. Two
complementary characterizations of the same 6 buoys, not two unrelated
results.

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

**Resolved (partially) - see Priority 3.** Mann-Kendall on annual mean/
p95 Hs found NO significant trend at Westhinder. The Stage 13 `[A]`
mixture-artifact signal is very likely a coverage-imbalance artifact of
uneven window sample counts, not genuine Hs-level drift - see Priority
3 for the full reasoning. **Still open**: whether the persistence
TIMESCALE itself (not Hs level) trends over the 36 years - the
26-253h per-segment tau spread was never directly tested for a
trend-over-time the way Hs itself now has been. If picked back up,
that's a distinct question from what Mann-Kendall already answered -
would need Mann-Kendall (or a simple regression) applied to per-segment
tau vs. segment midpoint date specifically, not Hs level.

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

## Priority 4 — Revisit the Zeebrugge tidal-frequency issue — BUILT, needs the real run

## Diurnal notch extension — BUILT, tested on real data, NOT ADOPTED

**Real result on 3 buoys: partial, inconsistent, doesn't clear the
threshold.** Diurnal ratio before->after: Westhinder 65.23->62.26
(~5% reduction), BolVanHeist 236.71->46.65 (80%), A2Buoy 252.91->104.33
(59%) - none dropped below the established >3 "clean" threshold, far
short of the synthetic validation's >99.99% reduction. M2's own cleanup
was unaffected by the joint fit (confirms no interference), so the
shortfall is specific to the diurnal term.

**Decision: do not adopt.** A partial, buoy-inconsistent improvement
doesn't justify the downstream re-validation cascade this would trigger
(11b, 04, 12 all depend on the detided series). Likely cause: land-sea
breeze phase shifts seasonally (sunrise/sunset timing, storm
interference), so a single fixed-phase sinusoid across a multi-year,
heavily-fragmented record only captures whatever fraction of cycles
happen to stay in phase - a genuinely different approach (time-varying
envelope, seasonal deseasonalizing) would be needed to do meaningfully
better, which is a bigger undertaking than this investigation's scope.
Parked, not pursued further for now.

**Fixed a real bug found from this same result**: the "power still
elevated" warning only ever checked M2's ratio, so when diurnal was the
actual problem (as in this real run), it printed a misleading M2-focused
suggestion ("try --fit-frequency") that had nothing to do with the real
issue. Fixed: M2 and diurnal now get independently checked and
diagnosed, with the diurnal warning specifically explaining why more
harmonics won't help (confirmed empirically, not just theorized) rather
than reusing M2's generic advice.

**`--diurnal-harmonics` flag stays in the codebase** (useful
infrastructure, correctly validated on synthetic data, and the finding
that it doesn't work well on real multi-year records is itself a
documented, useful negative result) - just not used in the default
pipeline path or recommended for routine use.



**Confirmed the resolution check first, as planned**: at Zeebrugge's real
14.333-year record, achievable frequency resolution is ~355x finer than
the M2/S2 spacing needs - comfortably sufficient (this would NOT have
held on the original 2-month record).

**Built `--fit-frequency` in `03b_tidal_notch.py`**: instead of assuming
the dominant tidal period is exactly M2 (12.4206h), finds the actual
spectral peak in a search band (default 11.5-13.5h) on the longest
contiguous segment, and uses that fitted period for the harmonic
regression instead of the hardcoded constant.

**Validated on synthetic data, dramatic result**: injected a known
tidal signal at exactly S2 (12.0000h, not M2) into 15 years of
synthetic data. Fixed-M2 assumption completely failed - post-notch
power ratio actually got slightly WORSE (669162 -> 669434), since it
was fitting harmonics at the wrong frequency entirely and effectively
adding noise rather than removing signal. `--fit-frequency` recovered
the true injected period exactly (12.0000h, zero error) and cleaned the
signal ~59,000x (668932 -> 11.26). Control test (inject genuine M2,
confirm no spurious drift): fitted 12.4208h vs. true 12.4206h -
essentially exact, correctly stayed at M2 rather than wandering.

**Minor open observation, not blocking**: the M2 control's post-notch
ratio didn't clean up as dramatically as the S2 test's despite both
correctly identifying their true frequency (36184 vs 11.26) - likely a
periodogram bin-alignment/spectral-leakage artifact specific to that
frequency/record-length combination, not a flaw in the frequency
detection itself (both tests recovered their true period to within
0.0002h). Worth watching whether this recurs on the real run.

**Ran on real Zeebrugge data - hypothesis REJECTED, not confirmed.**
Fitted period: 12.4163h, delta from M2 = -0.0043h (essentially exact
M2); delta from S2 = +0.4163h (nowhere close). The dominant frequency
genuinely IS M2 at this site - the "maybe it's S2 or a shoaling-shifted
frequency" hypothesis standing since the very first Zeebrugge
investigation is now closed, rejected rather than confirmed.

**Surprising part: even with the correct frequency, the notch made
things WORSE** (power ratio 1142.37 -> 1591.55, same failure signature
as fitting the *wrong* frequency in the synthetic S2 validation test).
This rules out "wrong fundamental frequency" as the explanation
entirely - something else is broken.

**Caveat on the frequency estimate itself**: computed on the longest
single segment only (60,343 samples ~= 1.72 years - Zeebrugge
fragmented into 105 segments from 20.8% gap coverage, far more
fragmented than the "355x resolution" sanity check assumed, since that
number was calculated against the full 14.3-year record, not the
1.72-year effective baseline actually used). Resolution used (42.6x)
still clears the sufficiency bar, but on a shorter effective baseline
than the plan's check implied - the M2 conclusion is probably still
right, just resting on less data than intended.

**Priority 4 - CLOSED as "investigated, genuinely hard, not a
parameter-tuning problem."** Three separate attempts (harmonics=2,
harmonics=3, fit-frequency) - none improved things, the last one made
it worse. Two remaining hypotheses, neither pursued further for now:
1. **Compound/shallow-water tide** (e.g. MS4, a nonlinear M2/S2
   interaction product) - not an integer multiple of M2's frequency,
   so the current harmonic-multiples-of-one-fundamental basis
   structurally can't represent it regardless of which fundamental is
   used.
2. **Time-varying tidal parameters** - a single fixed-amplitude,
   fixed-phase sinusoid fit across 14 years can't capture nodal-cycle
   modulation or genuine physical change at a harbor site over that
   span (dredging, mooring relocation - plausible at exactly this kind
   of site).

**Action taken**: reverted Zeebrugge's working pipeline output back to
the harmonics=2, non-fitted version - imperfect, but doesn't actively
make things worse the way the fitted version did. Zeebrugge's
downstream numbers (11b persistence, EVA, etc.) stay as they were
before this investigation - no re-run needed.

## Priority 5 — Pair ERA5 with the buoys — BUILT, tested, ready to run

**`utils.load_era5_for_buoy()`** built - nearest-grid-cell extraction,
concatenates across every monthly ERA5 file, derives wind_speed and
wind_dir_from_deg (meteorological "FROM" convention, hand-verified
against all 4 cardinal directions before trusting it) and
air_sea_temp_diff_c. Tested against synthetic multi-file fixtures
matching the real download naming - correctly concatenates months in
order, correct row count.

**`16_wind_wave_coupling.py`** built - wind/Hs and MSLP/Hs
cross-correlation, Hs~wind_speed(²) regression, directional alignment
by regime (CadzandBoei/Deurlo only).

**Caught and fixed a real bug before it ran on real data**: the CCF
lag-sign convention was backwards in the first draft (docstring/prints
claimed negative lag = x leads y). Verified empirically with a
known-lag synthetic case rather than trusting the index algebra by
eye - confirmed the actual convention is **positive lag = x leads y**,
fixed the docstring and every downstream interpretation (wind-leads/
MSLP-leads messages, plot axis labels) to match. This is exactly the
kind of easy-to-get-backwards detail flagged as a risk when this stage
was outlined - worth remembering as a general lesson: verify lag/sign
conventions empirically, don't just read the code and assume it's
right.

**End-to-end validated** with a physically realistic synthetic case
(wind genuinely leading Hs by 6h, MSLP leading with the expected
pressure-drop-precedes-Hs-rise sign): recovered exactly `lag=6h, wind
leads`, r=0.95; MSLP correctly identified as leading with the right
sign relationship flagged. First test attempt actually caught a
separate, unrelated bug in the TEST setup itself (stale ERA5 fixture
files left over from testing the loader in isolation, not matching the
new test's intended wind signal) - worth remembering that a test
giving a nonsense result (R²=0.000 on a strong injected relationship)
is itself informative and worth chasing down, not just re-running.

**Run it**:
```bash
python 16_wind_wave_coupling.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
```
Repeat per buoy of interest. Directional alignment only produces output
for CadzandBoei/Deurlo (the only buoys with VMDR) - expected, not a bug,
for the other 17.

**Run it**:
```bash
python 16_wind_wave_coupling.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
```
Repeat per buoy, or use `run_all_wind_coupling.py` (built alongside -
loops every buoy with Stage 0 output, continues past individual
failures, aggregates into one comparison CSV) for the network-wide run.

## Priority 5 — RESULTS from the full 19-buoy network run — CLOSED

19/19 succeeded. **Remarkably consistent signature across 18 of 19
buoys**: wind lag 0-3h (median 3h, two buoys - Deurlo, KwintebankBuoy -
at 0h/contemporaneous, both somewhat more sheltered/inshore), MSLP
consistently leading with the correct pressure-drop-precedes-Hs-rise
sign at every one of those 18, R² (linear) 0.41-0.84 (network median
0.565). This is a genuine, strong, network-wide confirmation - wind
drives Hs almost instantaneously across the whole Belgian coast, not
just at the one buoy tested first.

**Zeebrugge is a dramatic, isolated outlier - and this is a real
finding, not noise.** R²=0.093 (vs. network median 0.565), wind/Hs
correlation collapses to 0.32 (vs. 0.69-0.84 everywhere else), and
**MSLP relationship inverts entirely** - lags instead of leads (-18h),
positive correlation instead of negative. This is a different physical
regime, not a weaker version of the same one. **Fourth independent
signal now pointing at Zeebrugge being genuinely decoupled from
open-sea dynamics** (previous three: worst tidal notch performance
across all attempts including the fitted-frequency one, its own
singleton spatial cluster, unstable Stage 13 distribution verdict).
Plausible mechanism: as the most sheltered, harbor-interior site in
the network, local wave conditions there are likely dominated by
harbor geometry/vessel wake/reflected energy rather than direct
open-water wind forcing - a genuinely different causal pathway, which
would explain both the collapsed R² and the reversed MSLP sign.

**Recommendation, not yet acted on**: stop trying to fit Zeebrugge into
the same open-sea characterization framework as the other 18 buoys.
Four independent analyses (tidal notch, spatial clustering, Stage 13
stability, wind coupling) now agree it's structurally different, not
just noisier. Worth explicitly labeling it as a harbor/sheltered-site
case study rather than continuing to treat deviations from it as
problems to fix.

Deurlo's directional data (only VMDR buoy that succeeded here) showed a
real regime-dependent pattern worth remembering: 78.5° wind/wave
misalignment in the moderate regime vs. ~30° in energetic/storm regimes
- physically sensible (storms = local wind-sea dominates and aligns;
calmer conditions = older decoupled swell has more relative influence).

## Priority 6 — Forecasting pipeline — implementation outline, NOW UNBLOCKED

**Priority 5 is done - the "don't start until it exists" condition is
satisfied.** And it did more than just unblock this: it gave a concrete
empirical number to design around. Wind alone explains 41-84% of Hs
variance (network median R²=0.565) at 0-3h lag across 18 of 19 buoys -
that's real headroom for an ARMAX model over plain ARMA, stronger
justification than the outline originally anticipated when this was
just "worth checking whether it beats persistence by much."

**Recommended next step specifically**: start with Stage A (persistence
baseline) + Stage B (backtest harness) on Westhinder, per the outline
below - build and validate the harness itself first, since every later
model (ARMA, ARMA-GARCH, ARMAX) reuses it. Don't jump straight to
ARMAX just because Priority 5's R² is enticing - the backtest harness
needs to exist and be trustworthy before any model result means
anything, same discipline as building 11b/Stage 12 before trusting any
uncertainty number this whole project.

**One thing worth deciding before Stage E (ARMAX) specifically**:
Zeebrugge's wind-coupling collapse (Priority 5) means an ARMAX model
using wind as an exogenous input would need very different expectations
there than at the other 18 buoys - probably worth excluding Zeebrugge
from the initial ARMAX buildout entirely rather than expecting one
model design to work network-wide, consistent with treating it as a
structurally different site per the Priority 5 recommendation above.

**Stage A + B — BUILT and tested.** `forecast_utils.py` (shared harness:
`persistence_forecast`, `rolling_origin_backtest`, `summarize_backtest`,
`skill_score`) and `17_forecast_baseline.py` (runnable Stage A script).

**Caught a real off-by-one bug before trusting it, same discipline as
Priority 5's lag-sign fix**: verified the harness analytically first,
not just by eye - ran persistence forecasting on a deterministic linear
ramp, where the exact expected error at horizon h is known in closed
form (`h * slope`). First version was off: `history = series.iloc[:origin]`
excluded the origin index itself, so the "last known value" was
actually at `origin-1`, making every horizon effectively one sample
longer than labeled. Fixed by making `origin` explicitly mean "the last
OBSERVED sample" (history inclusive of origin), re-verified the
analytical test matches exactly, plus a self-comparison sanity check
(persistence scored against itself gives exactly `skill_score=0.0` at
every horizon - the correctness property that has to hold before this
harness can be trusted to score a *different* model against persistence
later).

**Gap handling note - different category than the lag-based fixes
earlier**: point forecasting doesn't need `longest_contiguous_segment`/
`all_contiguous_segments` the way ACF-based stages did (a forecast only
needs one valid origin point and one valid target point, not a whole
contiguous stretch) - but it DOES need the full regularized grid with
gaps preserved as NaN (not dropna'd), since horizons are converted
between hours and samples assuming uniform spacing. Using a dropna'd
series here would silently make every horizon wrong in a new way.
Verified on a synthetic autocorrelated series with an injected gap:
RMSE increases monotonically with horizon as expected (0.071m at 1h -
0.321m at 24h), gap correctly causes a small, sensible dropoff in
usable origin count near the gap rather than a crash or corrupted
horizon.

**Fixed a real O(n²) performance bug found on the actual Westhinder
run (630,262 samples), not caught by the synthetic tests (only 20,000
samples - not large enough to expose it).** `rolling_origin_backtest`
was passing `persistence_forecast` the FULL history from the start of
the series at every origin (`series.iloc[:origin+1]`), even though
persistence only ever needs the last few samples before the origin
(bounded by `max_lookback_samples`). Calling `.dropna()` on an
ever-growing slice at ~52,000 origins turned the whole backtest into
O(n²) work - user reported it still running after 5+ minutes. Fixed by
adding a `history_window` parameter that bounds the slice passed to
`forecast_fn`; Stage A passes `max_lookback_hours + small buffer`,
since that's structurally all persistence can ever use. Re-verified all
three correctness tests still pass with the bounded window (identical
results, not just faster), then measured the actual fix at real scale:
**17 seconds at 630,262 samples**, down from 5+ minutes and still not
finished. Worth remembering for Stage C (ARMA) and beyond: those models
genuinely need more history than persistence does, but "how much" is a
real design question, not "just pass everything" - refitting from
scratch at every one of tens of thousands of origins is its own
performance problem to solve deliberately when that stage is built, not
something to discover the hard way again.

**Run it**:
```bash
python 17_forecast_baseline.py --buoy WesthinderBuoy --var VHM0
```
This establishes the actual baseline numbers Stage C (ARMA) will need
to beat - run this before building anything further.

Outline below unchanged from before for Stages C-F, still the plan for
how to actually build those when picked up. Reminder of the primary
metric these will use: **skill score = 1 - (model_MSE /
persistence_MSE) per horizon**, not raw RMSE - matches the standing
plan note that raw RMSE isn't the right comparison metric here, and
`forecast_utils.skill_score()` already implements exactly this.

**Stage C — BUILT and tested.** `forecast_utils.make_arma_forecast_fn`
+ `18_forecast_arma.py`. Design deviated from the original outline in
one deliberate way: forecasts **raw Hs directly**, not the Box-Cox/
detided series - avoids needing to invert the transform and reconstruct
the tidal component (Stage 03b doesn't currently save the harmonic
coefficients for that), and since M2 is itself deterministic and
predictable at 30-min sampling, an ARMA model picking up on it is
legitimate short-horizon skill, not something to filter out first.
Fits ONCE (order search on a modest chunk, final fit on a larger one),
reused via `.apply(refit=False)` at each backtest origin with a
BOUNDED state window - verified in isolation first that this API
preserves fixed parameters exactly (not silently re-estimating) before
trusting it in the full pipeline.

**Two real bugs caught by testing before this was trusted:**
1. **Editing accident**: an earlier `str_replace` while adding
   `make_arma_forecast_fn` accidentally orphaned `skill_score`'s body as
   dead code inside the new function (still syntactically valid Python,
   so `py_compile` didn't catch it - only an actual import check did).
   Fixed, then re-verified all of Stage B's original correctness tests
   still pass, not just that it imports.
2. **Real design bug, not just a typo**: hardcoded `d=1` by copying
   Stage 04's differencing choice, without checking whether it actually
   applies to raw Hs (Stage 04's `d=1` was tuned for the detided
   series, a different object). Caught via a synthetic validation test:
   built a strongly mean-reverting (genuinely stationary, `d=0`) AR(1)
   process, and the forced-`d=1` version over-differenced it, destroyed
   the mean-reversion signal, and lost to persistence at every horizon
   beyond 1h - exactly backwards from what a correctly-specified model
   should do here. Fixed by searching `d` instead of assuming it;
   re-ran the same test and the model correctly identified `(1,0,0)` -
   exactly the true generating process - with skill scores now positive
   and GROWING with horizon (0.08 at 1h -> 0.29 at 24h), the correct
   physical signature since persistence never reverts to the mean while
   a correctly-fit AR(1) does. This is worth remembering as a general
   principle, not just a one-off fix: don't carry a modeling assumption
   from one series to a structurally different one without re-checking
   it, even within the same project.

**Fixed - repeated ValueWarning flood on real data**: `.apply()` inside
`make_arma_forecast_fn` fired a `ValueWarning` on every single backtest
origin (~2000+ times on the real run) because slicing a DatetimeIndex
drops its `.freq` attribute even though the underlying spacing is still
regular - statsmodels then re-infers frequency on every call and warns
each time. Fixed by explicitly setting the frequency once (derived from
`dt_hours`, already known) instead of leaving it to be re-inferred.
Verified: 0 warnings after the fix, byte-identical results before/after
- confirms this only cleaned up output, didn't change any number.

**Fixed a real, more important bug found on the actual Westhinder
run**: ARMA's backtest only covered 6.9% of Westhinder's record (the
36,170-sample longest contiguous segment, ~2012-2014) - but the
reported skill score was comparing against Stage A's persistence
baseline computed over the FULL 36-year record. Not a fair comparison -
if that one segment happens to be atypically easy or hard to forecast
relative to the rest of the record, the skill score would be biased by
that mismatch alone, independent of whether ARMA is actually good.
Fixed: `18_forecast_arma.py` now computes its own LOCAL persistence
baseline on the exact same segment and origins ARMA was tested on, and
that's what the skill score is computed against. The full-record
baseline is still shown, but explicitly labeled "for reference only -
NOT a fair comparison, different coverage" rather than silently
implied to be comparable.

**One genuinely good sign from the same real run, worth keeping**: the
order search picked `d=0` (no differencing) on the real Westhinder
segment - matching Stage 03's own ADF result for this buoy (already
said "stationary"). Two independently-built stages agreeing is a
positive cross-check that the model isn't badly misspecified, even
though the skill-score comparison itself needed fixing.

**Run it** (after Stage A has been run for the same buoy - only needed
for the "reference only" full-record context now, since the actual
skill score uses a local baseline computed automatically):
```bash
python 17_forecast_baseline.py --buoy WesthinderBuoy --var VHM0
python 18_forecast_arma.py --buoy WesthinderBuoy --var VHM0
```
Re-run on Westhinder with the coverage fix to get the corrected,
fair skill-score numbers - the previous run's 0.057-0.193 range should
be treated as provisional until this re-run confirms or revises it.

**Ran on 3 more buoys (BolVanHeistBuoy, ZeebruggeZandopvangkadeBuoy,
A2Buoy) - cross-buoy picture, and two real anomalies investigated and
resolved.**

- **BolVanHeistBuoy**: cleanest result of all 4 buoys tested. Same
  order (2,0,2) as Westhinder, skill positive and growing most cleanly
  with horizon (0.063 -> 0.233) - strongest confirmation the
  growing-skill pattern is a real network property, not Westhinder-
  specific.
- **A2Buoy - investigated and FIXED, general consequence for the
  pipeline, not just this buoy.** Original run picked `d=1` (the only
  one of 4 buoys to do so) and showed NEGATIVE skill at every horizon
  beyond 1h (-0.007 to -0.019). Hypothesis: an ARIMA(p,1,q) model
  structurally discards mean-reversion and degenerates toward
  persistence at long horizons - AIC (in-sample fit) can't penalize
  this for a forecasting use case. **Confirmed by forcing `--d-range 0`**:
  every horizon flipped positive (0.074 to 0.216), matching Westhinder/
  BolVanHeist's pattern almost exactly. **Fixed the default**:
  `18_forecast_arma.py --d-range` now defaults to `0` only (was
  `0,1`) - d=1 requires explicit opt-in. This is a real, generalizable
  finding: AIC-based `d` selection can actively hurt forecast skill
  even while improving in-sample fit; don't trust it blindly for this
  purpose again.
- **ZeebruggeZandopvangkadeBuoy - investigated, sharpened the
  hypothesis, not yet fully resolved.** Persistence RMSE was
  non-monotonic with a dip at 12h, and ARMA's skill went negative at
  exactly that horizon (only negative case among all 4 buoys at any
  horizon besides A2's d=1 issue). Fine-grained horizon sweep (8-16h in
  ~1h steps) confirmed a real local minimum at 12h - but 12.5h came out
  nearly identical to 12h (0.1499 vs 0.1493), meaning the dip is NOT
  precisely centered on M2's exact half-period (6.21h... i.e. the 12h
  point relative to M2's 12.4206h full period) the way pure residual M2
  leakage would look. **Revised hypothesis: a DIURNAL (24h) signal,
  not semi-diurnal M2** - 12h would be exactly its half-period. Added a
  general diurnal-band check to `02_eda_diagnostics.py` (not
  Zeebrugge-specific - useful network-wide, since Stage 03b's notch has
  only ever targeted M2, so a real 24h signal would currently pass
  through completely untouched at every buoy). Validated on a synthetic
  case with a known injected 24h signal before trusting it (correctly
  detected, ratio scaled with injected amplitude as expected).
  **Not yet run on real Zeebrugge data** - next step:
  ```bash
  python 02_eda_diagnostics.py --buoy ZeebruggeZandopvangkadeBuoy --var VHM0
  ```
  Check the new `diurnal_ratio_raw` value in the output - a high ratio
  (like the M2 ratio has always been for this buoy) would confirm the
  revised hypothesis and point at a genuinely new characterization gap
  (Stage 03b's notch needs a diurnal term, not just M2 harmonics) worth
  its own fix, separate from the M2-frequency work already closed in
  Priority 4.

**Stage D - ARMA-GARCH — BUILT and tested.** `19_forecast_arma_garch.py`.
Two-step approach (ARMA mean via the same order-search+fit as Stage C,
GARCH(1,1) fit separately on ARMA's in-sample residuals) rather than
joint estimation - simpler, and the `arch` package doesn't support
general ARMA with MA terms as its mean model anyway. `select_order`
moved into `forecast_utils.py` so Stage D could reuse it without
duplicating Stage C's logic (Python module names can't start with a
digit, so direct import from "18_forecast_arma" isn't possible) -
verified Stage C still produces identical results after the refactor
before trusting it.

**Real bug found and fixed, not superficial - this was the actual
substance of building this stage correctly.** First implementation used
GARCH's per-step innovation-variance forecast directly as the
prediction-interval width. Tested on a synthetic AR(1)+GARCH(1,1)
process with genuine persistence (phi=0.98, matching this pipeline's
real situation per 11b's confirmed long timescales) before trusting it
on real data - and it failed clearly: PI coverage collapsed from 0.845
(1h) to 0.358 (24h) against a 95% target, while interval width stayed
nearly flat (0.299 to 0.303) even as empirical RMSE tripled (0.111 to
0.360). Root cause: GARCH's per-step variance is the variance of the
innovation AT step h, not the variance of the ACCUMULATED forecast
error by step h - these coincide only at h=1. For a persistent process,
GARCH's per-step variance converges to a flat unconditional value
quickly while true forecast-error variance keeps growing as long as AR
memory hasn't decayed - exactly the observed symptom. Fixed by properly
weighting and summing GARCH's per-step variances through the ARMA
model's own impulse-response (MA(infinity)) coefficients - confirmed
`statsmodels`'s `.impulse_responses()` matches the analytical AR(1)
formula before using it, and that these weights depend only on the
fixed ARMA coefficients (computed once before the backtest loop, not
per origin - same performance discipline as everywhere else in this
forecasting work). After the fix: coverage 0.933-0.953 at every horizon
(right at the 95% target, and stable, not degrading), interval width
now properly grows with horizon (0.418 to 1.338) tracking RMSE's own
growth shape. This is standard practice for ARMA-GARCH multi-step
forecasting, documented here specifically because it's easy to get
wrong by treating GARCH's raw output as already being the answer -
worth remembering if this pattern comes up again in any future model.

**Run it**:
```bash
python 19_forecast_arma_garch.py --buoy WesthinderBuoy --var VHM0
```
Real per-horizon coverage and interval-width numbers will show whether
this calibration quality holds on actual Hs data, not just the
synthetic validation case - the true test, same as every other stage
built today.

**Real result on Westhinder: coverage 0.939-0.968 across all 5
horizons** - excellent, matching (even slightly exceeding) synthetic
validation quality. But it surfaced a genuine numerical edge case,
investigated and fixed properly rather than either dismissed or
over-alarmed about.

`alpha+beta` fit to `1.0000000000272715` - `1-alpha-beta = -2.727e-11`,
floating-point noise around exactly 1.0 (IGARCH), not a meaningfully
explosive process. The standard closed-form formula technically divides
by ~zero here (`uncond_var` undefined). Investigated why the results
still came out fine: with `alpha+beta` this close to 1,
`(alpha+beta)^(h-1)` stays ~1.0 at every tested horizon, which
algebraically cancels the degenerate `uncond_var` out of the formula -
it happened to gracefully degenerate to the correct IGARCH limit
(constant variance forecast, since shocks are permanent) through
cancellation, not by design. Confirmed this wasn't a coincidence
specific to bad luck: `sigma2`'s own starting point also relied on
`beta` decaying away a similarly-degenerate initial condition within
the 400-sample state window (true here: beta=0.83, 0.83^400 ~ 3e-33,
but not guaranteed for a higher beta).

**First fix attempt was itself WRONG - caught on the real re-run, corrected
a second time.** `garch_state_and_forecast` was changed to (1) initialize
`sigma2` from the window's own sample variance instead of the possibly-
degenerate `uncond_var` (this part was correct and stayed); (2) return a
FLAT variance at every horizon for the near-unit-root case, reasoning
"shocks don't decay, so the forecast shouldn't change." **This reasoning
was backwards.** Re-ran on real Westhinder data after this "fix" and
coverage got WORSE, not better - collapsed to 0.803 at 24h (vs. 0.968
before touching anything), a real regression, not an improvement.

Re-derived properly: the general GARCH(1,1) multi-step recursion is
`f(h) = omega + (alpha+beta)*f(h-1)`. At `alpha+beta=1` exactly, this
becomes `f(h) = omega + f(h-1)` - LINEAR growth
(`f(h) = f(1) + (h-1)*omega`), not a constant. "Shocks never decay"
means uncertainty about future variance keeps ACCUMULATING without
bound as horizon grows - the opposite of what a flat forecast implies.
A flat forecast understates genuinely growing uncertainty, which is
exactly why coverage collapsed at long horizons (interval too narrow).
This also means the ORIGINAL numbers (before touching any of this) had
been correct all along via legitimate algebraic cancellation in the
general closed-form formula - not a coincidence, real math working out
- and my first "fix" broke something that wasn't actually broken.

**Second fix, validated more rigorously than the first one was**: linear
formula `sigma2 + h*omega`. This time: (1) direct numerical iteration of
the TRUE recursion at alpha+beta=1.0 exactly, confirmed EXACT match
(not just close) against the closed-form linear formula; (2) unit test
on Westhinder's real fitted parameters confirms proper linear growth
(step size = omega exactly, `0.002049 -> 0.006279` over 48 steps);
(3) stationary synthetic test (alpha+beta=0.90, unaffected by this
branch) still byte-identical - no regression. **Lesson for next time**:
the first fix was shipped after only checking it didn't crash and
produced "reasonable-looking" flat numbers, not after checking it
against a known-correct reference the way the ORIGINAL closed-form
formula was checked against the `arch` library early on. A result that
runs without error and looks plausible is not the same as being
verified - worth remembering given how much of today's actual value
came from exactly this kind of check catching real bugs before they
shipped quietly.

**Standalone finding worth keeping**: `alpha+beta >= 0.98` (IGARCH-like
persistence) at Westhinder is itself real and substantive, not just a
numerical curiosity - variance shocks barely decay within any
practically forecastable horizon here, consistent with (and actually
stronger than) Stage 07's universal ARCH finding and 11b's long
persistence timescale for this specific buoy.

**Stage E - ARMAX — BUILT and validated.** `20_forecast_armax.py`. Uses
`SARIMAX` with lagged ERA5 wind speed as exogenous regressor, lag
default 3h matching Priority 5's Westhinder finding.

**Critical design issue handled explicitly, not glossed over**: genuine
multi-step forecasting needs FUTURE exog values at the forecast
horizon, but real future wind isn't known any more than real future Hs
is - using the true future ERA5 wind in a backtest would be look-ahead
bias (the model would appear to work by secretly being handed the
answer). Handled honestly: exog is only "real" up to (origin +
lag_samples); beyond that, the last genuinely-known wind observation is
persisted forward (same honest simplification Stage A uses for Hs
itself), not real future data. This makes a specific, testable
prediction: ARMAX's advantage over plain ARMA should be concentrated
within the lag window and should shrink or reverse beyond it.

**Validated on synthetic data built specifically to have wind-driven
predictability ARMA alone can't see** (Hs driven by wind lagged 3h,
deliberately WEAK own autocorrelation so plain ARMA has little to work
with). Result matched the design prediction almost exactly: skill vs.
plain ARMA (Stage C) = **0.79 at 1h, 0.95 at 3h** (dramatic, right at
the lag boundary) -> 0.42 at 6h -> 0.08 at 12h -> **NEGATIVE -0.34 at
24h** (worse than plain ARMA once relying on stale persisted wind).
This is strong evidence the honest-exog design is actually working as
intended, not secretly cheating - if it were using real future wind,
skill would stay high or improve at long horizons instead of collapsing
and reversing sign.

**Practical implication, worth remembering when this gets used for
real**: this isn't "ARMAX is better than ARMA," it's "ARMAX is much
better within the wind lag window and can be actively worse beyond it"
- an operational deployment should probably switch between ARMAX (short
horizons) and plain ARMA or persistence (longer horizons) rather than
use one model uniformly across all horizons.

**Run it on real Westhinder data**:
```bash
python 20_forecast_armax.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0 --exog-lag-hours 3.0
```
Requires Stage 17 and Stage 18 to have already run for this buoy (for
the persistence and plain-ARMA comparisons). Given real Hs has much
more complex dynamics than the clean synthetic test, don't expect
numbers this dramatic - but the qualitative shape (strong short-horizon
gain, decaying or reversing beyond the ~3h lag window) is the real
thing to check for.

**Real Westhinder result - diverged from the synthetic prediction in an
informative way, not a red flag.** Skill vs. plain ARMA came back
positive at EVERY horizon tested (0.046 to 0.253), peaking at 6h
(0.253) rather than reversing negative by 24h the way the synthetic
test predicted. Likely explanation: the synthetic test used a fast
~20h-period oscillating wind signal specifically to isolate the lag
effect cleanly - real wind has genuine persistence of its own (storms
evolve over many hours, not instantly), so "persist the last known
wind value" beyond the 3h lag - the honest fallback - is a much less
damaging approximation for real wind than for the synthetic driver.
This is physically sensible, not a sign the honest-exog design is
broken.

**Follow-up check built to test this explanation directly, not just
assert it**: added a wind persistence timescale check to
`16_wind_wave_coupling.py`, reusing 11b's exact ACF/significance-band
method (moved `integral_timescale()` into `utils.py` so it could be
shared - same reasoning as moving `select_order` for Stage D/E, Python
module names can't start with a digit so direct import from
"11b_dependence_structure" isn't possible). Verified the refactor
didn't change 11b's own behavior (regression-tested, identical output)
before using it elsewhere.

Validated on a synthetic AR(1) wind process with a KNOWN true
persistence timescale (phi=0.9, theoretical tau=57.0h) - recovered
29.7h, a real ~2x undershoot. Not treated as a new bug: same
significance-band method already trusted throughout this project
(11b), and finite-sample noise in a slowly-decaying ACF can trigger the
"5 consecutive lags inside the band" stopping criterion earlier than
the true population value justifies - a known character of this
method at moderate sample sizes (this test used 3000 samples). Real
ERA5 has far more (48,260 for Westhinder, 2010-2026) so should be less
affected, but worth remembering any single wind-persistence estimate
from this check is directionally informative, not exact - same caveat
that already applies to every other integral-timescale number in this
project.

**Run it on real Westhinder data**:
```bash
python 16_wind_wave_coupling.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
```
Check whether the recovered wind persistence timescale is in the same
ballpark as where ARMAX skill actually peaks (6h) and stays positive
through (24h) - if wind's own persistence extends out that far or
beyond, that directly confirms the explanation for why ARMAX didn't
degrade the way the synthetic test predicted.

**Real result: raw wind gave 343.0h (~14 days) - too large to be
plausible, and NOT confirming the explanation as cleanly as hoped.**
343h is far beyond any real synoptic/storm timescale. Diagnosed as the
SAME category of mistake 11b itself originally made (before the
raw-vs-detided fix): an ACF-based integral timescale can't distinguish
genuine short-term memory decay from real-but-slow deterministic
seasonal cycling (wind is windier in winter, calmer in summer - very
plausibly the primary physical driver of the Hs seasonality Stage 15
already confirmed) - both inflate the apparent timescale. Just like 11b
needed the DETIDED series, this needed a DESEASONALIZED wind series,
not raw ERA5 wind directly.

**Fixed**: subtract a 30-day rolling mean before computing the
timescale (removes seasonal-scale drift, preserves genuine day-to-day
storm-scale persistence - a simple approach, not a full STL
decomposition, stated as such rather than overclaiming rigor). Both raw
and deseasonalized values are now printed, so the contamination is
visible rather than hidden. Validated on a synthetic case with BOTH a
known short-term AR(1) persistence (phi=0.9, theoretical tau=57.0h) AND
a real injected annual seasonal cycle: raw computation gave 1132.7h
(~47 days - same order-of-magnitude inflation pattern as the real
343h result), deseasonalized gave 41.4h - much closer to the true 57h,
consistent with the ~2x-undershoot finite-sample bias already
characterized for this method. Confirms the fix actually separates
seasonal contamination from genuine short-term persistence, not just
theoretically but on a controlled test with a known answer.

**Second bug found on the real re-run - the rolling-mean fix was ALSO
wrong, in a different way.** Real Westhinder result came back
`DESEASONALIZED: -11.0h` - a mathematically impossible negative
timescale. Root cause: subtracting a rolling mean is a known source of
spurious NEGATIVE autocorrelation in the residual - a mechanical
artifact of local-average filtering itself (nearby points tend to
deviate from their rolling average in opposite directions), not a real
property of wind. This directly corrupted the `1 + 2*sum(rho)` sum the
timescale formula depends on.

**Fixed by reusing a method already TRUSTED in this exact pipeline
instead of a new, less-tested one**: harmonic regression at the annual
period (365.25 days) - literally the same approach Stage 03b already
uses for the M2 tidal notch, just a different fundamental frequency. A
fitted deterministic sinusoid doesn't mechanically induce the
rolling-filter artifact. Also added an explicit check: if the
timescale still comes out negative, print a clear "mathematically
impossible, do not use this number" warning rather than silently
reporting a nonsensical value, as a safety net regardless of the
deseasonalizing method's own correctness. Re-validated on the SAME
synthetic AR(1)+seasonal test: raw still 1132.7h (unchanged, as
expected), deseasonalized now **78.1h** - positive, same order of
magnitude as the true 57h (this time slightly over rather than under,
a different direction than the rolling-mean version's 41.4h - both
within the range of expected estimator variability for this method,
neither pathological).

**Re-run needed on real Westhinder data with the SECOND fix** - the
343h AND -11.0h numbers should both no longer be trusted:
```bash
python 16_wind_wave_coupling.py --data-dir data_multiyear --buoy WesthinderBuoy --var VHM0
```
Check the NEW deseasonalized value specifically (not the raw one) as
the actual comparison point against where ARMAX skill peaks/decays
(6-24h) - this is the number that should actually answer the original
question, which the raw (seasonally-contaminated) value could not.

**RESOLVED - real result confirms the hypothesis cleanly.** Second fix
gave 57.3h on real Westhinder data - positive, sensible, no artifact,
and close to the synthetic validation's recovered value (78.1h vs. true
57h - same order of magnitude, consistent with the method's known
precision rather than a fluke). This directly explains the original
finding: the entire ARMAX-tested horizon range (1-24h) sits well inside
wind's own 57.3h decorrelation timescale - at 24h, wind has moved
through less than half its own persistence timescale, so it hasn't
meaningfully decorrelated by the longest tested horizon. That's why
"persist last known wind" stayed a reasonable approximation throughout
(ARMAX skill vs. plain ARMA stayed positive at every horizon,
0.046-0.253), unlike the synthetic validation test's fast ~20h-period
oscillating wind, which decorrelates on a timescale close to its own
period and made the persistence approximation break down quickly. This
whole investigation thread (real ARMAX result -> wind persistence
hypothesis -> seasonal contamination bug -> rolling-mean artifact bug ->
confirmed answer) is now closed with a quantified, not just plausible,
explanation.

**Stage F - exceedance forecasting — BUILT and validated. LAST STAGE IN
THE A-F FORECASTING SEQUENCE - all six now built and tested.**
`21_forecast_exceedance.py`. Reframes as "will Hs exceed threshold X at
ANY point in the next H hours" - a window/max-based label (the
operationally useful framing), not a single-point check. Logistic
regression on 3 features computed strictly from data at or before the
origin (no look-ahead): current Hs level, 3h trend, 24h rolling
volatility (motivated directly by Stage 07's universal ARCH finding -
volatility itself should carry exceedance-risk information, not just
level and trend). Thresholds default to the buoy's own 90th/95th Hs
percentiles (consistent with EVA's percentile-based approach elsewhere
in this pipeline). Evaluated via Brier score and ROC-AUC against a
naive constant-probability (historical base rate) baseline, giving a
Brier Skill Score directly analogous to the point-forecast skill score
used in Stages C/D/E.

Simpler than Stages C/D/E in one respect: logistic regression
prediction is O(1) per call (no Kalman-filter-style state to manage),
so this doesn't need the fit-once-apply-cheaply performance pattern
those stages required - fit once on a training window, predict
directly at each backtest origin.

**Validated on synthetic data with genuine level/trend/volatility-driven
exceedance risk** (AR(1) with occasional stochastic volatility bursts):
Brier skill score +0.456 to +0.683 across all threshold/horizon
combinations (all positive, meaningfully beating the naive baseline),
ROC-AUC 0.871-0.985 (far better than random). Confirms the full
machinery - feature construction, window-based labeling, train/test
split, scoring - works correctly and can detect real signal when it's
present.

**Connects to the parked industrial-monitoring idea from an earlier
external review** (CUSUM/EWMA) - this is where that would actually plug
in, as a monitoring layer consuming this stage's exceedance
probabilities, not a standalone addition. Not built now, still parked,
but the natural foundation for it now exists.

**Run it on real data**:
```bash
python 21_forecast_exceedance.py --buoy WesthinderBuoy --var VHM0
```
Real result: Brier skill +0.424 to +0.680, ROC-AUC 0.827-0.967 across
all threshold/horizon combos - strong, and closely matching the
synthetic validation's numbers. Read carefully rather than just
celebrated: this buoy already shows extreme persistence across
multiple independent analyses today (11b's ~100-116h timescale,
near-IGARCH variance persistence in Stage D, ARMAX skill staying
positive to 24h because wind barely decorrelates that fast) - so
"current level+trend+volatility predicts near-term exceedance very
well" is largely a restatement of that same underlying persistence in
a new frame, not brand-new information about the physics. Genuinely
useful operationally either way (a DEME-style user cares about
threshold crossings, not about the underlying statistics).

**Wind-augmented version — BUILT and validated.** Added `--use-wind` to
`21_forecast_exceedance.py`, current ERA5 wind_speed as a 4th feature.
Simpler than Stage E in one respect: exceedance classification only
ever needs wind AT OR BEFORE the origin as a "current conditions"
snapshot - unlike ARMAX, it's not forecasting the exog variable
forward, so none of Stage E's honest-persistence-beyond-the-lag
complexity applies here.

Validated on a synthetic case designed so wind directly drives
exceedance while Hs itself has WEAK own memory (the opposite regime
from the earlier Hs-only validation, which had strong Hs persistence) -
confirms wind-augmentation correctly detects real signal when Hs's own
history can't see it: Brier skill delta and AUC delta both positive at
every threshold/horizon combination tested, strongest at the shortest
horizon (6h, where wind's contemporaneous relationship with Hs matters
most) and decaying at longer horizons - the expected pattern for a
"current state" snapshot feature.

**Run on real Westhinder data**:
```bash
python 21_forecast_exceedance.py --buoy WesthinderBuoy --var VHM0 --use-wind --data-dir data_multiyear --era5-dir meteo_era5
```
**Real result: wind delta positive at every threshold/horizon (6/6),
roughly flat with horizon (p90: 0.028/0.034/0.030 across 6h/12h/24h).**
Prediction that this would be small (given Westhinder's already-strong
Hs-only persistence) was wrong - the delta is real and non-negligible
(~5-10% relative improvement on top of an already-strong base).

**Ran A2Buoy as the comparison test - hypothesis about WHY the delta
was flat at Westhinder was refuted, but a better, more general
explanation emerged that survives both buoys.** Original hypothesis:
Westhinder's flat delta was specifically due to its unusually long wind
persistence (57.3h, confirmed earlier); predicted A2 (presumably less
wind-persistent) would show the delta DECAYING with horizon instead.
**Wrong direction**: A2's delta actually GROWS with horizon, more than
Westhinder's (p90: 0.043/0.053/0.057 across 6h/12h/24h) - the opposite
of the prediction.

**Revised explanation, now supported by both results instead of just
one**: working through the actual numbers, at BOTH buoys the wind-
augmented skill decays only slightly less steeply than the Hs-only
skill as horizon grows (Westhinder: 0.250 vs 0.247 drop 6h->24h; A2:
0.247 vs 0.233 drop). The growing delta isn't really about wind
"lasting longer" at one buoy - it's that Hs's own level/trend/
volatility features become progressively less informative as horizon
extends, while wind (reflecting broader synoptic conditions, not just
the buoy's own recent state) retains more RELATIVE relevance over the
same window. This is a general property of the method, not a buoy-
specific wind-persistence-timescale story - the original mechanistic
explanation was too narrow and didn't survive a second data point.

**Practical bottom line, now on firmer ground (2/2 buoys)**: wind-
augmentation helps consistently and doesn't fade with horizon across
both buoys tested - a real, generalizable finding for exceedance
forecasting, even though the first proposed mechanism for it was wrong.

**Suggested script breakdown** (numbering TBD when actually built,
don't lock it in now): a shared backtest harness module, then one
script per model (baseline, ARMA, ARMA-GARCH, ARMAX) that all call the
same harness, rather than duplicating the walk-forward logic four
times. Start with Westhinder only (cleanest, longest, most-validated
record in the whole network) before generalizing to other buoys.

## GPD shape parameter — external literature cross-check, NEW action item

From an external review: internal cross-validation (three-plus analyses
agreeing on something) is good practice but stays self-referential -
worth anchoring at least one finding against a published external
source rather than only against this project's own other stages.

**Done**: found a directly relevant published reference - Caires
(2011), JCOMM Technical Report No. 57 (WMO/IOC), a shallow-water North
Sea GPD analysis (Schiermonnikoog noord buoy, 19m depth, comparable
shelf conditions to the BCZ). Reported xi = -0.12 to -0.13 (95% CI:
-0.37, 0.09), vs. -0.07 to -0.08 for a deep-water comparison site
(NDBC 46005, Pacific, 2780m). The paper's own explanation for the
shallow/deep difference: shallow-water waves are depth-limited -
**the exact physical mechanism already hypothesized for this
project's own bounded-tail finding**, now with independent published
confirmation of the sign and general shallow-water magnitude, not just
internal agreement.

**Real discrepancy worth checking, not glossing over**: this network's
GPD xi range (-1.31 to -0.04, per the README) extends far more negative
than the published shallow-water reference point (-0.13). Two
non-exclusive explanations:
1. Several BCZ buoys are plausibly shallower than the 19m reference
   site (Nieuwpoort, Blankenberge, Zeebrugge-area) - genuinely stronger
   depth-limiting could be real, not an artifact.
2. The most extreme xi values may be the same small-peak-count buoys
   already known to have unreliable, wide confidence intervals (e.g.
   A2Buoy's CI crossed zero even at a properly window-corrected EVA
   run).

**RESOLVED - real network run, and it's a bigger finding than
expected.** The old -1.31 extreme was itself stale, from BEFORE
today's Stage 08/11b fix (dynamic declustering window per buoy). The
corrected, properly-declustered network run gives a dramatically
tighter range: **-0.544 to +0.327** - the fix itself already resolved
most of this concern, not just a footnote to it.

**Real outlier is Blankenberge, not the most-negative buoy**: the ONLY
positive xi in the network (+0.327, CI crosses zero), with both the
fewest storm peaks (43) AND the shortest record (2.77 years) of any
buoy. Positive xi implies an unbounded tail - physically implausible
for depth-limited shallow-sea waves. Near-certainly a small-sample
artifact. **The automated cross-check missed this** (tuned for
|xi|>0.5, Blankenberge's magnitude is 0.327) - worth extending the
check to also flag positive xi specifically, since that's the
physically implausible direction regardless of magnitude.

**Raversijde1Buoy (xi=-0.544, the most negative) - automated flag was
slightly overstated, worth correcting**: CI is [-0.929, -0.271] - wide
(0.658), but entirely NEGATIVE, doesn't cross zero. That's not "likely
artifact," it's "confidently bounded, imprecisely estimated" (67
peaks, second-fewest in network) - a real distinction the check's
OR-based logic (`wide CI OR crosses zero`) blurred together. Consider
splitting these into separate messages if this check gets used again.

**The reliable core of the network** (many peaks, narrow CIs -
Trapegeer, ScheurWielingen, AkkaertSouthwest, Wandelaar) clusters
tightly around **-0.27 to -0.32** - still ~2x more negative than
Caires's shallow-water reference (-0.12 to -0.13), down from the old
~10x gap, and now explicable as BCZ buoys being shallower/more
depth-limited than that reference site's 19m, rather than a
methodology concern.

**Action item, still open**: update `README.md`'s findings summary,
which currently states the stale "-1.31 to -0.04" range - replace with
the corrected range and the Blankenberge/Raversijde1 nuance above.

## Priority 7 — Provenance metadata — BUILT and tested, one stage wired in as the pattern

`utils.get_provenance(input_path=None, random_seed=None)` added:
timestamp, Python version, key package versions (numpy/pandas/scipy/
statsmodels/xarray/matplotlib/arch/scikit-learn/copernicusmarine/
cdsapi - gracefully `None` if not installed, not an error), git commit
hash + working-tree-dirty flag, and an input-file fingerprint
(mtime+size, not a full content hash - deliberate tradeoff, stated in
the docstring: hashing a 500k+ row CSV on every run would cost real
time for a "did this change" question mtime+size already answers
correctly for every realistic scenario here).

**Tested properly, not just run once**: git hash retrieval and dirty-
tree detection verified in an isolated real git repo (not just this
sandbox, which isn't a git clone and correctly returns `None`) - clean
tree correctly detected as clean immediately after a commit, correctly
detected as dirty after modifying a tracked file. Caught my own test
design flaw along the way (first attempt copied `utils.py` into the
test repo, which made an untracked file that correctly counted as
"dirty" for the wrong reason) and fixed the test, not the code - worth
remembering that a surprising test result isn't always a code bug.

**Wired into Stage 01 (`01_load_clean.py`) as the worked example**, not
all ~25 scripts at once - deliberate scope limit given this needs real
validation on the actual HPC repo (git hash population specifically)
before blindly copying the same one-line pattern everywhere.
Verified end-to-end: runs cleanly, provenance block appears correctly
in the real output JSON alongside everything else Stage 01 already
reports.

**CONFIRMED on the real HPC environment - Priority 7 fully closed.**
`git_commit_hash` populated correctly (real 40-char SHA-1, verified
programmatically not just eyeballed), `git_working_tree_dirty=true`
correctly reflecting real uncommitted state, `copernicusmarine`/`cdsapi`
versions populated (2.4.1/0.7.7) where the sandbox had shown `null`
(never installed there). Real environment differences already visible
in the first real run - Python 3.10.20 (HPC) vs 3.12.3 (dev sandbox),
and a notable `pandas` version gap (2.3.3 vs 3.0.2, a real major-version
jump with actual breaking changes in places) - a concrete, immediate
demonstration of why this was worth building, not just a nice-to-have.

**Next natural extension, not done yet**: the actual copernicusmarine
API-signature surprises that motivated this whole priority happened in
the DOWNLOAD scripts (`download_belgian_wave_buoys_history.py`,
`download_era5_meteo.py`), not in the numbered analysis stages -
worth prioritizing provenance tracking there specifically over further
analysis-stage rollout, since that's literally where the pain occurred.
Otherwise, the same 2-line pattern can be copied to other stages
incrementally whenever convenient, per the original "cheap enough to
slot in anytime" framing.

## Priority 8 — Physical interpretation stage — BUILT and tested

Added as extra printed blocks + small JSON outputs in
`08_extreme_value_analysis.py` (`interpret_gpd_shape()`) and
`06_distribution_fit.py` (`interpret_distribution()`), not a new stage
- exactly as scoped.

**Stage 08's GPD interpretation is directly calibrated by today's real
Blankenberge finding**, not generic textbook text: heavy tail (xi>0.05)
combined with few peaks (<100) explicitly triggers a small-sample-
artifact warning recommending the Stage 12 CI check - tested against
all three real buoy shapes from today's actual network run
(Westhinder-like bounded, Blankenberge-like heavy-tail-few-peaks,
A2Buoy-like near-zero) plus a positive-control case (heavy tail with
many peaks, which correctly gets a softer "worth double-checking"
message instead of the artifact warning) - all four produced the
right, distinctly-worded interpretation. Verified end-to-end through
the real script (not just the isolated function) on a synthetic
fixture - appears correctly in both console output and the saved JSON.

**Stage 06's lognormal interpretation is similarly earned, not
generic**: explicitly cites Stage 13's real Westhinder finding (every
calendar-era window individually preferred Weibull, yet the pooled
record came out lognormal) as the reason to check Stage 13 before
trusting a lognormal win as a clean physical result. Weibull/Rayleigh
get straightforward physical reads. Tested end-to-end on genuinely
lognormal-shaped synthetic data - correctly triggers the caveat, saved
to a new dedicated `_fit_interpretation.json` rather than disturbing
the existing `fit_summary.csv` schema `summarize_results.py` already
depends on.

**Field-tested on real data, both branches of the logic confirmed**:
Westhinder's Stage 06 run correctly triggered the lognormal/Stage 13
caveat for real (not just synthetic). BlankenbergeBuoy's Stage 08 run
correctly triggered the small-sample-artifact warning for real
(xi=0.2955, 69 peaks) - independently corroborated by its own return
levels (17.38m at 50 years, physically absurd for this coast) and the
`<2x record length` illustrative-only flags on most of its return
periods.

**Clarification worth remembering, found during field-testing, not a
bug**: running a stage script standalone (directly from the command
line) does NOT inherit `run_all_buoys.py`'s dynamic per-buoy argument
injection (Stage 08's persistence-based `--min-separation-hours`,
Stage 10's `--include-period`) unless passed explicitly. A standalone
Westhinder EVA run defaulted to 48h separation (1031 peaks, xi=-0.2315)
instead of the network run's actual 231.4h-declustered result
(483 peaks, xi=-0.3629) - both are "correct" outputs for what was
actually asked, but only one matches the real network's numbers.
Similarly, BlankenbergeBuoy's standalone run reported 69 peaks at the
default 48h separation, not the 43 peaks the earlier network batch run
found at its own dynamically-injected window - same cause, don't be
surprised if a standalone re-run of any buoy doesn't exactly match its
number from a full `run_all_buoys.py` pass unless the same arguments
are passed explicitly.

## Priority 9 — Stage X: four-question diagnostics report — BUILT and tested

`22_diagnostics_report.py` synthesizes across ~12 stage outputs
(01, 03, 03b, 05, 06+interpretation, 07, 08+interpretation, 10, 11,
11b, 12, 13) into four narrative paragraphs per buoy, reading gracefully
- doesn't compute anything new, doesn't require every stage present.

**Tested for the more important property first**: zero stage outputs
present (a buoy nothing has been run for yet) - correctly reports "not
run" for every question instead of crashing. Then tested with a
realistic partial fixture (Stage 01/03/06/08/11b present, Stage 12/13
deliberately absent) - correctly synthesizes the present stages into
readable prose AND correctly flags the missing uncertainty
quantification with an explicit caveat rather than silently omitting
it or leaving a blank section. Saved markdown output verified
well-formed.

**Field-tested on Westhinder (full stage history) - excellent result,
with two real bugs found and fixed along the way.** Every number in
the real report cross-checked exactly against prior real runs (16.3%
missing, 345 segments/61.4%/50 used, M2 ratio 642.4->24.3, xi=-0.363/
483 peaks, CI [-0.518,-0.319], persistence 115.7h - all matched
precisely). But two sections were silently empty that shouldn't have
been: Stage 13 (moving-window stability) and Stage 11 (spatial
cluster), given both had genuinely been run for Westhinder. Checked
rather than assumed "not run" was the real reason - it wasn't:

1. **Stage 13's actual field is nested and a fraction, not flat and a
   percentage**: script assumed `stab.get("pct_windows_agree_with_overall")`,
   the real structure is `stab["window_stability"]["fraction_windows_agreeing"]`
   (0-1 scale). Fixed - now correctly navigates the nested key and
   converts to percentage for display.
2. **Stage 11's actual filename has a `{var}_` prefix the script
   didn't account for**: assumed `buoy_clusters.csv`, real file is
   `{var}_buoy_clusters.csv` (e.g. `VHM0_buoy_clusters.csv`). Fixed.

Both re-tested against fixtures matching the CONFIRMED real schema
(not guessed a second time) before trusting the fix - Stage 13 test
correctly showed 25% agreement (plausibly matching Westhinder's
original, pre-calendar-windowing-fix Stage 13 result from earlier in
the project), Stage 11 test correctly identified singleton cluster
status. **General lesson worth remembering**: graceful degradation
("if the file/field isn't there, say so cleanly") is the right design,
but it also means a wrong assumed schema fails SILENTLY as "stage not
run" instead of loudly as a crash - worth spot-checking a full-history
buoy specifically to catch this category of bug, not just partial
synthetic fixtures where "empty section" and "wrong schema" look
identical from the outside.

## Priority 10 — Statistical fingerprint — BUILT and tested

`23_statistical_fingerprint.py`, two modes: single-buoy (one-page visual
- text summary, regime-fraction bar chart, xi+CI error bar) and
`--network` (fingerprint table across all discovered buoys + PCA 2D
visualization + KMeans clustering by STATISTICAL similarity - a
genuinely different lens than Stage 11's correlation-based spatial
clustering: two buoys could be geographically uncorrelated but
statistically similar, or vice versa). All fields read from schemas
CONFIRMED against actual stage source code before use, not guessed -
applying Priority 9's lesson from the start rather than repeating the
same mistake.

**Real bug found and fixed via testing, not assumed correct**: initial
version hardcoded `n_clusters=4` regardless of dataset size. Tested on
a deliberately-differentiated synthetic set (6 buoys, 2 clean groups -
3 "calm/low-persistence", 3 "stormy/high-persistence") specifically to
check whether clustering recovers real structure, not just "runs
without crashing" - it over-split into 4 clusters instead of the true
2 (no cross-contamination between the true groups, but not the clean
split that existed). Fixed by selecting k via silhouette score across
a range instead of a fixed guess - re-tested, now correctly and
exactly recovers the true 2-group structure (silhouette=0.743, zero
misclassification). This k-selection issue would likely have been
invisible at the real n=19 scale (where a fixed k=4 happens to be
plausible, similar to Stage 11's own 4 geographic clusters) - worth
remembering that a bug can hide at production scale and only surface
at a smaller, deliberately-designed test size, which is exactly why
that test was worth building rather than only testing at "realistic"
scale.

**Real network run (19 buoys, k=5, silhouette=0.346)**: three findings,
each independently interesting - Zeebrugge its own cluster, Blankenberge
its own cluster, CadzandBoei+Deurlo paired separately from everyone
else. Two large clusters (6 vs. 9 buoys) suspiciously tracked
deployment era (shorter-record buoys vs. longer-record ones) rather
than obviously-behavioral properties - `record_years`/`pct_missing`
were included as raw features alongside genuinely behavioral ones,
risking exactly this confound.

**`--exclude-record-length` re-run - a real correction, not just
confirmation.** Prediction was partially right, partially wrong in an
informative way:
- The 6-vs-9 deployment-era split **dissolved completely** as
  predicted - confirms it was a data-availability artifact.
- **CadzandBoei+Deurlo's pairing survived, and got STRONGER**
  (silhouette improved 0.346->0.391 after exclusion) - now the more
  defensible finding of the two, not a weaker one. Consistent with the
  original instrument-class hypothesis from the very start of this
  project (10-min sampling, VTPK/VMDR sensor set).
- **Zeebrugge and Blankenberge's singleton status did NOT survive** -
  both absorbed into one 17-buoy cluster. Not predicted, and requires
  walking back part of the earlier framing, not just noting a null
  result:
  - **Zeebrugge**: this fingerprint-clustering result should be
    downgraded from "confirms" to "was ambiguous, driven by
    record-length/missingness metadata (14.3yr, 20.8% missing - both
    real outliers among the 19 buoys) rather than core behavioral
    properties this method can detect." The other independent lines of
    evidence (tidal notch failure, correlation-based spatial cluster,
    Stage 13 instability, wind-coupling collapse) are unaffected - none
    of them depend on this feature set, so Zeebrugge's overall case
    doesn't weaken, but this specific piece of evidence for it does.
  - **Blankenberge**: disappearing here is actually CONFIRMATORY, not
    contradictory - consistent with the GPD xi cross-check's own
    conclusion that Blankenberge is a small-sample artifact (shortest
    record, fewest peaks, widest CI), not a genuine physical outlier.
    Losing its distinctiveness once the most direct record-length
    proxies are removed is exactly what that conclusion predicts.

**Lesson worth keeping**: a check that changes a conclusion (Zeebrugge
downgraded here) is more valuable than one that only confirms it -
worth remembering when deciding whether a "the result didn't replicate
under a stricter test" finding is a failure or the check doing its job
correctly.

## Priority 11 — Assumptions-per-test summary — SUPERSEDED, not built separately

Per this priority's own stated condition ("only worth building if it
actually consolidates... otherwise it's restating existing output with
extra steps"): **Priority 9's diagnostics report already does exactly
this consolidation**. Its "Can I trust the ASSUMPTIONS?" section
already pulls together Stage 03's ADF/KPSS agreement, Stage 06's
interpretation (which carries the KS-validity caveat), Stage 03b's M2
notch quality, Stage 05's Ljung-Box result, and Stage 07's ARCH
detection - built and field-tested on real data. Building a second,
separate script for the same purpose would violate the very condition
that justified building this one at all. Closed without building,
deliberately - the right call here is recognizing redundancy, not
padding out a checklist.

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

## Post-priority extension: HMM regime identification (Stage 24) — BUILT and tested

From a signal-processing-methods review: of the suggested additions
(Kalman/state-space, particle filters, wavelets, matched filtering,
Wiener filtering, HMM, Bayesian change-point detection, robust
statistics, SSA, ICA, EMD, total variation denoising, Gaussian
processes, sequential Monte Carlo), most don't fit this project's
identity as a characterization framework rather than a real-time
estimation one. Three were judged genuinely well-motivated by specific
open questions already in this project (not generic domain popularity):
HMM (extends Stage 10, cheapest), Bayesian change-point detection
(directly addresses the Westhinder drift question - Mann-Kendall tests
for smooth trend, a discrete step-change is a different hypothesis it
can't detect), and SSA (could address Zeebrugge's unresolved compound-
tide/time-varying-parameter hypotheses, harder to implement - Hankel
embedding at real record scale needs truncated/randomized SVD, and
still needs the longest-segment restriction like everything else at
Zeebrugge). Kalman/state-space, wavelets, and the rest were judged
either off-mission or not motivated by a specific question here, per
the review discussion.

**HMM built first (cheapest, easiest reuse of existing infrastructure)**.
`24_regime_hmm.py`, using `hmmlearn.GaussianHMM`. Extends Stage 10's
static regime fractions with TRANSITION PROBABILITY structure - given
you're in a storm regime, what's the probability you're still in one
next timestep - which neither Stage 10 (no memory) nor 11b (one
aggregate Hs timescale, not per-regime) can answer.

**Critical design point, verified empirically before trusting it, not
assumed from documentation**: `hmmlearn`'s `lengths` parameter is
required to fit on contiguous segments separately, for the same reason
nearly every lag-based stage in this pipeline needed gap-awareness.
Tested directly: on a 50-segment synthetic case with a deliberate
state-mismatch at every segment boundary, omitting `lengths` produced
spurious cross-state transition probability of ~0.025 (material
contamination); with `lengths`, the same spurious transitions were
numerically zero (~1e-134). Confirmed the effect scales with segment
count via a negative control - matters MORE on a heavily-fragmented
record (Westhinder: 1905 segments) than a lightly-fragmented one,
opposite of what an unaware implementation might assume.

**Validated end-to-end with a known transition structure**: synthetic
4-state Markov chain with specified transition probabilities (true
dwell times 200/100/66.7/50 samples) and realistic injected gap
fragmentation (28 segments after cleaning). Recovered regime means
essentially exact (true [0.5,1.2,2.2,3.5] vs. recovered
[0.500,1.203,2.199,3.493]); recovered dwell times within 5-18% of true
values (largest error on the stickiest regime, which has the fewest
independent excursions to estimate from - expected finite-sample
behavior, not a flaw); the tridiagonal transition structure (only
adjacent regimes connect in the true generating process) was correctly
recovered with near-zero probability on every non-adjacent "skip"
transition.

**Run it**:
```bash
python 24_regime_hmm.py --buoy WesthinderBuoy --var VHM0
```

**Real result on Westhinder - prediction above was WRONG, and understanding
why is the actually valuable part.** Ran: 697 segments, 96.6% of valid
data used (notably better coverage than 11b's own 61.4%, since this
stage has no analogue of 11b's `--max-segments` cap). Regime means
cleanly separated and monotonic (0.40/0.73/1.17/2.05m, calm to storm).

Dwell times: calm=26.1h, moderate=11.8h, energetic=11.5h,
**storm=21.2h** - a U-shape, not the predicted monotonic decrease, and
**11b's aggregate persistence (115.7h) is 4-10x LARGER than every
single regime's dwell time**, not sitting between them as predicted.

**Why the prediction was wrong - a real methodological distinction,
not just an error to note and move on from**: 11b measures
autocorrelation decay of the continuous Hs VALUE (how long until
Hs(t) and Hs(t+k) stop being related); HMM dwell time measures how
long Hs stays inside one DISCRETIZED bin. A smoothly, continuously
evolving storm (Hs climbing steadily over ~2 days, then declining)
keeps the raw signal highly autocorrelated throughout - driving 11b's
long timescale - while that same smooth climb crosses several regime
boundaries along the way, registering as multiple discrete HMM
transitions despite nothing discontinuous happening physically. These
are genuinely different quantities that happen to both get called
"persistence" loosely; there was no real basis for expecting them to
land in the same numeric range, and stating that expectation before
thinking it through (rather than after) was the actual mistake, not
just getting a number wrong.

**The U-shape itself looks like a real, physically coherent finding,
not noise**: calm and storm are the two "sticky" states (26.1h,
21.2h); moderate and energetic are faster-transiting states Hs passes
through on the way between them (11.8h, 11.5h) - consistent with
weather systems having persistent baseline/peak states connected by
faster ramp-up/ramp-down phases. ~21h for the PEAK portion of a storm
(not the whole event) is physically plausible for this coast.

**What held up correctly**: every non-adjacent "skip" transition came
out at exactly 0.0000 - confirming the tridiagonal structure (Hs only
ever moves through adjacent regimes, never jumps directly from calm to
storm) is a genuine real property, not an artifact specific to how the
synthetic validation test was constructed.

**Not yet built**: Bayesian change-point detection (next, given its
direct connection to the Westhinder drift question) and SSA (hardest,
only worth the 1-2 session investment if there's real appetite to
revisit Zeebrugge specifically).

## Post-priority extension: Bayesian change-point detection (Stage 25) — BUILT and tested

Directly targets the Westhinder drift question (README/METHODS.md
Section 6, provisionally resolved via Mann-Kendall finding no trend) -
but Mann-Kendall tests for a smooth monotonic trend, a different
hypothesis than a discrete step-change (mooring relocation, sensor
upgrade, processing-method change at a specific era boundary). A real
step-change could exist even when Mann-Kendall correctly reports "no
trend," since these are different questions, not the same question at
different sensitivity.

**Deliberate scope decision, not a shortcut**: uses PELT (via
`ruptures`) with an L2 mean-shift cost model, not full MCMC-based
Bayesian inference. PELT answers "does a change point exist, and
roughly where" - the actual question here - without the added
complexity of prior specification and convergence diagnostics a fully
Bayesian version would need for the same practical answer, per the
scoping decision made when this was first estimated.

Penalty (controls how many change points get detected - higher =
fewer, more conservative) is a real, acknowledged design choice, not
hidden behind one arbitrary number: reported across a sweep of 4
multipliers (0.5x/1x/2x/4x a BIC-style default) for both annual mean
and annual p95 Hs, so sensitivity to this choice is visible.

**Validated two ways, both informative**:
1. Known injected step-change (level 1.0 -> 1.35 at year 15 of a
   30-year synthetic record, true AR(1) noise on top): detected at
   EXACTLY the true year, consistently across all 4 penalty levels, for
   both mean and p95 series - no ambiguity, no near-miss.
2. No-change control (flat level, only noise): correctly zero change
   points at the default penalty and above. But at the most permissive
   setting (0.5x), a SPURIOUS change point appeared at the record's
   midpoint - a genuine false positive that the sweep design is
   specifically meant to catch. This is a real demonstration that the
   sweep isn't superfluous: if only the permissive setting had been
   run and reported, this false positive would have looked like a
   finding; because it only appears at the least conservative penalty
   (not persisting across the sweep), the tool's own printed guidance
   correctly flags it as the less-robust kind of result.

**Run it on real Westhinder data**:
```bash
python 25_changepoint_detection.py --buoy WesthinderBuoy --var VHM0
```
The real test: does a change point appear anywhere in the annual mean
or p95 series, and if so, does it persist across most/all penalty
levels (robust) or only the most permissive one (likely noise, per the
validated false-positive behavior above)? If nothing robust is found,
that's a real, additional piece of evidence the Westhinder drift
question is genuinely resolved - a second, differently-motivated method
agreeing with Mann-Kendall's "no trend" conclusion, rather than the
same test run twice. If something robust IS found, that would be a
genuinely new finding Mann-Kendall structurally could not have
surfaced, worth investigating against real deployment history metadata
(the external-review suggestion from earlier - checking Westhinder's
actual mooring/sensor history - would be the natural next step if so).
