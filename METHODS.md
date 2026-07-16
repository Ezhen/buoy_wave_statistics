# Methods

## 1. Study area and data

### 1.1 Wave observation network

Statistical characterization was performed for 19 in-situ wave buoys operated within the Belgian Coastal Zone (BCZ), obtained from the Copernicus Marine Environment Monitoring Service (CMEMS) North West Shelf regional in-situ product (`INSITU_NWS_PHYBGCWAV_DISCRETE_MYNRT_013_036`). Significant wave height (Hs, variable `VHM0`) was the primary analysis variable at all 19 buoys. Two buoys (CadzandBoei, Deurlo) additionally report peak period (`VTPK`) and mean wave direction (`VMDR`) at finer temporal resolution (10–15 min) than the remaining 17 buoys (30 min); this variable and sampling-rate heterogeneity is consistent with a distinct sensor class at these two sites and was confirmed independently via statistical-fingerprint clustering (Section 3.13).

Two acquisition pathways were used. Near-real-time (NRT) data were obtained via the `copernicusmarine` `subset()` service, which for in-situ products only returns the rolling ~30-day `latest` dataset part regardless of the requested date range. Full multi-year archives were obtained via the Files service (`copernicusmarine.get()`, `dataset_part="history"`), which serves one NetCDF file per platform under the naming convention `history/<type>/<network>_TS_<platform_type>_<name>.nc`. Per-buoy record length under this multi-year acquisition is markedly uneven: six buoys (Westhinder, TrapegeerBuoy, BolVanHeistBuoy, ScheurWielingenBuoy, WandelaarBuoy, OstendEasternPalisadeBuoy) have continuous deployment history extending back to 1990–1997; the remaining 13 buoys were deployed between 2009 and 2021. All records extend to the most recent download date (2026-06-30).

### 1.2 Atmospheric reanalysis

10 m wind components (u10, v10), mean sea-level pressure (MSLP), 2 m air temperature, and sea-surface temperature were obtained from ERA5 (Copernicus Climate Data Store, `cdsapi`) at 3-hourly resolution for 2010–2026, covering the geographic extent of the buoy network. For each buoy, the nearest ERA5 grid cell (0.25° resolution) was extracted; wind speed and wind direction (meteorological "from" convention) were derived from the u/v components.

### 1.3 Data regularization and quality control

Raw buoy time series were regularized onto a uniform temporal grid at each buoy's native sampling interval. Gaps of short duration (≤3 consecutive missing samples) were linearly interpolated; longer gaps were retained as missing (NaN) rather than bridged, since interpolation across a gap of unknown duration would fabricate structure not present in the observations. Sanity bounds were applied to flag physically implausible values. Duplicate timestamps in the raw record were identified and resolved prior to regularization.

Given the multi-year records are substantially fragmented by real sensor gaps (e.g. Westhinder: 1905 contiguous segments after gap removal, longest single segment covering 6.9% of valid data), a general principle was applied throughout: any statistic computed via a lag- or order-based method (autocorrelation, differencing, block bootstrap) that is sensitive to genuine temporal adjacency between consecutive samples was restricted to contiguous segments, aggregated across all sufficiently long segments (weighted by segment length) rather than computed naively on the full record. A naive concatenation of pre-gap and post-gap data — the default behaviour of `dropna()` on a series with internal NaN gaps — silently treats temporally distant periods as adjacent and corrupts any statistic that sums or differences across lags. Distribution- and threshold-based statistics that do not depend on sample ordering (e.g. distribution fitting, peak detection for extreme value analysis) were computed on the full valid record.

## 2. Periodicity removal

### 2.1 Tidal signal identification and removal

Periodogram analysis of the raw Hs series at every buoy revealed a dominant spectral peak near the M2 tidal period (12.4206 h), with peak-to-baseline power ratios ranging from 6× to over 1000× depending on buoy and record segment, reflecting genuine spatial variation in tidal influence across the network. The M2 constituent was removed via least-squares harmonic regression (fundamental frequency plus up to 2 harmonics) fitted to the Box–Cox-transformed level series. The regression design matrix is evaluated on the full regularized time grid (a deterministic function of elapsed time), while the fit itself uses only valid (non-missing) observations; gap positions are preserved as missing in the detided output rather than being collapsed, consistent with the general gap-handling principle in Section 1.3.

At one site (Zeebrugge), the notch consistently failed to fully suppress the M2 signature regardless of harmonic count. A frequency-fitting variant was tested in which the dominant period was estimated directly from the periodogram rather than assumed to be exactly 12.4206 h; the fitted period (12.4163 h) was found to be effectively identical to the nominal M2 frequency, ruling out constituent misidentification (e.g. S2, 12.0000 h) as the explanation. Application of the correctly-identified frequency did not improve — and marginally worsened — the notch performance, indicating the residual contamination is not attributable to a fundamental-frequency error. This buoy is treated as a structurally distinct site in subsequent analysis (Section 6).

### 2.2 Diurnal signal

A secondary, network-wide periodogram peak near 24 h was identified in Hs at every buoy tested (power ratios 85×–289× baseline). Cross-referencing against ERA5 wind speed at the same locations showed a comparable diurnal signature in wind itself (75×–1986× baseline), consistent with a land–sea breeze mechanism rather than an unmodelled astronomical tidal constituent. A joint M2 + diurnal harmonic notch was evaluated but gave only partial, buoy-inconsistent improvement (5–80% power reduction across test buoys, never reaching the same clean-notch threshold achieved for M2 alone) and was not adopted, consistent with land–sea breeze being a stochastic, weather-conditional process poorly represented by a single fixed-phase sinusoid over a multi-year record.

## 3. Statistical characterization

### 3.1 Stationarity

The Augmented Dickey–Fuller (ADF) and Kwiatkowski–Phillips–Schmidt–Shin (KPSS) tests were applied to the raw Hs series at each buoy. The two tests frequently disagree, which is expected given their opposing null hypotheses in the presence of strong periodic structure and long-range persistence rather than a genuine trend, and is not itself diagnostic of a data-quality problem.

### 3.2 Distribution fitting

Rayleigh, Weibull, and log-normal distributions were fitted to the raw (untransformed) Hs level series via maximum likelihood, and compared via the Kolmogorov–Smirnov (KS) statistic. Two independent limitations of the KS p-value in this context are noted: the classical test assumes i.i.d. samples, which is violated given the confirmed autocorrelation in Hs (Section 3.4); and at the sample sizes involved (up to ~10⁶), the KS p-value floors near zero regardless of fit quality (a large-*n* hypersensitivity effect). Point estimates of shape and scale are not affected by either limitation; p-values are reported for completeness but not used to reject candidate distributions.

Weibull provided the best fit at the majority of buoys under this criterion; log-normal was preferred at a minority, including Westhinder. A moving-window stability check (Section 3.7) established that a log-normal "win" on the pooled full record can arise from pooling multiple genuinely-Weibull sub-periods with slowly varying parameters, rather than reflecting a true log-normal generating process; this was confirmed directly at Westhinder, where every individual calendar-era window independently preferred Weibull despite the pooled fit preferring log-normal.

### 3.3 Volatility clustering

The Engle ARCH-LM test was applied to the residual series (Section 3.5) at each buoy. Statistically significant ARCH effects were found at every buoy in the network, indicating universal volatility clustering in Hs consistent with storm-driven variance dynamics.

### 3.4 Dependence structure and persistence timescale

The integral (decorrelation) timescale was estimated from the sample autocorrelation function (ACF) of the detided series, using a significance-band stopping criterion (five consecutive lags within ±1.96/√n) rather than a single zero-crossing, which is sensitive to spurious early crossings from sampling noise. For heavily fragmented multi-year records, the timescale was estimated independently on every contiguous segment exceeding a minimum length and combined as a length-weighted mean across segments, rather than computed on the single longest segment alone (which for the most fragmented records covers under 10% of valid data). At Westhinder, this yields a persistence timescale of 115.7 h (range 26–253 h across segments used), with a genuine open question as to whether this variation reflects estimation noise in shorter segments or real multi-decadal drift in persistence.

This timescale directly informs two downstream analyses: the block length for block-bootstrap resampling (Section 3.6), and the minimum peak-separation window for extreme value declustering (Section 3.5).

### 3.5 Extreme value analysis

Storm peaks were identified via peaks-over-threshold (POT) at the 95th percentile of Hs, declustered using a minimum separation window derived from twice the per-buoy integral timescale (Section 3.4) rather than a fixed round-number default; this connection was found to matter materially — at Westhinder, increasing the declustering window from a naive 48 h default to the persistence-justified 231 h reduced the independent peak count from 1031 to 483 and shifted the fitted shape parameter (ξ) from −0.23 to −0.36, while leaving return-level estimates at all tested horizons changing by under 0.15 m even at 50-year return periods — indicating the practically material output (return levels) is robust to this declustering choice even though the underlying parameters are not. A generalized Pareto distribution (GPD) was fitted to the exceedances above threshold via maximum likelihood.

Network-wide, GPD shape parameters ranged from −0.544 to +0.327, with the majority of well-instrumented buoys (many storm peaks, narrow bootstrap confidence intervals) clustering between −0.27 and −0.32 (bounded upper tail, consistent with a fetch- and/or depth-limited shelf sea). This range brackets a published shallow-water North Sea reference (Caires, 2011, JCOMM Technical Report No. 57: ξ ≈ −0.12 to −0.13 at a comparable 19 m-depth site), with the BCZ network running systematically more negative — plausibly reflecting shallower average depth across the BCZ relative to the reference site. One buoy (BlankenbergeBuoy) returned a positive shape parameter (ξ = +0.327), physically implausible for a depth-limited system; this buoy has both the shortest record (2.8 years) and fewest storm peaks (43) in the network, and its bootstrap confidence interval crosses zero, consistent with a small-sample estimation artifact rather than a genuine unbounded tail.

### 3.6 Uncertainty quantification

Confidence intervals for the Hs mean, selected quantiles, and the Weibull fit were computed via moving block bootstrap, with block length set from the persistence timescale (Section 3.4) rather than a generic heuristic. Confidence intervals for the GPD shape parameter were computed via bootstrap resampling of the declustered peaks. Pairwise spatial correlations (Section 3.8) were assigned confidence intervals via the Fisher *z* transform with an effective-sample-size correction (Dawdy–Matalas, using the lag-1 autocorrelation from Section 3.4) to account for the substantial loss of independent information under temporal autocorrelation; across the network this correction reduces the effective sample size by approximately 99% relative to the raw sample count, confirming that an uncorrected correlation confidence interval would be materially overconfident.

For records exceeding approximately 500,000 samples, the standard maximum-likelihood Weibull refit inside the bootstrap loop becomes computationally prohibitive (approximately 2 s per refit; 30–40 minutes for 1000 resamples). This was replaced, for the bootstrap loop specifically, with a closed-form method-of-moments estimator based on the coefficient-of-variation relation for the Weibull distribution, validated against full maximum likelihood on independent test data (0.02% relative parameter error, ~630× speed-up); the primary point estimate reported in Section 3.2 continues to use full maximum likelihood.

### 3.7 Stability analysis

Three complementary stability checks were performed. First, the distribution fit (Section 3.2) was recomputed within four equal-duration calendar-time windows per buoy and compared against the pooled full-record fit, to distinguish genuine record-wide distributional properties from artifacts of non-stationary pooling (see Section 3.2 for the Westhinder finding). Second, a jackknife sensitivity check recomputed the Weibull and GPD fits after removing the single largest storm from the record, to assess whether headline results are unduly influenced by one event; at Westhinder the resulting parameter shifts were small relative to the independently-estimated bootstrap confidence interval width. Third, regime-fraction confidence intervals (Section 3.9) were computed via block bootstrap restricted to contiguous segments detected by elapsed-time gaps in the regime-label sequence (which, unlike the raw Hs series, has gap positions entirely absent rather than NaN-marked, requiring a different segment-detection method).

### 3.8 Spatial structure

Pairwise Pearson correlation of Hs between all buoy pairs was computed using pairwise-complete observations (i.e., each pair's correlation uses only their genuinely overlapping time period, correctly handling the substantial disparity in deployment history across the network without requiring an explicit overlap-window calculation). Correlation was found to decrease with inter-buoy distance (Spearman ρ ≈ −0.47 to −0.51, reproduced independently on both the short NRT window and the full multi-year record). Hierarchical clustering on the correlation matrix identified four geographically coherent groups, with one buoy (Zeebrugge) forming a persistent singleton cluster across both data windows.

### 3.9 Regime identification

A four-component Gaussian mixture model was fitted to Hs (and, where available, peak period) to identify calm, moderate, energetic, and storm regimes. Regime time-fractions and their bootstrap confidence intervals (Section 3.6) were computed per buoy.

### 3.10 Trend and seasonal analysis

Restricted to the six buoys with 30+ year continuous records. A Mann–Kendall trend test with the Hamed–Rao variance correction for autocorrelation (rather than the uncorrected textbook test, which the confirmed network-wide serial correlation would otherwise inflate the false-positive rate of) was applied to annual mean and annual 95th-percentile Hs. Five of six buoys showed no significant trend in either series; one buoy showed a marginal increase in storm-intensity (p95) Hs that does not survive a Bonferroni correction for the twelve tests performed across the six buoys, and is treated as a lead for future investigation rather than an established finding.

Seasonal STL (Seasonal-Trend decomposition using LOESS) decomposition was applied to monthly mean and p95 Hs. A consistent winter-storm seasonal cycle was found at all six buoys (peak month November–January at every buoy, amplitude 0.49–1.23 m for the mean and 1.3–1.8 m for p95), with storm-intensity seasonal amplitude consistently 2–3× the mean-Hs seasonal amplitude — indicating winter conditions specifically inflate the extreme end of the distribution rather than uniformly raising typical conditions.

### 3.11 Statistical fingerprinting and network-wide similarity

A per-buoy feature vector (distribution type and goodness-of-fit, GPD shape parameter, storm-peak count, persistence timescale, confidence-interval widths, distributional stability fraction, regime time-fractions) was constructed from the outputs of Sections 3.2–3.9, standardized, and clustered via *k*-means with the number of clusters selected by silhouette score (rather than fixed a priori). This provides a statistical-similarity clustering distinct from the correlation-based spatial clustering of Section 3.8 — two buoys can be geographically uncorrelated but statistically similar, or vice versa.

On the full 19-buoy network, this identified CadzandBoei and Deurlo as a distinct pair separate from the remaining 17 buoys, consistent with their distinct sensor class (Section 1.1). Two additional apparent singleton clusters (Zeebrugge, BlankenbergeBuoy) were found to be substantially attributable to record-length and missing-data-fraction differences rather than core behavioural properties: re-clustering with record-length and completeness metadata excluded from the feature set dissolved both singleton clusters, while the CadzandBoei–Deurlo pairing not only persisted but produced a higher silhouette score, indicating it reflects genuine behavioural similarity rather than shared data characteristics. Zeebrugge's overall case for structural distinctiveness (Section 6) rests on the independent lines of evidence in Sections 2.1, 3.8, 3.7, and 4, none of which depend on this feature set.

## 4. Wind–wave coupling

ERA5 wind speed and MSLP (Section 1.2) were cross-correlated against Hs at each buoy. Across 18 of 19 buoys, wind speed leads Hs by 0–3 h (R² up to 0.71–0.77 at Westhinder), and MSLP leads with the physically expected sign (pressure decrease preceding Hs increase). One buoy (Zeebrugge) shows collapsed wind coupling (R² = 0.09) and an inverted MSLP lead–lag relationship, a further independent indicator of structural distinctiveness at this site (Section 6).

## 5. Forecasting

A rolling-origin backtest harness was implemented once and reused across all forecasting models: at each origin in a held-out test period, a forecast is generated from data available at or before that origin only, and scored against the subsequently observed value at each of a set of forecast horizons (1, 3, 6, 12, 24 h for point forecasts; 6, 12, 24 h windows for exceedance forecasting). All models are fitted once on a training window and applied with fixed parameters across the backtest, updating only the model's internal state (not its coefficients) at each origin, for computational tractability at the sample sizes involved.

### 5.1 Persistence baseline

The naive forecast Hs(t+h) = Hs(t) (most recent valid observation, with a maximum staleness limit beyond which no forecast is issued) establishes the baseline every subsequent model is compared against.

### 5.2 ARMA

An ARIMA(p,d,q) model was fitted to the raw Hs series (not the detided series, since the M2 signal is itself deterministic and predictable at sub-daily forecast horizons, and excluding it would discard genuine short-horizon predictive information). Order was selected via AIC grid search on a training subset. The differencing order *d* is selected via the same grid search rather than assumed equal to the value used for the characterization pipeline's stationarity work (Section 3.1), following the identification that a fixed *d*=1 assumption, appropriate for the detided series that Section 3.1 concerns, degrades forecast skill when applied to the raw series used here. At Westhinder, ARMA forecasts beat the persistence baseline at every tested horizon, with skill increasing with horizon (0.06 at 1 h to 0.19 at 24 h).

### 5.3 ARMA–GARCH

A GARCH(1,1) model was fitted to the ARMA residuals to obtain a calibrated forecast interval alongside the point forecast, motivated directly by the confirmed universal ARCH effect (Section 3.3). Multi-step forecast-error variance was obtained by combining the per-step GARCH conditional-variance forecast with the ARMA model's own impulse-response weights (the forecast-error variance at horizon *h* is a weighted sum of all intervening innovation variances, not the single-step variance alone). At Westhinder, the fitted persistence parameter (α+β) was found to lie within numerical precision of the IGARCH boundary (α+β = 1.0000000000273), for which the standard closed-form multi-step formula is undefined; the correct limiting behaviour (linear variance growth with horizon, rather than a constant forecast) was derived and validated against direct numerical iteration of the underlying recursion. The resulting 95% prediction intervals achieve empirical coverage of 0.94–0.97 across all tested horizons.

### 5.4 ARMAX

Lagged ERA5 wind speed was added as an exogenous regressor to the ARMA mean equation, using the lag identified in Section 4. For genuine multi-step forecasting this requires addressing that future wind is no more knowable at forecast time than future Hs: exogenous input is treated as genuinely known only within the empirically-identified lead time, and is persisted forward (not fabricated) beyond it. At both buoys tested (Westhinder, A2Buoy), ARMAX improved on plain ARMA at every forecast horizon; the improvement did not decay with horizon at either buoy, which was traced to wind's own persistence timescale (57.3 h at Westhinder, independently estimated via the same method as Section 3.4 applied to deseasonalized wind speed) exceeding the longest tested forecast horizon.

### 5.5 Exceedance forecasting

Point forecasting was supplemented with a probabilistic reframing: the probability that Hs exceeds a given threshold (90th or 95th percentile of the buoy's own record) at any point within a forecast window (6, 12, or 24 h), via logistic regression on three features computed strictly from data at or before the forecast origin (current level, short-term trend, recent volatility — the latter motivated by the confirmed universal ARCH effect). Performance was assessed via Brier score (against a constant base-rate baseline) and ROC-AUC. At Westhinder, Brier skill scores of 0.42–0.68 and ROC-AUC of 0.83–0.97 were obtained across tested thresholds and horizons; a substantial fraction of this skill reflects the buoy's already-established strong persistence (Section 3.4) rather than new information. Addition of current wind speed as a fourth feature improved skill consistently at both buoys tested, without the horizon-decay predicted by a synthetic validation case, for the same reason identified in Section 5.4.

## 6. Known limitations

**Zeebrugge.** Five independent analyses (tidal-notch failure persisting across all tested harmonic and frequency configurations, Section 2.1; singleton spatial cluster reproduced on independent data windows, Section 3.8; unstable distribution-fit verdict across moving windows, Section 3.7; collapsed and sign-inverted wind–wave coupling, Section 4; harbor-interior, most-sheltered location in the network) indicate this buoy is dynamically distinct from the rest of the network, plausibly reflecting local wave conditions dominated by harbor geometry and vessel wake rather than open-water wind forcing. A statistical-fingerprint-based sixth line of evidence (Section 3.11) did not survive controlling for record-length/completeness confounds and is not counted toward this conclusion. Two mechanisms for the persistent tidal-notch failure remain untested: a compound shallow-water tidal constituent (e.g. MS4) not representable by a fixed-fundamental-frequency harmonic basis, and genuine multi-year drift in the effective tidal parameters.

**Westhinder persistence drift.** Two independent signals (wide per-segment persistence-timescale spread, 26–253 h across qualifying segments, Section 3.4; a moving-window mixture artifact in the distribution fit, Section 3.2) motivated investigation of possible multi-decadal drift at this site. The corrected Mann–Kendall test (Section 3.10) found no significant trend in either mean or storm-intensity Hs, and the moving-window mixture artifact is now attributed to uneven gap-coverage across the four windowing periods rather than genuine physical drift; this question is considered resolved.

**Record heterogeneity.** Deployment history is highly uneven across the network (1990–1997 start for 6 buoys, 2009–2021 for the remaining 13), and sensor configuration differs at two buoys (Section 1.1). Any cross-buoy comparison should account for the specific overlap period and sensor class involved rather than assuming uniform data characteristics across the network.

**Text-mined external validation.** GPD shape parameter ranges were compared against one published external source (Caires, 2011); confidence in this comparison would benefit from a broader literature survey beyond the single reference used here.
