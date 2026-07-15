"""
Stage D - ARMA-GARCH: point forecast (same ARMA mean as Stage C) PLUS a
calibrated prediction interval, justified directly by Stage 07's
universal ARCH finding (not a generic upgrade - every buoy in this
network showed significant volatility clustering).

Two-step approach, not jointly estimated: fit ARMA on the mean (same
order-search + fit as Stage C), extract its in-sample residuals, fit
GARCH(1,1) on those residuals separately. Simpler and more transparent
than joint estimation, and the `arch` package doesn't support general
ARMA (with MA terms) as its mean model anyway - only pure AR.

Performance: same lesson as Stage C, applied again - fit GARCH ONCE,
not at every origin. At each origin, the GARCH state (current
conditional variance) is updated via the closed-form recursion over a
BOUNDED window of ARMA residuals using the FIXED fitted (omega, alpha,
beta) - no re-fitting. The ARMA impulse-response (MA(infinity))
weights, needed to correctly turn GARCH's per-step innovation-variance
forecasts into actual multi-step forecast-error variance, depend only
on the fixed ARMA coefficients - computed ONCE before the backtest
loop, not per origin.

IMPORTANT, found via testing (not assumed correct on the first try):
GARCH's per-step variance forecast is the variance of the innovation AT
step h, not the variance of the ACCUMULATED forecast error by step h -
these are the same only for h=1. For any series with real persistence
(phi close to 1, exactly this pipeline's situation per 11b's confirmed
long timescales), naively using the per-step variance as the interval
width badly underestimates uncertainty at longer horizons, since GARCH
variance converges to a flat unconditional value quickly while the true
multi-step forecast error keeps growing as long as AR memory hasn't
decayed. Confirmed on a synthetic AR(1)+GARCH(1,1) test (phi=0.98):
naive implementation gave PI coverage collapsing from 0.845 (1h) to
0.358 (24h) against a 95% target - badly overconfident, worsening with
horizon. Fixed by properly accumulating GARCH's per-step variances
through the ARMA's impulse-response weights (see
accumulate_forecast_error_variance) - this is standard practice for
ARMA-GARCH multi-step forecasting, just easy to get wrong by treating
GARCH's output as already being the answer.

Evaluation is different from Stage C's, not just point RMSE: reports
(1) point RMSE - should closely match Stage C's, since the mean model
is unchanged; a sanity check, not the interesting number here; (2) 95%
prediction interval coverage - the fraction of actual values that fall
within the GARCH-implied interval; should be close to 95% if the
model's uncertainty is honestly calibrated, not just accurate on
average; (3) mean interval width per horizon - narrower is more
informative, provided coverage holds.

Requires the `arch` package (pip install arch) - not the same as
statsmodels' limited built-in ARCH support.

Usage:
    python 19_forecast_arma_garch.py --buoy WesthinderBuoy --var VHM0
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.arima.model import ARIMA

from utils import default_paths, longest_contiguous_segment
from forecast_utils import select_order

try:
    from arch import arch_model
except ImportError:
    arch_model = None


def garch_state_and_forecast(resid_window: np.ndarray, omega: float, alpha: float,
                               beta: float, max_horizon: int):
    """Update GARCH(1,1) conditional variance over a bounded window of
    residuals (fixed params, no re-fitting), then produce closed-form
    multi-step variance forecasts for the INNOVATION at each future step
    (NOT yet the forecast error variance - see
    accumulate_forecast_error_variance for why these are different).
    Cross-validated against arch's own .forecast() before trusting this -
    see module docstring.

    Handles the near-unit-root (IGARCH, alpha+beta ~ 1) case explicitly,
    found necessary on real Westhinder data: alpha+beta fit to
    1.0000000000273 - 1-alpha-beta = -2.7e-11, i.e. floating-point noise
    around exactly 1.0, not a meaningfully explosive process. The
    standard closed-form formula (which needs uncond_var =
    omega/(1-alpha-beta)) technically divides by ~zero here. It turned
    out to still give correct results THROUGH ALGEBRAIC CANCELLATION -
    the (alpha+beta)^(h-1) factor multiplying uncond_var also -> 1,
    which cancels the degenerate uncond_var out of the formula - but
    relying on that cancellation is fragile (it depends on how many
    horizons are requested and floating-point rounding), not a
    deliberate design. Two fixes: (1) sigma2's own starting point now
    uses the window's sample variance (well-defined regardless of
    stationarity) instead of uncond_var, which itself was previously
    also relying on beta decaying it away within the window - true here
    (beta=0.83, window=400 samples, 0.83^400 ~ 3e-33) but not guaranteed
    for a higher beta; (2) explicitly detect |alpha+beta-1| below a
    threshold and use the exact IGARCH limit (flat variance = sigma2 at
    every horizon, since shocks are permanent) instead of computing
    uncond_var at all.
    """
    is_near_unit_root = abs(1 - alpha - beta) < 1e-6

    sigma2 = np.var(resid_window) if len(resid_window) > 1 else omega
    for e in resid_window:
        sigma2 = omega + alpha * e ** 2 + beta * sigma2

    if is_near_unit_root:
        # CORRECTED - a real bug in the previous version of this fix, found
        # by comparing against real Westhinder results: "shocks don't
        # decay" does NOT mean "forecast stays flat" - it means the
        # opposite. The general recursion is f(h) = omega + (alpha+beta)*
        # f(h-1); at alpha+beta=1 exactly this becomes f(h) = omega +
        # f(h-1), i.e. LINEAR growth (f(h) = f(1) + (h-1)*omega), not a
        # constant. A flat forecast understates the genuinely growing
        # uncertainty when variance shocks are permanent - confirmed on
        # real data: the flat version gave coverage collapsing to 0.803 at
        # 24h (undercovering, interval too narrow), while the ORIGINAL
        # cancellation-based formula (before either "fix") had actually
        # been correct all along at 0.968 - it was computing the right
        # mathematical limit through legitimate cancellation, not luck.
        # This linear formula should closely reproduce those original
        # numbers, via an explicit route that doesn't depend on that
        # cancellation - verify this before trusting it again.
        return [sigma2 + h * omega for h in range(max_horizon)]

    uncond_var = omega / (1 - alpha - beta)
    variances = [sigma2]
    for h in range(2, max_horizon + 1):
        variances.append(uncond_var + (alpha + beta) ** (h - 1) * (sigma2 - uncond_var))
    return variances  # index k-1 = variance of the innovation k steps ahead


def accumulate_forecast_error_variance(step_variances: list, psi: np.ndarray, max_horizon: int):
    """The h-step-ahead ARMA forecast error is a weighted sum of ALL
    innovations between now and h (weighted by the model's impulse-
    response/MA(infinity) coefficients psi), not just the innovation at
    step h alone: error(h) = sum_{j=0}^{h-1} psi[j] * eps(t+h-j), so
    Var(error(h)) = sum_{j=0}^{h-1} psi[j]^2 * Var(eps(t+h-j)).

    Found necessary after a real bug: using step_variances[h-1] directly
    as the forecast-error variance (what a naive implementation does)
    badly underestimates uncertainty for any series with real
    persistence (phi close to 1) - GARCH's per-step variance forecast
    converges quickly to the unconditional variance and then stays flat,
    while true multi-step forecast error keeps growing as long as the
    AR dynamics haven't fully decayed. Confirmed empirically on a
    synthetic AR(1)+GARCH(1,1) test with phi=0.98: interval width was
    nearly flat (0.299 at 1h vs 0.303 at 24h) while empirical RMSE grew
    3x (0.111 to 0.360) over the same range - PI coverage collapsed from
    0.845 to 0.358, badly overconfident and getting worse with horizon.
    """
    psi = np.asarray(psi)
    acc = np.zeros(max_horizon)
    for h in range(1, max_horizon + 1):
        total = 0.0
        for j in range(h):
            k = h - j
            total += psi[j] ** 2 * step_variances[k - 1]
        acc[h - 1] = total
    return acc


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--buoy", default="WesthinderBuoy")
    parser.add_argument("--var", default="VHM0")
    parser.add_argument("--horizons-hours", default="1,3,6,12,24")
    parser.add_argument("--origin-step-hours", default=6.0, type=float)
    parser.add_argument("--min-history-hours", default=24.0, type=float)
    parser.add_argument("--order-search-samples", default=4000, type=int)
    parser.add_argument("--fit-samples", default=8000, type=int)
    parser.add_argument("--state-window-hours", default=200.0, type=float)
    parser.add_argument("--d-range", default="0", help="see 18_forecast_arma.py - "
                         "defaults to d=0 only, same reasoning")
    args = parser.parse_args()

    if arch_model is None:
        print("The 'arch' package is not installed - required for this stage. "
              "Run: pip install arch")
        return

    clean_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_clean.csv"
    load_summary_path = Path("pipeline_out/01_load_clean") / f"{args.buoy}_{args.var}_load_summary.json"

    hs_full = pd.read_csv(clean_path, index_col=0, parse_dates=True)[args.var]
    dt_hours = 0.5
    if load_summary_path.exists():
        with open(load_summary_path) as f:
            dt_hours = json.load(f).get("sampling_interval_hours", 0.5)

    print(f"--- {args.buoy} / {args.var} ARMA-GARCH forecasting ---")

    segment, seg_meta = longest_contiguous_segment(hs_full)
    print(f"Longest contiguous segment: {len(segment)} samples "
          f"({seg_meta['pct_of_valid_used']}% of valid data)")
    if len(segment) < args.fit_samples + 1000:
        print(f"Segment too short (need > {args.fit_samples + 1000}, "
              f"have {len(segment)}) - stopping.")
        return

    # --- Fit ARMA mean model, same as Stage C ---
    d_range = [int(d) for d in args.d_range.split(",")]
    order_train = segment.iloc[:args.order_search_samples]
    print(f"\nSearching ARMA order on first {len(order_train)} samples "
          f"(p,q in 0..3, d in {d_range})...")
    order, aic = select_order(order_train, p_range=range(0, 4), q_range=range(0, 4), d_range=d_range)
    if order is None:
        print("Order search failed - stopping.")
        return
    print(f"Selected order: {order} (AIC={aic:.1f})")

    fit_train = segment.iloc[:args.fit_samples]
    print(f"Fitting ARMA on first {len(fit_train)} samples...")
    arma_fitted = ARIMA(fit_train.values, order=order).fit()

    # --- Fit GARCH(1,1) on ARMA's in-sample residuals, ONCE ---
    arma_resid = arma_fitted.resid
    print(f"Fitting GARCH(1,1) on {len(arma_resid)} ARMA residuals...")
    garch_fitted = arch_model(arma_resid, mean="Zero", vol="GARCH", p=1, q=1).fit(disp="off")
    omega = garch_fitted.params["omega"]
    alpha = garch_fitted.params["alpha[1]"]
    beta = garch_fitted.params["beta[1]"]
    print(f"GARCH params: omega={omega:.5f}, alpha={alpha:.4f}, beta={beta:.4f} "
          f"(persistence alpha+beta={alpha + beta:.4f})")
    # Full precision here deliberately - 4 decimals can round e.g. 0.99997
    # or 1.00003 to the same displayed "1.0000", but those are qualitatively
    # different cases for uncond_var = omega/(1-alpha-beta): the first gives
    # a large-but-finite positive number, the second gives a NEGATIVE
    # number (which the interval-width code silently clips to 0 via
    # max(var_h, 0), degrading calibration without any visible sign of it).
    # Don't trust the rounded display to tell you which side of the
    # boundary you're actually on.
    denom = 1 - alpha - beta
    print(f"  Full precision: alpha+beta={alpha + beta!r}, "
          f"1-alpha-beta={denom!r}")
    if abs(denom) < 1e-6:
        print(f"  Near-unit-root (IGARCH-like): |1-alpha-beta| < 1e-6 - "
              f"handled explicitly (LINEAR variance growth with horizon, "
              f"+omega per step - the correct IGARCH limit; an earlier "
              f"version of this code incorrectly used a flat forecast here, "
              f"which undercovered badly at long horizons - fixed), not via "
              f"the standard closed-form formula, which would divide by "
              f"~zero here.")
    elif denom > 0:
        uncond_var = omega / denom
        print(f"  Unconditional variance = {uncond_var:.6f} "
              f"(finite and positive - standard closed-form formula applies)")
    else:
        print(f"  WARNING: 1-alpha-beta is meaningfully negative (not just "
              f"floating-point noise near zero) - genuinely explosive "
              f"variance process. Results are likely unreliable; this case "
              f"isn't handled by the IGARCH fallback either, since that's "
              f"only for the near-exactly-1 case.")
    if alpha + beta >= 0.98:
        print(f"  NOTE: persistence >= 0.98 is itself a real, substantive "
              f"finding, not just a numerical edge case - variance shocks "
              f"barely decay within any practically forecastable horizon. "
              f"Consistent with (and stronger than) Stage 07's universal "
              f"ARCH finding and 11b's long persistence timescales for this "
              f"buoy - real storms leave a very long footprint in volatility.")

    # --- Walk-forward: point forecast (same as Stage C) + GARCH interval ---
    horizons_hours = [float(h) for h in args.horizons_hours.split(",")]
    horizons_samples = [max(1, round(h / dt_hours)) for h in horizons_hours]
    max_h = max(horizons_samples)
    origin_step = max(1, round(args.origin_step_hours / dt_hours))
    min_history = max(args.fit_samples, round(args.min_history_hours / dt_hours))
    state_window = max(1, round(args.state_window_hours / dt_hours))
    freq_str = f"{int(round(dt_hours * 60))}min"

    print(f"\nBacktesting (origin step {args.origin_step_hours}h, "
          f"state window {args.state_window_hours}h)...")

    # ARMA impulse-response weights - depend only on the fixed fitted
    # coefficients, computed ONCE, reused at every origin (see module
    # docstring for why this matters, not just an optimization).
    psi = arma_fitted.impulse_responses(steps=max_h)
    print(f"Impulse response weights (first few): {np.round(psi[:5], 3)} "
          f"(psi[0]=1 by construction, decay rate reflects the ARMA's own memory)")

    n = len(segment)
    values = segment.values
    rows = []
    for origin in range(min_history, n, origin_step):
        window = segment.iloc[max(0, origin + 1 - state_window):origin + 1].dropna()
        if len(window) < 5:
            continue
        if window.index.freq is None:
            window = window.asfreq(freq_str)
        try:
            applied = arma_fitted.apply(window, refit=False)
            point_fc = applied.forecast(steps=max_h)
            resid_window = applied.resid
        except Exception:
            continue

        step_variances = garch_state_and_forecast(resid_window, omega, alpha, beta, max_h)
        forecast_error_variances = accumulate_forecast_error_variance(step_variances, psi, max_h)

        for h in horizons_samples:
            target_idx = origin + h
            if target_idx >= n:
                continue
            actual = values[target_idx]
            if np.isnan(actual):
                continue
            mean_h = float(point_fc.iloc[h - 1])
            var_h = forecast_error_variances[h - 1]
            half_width = 1.96 * np.sqrt(max(var_h, 0))
            covered = (actual >= mean_h - half_width) and (actual <= mean_h + half_width)
            rows.append({
                "origin_idx": origin, "horizon_samples": h,
                "actual": actual, "predicted": mean_h,
                "error": actual - mean_h, "half_width": half_width,
                "covered": covered,
            })

    results = pd.DataFrame(rows)
    if len(results) == 0:
        print("No valid forecasts produced - stopping.")
        return

    summary_rows = []
    for h, group in results.groupby("horizon_samples"):
        rmse = float(np.sqrt(np.mean(group["error"] ** 2)))
        coverage = float(group["covered"].mean())
        mean_width = float((2 * group["half_width"]).mean())
        summary_rows.append({
            "horizon_samples": h, "horizon_hours": h * dt_hours, "n": len(group),
            "rmse": rmse, "pi_coverage_95": coverage, "mean_pi_width": mean_width,
        })
    summary = pd.DataFrame(summary_rows).sort_values("horizon_samples")

    print(f"\nUsed {results['origin_idx'].nunique()} distinct origins, "
          f"{len(results)} scored pairs.")
    print("\nARMA-GARCH per-horizon result:")
    print(summary.to_string(index=False))
    print("\n(pi_coverage_95 should be close to 0.95 if well-calibrated - "
          "notably lower means the interval is too narrow/overconfident, "
          "notably higher means it's too wide/conservative)")

    out_dir = default_paths("19_forecast_arma_garch")
    summary.to_csv(out_dir / f"{args.buoy}_{args.var}_arma_garch_summary.csv", index=False)
    with open(out_dir / f"{args.buoy}_{args.var}_arma_garch_meta.json", "w") as f:
        json.dump({
            "arma_order": list(order), "garch_omega": float(omega),
            "garch_alpha": float(alpha), "garch_beta": float(beta),
            "garch_persistence": float(alpha + beta),
            "n_origins": int(results["origin_idx"].nunique()),
        }, f, indent=2)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    axes[0].plot(summary["horizon_hours"], summary["pi_coverage_95"], marker="o")
    axes[0].axhline(0.95, color="firebrick", ls="--", label="nominal 95%")
    axes[0].set_xlabel("horizon (hours)")
    axes[0].set_ylabel("empirical PI coverage")
    axes[0].set_title(f"{args.buoy} — 95% interval coverage")
    axes[0].legend()
    axes[0].grid(alpha=0.3)

    axes[1].plot(summary["horizon_hours"], summary["mean_pi_width"], marker="s", color="darkorange")
    axes[1].set_xlabel("horizon (hours)")
    axes[1].set_ylabel(f"mean 95% PI width ({args.var}, m)")
    axes[1].set_title(f"{args.buoy} — interval width vs. horizon")
    axes[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_dir / f"{args.buoy}_{args.var}_arma_garch_coverage.png", dpi=150)
    plt.close(fig)

    print(f"\nSaved: {out_dir}")


if __name__ == "__main__":
    main()
