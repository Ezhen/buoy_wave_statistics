"""
Shared rolling-origin backtest harness for the forecasting stages.

Built once here, reused by every model (persistence, ARMA, ARMA-GARCH,
ARMAX) rather than duplicating the walk-forward loop per model - the
harness itself doesn't know or care what forecast_fn does internally,
only that it takes (history, horizons_samples) and returns
{horizon_samples: predicted_value}.

Gap-aware in a way appropriate to POINT forecasting (different from the
lag-based gap handling elsewhere in this pipeline - a forecast only
needs a valid origin point and a valid target point, not a whole
contiguous stretch, so this doesn't need longest_contiguous_segment/
all_contiguous_segments the way ACF-based stages did):
  - An origin is only used if there's a genuinely RECENT valid sample to
    forecast from (within max_lookback_samples) - "persisting" a value
    from days ago across a gap isn't a meaningful forecast.
  - A (origin, horizon) pair is skipped if the target falls in a gap -
    can't score a forecast against a missing actual.
"""

import numpy as np
import pandas as pd


def persistence_forecast(history: pd.Series, horizons_samples: list, max_lookback_samples: int):
    """Naive Hs(t+h) = Hs(t) for every horizon. Looks backward from the
    end of history for the most recent valid sample, refusing to persist
    across a gap wider than max_lookback_samples."""
    valid = history.dropna()
    if len(valid) == 0:
        return None
    last_valid_pos = history.index.get_loc(valid.index[-1])
    lookback = len(history) - 1 - last_valid_pos
    if lookback > max_lookback_samples:
        return None  # last valid value too stale to persist from
    last_value = valid.iloc[-1]
    return {h: last_value for h in horizons_samples}


def rolling_origin_backtest(series: pd.Series, forecast_fn, horizons_samples: list,
                              origin_step: int, min_history: int,
                              max_lookback_samples: int, history_window: int = None,
                              **forecast_fn_kwargs):
    """Walk forward through the series, calling forecast_fn at each
    origin, scoring against the actual value at origin+horizon.

    Convention (verified with an analytical test before trusting it -
    this is exactly the kind of off-by-one that's easy to get wrong):
    `origin` is the index of the LAST OBSERVED sample - "now" - not the
    first unobserved one. So history = series.iloc[:origin+1] (inclusive
    of origin), and a horizon of h genuinely means h samples strictly
    after "now", i.e. target_idx = origin + h.

    `history_window`: if set, only the last `history_window` samples
    before/including origin are passed to forecast_fn, instead of the
    full history from the start of the series. Found necessary after a
    real performance bug: passing the full ever-growing history to a
    model that only needs a small trailing window (like persistence,
    which never looks further back than max_lookback_samples) turns
    each of ~50k+ origins into O(origin) work instead of O(1) -
    effectively O(n^2) over the whole backtest. Models that genuinely
    need more history (e.g. ARMA needing enough data to fit) should
    pass a larger history_window or leave it None, but should also
    consider whether they truly need to refit from scratch at every
    single origin - a separate performance question for those stages
    when built, not solved by this parameter alone.

    Returns a long-format DataFrame: one row per (origin, horizon)."""
    n = len(series)
    values = series.values
    rows = []

    for origin in range(min_history, n, origin_step):
        if history_window is not None:
            start = max(0, origin + 1 - history_window)
            history = series.iloc[start:origin + 1]
        else:
            history = series.iloc[:origin + 1]  # inclusive of origin = "now"
        preds = forecast_fn(history, horizons_samples, max_lookback_samples, **forecast_fn_kwargs)
        if preds is None:
            continue

        for h in horizons_samples:
            target_idx = origin + h
            if target_idx >= n:
                continue
            actual = values[target_idx]
            if np.isnan(actual):
                continue
            predicted = preds[h]
            rows.append({
                "origin_idx": origin,
                "origin_time": series.index[origin],
                "horizon_samples": h,
                "actual": actual,
                "predicted": predicted,
                "error": actual - predicted,
            })

    return pd.DataFrame(rows)


def summarize_backtest(results: pd.DataFrame, dt_hours: float):
    """Per-horizon RMSE/MAE/n, with horizon converted to hours for
    readability."""
    if len(results) == 0:
        return pd.DataFrame()
    rows = []
    for h, group in results.groupby("horizon_samples"):
        rmse = float(np.sqrt(np.mean(group["error"] ** 2)))
        mae = float(np.mean(np.abs(group["error"])))
        rows.append({
            "horizon_samples": h,
            "horizon_hours": h * dt_hours,
            "n": len(group),
            "rmse": rmse,
            "mae": mae,
            "mse": float(np.mean(group["error"] ** 2)),
        })
    return pd.DataFrame(rows).sort_values("horizon_samples")


def make_arma_forecast_fn(fitted_result, max_horizon_samples: int, freq_str: str = None):
    """Wraps a statsmodels SARIMAX/ARIMA fit into a forecast_fn compatible
    with rolling_origin_backtest.

    Deliberately does NOT refit at every origin - full MLE refitting at
    each of tens of thousands of origins would repeat the exact
    performance mistake found and fixed in the persistence baseline,
    at far higher per-call cost (numerical optimization, not just a
    dropna() scan). Parameters are fixed from one initial fit; each
    origin only updates the model's internal state on its (bounded)
    history window via .apply(refit=False), which reuses the fixed
    coefficients and just re-runs the Kalman filter on the given data -
    still O(window), but a much smaller constant.

    freq_str: explicit pandas offset alias (e.g. "30min"). Slicing a
    DatetimeIndex drops its .freq attribute even when the underlying
    spacing is still perfectly regular - without this, statsmodels has
    to re-infer the frequency on every single .apply() call and emits a
    ValueWarning each time, flooding the output at ~2000+ backtest
    origins. Setting it explicitly (we already know it - it's dt_hours
    from Stage 0) fixes this at the source instead of suppressing the
    symptom.
    """
    def forecast_fn(history: pd.Series, horizons_samples: list, max_lookback_samples: int):
        valid = history.dropna()
        if len(valid) < 5:
            return None
        if freq_str is not None and valid.index.freq is None:
            valid = valid.asfreq(freq_str)
        try:
            applied = fitted_result.apply(valid, refit=False)
            fc = applied.forecast(steps=max_horizon_samples)
            return {h: float(fc.iloc[h - 1]) for h in horizons_samples}
        except Exception:
            return None
    return forecast_fn


def select_order(train: pd.Series, p_range, q_range, d_range):
    """AIC grid search for ARIMA order. Moved here from 18_forecast_arma.py
    so Stage D (ARMA-GARCH) can reuse it - Python module names can't start
    with a digit, so `import` from "18_forecast_arma" isn't possible;
    shared logic needs to live in this module instead."""
    from statsmodels.tsa.arima.model import ARIMA
    best_aic = np.inf
    best_order = None
    for d in d_range:
        for p in p_range:
            for q in q_range:
                if p == 0 and q == 0:
                    continue
                try:
                    fit = ARIMA(train.values, order=(p, d, q)).fit()
                    if fit.aic < best_aic:
                        best_aic = fit.aic
                        best_order = (p, d, q)
                except Exception:
                    continue
    return best_order, best_aic


def skill_score(model_summary: pd.DataFrame, baseline_summary: pd.DataFrame):
    """1 - (model_MSE / baseline_MSE) per horizon - positive means the
    model beats the baseline, 0 means identical, negative means worse
    than the baseline. Joins on horizon_samples."""
    merged = model_summary.merge(baseline_summary, on="horizon_samples",
                                  suffixes=("_model", "_baseline"))
    merged["skill_score"] = 1 - (merged["mse_model"] / merged["mse_baseline"])
    return merged[["horizon_samples", "horizon_hours_model", "mse_model",
                    "mse_baseline", "skill_score"]].rename(
        columns={"horizon_hours_model": "horizon_hours"})
