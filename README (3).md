# Forecasting Commodity Volatility: GARCH vs. Smooth Transition Exponential Smoothing

**One-step-ahead volatility forecasting across five commodities, comparing the GARCH family against a from-scratch implementation of Smooth Transition Exponential Smoothing (STES), with statistically rigorous model selection.**

---

## TL;DR

I built and benchmarked **10 forecasting models** to predict next-day volatility for five commodities (WTI crude oil, natural gas, gold, corn, palladium). The project answers a practical question: *when you have several competing forecasting models, which one actually wins — and can you prove the difference is real rather than luck?*

- **Implemented STES from scratch** (four transition-variable specifications) — not a library call, a full implementation of the method including the adaptive smoothing parameter and Nelder-Mead optimisation.
- **Benchmarked against six GARCH-family models** (sGARCH, eGARCH, GJR-GARCH, iGARCH, apARCH, TGARCH) with Student-t errors.
- **Evaluated out-of-sample** using rolling 1,000-day windows, parallelised across CPU cores.
- **Tested statistical significance** of the results using Hansen's **Model Confidence Set** (2,000 bootstrap replications) — so the conclusions are backed by significance testing, not just a leaderboard.

**Headline result:** there is no universally best model. GARCH wins under QLIKE loss (which heavily penalises under-predicting volatility — the metric that matters for risk); STES wins under absolute-error metrics (better at tracking the *magnitude* of moves). The right model depends on the asset's volatility structure and on what error you most want to avoid.

---

## Why this project is relevant beyond finance

The finance framing is the domain, but the transferable skill is **time-series forecasting and rigorous model comparison** — the same machinery used for churn, retention, LTV, and demand forecasting.

| What this project does | The equivalent analytics problem |
|---|---|
| Forecast next-day volatility from past returns | Forecast next-period churn / conversion from past behaviour |
| Compare 10 models on out-of-sample data | Compare competing predictive models before shipping one |
| Choose the loss function that matches the business cost | Decide whether false positives or false negatives cost more |
| Prove the best model is *statistically* better (MCS) | Avoid deploying a model that only looked better by chance |

The core competency on display: **take a forecasting problem, implement competing models properly, evaluate them out-of-sample with the right loss function, and prove statistically which one genuinely wins.** That final step — testing whether a model is *actually* better rather than better-on-this-sample — is what separates a considered analysis from running defaults.

---

## Methodology

**Data.** Daily closing prices for five commodities, sourced from the public [Historical Prices of Major Natural Resources](https://www.kaggle.com/datasets/albertobircoci/historical-prices-of-major-natural-resource) dataset on Kaggle (Albert5913, 2025). Log returns computed; squared daily returns used as the realised-variance proxy (a standard, if noisy, choice — see Limitations).

**Models.**
- *GARCH family (6):* sGARCH, eGARCH, GJR-GARCH, iGARCH, apARCH, TGARCH — all GARCH(1,1) with Student-t innovations, via `rugarch`.
- *STES (4):* Smooth Transition Exponential Smoothing with Absolute Error (AE), Squared Error (SE), and extended sign+magnitude variants (E-AE, E-SE) as transition variables. Implemented from scratch — logistic transition function, adaptive smoothing parameter α_t, parameters estimated by minimising the sum of squared one-step errors (Nelder-Mead).

**Evaluation.** Rolling-window one-step-ahead forecasts (1,000-observation windows), parallelised. Four loss functions: RMSE, MAE, MedAE, and QLIKE. Models ranked per-commodity and averaged.

**Statistical model selection.** Hansen, Lunde & Nason's (2011) **Model Confidence Set** at the 10% level, 2,000 bootstrap replications, using both Tmax and TR statistics, to identify the set of models that are not statistically distinguishable from the best.

---

## Key findings

- **No universal winner.** Performance depends on the asset's volatility structure *and* the evaluation metric.
- **QLIKE favours GARCH.** Because QLIKE punishes under-estimating volatility, GARCH — which folds past shocks directly into the variance equation — dominated the Model Confidence Set under this loss across four of five commodities. This is the loss to care about for value-at-risk and risk-capital use cases.
- **Absolute-error metrics favour STES.** STES-AE and STES-EAE entered the Superior Set under MAE/MedAE across three commodities — they track the magnitude of gradual volatility changes well, but respond to sudden spikes with a one-step lag.
- **Asset structure drives model choice.** Symmetric volatility clusters (WTI) favoured standard GARCH; leverage/asymmetry (palladium, corn) favoured GJR-GARCH/TGARCH/eGARCH; relatively stable series with occasional bursts (gold) split by metric.
- **Simpler STES beat the extended variants.** The squared-error transition specifications added instability without accuracy — a reminder that added model complexity has to earn its place.

---

## Repository structure

```
├── data/
│   └── Resources_Dataset.csv        # Public Kaggle commodity price data
├── src/
│   ├── GARCH_V10.R                  # Six GARCH-family specs, rolling forecasts
│   ├── STES_V10.R                   # From-scratch STES (4 variants) + adaptive alpha
│   ├── Errors_V10.R                 # RMSE / MAE / MedAE / QLIKE + model ranking
│   ├── MCS_V10.R                    # Model Confidence Set (bootstrap significance)
│   ├── Oil_V10.R / Gas_V10.R / Gold_V10.R / Corn_V10.R / Palladium_V10.R
│                                    # Per-commodity data prep & driver scripts
└── README.md
```

## How to run

Requires R with the following packages:

```r
install.packages(c("rugarch", "MCS", "foreach", "doParallel", "xts", "PerformanceAnalytics"))
```

Each per-commodity script prepares the return series and sets `train_size`, `n_test`, `eps`, and `squared_returns`; `GARCH_V10.R` and `STES_V10.R` then produce the forecasts, and `Errors_V10.R` / `MCS_V10.R` evaluate them. Rolling-window estimation is parallelised across available cores.

---

## Limitations & honest scope

- **Noisy variance proxy.** Squared daily returns proxy realised variance; intraday high-frequency data would give a cleaner target.
- **No cross-commodity spillovers.** Each series is modelled univariately, though oil demonstrably influences gas, corn, etc.
- **GARCH restricted to (1,1)** while multiple STES specifications were explored — a slight asymmetry in the comparison.
- **No exogenous drivers.** Geopolitical, weather, and macro shocks that move commodity volatility are outside the models.

These are stated plainly because knowing what a model *doesn't* capture is part of using it responsibly.

---

## Context

Originally developed as my MSc Business Analytics dissertation at Lancaster University (2025), supervised by Dr Alisa Yusupova. Reworked here as a standalone project. The dissertation is available on request.

**Data:** Albert5913 (2025), *Historical Prices of Major Natural Resource*, Kaggle. **Method references:** Taylor (2004) on STES; Hansen, Lunde & Nason (2011) on the Model Confidence Set.
