# Forecasting Commodity Volatility: GARCH vs. Smooth Transition Exponential Smoothing

**One-step-ahead volatility forecasting across five commodities, comparing the GARCH family against a from-scratch implementation of Smooth Transition Exponential Smoothing (STES)**

---



Built and benchmarked **10 forecasting models** to predict next-day volatility for five commodities (WTI crude oil, natural gas, gold, corn, palladium).

- **Implemented STES from scratch** (four transition-variable specifications) - a full implementation of the method including the adaptive smoothing parameter and Nelder-Mead optimisation.
- **Benchmarked against six GARCH-family models** (sGARCH, eGARCH, GJR-GARCH, iGARCH, apARCH, TGARCH) with Student-t errors.
- **Evaluated out-of-sample** using rolling 1,000-day windows, parallelised across CPU cores.
- **Tested statistical significance** of the results using Hansen's **Model Confidence Set** (2,000 bootstrap replications)
  




## Methodology

**Data.** Daily closing prices for five commodities, sourced from the public [Historical Prices of Major Natural Resources](https://www.kaggle.com/datasets/albertobircoci/historical-prices-of-major-natural-resource) dataset on Kaggle (Albert5913, 2025). Log returns computed; squared daily returns used as the realised-variance proxy (a standard, if noisy, choice - see Limitations).

**Models.**
- *GARCH family (6):* sGARCH, eGARCH, GJR-GARCH, iGARCH, apARCH, TGARCH - all GARCH(1,1) with Student-t innovations, via `rugarch`.
- *STES (4):* Smooth Transition Exponential Smoothing with Absolute Error (AE), Squared Error (SE), and extended sign+magnitude variants (E-AE, E-SE) as transition variables. Implemented from scratch - logistic transition function, adaptive smoothing parameter α_t, parameters estimated by minimising the sum of squared one-step errors (Nelder-Mead).

**Evaluation.** Rolling-window one-step-ahead forecasts (1,000-observation windows), parallelised. Four loss functions: RMSE, MAE, MedAE, and QLIKE. Models ranked per-commodity and averaged.

**Statistical model selection.** Hansen, Lunde & Nason's (2011) **Model Confidence Set** at the 10% level, 2,000 bootstrap replications, using both Tmax and TR statistics, to identify the set of models that are not statistically distinguishable from the best.

---



## Repository structure

```
├── data/
│   └── Resources_Dataset.csv        # Public Kaggle commodity price data
├── src/
│   ├── GARCH.R                  # Six GARCH-family specs, rolling forecasts
│   ├── STES.R                   # From-scratch STES (4 variants) + adaptive alpha
│   ├── Errors.R                 # RMSE / MAE / MedAE / QLIKE + model ranking
│   ├── MCS.R                    # Model Confidence Set (bootstrap significance)
│   ├── Oil.R / Gas.R / Gold.R / Corn.R / Palladium.R
│                                    # Per-commodity data prep & driver scripts
└── README.md
```

## How to run

Requires R with the following packages:

```r
install.packages(c("rugarch", "MCS", "foreach", "doParallel", "xts", "PerformanceAnalytics"))
```

Each per-commodity script prepares the return series and sets `train_size`, `n_test`, `eps`, and `squared_returns`; `GARCH_V10.R` and `STES_V10.R` then produce the forecasts, and `Errors_V10.R` / `MCS_V10.R` evaluate them. Rolling-window estimation is parallelised across available cores.

