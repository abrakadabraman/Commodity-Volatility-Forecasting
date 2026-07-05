# Forecasting Commodity Volatility: GARCH vs. STES

One-step-ahead volatility forecasting for five commodities (WTI crude oil, natural gas, gold, corn, palladium), comparing six GARCH-family models against a from-scratch implementation of Smooth Transition Exponential Smoothing (STES), with model selection via Hansen's Model Confidence Set.



## Data

Daily commodity closing prices from the public [Historical Prices of Major Natural Resources](https://www.kaggle.com/datasets/albertobircoci/historical-prices-of-major-natural-resource) dataset on Kaggle (Albert5913, 2025). Included as `Resources_Dataset.csv`.

## Requirements

R, with:

```r
install.packages(c("rugarch", "MCS", "foreach", "doParallel", "xts", "PerformanceAnalytics"))
```

## Files

```
Resources_Dataset.csv     # Price data
oil.R / gas.R / gold.R / corn.R / palladium.R   # Per-commodity data prep
garch.R                   # Six GARCH-family specifications
stes.R                    # STES (four variants), implemented from scratch
error_metrics.R           # RMSE / MAE / MedAE / QLIKE + rankings
model_confidence_set.R    # Model Confidence Set (bootstrap significance test)
```

## How to run

Keep `Resources_Dataset.csv` in the same folder as the scripts, then run in this order:

1. **One commodity prep script** (e.g. `oil.R`) - loads the data and sets `train_size`, `n_test`, `eps`, `squared_returns`.
2. **`garch.R`** — produces GARCH forecasts (`results`).
3. **`stes.R`** — produces STES forecasts (`pred_ae`, `pred_se`, `pred_e_ae`, `pred_e_se`).
4. **`error_metrics.R`** — computes error metrics and rankings.
5. **`model_confidence_set.R`** — runs the Model Confidence Set.

Repeat with a different prep script (`gas.R`, `gold.R`, etc.) to evaluate another commodity. Rolling-window estimation is parallelised across available cores.

## References

Data: Albert5913 (2025), *Historical Prices of Major Natural Resource*, Kaggle. Methods: Taylor (2004), STES; Hansen, Lunde & Nason (2011), Model Confidence Set.
