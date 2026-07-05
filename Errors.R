#### Error Metrics #######

# Loss functions (with positivity clip for QLIKE)
floor_val <- 1e-6
qlike_safe <- function(forecast, actual) {
  vhat <- pmax(forecast, floor_val)
  vtrue <- pmax(actual,   floor_val)
  mean(log(vhat) + vtrue / vhat, na.rm = TRUE)
}
rmse  <- function(forecast, actual) sqrt(mean((forecast - actual)^2, na.rm=TRUE))
mae   <- function(forecast, actual) mean(abs(forecast - actual), na.rm=TRUE)
medae <- function(forecast, actual) median(abs(forecast - actual), na.rm=TRUE)


# Realized variance in test set
realized <- squared_returns[(train_size+1):(train_size+n_test)]
stopifnot(length(realized) == n_test)

# Forecasts from garch 
garch_forecasts <- as.data.frame(results)

if (is.null(colnames(garch_forecasts))) {
  colnames(garch_forecasts) <- c("sGARCH","eGARCH","gjrGARCH","iGARCH","apARCH","TGARCH")
}

#Forecasts from STES
stes_forecasts <- data.frame(
  STES_AE   = pred_ae,
  STES_SE   = pred_se,
  STES_E_AE = pred_e_ae,
  STES_E_SE = pred_e_se
)

# Combine all models
all_models <- cbind(garch_forecasts, stes_forecasts)
common_idx <- is.finite(realized)
for (j in seq_len(ncol(all_models))) common_idx <- common_idx & is.finite(all_models[[j]])
realized_c   <- realized[common_idx]
all_models_c <- all_models[common_idx, , drop = FALSE]


# Compute metrics -
metric_names <- c("RMSE","MAE","MedAE","QLIKE")
error_table <- data.frame(
  Model = colnames(all_models),
  RMSE  = NA_real_,
  MAE   = NA_real_,
  MedAE = NA_real_,
  QLIKE = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(all_models)) {
  f <- all_models[[i]]
  error_table$RMSE[i]  <- rmse(f, realized)
  error_table$MAE[i]   <- mae(f, realized)
  error_table$MedAE[i] <- medae(f, realized)
  error_table$QLIKE[i] <- qlike_safe(f, realized)
}

#Rank
error_table$RMSE_rank  <- rank(error_table$RMSE,  ties.method="min")
error_table$MAE_rank   <- rank(error_table$MAE,   ties.method="min")
error_table$MedAE_rank <- rank(error_table$MedAE, ties.method="min")
error_table$QLIKE_rank <- rank(error_table$QLIKE, ties.method="min")
error_table$AvgRank <- rowMeans(error_table[, c("RMSE_rank","MAE_rank","MedAE_rank","QLIKE_rank")])

# Sort by AvgRank
error_table_sorted <- error_table[order(error_table$AvgRank, error_table$QLIKE), ]

print(error_table_sorted)

# Winners by metric
cat("Best models by metric:\n")
for (m in metric_names) {
  best_idx <- which.min(error_table_sorted[[m]])
  cat(sprintf("%-6s: %s\n", m, error_table_sorted$Model[best_idx]))
}



