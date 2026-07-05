
##### GARCH ####

# Window 
window_size <- 1000
ncores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(ncores)
registerDoParallel(cl)

#  Helper
safe_sigma2 <- function(fit) {
  out <- tryCatch({
    fc <- ugarchforecast(fit, n.ahead = 1)
    as.numeric(sigma(fc))^2
  }, error = function(e) NA_real_)
  out
}

# GARCH Student t
spec_sGARCH <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
spec_eGARCH <- ugarchspec(
  variance.model = list(model = "eGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
spec_gjr <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
spec_iGARCH <- ugarchspec(
  variance.model = list(model = "iGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
spec_apARCH <- ugarchspec(
  variance.model = list(model = "apARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)
spec_tgarch <- ugarchspec(
  variance.model = list(model = "fGARCH", submodel = "TGARCH", garchOrder = c(1,1)),
  mean.model     = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)

# cluster
clusterExport(cl, c("spec_sGARCH","spec_eGARCH","spec_gjr","spec_iGARCH","spec_apARCH","spec_tgarch",
                    "safe_sigma2"), envir = environment())

## Rolling 1-step-ahead forecasts
results <- foreach(i = seq_len(n_test), .combine = rbind, .packages = "rugarch") %dopar% {
  idx_end   <- train_size + i - 1         
  idx_start <- max(1, idx_end - window_size + 1)
  wdata <- eps[idx_start:idx_end]         
  
  out <- rep(NA_real_, 6)
  
  fit1 <- tryCatch(ugarchfit(spec_sGARCH, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[1] <- if (is.null(fit1)) NA_real_ else safe_sigma2(fit1)
  
  fit2 <- tryCatch(ugarchfit(spec_eGARCH, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[2] <- if (is.null(fit2)) NA_real_ else safe_sigma2(fit2)
  
  fit3 <- tryCatch(ugarchfit(spec_gjr, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[3] <- if (is.null(fit3)) NA_real_ else safe_sigma2(fit3)
  
  fit4 <- tryCatch(ugarchfit(spec_iGARCH, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[4] <- if (is.null(fit4)) NA_real_ else safe_sigma2(fit4)
  
  fit5 <- tryCatch(ugarchfit(spec_apARCH, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[5] <- if (is.null(fit5)) NA_real_ else safe_sigma2(fit5)
  
  fit6 <- tryCatch(ugarchfit(spec_tgarch, data = wdata, solver = "hybrid"), error = function(e) NULL)
  out[6] <- if (is.null(fit6)) NA_real_ else safe_sigma2(fit6)
  
  out
}

stopCluster(cl)

#  Label columns 
colnames(results) <- c("sGARCH","eGARCH","gjrGARCH","iGARCH","apARCH","TGARCH")

# Realized variance for test set
realized_test_var <- squared_returns[(train_size+1):(train_size+n_test)]
stopifnot(nrow(results) == length(realized_test_var))

# loss functions
qlike <- function(vhat, vtrue) mean(log(pmax(vhat, .Machine$double.eps)) + vtrue / pmax(vhat, .Machine$double.eps))
mse   <- function(vhat, vtrue) mean((vhat - vtrue)^2)
mae   <- function(vhat, vtrue) mean(abs(vhat - vtrue))

score_tbl <- data.frame(
  model = colnames(results),
  QLIKE = apply(results, 2, qlike, vtrue = realized_test_var),
  MSE   = apply(results, 2, mse,   vtrue = realized_test_var),
  MAE   = apply(results, 2, mae,   vtrue = realized_test_var)
)
print(score_tbl[order(score_tbl$QLIKE), ], row.names = FALSE)

# Plots 
plot_len <- length(realized_test_var)
matplot(1:plot_len, cbind(realized_test_var, results),
        type="l", lty=1,
        col=c("black","red","blue","forestgreen","orange","purple","brown"),
        lwd=c(2, rep(1.5, 6)),
        main="GARCH 1-step Variance Forecasts vs Actual (Palladium)",
        ylab= "Realised Variance", xlab="Time (test)")
legend("topright", legend=c("Actual", colnames(results)),
       col=c("black","red","blue","forestgreen","orange","purple","brown"),
       lwd=c(2, rep(1.5, 6)), bty="n")

errs <- abs(results - matrix(realized_test_var, nrow = plot_len, ncol = ncol(results)))
matplot(errs, type="l", lty=1, lwd=1.3,
        col=c("red","blue","forestgreen","orange","purple","brown"),
        main="Absolute Forecast Errors (GARCH)",
        ylab="|Error|", xlab="Time (test)")
legend("topright", legend=colnames(results),
       col=c("red","blue","forestgreen","orange","purple","brown"), lwd=1.3, bty="n")

plot(density(errs[,1], na.rm=TRUE), col="red", lwd=2,
     main="|Error| Density – GARCH models", xlab="|Error|")
for(i in 2:ncol(results)) lines(density(errs[,i], na.rm=TRUE),
                                col=c("blue","forestgreen","orange","purple","brown")[i-1], lwd=2)
legend("topright", legend=colnames(results),
       col=c("red","blue","forestgreen","orange","purple","brown"), lwd=2, bty="n")


