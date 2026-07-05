
#### STES
#### loss functions
qlike <- function(vhat, vtrue) mean(log(pmax(vhat, .Machine$double.eps)) + vtrue / pmax(vhat, .Machine$double.eps))
mse   <- function(vhat, vtrue) mean((vhat - vtrue)^2)
mae   <- function(vhat, vtrue) mean(abs(vhat - vtrue))


# y   = ε_t^2   (target variance)
# eps = ε_t     (centered return residual)

stes_taylor <- function(y, eps, x0, modeltype = c("AE","SE","E_AE","E_SE")) {
  modeltype <- match.arg(modeltype)
  n <- length(y); stopifnot(length(eps) == n)
  
  sigma2 <- numeric(n)   
  alpha  <- numeric(n)
  
  B0 <- x0[2]; B1 <- x0[3]; B2 <- x0[4]
  sigma2_0 <- x0[1]
  sigma2[1] <- sigma2_0
  alpha[1]  <- 1/(1 + exp(-B0))  
  
  get_alpha <- function(e_lag){
    sign_term <- e_lag
    mag_term  <- switch(modeltype,
                        AE   = abs(e_lag),
                        SE   = e_lag^2,
                        E_AE = abs(e_lag),
                        E_SE = e_lag^2)
    lin <- if (modeltype %in% c("AE","SE")) B0 + B1*mag_term else B0 + B1*sign_term + B2*mag_term
    
    lin <- max(min(lin, 20), -20)
    1/(1 + exp(-lin))
  }
  
  for (t in 2:n) {
    a_tm1    <- get_alpha(eps[t-1])     
    alpha[t] <- a_tm1
    sigma2[t] <- sigma2[t-1] + a_tm1 * (y[t-1] - sigma2[t-1])  
  }
  
  ### one-step forecast errors on variance
  err <- y - sigma2
  list(sigma2 = sigma2, alpha = alpha, epsilon = err)
}

### minimise sse 
stesEst2 <- function(y, eps, modeltype, x0){
  obj <- function(p){
    fit <- stes_taylor(y, eps, p, modeltype)
    sum(fit$epsilon^2)
  }
  opt <- optim(par = x0, fn = obj, method = "Nelder-Mead",
               control = list(maxit = 4000))
  fit <- stes_taylor(y, eps, opt$par, modeltype)
  c(list(xopt=opt$par, fopt=opt$value), fit)
}

# Rolling forecast
rolling_stes_forecast_parallel <- function(y_full, eps_full, train_size, test_size, window,
                                           modeltype, ncores = max(1, parallel::detectCores()-1)) {
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("stes_taylor", "stesEst2"), envir = environment())
  
  preds <- foreach::foreach(i = seq_len(test_size), .combine = c) %dopar% {
    idx2 <- train_size + i - 1
    idx1 <- max(1, idx2 - window + 1)
    
    y_tr   <- y_full[idx1:idx2]
    eps_tr <- eps_full[idx1:idx2]
    
    # starting values
    x0 <- c(mean(y_tr), qlogis(0.1), 0, 0)
    
    fit <- stesEst2(y_tr, eps_tr, modeltype, x0)
    
    last_sigma2 <- tail(fit$sigma2, 1)    
    e_last      <- tail(eps_tr, 1)        
    
    # α_t from ε_t
    B0 <- fit$xopt[2]; B1 <- fit$xopt[3]; B2 <- fit$xopt[4]
    sign_term <- e_last
    mag_term  <- if (modeltype %in% c("AE","E_AE")) abs(e_last) else e_last^2
    lin       <- if (modeltype %in% c("AE","SE")) B0 + B1*mag_term else B0 + B1*sign_term + B2*mag_term
    lin       <- max(min(lin, 20), -20)
    a_next    <- 1/(1 + exp(-lin))
    
    y_last <- tail(y_tr, 1)               
    # next-day variance forecast
    next_forecast <- last_sigma2 + a_next * (y_last - last_sigma2)
    next_forecast
  }
  preds
}

## all four STES variants
y_full   <- squared_returns
eps_full <- eps
test_size <- n_test

ncores <- max(1, parallel::detectCores() - 1)

pred_ae   <- rolling_stes_forecast_parallel(y_full, eps_full, train_size, test_size, window, "AE",   ncores)
pred_se   <- rolling_stes_forecast_parallel(y_full, eps_full, train_size, test_size, window, "SE",   ncores)
pred_e_ae <- rolling_stes_forecast_parallel(y_full, eps_full, train_size, test_size, window, "E_AE", ncores)
pred_e_se <- rolling_stes_forecast_parallel(y_full, eps_full, train_size, test_size, window, "E_SE", ncores)

#test
actual_test <- y_full[(train_size+1):(train_size+length(pred_ae))]

# metrics
scores <- data.frame(
  model = c("STES_AE","STES_SE","STES_E_AE","STES_E_SE"),
  QLIKE = c(qlike(pred_ae, actual_test),
            qlike(pred_se, actual_test),
            qlike(pred_e_ae, actual_test),
            qlike(pred_e_se, actual_test)),
  MSE   = c(mse(pred_ae, actual_test),
            mse(pred_se, actual_test),
            mse(pred_e_ae, actual_test),
            mse(pred_e_se, actual_test)),
  MAE   = c(mae(pred_ae, actual_test),
            mae(pred_se, actual_test),
            mae(pred_e_ae, actual_test),
            mae(pred_e_se, actual_test))
)
print(scores)

# Plots
plot_len <- length(actual_test)
ylim_all <- range(c(actual_test, pred_ae, pred_se, pred_e_ae, pred_e_se), na.rm=TRUE)

plot(actual_test, type="l", col="black", lwd=2,
     main="STES Forecasts vs Actual (Palladium)",
     ylab="Realised Variance", xlab="Time (Test set)",
     ylim=ylim_all)
lines(pred_ae,   col="red",         lwd=1.5)
lines(pred_se,   col="blue",        lwd=1.5)
lines(pred_e_ae, col="forestgreen", lwd=1.5)
lines(pred_e_se, col="orange",      lwd=1.5)
legend("topright",
       legend=c("Actual","STES_AE","STES_SE","STES_E_AE","STES_E_SE"),
       col=c("black","red","blue","forestgreen","orange"), lwd=c(2,1.5,1.5,1.5,1.5), bty="n")

# Zoomed plot
range_idx <- 1:min(300, plot_len)
plot(actual_test[range_idx], type="l", col="black", lwd=2,
     main="Zoom: Actual vs STES Forecasts",
     ylab="ε_t^2", xlab="Time (Test subset)")
lines(pred_ae[range_idx],   col="red",         lwd=1.5)
lines(pred_se[range_idx],   col="blue",        lwd=1.5)
lines(pred_e_ae[range_idx], col="forestgreen", lwd=1.5)
lines(pred_e_se[range_idx], col="orange",      lwd=1.5)
legend("topright",
       legend=c("Actual","STES_AE","STES_SE","STES_E_AE","STES_E_SE"),
       col=c("black","red","blue","forestgreen","orange"),
       lwd=c(2,1.5,1.5,1.5,1.5), bty="n")

# Error distributions
err_ae   <- abs(actual_test - pred_ae)
err_se   <- abs(actual_test - pred_se)
err_e_ae <- abs(actual_test - pred_e_ae)
err_e_se <- abs(actual_test - pred_e_se)

plot(stats::density(err_ae,   na.rm=TRUE), col="red",         lwd=2,
     main="|Error| Density – STES Variants", xlab="|Error|")
lines(stats::density(err_se,   na.rm=TRUE), col="blue",        lwd=2)
lines(stats::density(err_e_ae, na.rm=TRUE), col="forestgreen", lwd=2)
lines(stats::density(err_e_se, na.rm=TRUE), col="orange",      lwd=2)
legend("topright",
       legend=c("STES_AE","STES_SE","STES_E_AE","STES_E_SE"),
       col=c("red","blue","forestgreen","orange"), lwd=2, bty="n")

### for alphas
extract_last_alpha <- function(modeltype){
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl))
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("stesEst2","stes_taylor","y_full","eps_full","train_size","test_size","window","modeltype"), envir = environment())
  
  foreach::foreach(i = seq_len(test_size), .combine = c) %dopar% {
    idx2 <- train_size + i - 1
    idx1 <- max(1, idx2 - window + 1)
    y_tr   <- y_full[idx1:idx2]
    eps_tr <- eps_full[idx1:idx2]
    x0 <- c(mean(y_tr), qlogis(0.1), 0, 0)
    fit <- stesEst2(y_tr, eps_tr, modeltype, x0)
    tail(fit$alpha, 1)
  }
}

#alpha plot
alphas_ae   <- extract_last_alpha("AE")
alphas_se   <- extract_last_alpha("SE")
alphas_e_ae <- extract_last_alpha("E_AE")
alphas_e_se <- extract_last_alpha("E_SE")
# Plot 
matplot(cbind(alphas_ae, alphas_se, alphas_e_ae, alphas_e_se),
        type="l", lwd=1.5, col=c("red","blue","forestgreen","orange"),
        main="Adaptive α and Volatility – Palladium", ylab="alpha", xlab="Time")

# secondary axis
par(new=TRUE)
plot(stes_test, type="l", col="black", lwd=1, axes=FALSE, xlab="", ylab="")
axis(side=4)   
mtext("Squared returns (test)", side=4, line=3)

legend("topright", legend=c("AE","SE","E_AE","E_SE","Test Series"),
       col=c("red","blue","forestgreen","orange","black"), lwd=1.5, bty="n")