
set.seed(123) # since bootsrapping

# realised
stopifnot(exists("squared_returns"), exists("train_size"), exists("n_test"))
realized <- as.numeric(squared_returns[(train_size + 1):(train_size + n_test)])
stopifnot(length(realized) == n_test)

# model forecast matrix 
# STES models 
stopifnot(exists("pred_ae"), exists("pred_se"), exists("pred_e_ae"), exists("pred_e_se"))
.align_to_n <- function(x, n) { x <- as.numeric(x); if (length(x) >= n) x[1:n] else c(x, rep(NA_real_, n - length(x))) }

stes_forecasts <- data.frame(
  STES_AE   = .align_to_n(pred_ae,   n_test),
  STES_SE   = .align_to_n(pred_se,   n_test),
  STES_E_AE = .align_to_n(pred_e_ae, n_test),
  STES_E_SE = .align_to_n(pred_e_se, n_test)
)

# GARCH models
if (exists("results")) {
  garch_forecasts <- as.data.frame(results)
  if (is.null(colnames(garch_forecasts)) || any(!nzchar(colnames(garch_forecasts)))) {
    colnames(garch_forecasts) <- c("sGARCH","eGARCH","gjrGARCH","iGARCH","apARCH","TGARCH")[seq_len(ncol(garch_forecasts))]
  }
  garch_forecasts[] <- lapply(garch_forecasts, as.numeric)
  stopifnot(nrow(garch_forecasts) == n_test)
  all_models <- cbind(garch_forecasts, stes_forecasts)
} else {
  message("Note: 'results' (GARCH forecasts) not found. MCS will run on STES models only.")
  all_models <- stes_forecasts
}
stopifnot(nrow(all_models) == n_test)

# Build loss matrices (AE, SE, QLIKE) 
eps_clip <- 1e-8 #for qlike
Loss_AE <- sapply(all_models, function(vhat) abs(as.numeric(vhat) - realized))
Loss_SE <- sapply(all_models, function(vhat) (as.numeric(vhat) - realized)^2)
Loss_QLIKE <- sapply(all_models, function(vhat) {
  vhat <- pmax(as.numeric(vhat), eps_clip)
  log(vhat) + realized / vhat
})

# Clean losses
clean_loss <- function(Loss) {
  Loss <- as.matrix(Loss)
  colnames(Loss) <- colnames(all_models)[seq_len(ncol(Loss))]
  drop_cols <- vapply(as.data.frame(Loss),
                      function(x) mean(!is.finite(x) | is.na(x)) > 0.20,
                      logical(1))
  if (any(drop_cols)) {
    message("Dropping models (>20% non-finite losses): ",
            paste(colnames(Loss)[drop_cols], collapse = ", "))
    Loss <- Loss[, !drop_cols, drop = FALSE]
  }
  for (j in seq_len(ncol(Loss))) {
    x <- Loss[, j]
    if (any(!is.finite(x) | is.na(x))) {
      med <- median(x[is.finite(x)], na.rm = TRUE)
      x[!is.finite(x) | is.na(x)] <- med
      Loss[, j] <- x
    }
  }
  Loss
}

Loss_list <- list(
  AE    = clean_loss(Loss_AE),
  SE    = clean_loss(Loss_SE),
  QLIKE = clean_loss(Loss_QLIKE)
)

# extract SSM from MCS object
get_ssm <- function(mcs_obj) {
  tbl <- as.data.frame(mcs_obj@show)  
  keep <- rep(FALSE, nrow(tbl))
  if ("MCS_M" %in% names(tbl)) keep <- keep | (tbl$MCS_M == 1)
  if ("MCS_R" %in% names(tbl)) keep <- keep | (tbl$MCS_R == 1)
  rownames(tbl)[keep]
}

#  Run MCS for each loss
B_boot <- 2000   
alpha  <- 0.10

cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

run_mcs <- function(Loss, loss_name) {
  cat("\n============================\n")
  cat("MCS on", loss_name, "loss\n")
  cat("============================\n")
  
  mcs_tmax <- MCSprocedure(Loss = Loss, alpha = alpha, B = B_boot,
                           cl = cl, statistic = "Tmax", verbose = TRUE)
  mcs_tr   <- MCSprocedure(Loss = Loss, alpha = alpha, B = B_boot,
                           cl = cl, statistic = "TR",   verbose = FALSE)
  
  cat("\n--- Tmax ---\n"); print(mcs_tmax)
  cat("\n--- TR ---\n");   print(mcs_tr)
  
  ssm_tmax <- get_ssm(mcs_tmax)
  ssm_tr   <- get_ssm(mcs_tr)
  ssm_int  <- intersect(ssm_tmax, ssm_tr)
  
  cat("SSM (Tmax):        ", paste(ssm_tmax, collapse = ", "), "\n")
  cat("SSM (TR):          ", paste(ssm_tr,   collapse = ", "), "\n")
  cat("SSM (Intersection):", paste(ssm_int,  collapse = ", "), "\n")
  
  invisible(list(tmax = mcs_tmax, tr = mcs_tr,
                 ssm_tmax = ssm_tmax, ssm_tr = ssm_tr, ssm_intersection = ssm_int))
}

mcs_results <- lapply(names(Loss_list), function(nm) run_mcs(Loss_list[[nm]], nm))
names(mcs_results) <- names(Loss_list)

# summary table 
summarise_ssm <- function(res_list) {
  out <- do.call(rbind, lapply(names(res_list), function(nm) {
    r <- res_list[[nm]]
    data.frame(
      Loss        = nm,
      SSM_Tmax    = paste(r$ssm_tmax,        collapse = "; "),
      SSM_TR      = paste(r$ssm_tr,          collapse = "; "),
      SSM_Both    = paste(r$ssm_intersection,collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

ssm_summary <- summarise_ssm(mcs_results)
print(ssm_summary)


