## ======================================================================
## Script  : Cox + FP +  multiples corrections
## Based on Liquet, Roux & Riou (2026)
## ======================================================================
library(survival)

## ----------------------------------------------------------------------
## Based Functions
## ----------------------------------------------------------------------

safe_log <- function(x, eps = 1e-10) {
  x[x <= 0] <- eps
  log(x)
}

validate_data <- function(data, time_var, status_var, X_var, 
                          covariates = NULL, min_n = 50) {
  errors <- character(0)
  
  if (is.null(data) || !is.data.frame(data) || nrow(data) == 0) {
    errors <- c(errors, "NA")
    return(list(valid = FALSE, errors = errors, n_complete = 0))
  }
  
  required_vars <- c(time_var, status_var, X_var)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    errors <- c(errors, paste("Missing variables:", paste(missing_vars, collapse = ", ")))
  }
  
  if (!is.null(covariates)) {
    missing_cov <- setdiff(covariates, names(data))
    if (length(missing_cov) > 0) {
      errors <- c(errors, paste("Missi,g Covariables:", paste(missing_cov, collapse = ", ")))
      covariates <- setdiff(covariates, missing_cov)
    }
  }
  
  if (!is.null(data[[time_var]])) {
    if (!is.numeric(data[[time_var]])) {
      errors <- c(errors, "The time variable must be numeric.")
    } else if (any(data[[time_var]] <= 0, na.rm = TRUE)) {
      errors <- c(errors, "The time variable must be strictly positive.")
    }
  }
  
  if (!is.null(data[[status_var]])) {
    sv <- data[[status_var]]
    if (is.factor(sv)) sv <- as.numeric(as.character(sv))
    if (!all(sv %in% c(0, 1) | is.na(sv))) {
      errors <- c(errors, "The status variable must be binary (0/1).")
    }
  }
  
  if (!is.null(data[[X_var]]) && !is.numeric(data[[X_var]])) {
    errors <- c(errors, "The variable X must be numeric.")
  }
  
  cols_to_check <- unique(c(required_vars, covariates))
  n_complete <- sum(complete.cases(data[, cols_to_check, drop = FALSE]))
  
  if (n_complete < min_n) {
    errors <- c(errors, paste("Sample too small:", n_complete, 
                              "complete observations (minimum:", min_n, ")"))
  }
  
  list(valid = length(errors) == 0, errors = errors, n_complete = n_complete)
}

make_X_positive <- function(X, offset = 0.1, verbose = FALSE) {
  if (is.null(X) || length(X) == 0) return(X)
  if (!is.numeric(X)) stop("X must be numeric")
  if (all(is.na(X))) return(X)
  
  minX <- min(X, na.rm = TRUE)
  if (!is.finite(minX)) return(X)
  
  if (minX <= 0) {
    shift <- abs(minX) + offset
    X <- X + shift
    if (verbose) message("Offset applied to X: ", round(shift, 4))
  }
  X
}

## ----------------------------------------------------------------------
## Terms of fractional polynomials
## ----------------------------------------------------------------------

get_fp_terms <- function(x, powers_list, eps = 1e-6) {
  if (!is.numeric(x)) stop("x must be numeric")
  
  if (is.list(powers_list) && length(powers_list) >= 1 && is.numeric(powers_list[[1]])) {
    powers <- powers_list[[1]]
  } else if (is.numeric(powers_list)) {
    powers <- powers_list
  } else {
    stop("powers_list must be a numeric vector or a list of vectors.")
  }
  
  if (any(x <= 0, na.rm = TRUE)) {
    x <- x + abs(min(x, na.rm = TRUE)) + eps
  }
  
  if (length(powers) == 1) {
    p <- powers[1]
    term1 <- if (abs(p) < eps) safe_log(x) else x^p
    if (any(!is.finite(term1) & !is.na(term1))) return(NULL)
    return(matrix(term1, ncol = 1))
  }
  
  if (length(powers) == 2) {
    p1 <- powers[1]; p2 <- powers[2]
    if (abs(p1 - p2) < eps) {
      if (abs(p1) < eps) {
        term1 <- safe_log(x)
        term2 <- safe_log(x)^2
      } else {
        term1 <- x^p1
        term2 <- x^p1 * safe_log(x)
      }
    } else {
      term1 <- if (abs(p1) < eps) safe_log(x) else x^p1
      term2 <- if (abs(p2) < eps) safe_log(x) else x^p2
    }
    if (any((!is.finite(term1) & !is.na(term1)) | (!is.finite(term2) & !is.na(term2)))) {
      return(NULL)
    }
    return(cbind(term1, term2))
  }
  
  stop("powers must contain 1 or 2 elements (FP°1 or FP°2)")
}

get_fp_pvalue <- function(fit, term_prefix = "fp") {
  tryCatch({
    coefs <- coef(fit)
    vcov_matrix <- vcov(fit)
    
    idx <- grep(paste0("^", term_prefix), names(coefs))
    if (length(idx) == 0) return(NA)
    
    if (length(idx) == 1) {
      z <- coefs[idx] / sqrt(vcov_matrix[idx, idx])
      return(2 * (1 - pnorm(abs(z))))
    }
    
    beta <- coefs[idx]
    V <- vcov_matrix[idx, idx, drop = FALSE]
    if (any(is.na(V)) || any(!is.finite(V)) || det(V) <= 1e-12) return(NA)
    
    w <- as.numeric(t(beta) %*% solve(V) %*% beta)
    if (!is.finite(w)) return(NA)
    1 - pchisq(w, df = length(idx))
  }, error = function(e) NA)
}

fit_fp_models <- function(data, time_var, status_var, covariates = NULL, 
                          X_var, powers_list, offset = 0.1, 
                          min_n_obs = 30, verbose = FALSE) {
  validation <- validate_data(data, time_var, status_var, X_var, covariates)
  if (!validation$valid) {
    return(list(error = paste(validation$errors, collapse = "; ")))
  }
  
  required_vars <- c(time_var, status_var, X_var, covariates)
  df <- data[, required_vars, drop = FALSE]
  df <- df[complete.cases(df), ]
  if (nrow(df) < min_n_obs) {
    return(list(error = paste("Sample too small:", nrow(df))))
  }
  
  n_events <- sum(df[[status_var]])
  if (n_events < 10) {
    return(list(error = paste("Too few events:", n_events)))
  }
  
  df$SurvObj <- with(df, Surv(df[[time_var]], df[[status_var]]))
  df$X_pos <- make_X_positive(df[[X_var]], offset = offset, verbose = verbose)
  
  fit_single_power <- function(powers) {
    fp_matrix <- get_fp_terms(df$X_pos, list(powers))
    if (is.null(fp_matrix)) return(NULL)
    if (any(apply(fp_matrix, 2, sd, na.rm = TRUE) < 1e-8)) return(NULL)
    
    df_model <- df
    fp_names <- paste0("fp", seq_len(ncol(fp_matrix)))
    for (j in seq_along(fp_names)) df_model[[fp_names[j]]] <- fp_matrix[, j]
    
    if (length(covariates) > 0) {
      formula_str <- paste("SurvObj ~", paste(c(covariates, fp_names), collapse = " + "))
    } else {
      formula_str <- paste("SurvObj ~", paste(fp_names, collapse = " + "))
    }
    
    fit <- tryCatch(coxph(as.formula(formula_str), data = df_model),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    
    p_val <- get_fp_pvalue(fit, term_prefix = "fp")
    aic_val <- AIC(fit)
    deviance_val <- fit$loglik[2]
    
    list(powers = powers, pvalue = p_val, aic = aic_val, 
         deviance = deviance_val, fit = fit)
  }
  
  all_results <- lapply(powers_list, fit_single_power)
  all_results <- Filter(Negate(is.null), all_results)
  if (length(all_results) == 0) {
    return(list(error = "No FP model could be adjusted."))
  }
  
  all_pvalues <- sapply(all_results, function(x) x$pvalue)
  best_idx <- which.min(all_pvalues)
  
  list(
    best_pval = all_pvalues[best_idx],
    best_power = all_results[[best_idx]]$powers,
    best_fit = all_results[[best_idx]]$fit,
    best_aic = all_results[[best_idx]]$aic,
    all_results = all_results,
    all_pvalues = all_pvalues,
    n_obs = nrow(df),
    data_used = df
  )
}

## ----------------------------------------------------------------------
## Conventional correction methods
## ----------------------------------------------------------------------

bonferroni_fp_correction <- function(all_pvalues) {
  K <- length(all_pvalues)
  min_p <- min(all_pvalues, na.rm = TRUE)
  p_corrected <- min(1, K * min_p)
  list(
    p_corrected = p_corrected,
    n_tests = K,
    method = "Bonferroni"
  )
}

royston_sauerbrei_procedure <- function(data, time_var, status_var, X_var, 
                                        covariates = NULL, powers_list, 
                                        offset = 0.1, alpha = 0.05, verbose = FALSE) {
  
  fp1_powers <- powers_list[sapply(powers_list, length) == 1]
  fp2_powers <- powers_list[sapply(powers_list, length) == 2]
  
  if (length(fp1_powers) == 0 || length(fp2_powers) == 0) {
    return(list(error = "The procedure requires FP1 and FP2."))
  }
  
  required_vars <- c(time_var, status_var, X_var, covariates)
  df <- data[, required_vars, drop = FALSE]
  df <- df[complete.cases(df), ]
  
  if (length(covariates) > 0) {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                               paste(covariates, collapse = " + ")))
  } else {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ 1"))
  }
  fit_null <- coxph(f_null, data = df)
  deviance_null <- -2 * fit_null$loglik[2]
  
  fp1_results <- fit_fp_models(data, time_var, status_var, covariates, X_var, fp1_powers, offset)
  fp2_results <- fit_fp_models(data, time_var, status_var, covariates, X_var, fp2_powers, offset)
  
  if (!is.null(fp1_results$error) || !is.null(fp2_results$error)) {
    return(list(error = "Failure to adjust FP models"))
  }
  
  best_fp1_idx <- which.max(sapply(fp1_results$all_results, function(x) x$deviance))
  best_fp2_idx <- which.max(sapply(fp2_results$all_results, function(x) x$deviance))
  
  best_fp1 <- fp1_results$all_results[[best_fp1_idx]]
  best_fp2 <- fp2_results$all_results[[best_fp2_idx]]
  
  deviance_fp2 <- -2 * best_fp2$fit$loglik[2]
  lrt_fp2_vs_null <- deviance_null - deviance_fp2
  df_fp2 <- length(best_fp2$powers)
  p_fp2_vs_null <- 1 - pchisq(lrt_fp2_vs_null, df = df_fp2)
  
  if (verbose) {
    cat("\n Royston & Sauerbrei procedure:\n")
    cat("  Test FP2 vs null: p =", format.pval(p_fp2_vs_null, digits = 4), "\n")
  }
  
  if (p_fp2_vs_null >= alpha) {
    return(list(
      final_model = "null",
      final_powers = NULL,
      p_corrected = p_fp2_vs_null,
      method = "Royston & Sauerbrei MFP",
      decision = "No significant effect",
      all_tests = list(fp2_vs_null = p_fp2_vs_null)
    ))
  }
  
  df_linear <- df
  df_linear$X_pos <- make_X_positive(df[[X_var]], offset = offset)
  if (length(covariates) > 0) {
    f_linear <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                                 paste(c(covariates, "X_pos"), collapse = " + ")))
  } else {
    f_linear <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ X_pos"))
  }
  fit_linear <- coxph(f_linear, data = df_linear)
  deviance_linear <- -2 * fit_linear$loglik[2]
  
  lrt_fp2_vs_linear <- deviance_linear - deviance_fp2
  df_diff <- df_fp2 - 1
  p_fp2_vs_linear <- 1 - pchisq(lrt_fp2_vs_linear, df = df_diff)
  
  if (verbose) {
    cat("  Test FP2 vs linear: p =", format.pval(p_fp2_vs_linear, digits = 4), "\n")
  }
  
  if (p_fp2_vs_linear >= alpha) {
    return(list(
      final_model = "linear",
      final_powers = 1,
      p_corrected = p_fp2_vs_null,
      method = "Royston & Sauerbrei MFP",
      decision = "Linear effect",
      all_tests = list(fp2_vs_null = p_fp2_vs_null, fp2_vs_linear = p_fp2_vs_linear)
    ))
  }
  
  deviance_fp1 <- -2 * best_fp1$fit$loglik[2]
  lrt_fp2_vs_fp1 <- deviance_fp1 - deviance_fp2
  df_diff2 <- df_fp2 - length(best_fp1$powers)
  p_fp2_vs_fp1 <- 1 - pchisq(lrt_fp2_vs_fp1, df = df_diff2)
  
  if (verbose) {
    cat("  Test FP2 vs FP1: p =", format.pval(p_fp2_vs_fp1, digits = 4), "\n")
  }
  
  if (p_fp2_vs_fp1 < alpha) {
    return(list(
      final_model = "FP2",
      final_powers = best_fp2$powers,
      p_corrected = p_fp2_vs_null,
      method = "Royston & Sauerbrei MFP",
      decision = "FP2 selected",
      all_tests = list(
        fp2_vs_null = p_fp2_vs_null,
        fp2_vs_linear = p_fp2_vs_linear,
        fp2_vs_fp1 = p_fp2_vs_fp1
      )
    ))
  } else {
    return(list(
      final_model = "FP1",
      final_powers = best_fp1$powers,
      p_corrected = p_fp2_vs_null,
      method = "Royston & Sauerbrei MFP",
      decision = "FP1 selected",
      all_tests = list(
        fp2_vs_null = p_fp2_vs_null,
        fp2_vs_linear = p_fp2_vs_linear,
        fp2_vs_fp1 = p_fp2_vs_fp1
      )
    ))
  }
}

## ----------------------------------------------------------------------
## Simulation under H0
## ----------------------------------------------------------------------

simulate_survival_H0 <- function(fit_null, df, time_var, status_var) {
  n <- nrow(df)
  lp <- predict(fit_null, newdata = df, type = "lp")
  
  bh <- basehaz(fit_null, centered = FALSE)
  times <- bh$time
  H0 <- bh$hazard
  
  if (times[1] > 0) {
    times <- c(0, times)
    H0 <- c(0, H0)
  }
  
  invert_cumhaz <- function(target) {
    idx <- findInterval(target, H0)
    if (idx == 0) {
      if (H0[1] > 0) return(times[1] * target / H0[1])
      return(0)
    }
    if (idx >= length(times)) return(times[length(times)])
    t1 <- times[idx]; t2 <- times[idx + 1]
    h1 <- H0[idx];    h2 <- H0[idx + 1]
    if (h2 == h1) return(t1)
    t1 + (t2 - t1) * (target - h1) / (h2 - h1)
  }
  
  U <- runif(n)
  target_H <- -log(U) / exp(lp)
  event_times <- sapply(target_H, invert_cumhaz)
  
  cens_obs <- df[[time_var]][df[[status_var]] == 0]
  if (length(cens_obs) < 5) cens_obs <- df[[time_var]]
  cens_star <- sample(cens_obs, n, replace = TRUE)
  
  obs_times <- pmin(event_times, cens_star)
  status <- as.numeric(event_times <= cens_star)
  
  list(time = obs_times, status = status)
}

## ----------------------------------------------------------------------
## Parametric Bootstrap
## ----------------------------------------------------------------------

bootstrap_parametric <- function(data, time_var, status_var, X_var, 
                                 covariates = NULL, powers_list, 
                                 offset = 0.1, B = 500) {
  required_vars <- c(time_var, status_var, X_var, covariates)
  df <- data[, required_vars, drop = FALSE]
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  if (n < 30) return(list(error = "Sample too small for bootstrap"))
  
  if (length(covariates) > 0) {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                               paste(covariates, collapse = " + ")))
  } else {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ 1"))
  }
  fit_null <- tryCatch(coxph(f_null, data = df), error = function(e) NULL)
  if (is.null(fit_null)) return(list(error = "Model failure"))
  
  fp_obs <- fit_fp_models(df, time_var, status_var, covariates, X_var, powers_list, offset)
  if (!is.null(fp_obs$error)) return(list(error = fp_obs$error))
  
  p_observed <- fp_obs$all_pvalues
  K <- length(p_observed)
  
  p_boot_matrix <- matrix(NA, nrow = B, ncol = K)
  n_success <- 0
  
  for (b in seq_len(B)) {
    sim <- tryCatch(
      simulate_survival_H0(fit_null, df, time_var, status_var),
      error = function(e) NULL
    )
    if (is.null(sim)) next
    
    df_boot <- df
    df_boot[[time_var]] <- sim$time
    df_boot[[status_var]] <- sim$status
    
    df_boot[[X_var]] <- sample(df[[X_var]])
    
    if (sum(df_boot[[status_var]]) < 5) next
    
    fp_boot <- tryCatch(
      fit_fp_models(df_boot, time_var, status_var, covariates, X_var, powers_list, offset, verbose = FALSE),
      error = function(e) list(error = e$message)
    )
    
    if (is.null(fp_boot$error) && !is.null(fp_boot$all_pvalues)) {
      if (length(fp_boot$all_pvalues) == K && all(!is.na(fp_boot$all_pvalues))) {
        p_boot_matrix[b, ] <- fp_boot$all_pvalues
        n_success <- n_success + 1
      }
    }
  }
  
  p_boot_valid <- p_boot_matrix[complete.cases(p_boot_matrix), , drop = FALSE]
  min_required <- max(15, 0.05 * B)
  if (nrow(p_boot_valid) < min_required) {
    return(list(
      error = paste("Too few valid bootstraps:", nrow(p_boot_valid), "/", B, 
                    "(minimum requirement:", min_required, ")"),
      n_bootstrap = nrow(p_boot_valid),
      success_rate = n_success / B
    ))
  }
  
  min_p_boot <- apply(p_boot_valid, 1, min)
  min_p_obs  <- min(p_observed)
  
  p_corr <- (1 + sum(min_p_boot <= min_p_obs)) / (nrow(p_boot_valid) + 1)
  
  list(
    p_corrected = p_corr,
    min_p_obs = min_p_obs,
    n_bootstrap = nrow(p_boot_valid),
    success_rate = n_success / B,
    method = "Parametric Bootstrap"
  )
}

## ----------------------------------------------------------------------
## Semi-parametric Bootstrap 
## ----------------------------------------------------------------------

bootstrap_semiparametric <- function(data, time_var, status_var, X_var,
                                     covariates = NULL, powers_list,
                                     offset = 0.1, B = 500) {
  df <- data[, c(time_var, status_var, X_var, covariates), drop = FALSE]
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  if (n < 30) return(list(error = "Sample too small"))
  
  if (length(covariates) > 0) {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                               paste(covariates, collapse = " + ")))
  } else {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ 1"))
  }
  fit_null <- coxph(f_null, data = df)
  
  fp_obs <- fit_fp_models(df, time_var, status_var, covariates, X_var, powers_list, offset)
  if (!is.null(fp_obs$error)) return(list(error = fp_obs$error))
  p_obs <- fp_obs$all_pvalues
  K <- length(p_obs)
  
  p_boot <- matrix(NA, nrow = B, ncol = K)
  n_success <- 0
  
  for (b in 1:B) {
    G <- rexp(n, rate = 1)
    
    fit_null_weighted <- tryCatch(
      coxph(f_null, data = df, weights = G),
      error = function(e) NULL
    )
    
    if (is.null(fit_null_weighted)) next
    
    lp_w <- predict(fit_null_weighted, newdata = df, type = "lp")
    
    bh_w <- basehaz(fit_null_weighted, centered = FALSE)
    times_w <- bh_w$time
    cumhaz_w <- bh_w$hazard
    
    if (times_w[1] > 0) {
      times_w <- c(0, times_w)
      cumhaz_w <- c(0, cumhaz_w)
    }
    
    invert_cumhaz <- function(target) {
      idx <- findInterval(target, cumhaz_w)
      if (idx == 0) {
        if (cumhaz_w[1] > 0) return(times_w[1] * target / cumhaz_w[1])
        return(0)
      }
      if (idx >= length(times_w)) return(times_w[length(times_w)])
      t1 <- times_w[idx]; t2 <- times_w[idx + 1]
      h1 <- cumhaz_w[idx]; h2 <- cumhaz_w[idx + 1]
      if (h2 == h1) return(t1)
      t1 + (t2 - t1) * (target - h1) / (h2 - h1)
    }
    
    U <- runif(n)
    target <- -log(U) / exp(lp_w)
    event_times <- sapply(target, invert_cumhaz)
    if (any(!is.finite(event_times))) next
    
    cens_obs <- df[[time_var]][df[[status_var]] == 0]
    if (length(cens_obs) < 5) cens_obs <- df[[time_var]]
    cens_star <- sample(cens_obs, n, replace = TRUE)
    obs_times <- pmin(event_times, cens_star)
    status_star <- as.numeric(event_times <= cens_star)
    if (sum(status_star) < 5) next
    
    df_boot <- df
    df_boot[[time_var]] <- obs_times
    df_boot[[status_var]] <- status_star
    
       df_boot[[X_var]] <- sample(df[[X_var]])
    
    fp_boot <- tryCatch(
      fit_fp_models(df_boot, time_var, status_var, covariates, X_var, powers_list, offset, verbose = FALSE),
      error = function(e) list(error = e$message)
    )
    
    if (is.null(fp_boot$error) && !is.null(fp_boot$all_pvalues)) {
      if (length(fp_boot$all_pvalues) == K && all(!is.na(fp_boot$all_pvalues))) {
        p_boot[b, ] <- fp_boot$all_pvalues
        n_success <- n_success + 1
      }
    }
  }
  
  p_boot_valid <- p_boot[complete.cases(p_boot), , drop = FALSE]
  min_required <- max(15, 0.05 * B)
  if (nrow(p_boot_valid) < min_required) {
    return(list(
      error = paste("Too few valid bootstraps:", nrow(p_boot_valid), "/", B, "(minimum requirement:", min_required, ")"),
      n_bootstrap = nrow(p_boot_valid),
      success_rate = n_success / B
    ))
  }
  
  min_p_boot <- apply(p_boot_valid, 1, min)
  min_p_obs <- min(p_obs)
  p_corr <- (1 + sum(min_p_boot <= min_p_obs)) / (nrow(p_boot_valid) + 1)
  
  list(
    p_corrected = p_corr,
    min_p_obs = min_p_obs,
    n_bootstrap = nrow(p_boot_valid),
    success_rate = n_success / B,
    method = "Semi-parametric Bootstrap"
  )
}

## ----------------------------------------------------------------------
##  Wild Bootstrap
## ----------------------------------------------------------------------

bootstrap_wild <- function(data, time_var, status_var, X_var,
                           covariates = NULL, powers_list,
                           offset = 0.1, B = 500) {
  df <- data[, c(time_var, status_var, X_var, covariates), drop = FALSE]
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  if (n < 30) return(list(error = "Sample too small"))
  
  if (length(covariates) > 0) {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                               paste(covariates, collapse = " + ")))
  } else {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ 1"))
  }
  fit_null <- coxph(f_null, data = df)
  
  mart_resid <- residuals(fit_null, type = "martingale")
  
  fp_obs <- fit_fp_models(df, time_var, status_var, covariates, X_var, powers_list, offset)
  if (!is.null(fp_obs$error)) return(list(error = fp_obs$error))
  p_obs <- fp_obs$all_pvalues
  K <- length(p_obs)
  
  p_boot <- matrix(NA, nrow = B, ncol = K)
  n_success <- 0
  
  for (b in 1:B) {
    W <- sample(c(-1, 1), n, replace = TRUE)
    perturbed_resid <- W * mart_resid
    
    bh <- basehaz(fit_null, centered = FALSE)
    lp <- predict(fit_null, type = "lp")
    
    cumhaz_at_obs <- approx(bh$time, bh$hazard, xout = df[[time_var]], 
                            method = "constant", f = 0, yleft = 0, 
                            yright = max(bh$hazard))$y
    
    expected_events <- cumhaz_at_obs * exp(lp)
    pseudo_delta <- perturbed_resid + expected_events
    
    threshold <- quantile(pseudo_delta, probs = sum(df[[status_var]]) / n)
    status_star <- as.numeric(pseudo_delta >= threshold)
    
    if (sum(status_star) < 5) next
    
    df_boot <- df
    df_boot[[status_var]] <- status_star
    
    
    df_boot[[X_var]] <- sample(df[[X_var]])
    
    fp_boot <- tryCatch(
      fit_fp_models(df_boot, time_var, status_var, covariates, X_var, powers_list, offset, verbose = FALSE),
      error = function(e) list(error = e$message)
    )
    
    if (is.null(fp_boot$error) && !is.null(fp_boot$all_pvalues)) {
      if (length(fp_boot$all_pvalues) == K && all(!is.na(fp_boot$all_pvalues))) {
        p_boot[b, ] <- fp_boot$all_pvalues
        n_success <- n_success + 1
      }
    }
  }
  
  p_boot_valid <- p_boot[complete.cases(p_boot), , drop = FALSE]
  min_required <- max(15, 0.05 * B)
  if (nrow(p_boot_valid) < min_required) {
    return(list(
      error = paste("Too few valid bootstraps:", nrow(p_boot_valid), "/", B, "(minimum requirement:", min_required, ")"),
      n_bootstrap = nrow(p_boot_valid),
      success_rate = n_success / B
    ))
  }
  
  min_p_boot <- apply(p_boot_valid, 1, min)
  min_p_obs <- min(p_obs)
  p_corr <- (1 + sum(min_p_boot <= min_p_obs)) / (nrow(p_boot_valid) + 1)
  
  list(
    p_corrected = p_corr,
    min_p_obs = min_p_obs,
    n_bootstrap = nrow(p_boot_valid),
    success_rate = n_success / B,
    method = "Wild Bootstrap (FIXED)"
  )
}

## ----------------------------------------------------------------------
## Martingale permutation 
## ----------------------------------------------------------------------

permutation_martingale <- function(data, time_var, status_var, X_var,
                                   covariates = NULL, powers_list,
                                   offset = 0.1, B = 500) {
  df <- data[, c(time_var, status_var, X_var, covariates), drop = FALSE]
  df <- df[complete.cases(df), ]
  n <- nrow(df)
  if (n < 30) return(list(error = "Sample too small"))
  
  if (length(covariates) > 0) {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~",
                               paste(covariates, collapse = " + ")))
  } else {
    f_null <- as.formula(paste("Surv(", time_var, ",", status_var, ") ~ 1"))
  }
  fit_null <- coxph(f_null, data = df)
  
  mart_resid <- residuals(fit_null, type = "martingale")
  
  fp_obs <- fit_fp_models(df, time_var, status_var, covariates, X_var, powers_list, offset)
  if (!is.null(fp_obs$error)) return(list(error = fp_obs$error))
  p_obs <- fp_obs$all_pvalues
  K <- length(p_obs)
  
  p_perm <- matrix(NA, nrow = B, ncol = K)
  n_success <- 0
  
  for (b in 1:B) {
    perm_idx <- sample(1:n)
    perturbed_resid <- mart_resid[perm_idx]
    
    bh <- basehaz(fit_null, centered = FALSE)
    lp <- predict(fit_null, type = "lp")
    
    cumhaz_at_obs <- approx(bh$time, bh$hazard, xout = df[[time_var]], 
                            method = "constant", f = 0, yleft = 0, 
                            yright = max(bh$hazard))$y
    
    expected_events <- cumhaz_at_obs * exp(lp)
    pseudo_delta <- perturbed_resid + expected_events
    
    threshold <- quantile(pseudo_delta, probs = sum(df[[status_var]]) / n)
    status_star <- as.numeric(pseudo_delta >= threshold)
    
    if (sum(status_star) < 5) next
    
    df_perm <- df
    df_perm[[status_var]] <- status_star
    

    df_perm[[X_var]] <- sample(df[[X_var]])
    
    fp_perm <- tryCatch(
      fit_fp_models(df_perm, time_var, status_var, covariates, X_var, powers_list, offset, verbose = FALSE),
      error = function(e) list(error = e$message)
    )
    
    if (is.null(fp_perm$error) && !is.null(fp_perm$all_pvalues)) {
      if (length(fp_perm$all_pvalues) == K && all(!is.na(fp_perm$all_pvalues))) {
        p_perm[b, ] <- fp_perm$all_pvalues
        n_success <- n_success + 1
      }
    }
  }
  
  p_perm_valid <- p_perm[complete.cases(p_perm), , drop = FALSE]
  if (nrow(p_perm_valid) < 20) {
    return(list(
      error = paste("Too few valid permutations:", nrow(p_perm_valid), "/", B),
      n_bootstrap = nrow(p_perm_valid),
      success_rate = n_success / B
    ))
  }
  
  min_p_perm <- apply(p_perm_valid, 1, min)
  min_p_obs <- min(p_obs)
  p_corr <- (1 + sum(min_p_perm <= min_p_obs)) / (nrow(p_perm_valid) + 1)
  
  list(
    p_corrected = p_corr,
    min_p_obs = min_p_obs,
    n_bootstrap = nrow(p_perm_valid),
    success_rate = n_success / B,
    method = "Permutation Martingale "
  )
}



cox_fp_analysis <- function(formula, data, var_interest, fp_powers,
                            method = c("naive", "bonferroni", "royston",
                                       "parametric_bootstrap", "semiparametric_bootstrap",
                                       "wild_bootstrap", "perm_martingale"),
                            B = 500, alpha = 0.05, offset = 0.1, verbose = TRUE) {
  
  method <- match.arg(method)
  
  response_var <- all.vars(formula)[1:2]
  time_var <- response_var[1]
  status_var <- response_var[2]
  covariates <- all.vars(formula)[-c(1:2)]
  if (length(covariates) == 0) covariates <- NULL
  
  validation <- validate_data(data, time_var, status_var, var_interest, covariates)
  if (!validation$valid) {
    stop(paste("Invalid data:", paste(validation$errors, collapse = "; ")))
  }
  
  if (verbose) message("Adjustment of FP models...")
  fp_results <- fit_fp_models(
    data = data, time_var = time_var, status_var = status_var,
    covariates = covariates, X_var = var_interest,
    powers_list = fp_powers, offset = offset, verbose = verbose
  )
  if (!is.null(fp_results$error)) stop(fp_results$error)
  
  p_raw <- fp_results$best_pval
  all_pvals <- fp_results$all_pvalues
  
  if (verbose) message("Application of the correction: ", method)
  
  correction_result <- switch(
    method,
    
    "naive" = list(
      p_corrected = p_raw,
      method = "Naive (no correction)",
      n_tests = length(all_pvals)
    ),
    
    "bonferroni" = bonferroni_fp_correction(all_pvals),
    
    "royston" = {
      royston_result <- royston_sauerbrei_procedure(
        data = data, time_var = time_var, status_var = status_var,
        X_var = var_interest, covariates = covariates,
        powers_list = fp_powers, offset = offset, alpha = alpha, verbose = verbose
      )
      if (!is.null(royston_result$error)) {
        warning(royston_result$error)
        list(p_corrected = p_raw, method = "Royston & Sauerbrei")
      } else royston_result
    },
    
    "parametric_bootstrap" = {
      boot_result <- bootstrap_parametric(
        data = data, time_var = time_var, status_var = status_var,
        X_var = var_interest, covariates = covariates,
        powers_list = fp_powers, offset = offset, B = B
      )
      if (!is.null(boot_result$error)) {
        warning(boot_result$error)
        list(p_corrected = p_raw, method = "Parametric Bootstrap")
      } else boot_result
    },
    
    "semiparametric_bootstrap" = {
      boot_result <- bootstrap_semiparametric(
        data = data, time_var = time_var, status_var = status_var,
        X_var = var_interest, covariates = covariates,
        powers_list = fp_powers, offset = offset, B = B
      )
      if (!is.null(boot_result$error)) {
        warning(boot_result$error)
        list(p_corrected = p_raw, method = "Semi-parametric Bootstrap ")
      } else boot_result
    },
    
    "wild_bootstrap" = {
      boot_result <- bootstrap_wild(
        data = data, time_var = time_var, status_var = status_var,
        X_var = var_interest, covariates = covariates,
        powers_list = fp_powers, offset = offset, B = B
      )
      if (!is.null(boot_result$error)) {
        warning(boot_result$error)
        list(p_corrected = p_raw, method = "Wild Bootstrap")
      } else boot_result
    },
    
    "perm_martingale" = {
      perm_result <- permutation_martingale(
        data = data, time_var = time_var, status_var = status_var,
        X_var = var_interest, covariates = covariates,
        powers_list = fp_powers, offset = offset, B = B
      )
      if (!is.null(perm_result$error)) {
        warning(perm_result$error)
        list(p_corrected = p_raw, method = "Permutation")
      } else perm_result
    },


  )
  
  final_results <- list(
    best_power = fp_results$best_power,
    best_pval_raw = p_raw,
    best_aic = fp_results$best_aic,
    best_fit = fp_results$best_fit,
    p_corrected = correction_result$p_corrected,
    correction_method = correction_result$method,
    all_powers = sapply(fp_results$all_results, function(x) x$powers),
    all_pvalues_raw = all_pvals,
    all_results = fp_results$all_results,
    n_obs = fp_results$n_obs,
    n_tests = length(all_pvals),
    alpha = alpha,
    additional_info = correction_result[setdiff(names(correction_result),
                                                c("p_corrected", "method"))]
  )
  class(final_results) <- c("cox_fp_analysis", "list")
  
  if (verbose) {
    cat("\n=== RESULTS ===\n")
    cat("Method :", final_results$correction_method, "\n")
    cat("n =", final_results$n_obs, " | K =", final_results$n_tests, "\n")
    cat("Best FP powers :", paste(final_results$best_power, collapse = ", "), "\n")
    cat("raw-p =", format.pval(final_results$best_pval_raw, digits = 4), "\n")
    cat("adjusted-p =", format.pval(final_results$p_corrected, digits = 4), "\n")
    cat("Significant (α =", alpha, ") :", 
        ifelse(final_results$p_corrected < alpha, "YES", "NO"), "\n")
  }
  
  final_results
}
