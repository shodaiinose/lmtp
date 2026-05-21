cf_density_ratios <- function(task, learners, mtp, control, pb) {
  ans <- vector("list", length = length(task$folds))
  
  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }
  
  # For production / parallel execution
  # for (fold in seq_along(task$folds)) {
  #   ans[[fold]] <- future::future({
  #     estimate_density_ratios(task, fold, learners, mtp, control, pb)
  #   },
  #   seed = TRUE)
  # }
  
  # NON-FUTURE FOR DEBUGGING
  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- estimate_density_ratios(task, fold, learners, mtp, control, pb)
  }
  
  ans <- future::value(ans)
  density_ratios= recombine(rbind_depth(ans, "ratios"), task$folds)
  ans <- list(density_ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]]))
  
  ans$density_ratios <- trim(ans$density_ratios, control$.trim)
  ans
}

estimate_density_ratios <- function(task, fold, learners, mtp, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  
  density_ratios <- matrix(nrow = nrow(natural$valid), ncol = task$time_horizon)
  fits <- vector("list", length = task$time_horizon)
  
  for (time in seq_len(task$time_horizon)) {
    i <- task$observed(natural$train, time - 1) %and% task$is_at_risk(natural$train, time)
    i <- rep(i, 2)
    A_t <- current_trt(task$vars$A, time)
    vars <- c("..i..lmtp_id", task$vars$history("A", time), task$vars$C[time], "..i..lmtp_stack_indicator") # remove A from vars
    
    vars <- na.omit(vars)
    stacked <- stack_data(natural$train, shifted$train, task$vars$A, task$vars$C, time)

    # ---- TREATMENT MODEL (no covariates, intercept-only) ----
    # treatment is randomized, so we fit A_t ~ 1 on the natural data only
    if (time == 1) {
      fit <- glm(as.formula(paste(A_t, "~ 1")),
        data = natural$train, # Run on the entire original sample
        family = binomial())
      
      
    # ---- CENSORING MODEL (W, Z, R, and an ensemble of learners) ----
    } else {
      fit <- run_ensemble(stacked[i, vars], "..i..lmtp_stack_indicator",
                               learners, "binomial", "..i..lmtp_id",
                               control$.learners_trt_folds, 
                               control$.discrete, 
                               control$.info
      )
      
    }
    
    # store fits
    if (control$.return_full_fits) {
      fits[[time]] <- fit
    } else {
      if (time == 1) {
        fits[[time]] <- list(glm_fit = fit)
      } else {
        fits[[time]] <- extract_sl_weights(fit)
      }
    }
    
    # ---- PREDICTION AND DENSITY RATIOS ----
    i <- task$observed(natural$valid, time - 1) %and% task$is_at_risk(natural$valid, time)
    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    
    if (time == 1) {
      pred[i, ] <- predict(fit, newdata = natural$valid[i, ], type = "response")
    } else {
      pred[i, ] <- predict(fit, natural$valid[i, ])
    }
    
    obs <- task$observed(natural$valid, time)
    at_risk <- task$is_at_risk(natural$valid, time)
    followed <- followed_rule_NEW(natural$valid, shifted$valid, A_t, mtp)
    
    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    density_ratios[, time] <- (pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))
    pb()
  }
  
  list(ratios = density_ratios, fits = fits)
}


stack_data <- function(natural, shifted, trt, cens, time) {
  shifted_half <- natural
  
  if (length(trt) > 1 || time == 1) {
    if (!is.na(trt[[time]])) {
      shifted_half[, trt[[time]]] <- shifted[, trt[[time]]]
    }
  }
  
  if (!is.null(cens)) {
    shifted_half[[cens[time]]] <- shifted[[cens[time]]]
  }
  
  out <- rbind(natural, shifted_half)
  out[["..i..lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
  out
}