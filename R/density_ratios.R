cf_density_ratios <- function(task, learners, mtp, control, pb) {
  ans <- vector("list", length = length(task$folds))
  
  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }
  
  # For production / parallel execution
  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_density_ratios(task, fold, learners, mtp, control, pb)
    },
    seed = TRUE)
  }
  
  # NON-FUTURE FOR DEBUGGING
  # for (fold in seq_along(task$folds)) {
  #   ans[[fold]] <- estimate_density_ratios(task, fold, learners, mtp, control, pb)
  # }
  
  ans <- future::value(ans)
  
  ans <- list(density_ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]]))
  
  ans$density_ratios <- trim(ans$density_ratios, control$.trim)
  ans
}

estimate_density_ratios <- function(task, fold, learners, mtp, control, pb) {
  # get natural and shifted data for this fold
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)
  
  density_ratios <- matrix(nrow = nrow(natural$valid), ncol = task$time_horizon)
  fits <- vector("list", length = task$time_horizon)
  
  for (time in seq_len(task$time_horizon)) {
    # indices: observed up to time - 1 and at risk at time
    i <- task$observed(natural$train, time - 1) %and% task$is_at_risk(natural$train, time)
    i <- rep(i, 2)
    
    # current treatment variable name
    A_t <- current_trt(task$vars$A, time)
    
    # stacked data for censoring model (A will not be used as a predictor)
    stacked <- stack_data(
      natural = natural$train,
      shifted = shifted$train,
      trt = task$vars$A,
      cens = task$vars$C,
      time = time
    )
    
    # ---- TREATMENT MODEL (separate, intercept-only) ----
    # treatment is randomized, so we fit A_t ~ 1 on the natural data only
    treat_fit <- glm(
      stats::as.formula(paste(A_t, "~ 1")),
      data = natural$train[task$observed(natural$train, time - 1) &
                             task$is_at_risk(natural$train, time), , drop = FALSE],
      family = binomial()
    )
    
    # ---- CENSORING MODEL (separate, SuperLearner) ----
    # covariates supplied for censoring model (baseline W, Z, R)
   
    vars_cens <- c("..i..lmtp_id", task$vars$C[time], "..i..lmtp_stack_indicator")
    vars_cens <- stats::na.omit(vars_cens)
    
    cens_fit <- run_ensemble(
      x = stacked[i, vars_cens, drop = FALSE],
      y = "..i..lmtp_stack_indicator",
      learners = learners,
      family = "binomial",
      id = "..i..lmtp_id",
      folds = control$.learners_trt_folds,
      discrete = control$.discrete,
      info = control$.info
    )
    
    # store fits
    if (control$.return_full_fits) {
      fits[[time]] <- list(
        treatment = treat_fit,
        censoring = cens_fit
      )
    } else {
      fits[[time]] <- list(
        treatment = treat_fit,
        censoring = if (is.null(cens_fit)) NULL else extract_sl_weights(cens_fit)
      )
    }
    
    # ---- PREDICTION AND DENSITY RATIOS (using censoring model only) ----
    i_valid <- task$observed(natural$valid, time - 1) %and%
      task$is_at_risk(natural$valid, time)
    
    pred <- matrix(-999L, nrow = nrow(natural$valid), ncol = 1)
    
    # run_ensemble returns a SuperLearner object; use its predict method
    pred[i_valid, ] <- predict(cens_fit, newdata = natural$valid[i_valid, , drop = FALSE])
    
    obs <- task$observed(natural$valid, time)
    at_risk <- task$is_at_risk(natural$valid, time)
    followed <- followed_rule_NEW(natural$valid, shifted$valid, A_t, mtp)
    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    density_ratios[, time] <- (pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))
    pb()
  }
  
  list(
    ratios = density_ratios,
    fits = fits
  )
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