# ---------- internal helpers ----------
.get_pretty_map <- function(fit) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  labs <- fit$param_labels %||% list(T_cov=character(0), T_lag=character(0),
                                     Y_cov=character(0), Y_lag=character(0))

  all_sum <- rstan::summary(fit$stan_fit)$summary
  rn <- rownames(all_sum)

  beta_tags  <- grep("^beta0\\[\\d+\\]$",  rn, value = TRUE)
  theta_tags <- grep("^theta0\\[\\d+\\]$", rn, value = TRUE)
  max_beta   <- if (length(beta_tags))  max(as.integer(sub("^beta0\\[(\\d+)\\]$","\\1",  beta_tags)))  else 0L
  max_theta  <- if (length(theta_tags)) max(as.integer(sub("^theta0\\[(\\d+)\\]$","\\1", theta_tags))) else 0L

  map <- character(0)

  k <- 0L
  for (nm in labs$T_cov) { k <- k+1L; map[paste0("beta0[", k, "]")] <- paste0("beta_T:", nm) }
  for (nm in labs$T_lag) { k <- k+1L; map[paste0("beta0[", k, "]")] <- paste0("theta_T_lag:", nm) }
  if (max_beta > k) {
    for (i in (k+1):max_beta) {
      map[paste0("beta0[", i, "]")] <- paste0("time_baseline_T[", i - k, "]")
    }
  }

  k <- 0L
  for (nm in labs$Y_cov) { k <- k+1L; map[paste0("theta0[", k, "]")] <- paste0("beta_Y:", nm) }
  for (nm in labs$Y_lag) { k <- k+1L; map[paste0("theta0[", k, "]")] <- paste0("theta_Y_lag:", nm) }
  if (max_theta > k) {
    for (i in (k+1):max_theta) {
      map[paste0("theta0[", i, "]")] <- paste0("time_baseline_Y[", i - k, "]")
    }
  }

  map["beta1"]  <- "treatment_effect_T"
  map["theta1"] <- "treatment_effect_Y"

  map["beta0_star"]  <- "beta_T_intercept"
  map["theta0_star"] <- "beta_Y_intercept"

  tl_tags <- grep("^thetaLag\\[\\d+\\]$", rn, value = TRUE)
  if (length(tl_tags)) {
    for (tg in tl_tags) {
      idx <- sub("^thetaLag\\[(\\d+)\\]$","\\1", tg)
      map[tg] <- paste0("theta_lag_extra[", idx, "]")
    }
  }
  map
}

.map_param_names <- function(fit, x) {
  map <- .get_pretty_map(fit)
  if (length(map) == 0L) return(x)
  repl <- unname(map[match(x, names(map))])
  x_new <- ifelse(is.na(repl), x, repl)
  x_new
}

Summary_lable <- function(fit) {
  stopifnot(inherits(fit$stan_fit, "stanfit"))
  sm <- rstan::summary(fit$stan_fit)$summary
  rn <- rownames(sm)
  rn <- .map_param_names(fit, rn)
  out <- as.data.frame(sm)
  out$Parameter <- rn
  rownames(out) <- NULL
  out <- out[, c("Parameter","mean","sd","2.5%","97.5%","n_eff","Rhat","se_mean"), drop = FALSE]
  names(out) <- c("Parameter","Mean","SD","2.5%","97.5%","n_eff","Rhat","MCSE")
  out$CI_width <- out$`97.5%` - out$`2.5%`
  out
}
