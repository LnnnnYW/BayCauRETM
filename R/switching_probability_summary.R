#' Compute switching-probability summary
#'
#' @description
#' Summarize how often subjects initiate treatment across follow-up intervals
#' on two complementary scales:
#' - mass: the cohort-level probability that initiation occurs at interval k
#'   (proportion of subjects whose first initiation is exactly at k);
#' - hazard: the cohort-level average of the per-subject initiation hazard at k
#'   among those at risk at k (initiation is defined as switching from 0 to 1).
#'
#' @param df A `data.frame` (e.g., from `preprocess_data()`) containing at least
#'   a subject id, a time index, and a binary treatment indicator.
#' @param covariates Optional character vector of column names in `df` to adjust for
#'  in the switching hazard model. Default `NULL` (unadjusted).
#' @param scale Character string choosing the summary/plot scale.
#'   "mass" (default) shows the marginal probability of initiating at interval k
#'   (proportion of the cohort that starts at k); "hazard" shows the conditional
#'   probability of initiating at k among those still at risk at k. This choice
#'   affects both the returned summary and what `plot()` displays.
#' @inheritParams base::print
#'
#' @return An object of class `switching_summary` with components:
#'   - `df_haz`: row-level data with indicators for at-risk status and the
#'     computed `switch_prob` on the chosen scale;
#'   - `by_k`: interval-level summary with columns
#'     `k_idx`, `n_at_risk`, `mean`, `sd_`, `se`, `lo95`, `hi95`, `p25`, `p50`, `p75`;
#'     here `mean` is computed as `sum(switch_prob, na.rm = TRUE) / N_ids`;
#'   - `model`: the fitted logistic regression used to estimate the initiation
#'     hazard on at-risk rows;
#'   - `id_col`, `time_col`, `treat_col`, `death_col`: resolved column names;
#'   - `scale`: the selected scale (`"mass"` or `"hazard"`).
#'
#' @details
#' The data are ordered by subject and time. Within each subject the function
#' constructs the lagged treatment `A_prev`, an at-risk indicator that requires
#' survival up to the start of interval k (if a terminal-event column is present)
#' and no prior initiation (`A_prev == 0`), and an initiation indicator at k.
#' A logistic regression is fit on the at-risk rows with response `init_k` and
#' predictors `factor(<time index column>)` (for example, `factor(k_idx)`) plus
#' any `covariates`. Predicted hazards are assigned to at-risk rows. On the mass
#' scale, the per-subject probability of first initiation at k is computed as the
#' subject's survival up to k-1 (the product over prior at-risk rows of one minus
#' the predicted hazard) multiplied by the predicted hazard at k. The interval-level
#' table `by_k` reports a cohort-level mean and Wald-type 95% limits based on a
#' normal approximation with standard error `sd_ / sqrt(N_ids)`, along with
#' selected quantiles of `switch_prob` for descriptive dispersion.
#'
#' @examples
#' \dontrun{
#' # default: marginal "mass" scale
#' sw <- switching_probability_summary(fit$data_preprocessed)
#' plot(sw, type = "boxplot")
#'
#' # conditional hazard view
#' sw_h <- switching_probability_summary(fit$data_preprocessed, scale = "hazard")
#' plot(sw_h, type = "line")
#' }
#'
#' @export
#'
#' @importFrom stats glm predict quantile
#' @importFrom dplyr arrange group_by mutate summarise ungroup lag
#' @importFrom utils modifyList
#' @export

switching_probability_summary <- function(df, covariates = NULL,
                                          scale = c("mass", "hazard")) {
  if(!scale %in% c("mass", "hazard"))
    stop("`scale` must be one of 'mass' or 'hazard'.")

  # helpers
  .match_any <- function(cands, pool) {
    pool_l <- tolower(pool)
    for (nm in cands) {
      hit <- which(pool_l == tolower(nm))
      if (length(hit)) return(pool[hit[1]])
    }
    NA_character_
  }
  .bt <- function(nm) { ok <- grepl("^[A-Za-z.][A-Za-z0-9_.]*$", nm); ifelse(ok, nm, sprintf("`%s`", nm)) }

  # columns (flexible matching)
  id_col    <- .match_any(c("pat_id","id"),                     names(df))
  time_col  <- .match_any(c("k_idx","k","time_idx","k_fac"),    names(df))
  treat_col <- .match_any(c("A","Ak","treat","treatment"),      names(df))
  death_col <- .match_any(c("T_obs","Tk","death","D"),          names(df))

  if (is.na(id_col) || is.na(treat_col))
    stop("Data must contain an ID column pat_id/id and a treatment column A/Ak.")

  if (is.na(time_col)) {
    df <- df |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::mutate(k_idx = dplyr::row_number()) |>
      dplyr::ungroup()
    time_col <- "k_idx"
  }

  # Order and keep only needed cols early
  df <- df[order(df[[id_col]], df[[time_col]]), , drop = FALSE]

  # build 'at risk' and initiation indicator
  # At risk at k : alive until k-1 AND not yet initiated (A_prev == 0)
  d <- df |>
    dplyr::group_by(.data[[id_col]]) |>
    dplyr::mutate(
      A_prev = dplyr::lag(.data[[treat_col]], default = 0),
      alive_to_km1 = if (!is.na(death_col)) {
        v <- .data[[death_col]]
        c(1L, cumprod(1L - v[-length(v)]))
      } else {
        rep(1L, dplyr::n())
      },
      at_risk = (alive_to_km1 == 1L) & (A_prev == 0)
    ) |>
    dplyr::ungroup()


  # initiation at k among those at risk at k (discrete-time hazard outcome)
  d <- d |>
    dplyr::mutate(init_k = as.integer(at_risk & (.data[[treat_col]] == 1)))

  # hazard model: init_k ~ factor(k) + covariates (on at-risk rows)
  at_risk_rows <- d$at_risk == 1L
  if (!any(at_risk_rows)) stop("No at-risk rows to fit the switching hazard model.")

  cov_map <- character(0)
  if (!is.null(covariates)) {
    cov_map <- covariates[covariates %in% names(d)]
    if (length(setdiff(covariates, cov_map)))
      warning("Dropped missing covariate(s): ",
              paste(setdiff(covariates, cov_map), collapse = ", "))
  }

  fac_time <- "__fac_time__"
  while (fac_time %in% names(d)) fac_time <- paste0(fac_time, "_")
  d[[fac_time]] <- factor(d[[time_col]])

  #max time in at risk rows
  max_time <- max(d[[time_col]][at_risk_rows])

  rhs_terms <- c(.bt(fac_time), .bt(cov_map))
  form <- stats::as.formula(paste0("init_k ~ ", paste(rhs_terms, collapse = " + ")))

  mf <- stats::model.frame(form, data = d[at_risk_rows, , drop = FALSE],
                           na.action = stats::na.omit)

  if (!nrow(mf)) stop("No rows remain after removing NA for hazard model.")
  rows_fit <- as.integer(rownames(mf))

  mod <- stats::glm(form, data = mf,
                    family = stats::binomial,
                    control = stats::glm.control(maxit = 100))

  #use only data that time <= max_time
  d <- d[d[[time_col]] <= max_time, , drop = FALSE]

  haz_hat <- rep(NA_real_, nrow(d))
  haz_hat <- stats::predict(mod, type = "response", newdata = d)

  d$hazard <- haz_hat


  if (scale == "mass") {
    d <- d |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::mutate(
        # one_minus_h = dplyr::if_else(at_risk & is.finite(hazard), 1 - hazard, 1),
        one_minus_h = 1 - hazard,
        S_km1       = dplyr::lag(cumprod(one_minus_h), default = 1),
        # switch_prob = dplyr::if_else(at_risk & is.finite(hazard), S_km1 * (1 - hazard), 0)
        switch_prob = S_km1 * (1 - hazard)
      ) |>
      dplyr::ungroup()
  } else {
    d$switch_prob <- d$hazard
  }

  d_plot <- d

  N <- length(unique(d[[id_col]]))

  by_k <- d |>
    dplyr::group_by(.data[[time_col]]) |>
    dplyr::summarise(
      n_at_risk = sum(at_risk),
      mean   = sum(switch_prob, na.rm = TRUE) / N,   # 总体质量 P(W=k)
      sd_    = stats::sd(switch_prob, na.rm = TRUE),     # 以“全体”为样本
      se     = sd_ / sqrt(N),
      lo95   = pmax(0, mean - 1.96 * se),
      hi95   = pmin(1, mean + 1.96 * se),
      p25    = stats::quantile(switch_prob, 0.25, na.rm = TRUE),
      p50    = stats::quantile(switch_prob, 0.50, na.rm = TRUE),
      p75    = stats::quantile(switch_prob, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::rename(k_idx = !!time_col)


  out <- list(
    df_haz    = d,
    by_k      = by_k,
    model     = mod,
    id_col    = id_col,
    time_col  = time_col,
    treat_col = treat_col,
    death_col = death_col,
    scale     = scale
  )
  class(out) <- "switching_summary"
  out
}


# print / summary / plot methods

#' @describeIn switching_probability_summary Print the interval-level summary table.
#' @export

print.switching_summary <- function(x, ...) {
  cat("Switching probability summary by interval (", x$scale, "):\n", sep = "")
  print(x$by_k)
  invisible(x$by_k)
}

#' @describeIn switching_probability_summary Summary the interval-level summary table
#' @param object A `switching_summary` object.
#' @method summary switching_summary
#' @export

summary.switching_summary <- function(object, ...) {
  print(object, ...)
}

#' @describeIn switching_probability_summary Plot switching summaries across intervals.
#' @param type One of `"boxplot"` (per-interval boxplots of per-subject
#'   indicators with mean -+ 95% Wald bars overlaid) or `"line"` (mean curve with
#'   95% ribbon). Defaults to `"boxplot"`.
#' @param ... Passed to `ggplot2::geom_boxplot()` when `type="boxplot"`.
#' @export

plot.switching_summary <- function(x, type = c("boxplot","line"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")

  by_k   <- x$by_k
  d_plot <- x$df_haz

  #create a data.frame with N * max_time rows, colunmns: id_col, time_col, switch_prob,
  #enter value from d_plot, fill missing with 0
  # figure out max_time
  max_time <- max(d_plot[[x$time_col]], na.rm = TRUE)
  id <- unique(d_plot[[x$id_col]])
  full_grid <- expand.grid(
    id = id,
    k_idx = seq_len(max_time)
  )
  names(full_grid)[1] <- x$id_col
  d_plot <- merge(full_grid, d_plot[, c(x$id_col, x$time_col, "switch_prob")],
                  by = c(x$id_col, x$time_col), all.x = TRUE)
  d_plot$switch_prob[is.na(d_plot$switch_prob)] <- 0


  if (type == "boxplot") {
    p <- ggplot2::ggplot() +
      ggplot2::geom_boxplot(
        data = d_plot,
        ggplot2::aes(x = factor(.data[[x$time_col]]), y = switch_prob),
        width = 0.5, outlier.size = 0.7, ...
      ) +
      ggplot2::geom_errorbar(
        data = by_k,
        ggplot2::aes(x = factor(k_idx), ymin = lo95, ymax = hi95),
        width = 0.15
      ) +
      ggplot2::geom_point(
        data = by_k,
        ggplot2::aes(x = factor(k_idx), y = mean),
        size = 1.2
      ) +
      ggplot2::geom_line(
        data = by_k,
        ggplot2::aes(x = as.numeric(factor(k_idx)), y = mean, group = 1)
      ) +
      ggplot2::labs(
        title = if (x$scale == "mass")
          "Switching probability over intervals (P(init at k))"
        else
          "Switching hazard over intervals (P(init at k | at risk))",
        x = "Interval (k_idx)", y = "Probability"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    p <- ggplot2::ggplot(by_k, ggplot2::aes(x = k_idx, y = mean)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lo95, ymax = hi95), alpha = 0.15) +
      ggplot2::geom_line() +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::scale_x_continuous(breaks = unique(by_k$k_idx)) +
      ggplot2::labs(
        title = if (x$scale == "mass")
          "Switching probability over intervals (mean -+ 95% bands)"
        else
          "Switching hazard over intervals (mean -+ 95% bands)",
        x = "Interval (k_idx)", y = "Probability"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12)
  }

  print(p)
  invisible(p)
}
