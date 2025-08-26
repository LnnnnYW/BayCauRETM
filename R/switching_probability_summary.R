#' Compute switching-probability summary
#'
#' @description
#' Summarize how often subjects switch (initiate) treatment across follow-up
#' intervals. Two complementary scales are supported:
#'
#' - mass (default): the marginal probability of initiating at interval \(k\),
#'  i.e., the proportion of the cohort that switches exactly at \(k\).
#' - hazard: the conditional probability of switching at interval \(k\)
#'  given the subject has not switched before \(k\).
#'
#' In discrete time these satisfy \(m(k) = S(k-1)\,h(k)\), where \(h(k)\) is the
#' hazard and \(S(k-1)\) the probability of not yet switched before \(k\).
#'
#' @param df A `data.frame` (e.g., from `preprocess_data()`) containing at least
#'   a subject id, a time index, and a binary treatment indicator.
#' @param scale Character, `"mass"` (default) or `"hazard"`. See Details.
#' @inheritParams base::print
#'
#' @return An object of class `switching_summary`, a list with components:
#' \describe{
#'   \item{by_k}{A data frame with one row per interval containing:
#'     k_idx, n (at-risk sample size), mean_p (mean
#'     probability on the chosen scale), Wald 95% limits lo95/hi95,
#'     and q25}/q50/q75 for descriptive dispersion.}
#'   \item{by_id_k}{Per-subject indicators used for boxplots:
#'     for scale="hazard", it is the 0/1 indicator of switching at \(k\)
#'     among the at-risk set; for scale="mass", it is the 0/1 indicator of
#'     “first switch occurs at \(k\)”.}
#'   \item{id_col,time_col,treat_col,death_col}{Resolved column names.}
#'   \item{scale}{The scale used ("mass" or "hazard").}
#' }
#'
#' @details
#' Let \(A_{i,k}\) be subject \(i\)'s treatment at interval \(k\) and define
#' \(R_{i,k}=\mathbf{1}\{\text{subject }i\text{ is alive and }k>1\}\). The event
#' of switching at \(k\) is \(E_{i,k}=\mathbf{1}\{A_{i,k}\neq A_{i,k-1}\}\).
#' \itemize{
#'   \item \strong{Hazard scale}: \(h(k)=\mathbb{E}(E_{i,k}\mid R_{i,k}=1)\),
#'         estimated by the at-risk mean of \(E_{i,k}\) at each \(k\).
#'   \item \strong{Mass scale}: \(m(k)=\mathbb{P}(\text{first switch at }k)\),
#'         estimated as the cohort mean of the first-switch indicator; this is
#'         numerically equivalent to \(S(k\!-\!1)\,h(k)\) with \(S(k\!-\!1)\) the
#'         product over \(j<k\) of \((1-h(j))\).
#' }
#' The summary is purely descriptive (nonparametric) and does not fit a model;
#' it is intended to contextualize causal estimates and reveal time-structured
#' policies (e.g., delayed initiation or systematic stopping).
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
  scale <- match.arg(scale)

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
    stop("Data must contain an ID column (pat_id/id) and a treatment column (A/Ak).")

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
      alive_to_km1 = if (!is.na(death_col))
        c(1L, cumprod(1L - .data[[death_col]][-length(.data[[death_col]])])) else 1L,
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

  rhs_terms <- c(.bt(fac_time), .bt(cov_map))
  form <- stats::as.formula(paste0("init_k ~ ", paste(rhs_terms, collapse = " + ")))

  mf <- stats::model.frame(form, data = d[at_risk_rows, , drop = FALSE],
                           na.action = stats::na.omit)
  if (!nrow(mf)) stop("No rows remain after removing NA for hazard model.")
  rows_fit <- as.integer(rownames(mf))

  mod <- stats::glm(form, data = mf,
                    family = stats::binomial,
                    control = stats::glm.control(maxit = 100))

  haz_hat <- rep(NA_real_, nrow(d))
  haz_hat[rows_fit] <- stats::predict(mod, type = "response")

  d$hazard <- haz_hat

  if (scale == "mass") {
    # S_{k-1} from predicted hazards (per id, cumulative product up to k-1)
    d <- d |>
      dplyr::group_by(.data[[id_col]]) |>
      dplyr::mutate(
        one_minus_h = dplyr::if_else(is.finite(hazard), 1 - hazard, 1),
        S_km1 = dplyr::lag(cumprod(one_minus_h), default = 1),
        switch_prob = S_km1 * hazard
      ) |>
      dplyr::ungroup()
  } else {
    d$switch_prob <- d$hazard
  }


  d_plot <- d[at_risk_rows, , drop = FALSE]

  by_k <- d_plot |>
    dplyr::group_by(.data[[time_col]]) |>
    dplyr::summarise(
      n      = dplyr::n(),
      mean   = mean(switch_prob, na.rm = TRUE),
      sd_    = stats::sd(switch_prob, na.rm = TRUE),
      se     = ifelse(n > 0, sd_ / sqrt(n), NA_real_),
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

#' @describeIn switching_probability_summary
#' Print the interval-level summary table.
#' @export

print.switching_summary <- function(x, ...) {
  cat("Switching probability summary by interval (", x$scale, "):\n", sep = "")
  print(x$by_k)
  invisible(x$by_k)
}

#' @describeIn switching_probability_summary
#' @param object A `switching_summary` object.
#' @method summary switching_summary
#' @export

summary.switching_summary <- function(object, ...) {
  print(object, ...)
}

#' @describeIn switching_probability_summary
#' Plot switching summaries across intervals.
#'
#' @param type One of `"boxplot"` (per-interval boxplots of per-subject
#'   indicators with mean ± 95% Wald bars overlaid) or `"line"` (mean curve with
#'   95% ribbon). Defaults to `"boxplot"`.
#' @param ... Passed to `ggplot2::geom_boxplot()` when `type="boxplot"`.
#' @export

plot.switching_summary <- function(x, type = c("boxplot","line"), ...) {
  type <- match.arg(type)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required for plotting.")

  by_k   <- x$by_k
  d_plot <- x$df_haz[x$df_haz$at_risk == 1L, , drop = FALSE]

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
          "Switching probability over intervals (mean ± 95% bands)"
        else
          "Switching hazard over intervals (mean ± 95% bands)",
        x = "Interval (k_idx)", y = "Probability"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::theme_minimal(base_size = 12)
  }

  print(p)
  invisible(p)
}
