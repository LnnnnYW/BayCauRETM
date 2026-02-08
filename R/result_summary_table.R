#' Summarize posterior parameters and g-computation contrasts
#'
#' @description
#' Extract posterior summaries from a Stan fit (`fit_out`) and causal contrasts
#' `delta(s, K+1)` from a `gcomp_out`, then merge them into a single result that
#' can be printed, exported, or rendered with knitr/gt.
#'
#' @param fit_out Output list from [fit_causal_recur()].
#' @param gcomp_out Output list from [g_computation()].
#' @param pars_to_report Character vector of Stan parameter base names to report.
#'   Vector-valued parameters include all indices; only names present in the
#'   fitted object are used. Default:
#'   `c("beta0", "beta1", "theta0", "theta1", "theta_lag")`.
#' @param s_vec Integer vector selecting which treatment-start intervals s
#'   (discrete times 1..K) to report, e.g. c(1, 3, 5). If NULL, all s values
#'   available in gcomp_out$delta are used; any values not present there are
#'   silently ignored.
#' @param format One of `"data.frame"`, `"kable"`, or `"gt"`. When the required
#'   package is unavailable, output falls back to `"data.frame"`.
#' @param export_file Optional file path. If supplied, writes CSV/XLSX export.
#'
#' @return An object of class `result_summary_table` (also inheriting
#'   `"baycar_results"`) with components:
#'   - `param_summary`: data frame with posterior summaries (`Parameter`, `Mean`,
#'     `2.5%`, `97.5%`, `Rhat`, `n_eff`, `MCSE`, `CI_width`). Parameter names are
#'     mapped to user-friendly labels when available.
#'   - `delta_summary`: data frame with `s`, `Mean`, `2.5%`, `97.5%`, `CI_width`
#'     for `delta(s, K+1)`, ordered by `s`.
#'   - `export_file`: path(s) of the file(s) written when `export_file` is given,
#'     otherwise `NULL`.
#'
#' @seealso [print.result_summary_table()], [g_computation()],
#'   [plot_posterior_causal_contrast_static()].
#'
#' @details
#' Posterior summaries are obtained from `rstan::summary()` for the requested
#' parameter blocks. Credible intervals reported in the tables are the equal-tailed
#' 2.5% and 97.5% quantiles returned by Stan. If `s_vec` is supplied, the delta
#' table is subset to those intervals (matching names like `"s=1"`, `"s=2"`, ...).
#'
#' @examples
#' \dontrun{
#' res <- result_summary_table(
#'   fit_out, gcomp_out,
#'   pars_to_report = c("beta_Y","beta_A"),
#'   s_vec = 1:5, format = "kable"
#' )
#' print(res)
#' }
#' @importFrom rstan summary
#' @importFrom knitr kable
#' @importFrom gt gt tab_header
#' @importFrom writexl write_xlsx
#' @importFrom utils write.csv
#' @export

result_summary_table <- function(fit_out,
                                 gcomp_out,
                                 pars_to_report = c("beta0","beta1","theta0","theta1","theta_lag"),
                                 s_vec         = NULL,
                                 format        = "data.frame",
                                 export_file   = NULL) {

  #check input types
  if (!inherits(fit_out$stan_fit, "stanfit"))
    stop("fit_out must contain be a 'stanfit' object")
  if (!is.list(gcomp_out) || is.null(gcomp_out$delta) || !is.list(gcomp_out$delta))
    stop("gcomp_out must be the output of g_computation")
  if (!is.character(pars_to_report) || length(pars_to_report) == 0)
    stop("'pars_to_report' must be a non-empty character vector")
  if (!is.null(s_vec) && (!is.numeric(s_vec) || any(s_vec <= 0)))
    stop("'s_vec' must be NULL or a numeric vector of positive integers")
  if (!is.character(format) || length(format) != 1)
    stop("'format' must be a single character string")
  if (!is.null(export_file) && (!is.character(export_file) || length(export_file) != 1))
    stop("'export_file' must be NULL or a single character string")

  df_par <- data.frame(
    Parameter = character(0),
    Mean      = numeric(0),
    `2.5%`    = numeric(0),
    `97.5%`   = numeric(0),
    Rhat      = numeric(0),
    n_eff     = numeric(0),
    MCSE      = numeric(0),
    CI_width  = numeric(0),
    stringsAsFactors = FALSE
  )

  if (!is.null(fit_out$stan_fit) && inherits(fit_out$stan_fit, "stanfit")) {
    all_rows <- rstan::summary(fit_out$stan_fit)$summary
    all_pars <- rownames(all_rows)
    keep     <- unlist(lapply(pars_to_report, function(p)
      grep(paste0("^", gsub("\\[.*\\]$", "", p), "(\\[.*\\])?$"), all_pars)))
    if (length(keep) > 0) {
      ss <- all_rows[keep, , drop = FALSE]
      df_par <- data.frame(
        Parameter = rownames(ss),
        Mean      = ss[, "mean"],
        `2.5%`    = ss[, "2.5%"],
        `97.5%`   = ss[, "97.5%"],
        Rhat      = ss[, "Rhat"],
        n_eff     = ss[, "n_eff"],
        MCSE      = ss[, "se_mean"],
        CI_width  = ss[, "97.5%"] - ss[, "2.5%"],
        row.names = NULL, stringsAsFactors = FALSE
      )
      df_par$Parameter <- .map_param_names(fit_out, df_par$Parameter)
    }
  }

  delta_list <- gcomp_out$delta
  if (!is.null(s_vec) && length(delta_list) > 0) {
    delta_list <- delta_list[paste0("s=", s_vec)]
    delta_list <- delta_list[!vapply(delta_list, is.null, logical(1))]
  }

  if (length(delta_list) > 0) {
    df_delta <- do.call(rbind, lapply(names(delta_list), function(nm) {
      x <- delta_list[[nm]]
      data.frame(
        s        = suppressWarnings(as.integer(sub("^s=", "", nm))),
        Mean     = as.numeric(x$mean),
        Lower    = as.numeric(x$CI_lower),
        Upper    = as.numeric(x$CI_upper),
        CI_width = as.numeric(x$CI_upper - x$CI_lower),
        stringsAsFactors = FALSE
      )
    }))
    df_delta <- df_delta[order(df_delta$s), , drop = FALSE]
  } else {
    df_delta <- data.frame(
      s        = integer(0),
      Mean     = numeric(0),
      Lower    = numeric(0),
      Upper    = numeric(0),
      CI_width = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  export_path <- NULL
  if (!is.null(export_file)) {
    ext <- tolower(tools::file_ext(export_file))
    if (ext %in% c("xlsx","xls")) {
      if (!requireNamespace("writexl", quietly = TRUE)) {
        warning("Package 'writexl' not installed; exporting CSVs instead.")
        write.csv(df_par, export_file, row.names = FALSE)
        dfile <- sub("(\\.[^.]+)?$", "_delta.csv", export_file)
        write.csv(df_delta, dfile, row.names = FALSE)
        export_path <- c(parameters = export_file, delta = dfile)
      } else {
        writexl::write_xlsx(list(parameters = df_par, delta = df_delta),
                            path = export_file)
        export_path <- export_file
      }
    } else {
      write.csv(df_par, export_file, row.names = FALSE)
      dfile <- sub("(\\.[^.]+)?$", "_delta.csv", export_file)
      write.csv(df_delta, dfile, row.names = FALSE)
      export_path <- c(parameters = export_file, delta = dfile)
    }
  }

  fmt <- match.arg(format, c("data.frame","kable","gt"))
  param_table <- delta_table <- NULL
  if (fmt == "kable") {
    if (!requireNamespace("knitr", quietly = TRUE))
      warning("Package 'knitr' not installed; falling back to data.frame.")
    else {
      param_table <- knitr::kable(df_par, caption = "Posterior Parameters")
      delta_table <- knitr::kable(df_delta, caption = "delta(s, K+1)")
    }
  } else if (fmt == "gt") {
    if (!requireNamespace("gt", quietly = TRUE))
      warning("Package 'gt' not installed; falling back to data.frame.")
    else {
      param_table <- gt::gt(df_par) |> gt::tab_header(title = "Posterior Parameters")
      delta_table <- gt::gt(df_delta) |> gt::tab_header(title = "delta(s, K+1)")
    }
  }

  out <- list(
    param_summary = df_par,
    delta_summary = df_delta,
    export_file   = export_path
  )
  class(out) <- c("result_summary_table", "baycar_results")
  invisible(out)

}


#' @describeIn result_summary_table Formatted console/knitr/gt output of the two summary tables.
#' @param x A `result_summary_table` object (output of `result_summary_table()`).
#' @param ... Ignored; present for S3 consistency.
#' @return `x`, invisibly.
#' @method print result_summary_table
#' @export

print.result_summary_table <- function(x, ...) {
  use_kable <- requireNamespace("knitr", quietly = TRUE)

  cat("----- Posterior Parameters -----\n")
  if (use_kable) {
    print(knitr::kable(x$param_summary, format = "pipe", digits = 6))
  } else {
    print(x$param_summary, row.names = FALSE)
  }

  cat("\n----- g-computation delta(s, K+1) -----\n")
  if (use_kable) {
    print(knitr::kable(x$delta_summary, format = "pipe", digits = 6))
  } else {
    print(x$delta_summary, row.names = FALSE)
  }
  invisible(x)
}

#' @describeIn result_summary_table Alias for `print()`.
#' @param object A `result_summary_table` object.
#' @param ... Ignored; present for S3 consistency.
#' @return `object`, invisibly.
#' @method summary result_summary_table
#' @export

summary.result_summary_table <- function(object, ...) {
  print(object)
  invisible(object)
}

