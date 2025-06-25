# File: R/result_summary_table.R

#' Summarise posterior parameters & g-computation contrasts
#'
#' @description
#' Extracts posterior summaries from a Stan fit (**fit_out**) and
#' causal contrasts Δ(s, K+1) from a **gcomp_out**, then merges them into
#' a single object that can be printed, exported, or rendered with
#' **knitr** / **gt**.
#'
#' @param fit_out   Output list of [fit_causal_recur()].
#' @param gcomp_out Output list of [g_computation()].
#' @param pars_to_report Character; Stan parameter names to keep.
#' @param s_vec     Integer vector of start intervals to keep.
#'                  `NULL` = all available.
#' @param filter_pars Optional tidy‐expression to filter the parameter table
#'                    (e.g. `Mean > 0 & Rhat < 1.1`).
#' @param sort_by   Column to sort the parameter table by.
#' @param sort_desc Logical; descending sort?
#' @param format    `"data.frame"`, `"kable"`, or `"gt"`.
#' @param export_file Optional path; if given, writes CSV / XLSX export.
#' @return A **result_summary_table** object (see `print()` method).
#'
#' @seealso [print.result_summary_table()], [g_computation()],
#'   [plot_posterior_causal_contrast_static()].
#'
#' @examples
#' \dontrun{
#' res <- result_summary_table(fit_out, gcomp_out,
#'                             pars_to_report = c("beta_Y","beta_A"),
#'                             s_vec = 1:5, format = "kable")
#' print(res)   # pretty table
#' }
#'
#' @importFrom rstan summary
#' @importFrom dplyr filter arrange desc
#' @importFrom knitr kable
#' @importFrom gt gt tab_header
#' @importFrom writexl write_xlsx
#' @importFrom utils write.csv
#' @export


result_summary_table <- function(fit_out,
                                 gcomp_out,
                                 pars_to_report = c("beta_Y","beta_A","gamma_Y","gamma_A"),
                                 s_vec         = NULL,
                                 filter_pars   = NULL,
                                 sort_by       = "Mean",
                                 sort_desc     = TRUE,
                                 format        = "data.frame",
                                 export_file   = NULL) {

  if (!is.null(fit_out$stan_fit) && inherits(fit_out$stan_fit, "stanfit")) {
    ss <- rstan::summary(fit_out$stan_fit, pars = pars_to_report)$summary
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
  } else {                               # empty fallback
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
  }

  delta_list <- gcomp_out$delta
  if (!is.null(s_vec)) delta_list <- delta_list[paste0("s=", s_vec)]

  if (length(delta_list) > 0) {
    df_delta <- do.call(rbind, lapply(names(delta_list), function(nm) {
      x <- delta_list[[nm]]
      data.frame(
        s        = as.integer(sub("^s=", "", nm)),
        Mean     = x$mean,
        `2.5%`   = x$CI_lower,
        `97.5%`  = x$CI_upper,
        CI_width = x$CI_upper - x$CI_lower,
        stringsAsFactors = FALSE
      )
    }))
    df_delta <- df_delta[order(df_delta$s), , drop = FALSE]
  } else {
    df_delta <- data.frame(
      s        = integer(0),
      Mean     = numeric(0),
      `2.5%`   = numeric(0),
      `97.5%`  = numeric(0),
      CI_width = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  if (!is.null(filter_pars) && nrow(df_par) > 0) {
    df_par <- dplyr::filter(df_par, !!rlang::enquo(filter_pars))
  }
  if (!is.null(sort_by) && sort_by %in% names(df_par) && nrow(df_par) > 0) {
    df_par <- if (sort_desc) dplyr::arrange(df_par, dplyr::desc(.data[[sort_by]]))
    else           dplyr::arrange(df_par, .data[[sort_by]])
  }

  export_path <- NULL
  if (!is.null(export_file)) {
    ext <- tolower(tools::file_ext(export_file))
    if (ext %in% c("xlsx","xls")) {
      writexl::write_xlsx(list(parameters = df_par, delta = df_delta),
                          path = export_file)
      export_path <- export_file
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
    param_table <- knitr::kable(df_par, caption = "Posterior Parameters")
    delta_table <- knitr::kable(df_delta, caption = "delta(s, K+1)")
  } else if (fmt == "gt") {
    param_table <- gt::gt(df_par) %>% gt::tab_header("Posterior Parameters")
    delta_table <- gt::gt(df_delta) %>% gt::tab_header("delta(s, K+1)")
  }

  out <- list(
    param_summary = df_par,
    delta_summary = df_delta,
    param_table   = if (is.null(param_table)) df_par else param_table,
    delta_table   = if (is.null(delta_table)) df_delta else delta_table,
    export_file   = export_path
  )
  class(out) <- c("result_summary_table", "baycar_results")

  invisible(out)
}

#' @title Print a result_summary_table
#' @description Formatted console/knitr/gt output of the two summary tables.
#' @param x A **result_summary_table** object.
#' @param ... Ignored.
#' @return *x*, invisibly.
#' @method print result_summary_table
#' @export
print.result_summary_table <- function(x, ...) {
  if (inherits(x$param_table, "knitr_kable") ||
      inherits(x$param_table, "gt_tbl")) {
    print(x$param_table)
  } else {
    cat("----- Posterior Parameters -----\n")
    print(x$param_summary)
  }
  cat("\n")
  if (inherits(x$delta_table, "knitr_kable") ||
      inherits(x$delta_table, "gt_tbl")) {
    print(x$delta_table)
  } else {
    cat("----- g-computation delta(s, K+1) -----\n")
    print(x$delta_summary)
  }
  invisible(x)
}

#' @title Summary for result_summary_table
#' @description Alias for [print.result_summary_table()].
#' @inheritParams print.result_summary_table
#' @param object A \code{result_summary_table} object.
#' @return *object*, invisibly.
#' @method summary result_summary_table
#' @export
summary.result_summary_table <- function(object, ...) {
  print(object)
  invisible(object)
}



