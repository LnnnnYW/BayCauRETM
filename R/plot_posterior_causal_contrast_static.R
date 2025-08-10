#' Posterior causal contrast delta(s, K+1) - static wrapper
#'
#' @description
#' Produce a static **ggplot2** figure of the posterior causal contrast
#' `delta(s, K+1)` against the treatment-start interval `s` for one or multiple
#' scenarios produced by [g_computation()].
#'
#' @param contrast_list Either a single `gcomp_out` object or a **named** list
#'   of such objects to compare multiple scenarios side by side.
#' @param s_vec Integer vector of start intervals to include. If `NULL`
#'   (default), all intervals contained in `$delta` are plotted.
#' @param theme_fn A ggplot2 theme function (default `ggplot2::theme_minimal`).
#' @param point_size Numeric; size of point markers (default `3`).
#' @param error_width Numeric; width of error bars (default `0.15`).
#' @param ref_line Numeric or `NULL`; if numeric, draws a horizontal dashed
#'   line at that y-value (typical use `ref_line = 0`).
#' @param ... Additional styling options forwarded to the underlying geoms and
#'   helpers (e.g., `line_size`, `ribbon_alpha`, `show_points`, `label_points`).
#'
#' @return A [ggplot2::ggplot] object.
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_posterior_causal_contrast_static(
#'   contrast_list = gcomp_out,
#'   s_vec = 1:K,
#'   ref_line = 0
#' )
#' print(p)
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_hline geom_ribbon geom_line geom_point
#'   geom_text labs theme position_dodge theme_minimal



plot_posterior_causal_contrast_static <- function(contrast_list,
                                                  s_vec        = NULL,
                                                  theme_fn     = ggplot2::theme_minimal,
                                                  point_size   = 3,
                                                  error_width  = .15,
                                                  ref_line     = NULL,
                                                  ...) {

  if (!is.list(contrast_list) ||
      (length(contrast_list) > 0 && !inherits(contrast_list[[1]], "gcomp_out"))) {
    contrast_list <- list(default = contrast_list)
  }
  if (is.null(names(contrast_list)) || any(names(contrast_list) == "")) {
    names(contrast_list) <- paste0("Scenario", seq_along(contrast_list))
  }

  df_plot <- dplyr::bind_rows(lapply(names(contrast_list), function(nm) {
    delta_list <- contrast_list[[nm]]$delta
    if (!is.null(s_vec)) delta_list <- delta_list[paste0("s=", s_vec)]
    do.call(rbind, lapply(names(delta_list), function(k) {
      x <- delta_list[[k]]
      data.frame(
        scenario = nm,
        s        = as.integer(sub("^s=", "", k)),
        mean     = x$mean,
        lower    = x$CI_lower,
        upper    = x$CI_upper,
        stringsAsFactors = FALSE
      )
    }))
  }), .id = NULL)
  df_plot$scenario <- factor(df_plot$scenario, levels = names(contrast_list))

  ggplot2::ggplot(df_plot,
                  ggplot2::aes(x = factor(s), y = mean,
                               colour = scenario, group = scenario)) +
    { if (!is.null(ref_line))
      ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed") } +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower, ymax = upper),
                           width = error_width,
                           position = ggplot2::position_dodge(width = .6)) +
    ggplot2::geom_point(size = point_size,
                        position = ggplot2::position_dodge(width = .6)) +
    ggplot2::labs(
      x = "Treatment-start interval s",
      y = expression(Delta(s, K+1)),
      title = expression("Posterior causal contrast " ~ Delta(s, K+1))
    ) +
    theme_fn() +
    ggplot2::theme(legend.position = "bottom")
}

