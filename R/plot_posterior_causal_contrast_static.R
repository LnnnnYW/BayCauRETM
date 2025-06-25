#' Posterior Causal Contrast delta(s, K + 1) - Static Wrapper
#'
#' @description
#' Produce a static **ggplot2** figure of the posterior causal contrast
#' \eqn{\Delta(s,\,K+1)} against treatment–start interval \eqn{s}, for one or
#' multiple scenarios produced by \code{\link{g_computation}}.
#'
#' @param contrast_list Either a single \code{gcomp_out} object or **named** list
#'   of such objects to compare several scenarios side-by-side.
#' @param s_vec Integer vector of start intervals to include.  If \code{NULL}
#'   (default) all intervals contained in \code{$delta} are plotted.
#' @param theme_fn A ggplot2 theme function (default
#'   \code{ggplot2::theme_minimal}).
#' @param line_size Numeric; width of mean line (default 1).
#' @param ribbon_alpha Numeric in 0 to 1; transparency of credible-interval
#'   ribbons (default 0.3).
#' @param show_points Logical; draw points at each \eqn{s}? (default \code{TRUE}).
#' @param label_points Logical; annotate mean value at points? (default \code{FALSE}).
#' @param ref_line Numeric or \code{NULL}; if numeric, draws a horizontal dashed
#'   line at that \emph{y}.  Typical use \code{ref_line = 0}.
#'
#' @return A \link[ggplot2]{ggplot} object.
#' @export
#'
#' @examples
#' \dontrun{
#' # single scenario
#' p <- plot_posterior_causal_contrast_static(gcomp_out, s_vec = 1:K, ref_line = 0)
#' print(p)
#'
#' # two scenarios
#' p2 <- plot_posterior_causal_contrast_static(
#'   list(adj = g_adj, unadj = g_unadj), s_vec = 1:K, show_points = TRUE
#' )
#' }
#'
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot aes geom_hline geom_ribbon geom_line geom_point
#'   geom_text labs theme position_dodge theme_minimal


plot_posterior_causal_contrast_static <- function(contrast_list,
                                                  s_vec        = NULL,
                                                  theme_fn     = ggplot2::theme_minimal,
                                                  ribbon_alpha = 0.3,
                                                  line_size    = 1,
                                                  show_points  = TRUE,
                                                  label_points = FALSE,
                                                  ref_line     = NULL) {

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
                  ggplot2::aes(x = s, y = mean,
                               colour = scenario, fill = scenario)) +
    { if (!is.null(ref_line))
      ggplot2::geom_hline(yintercept = ref_line, linetype = "dashed") } +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         alpha = ribbon_alpha, colour = NA,
                         position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::geom_line(linewidth = line_size,
                       position = ggplot2::position_dodge(width = 0.5)) +
    { if (show_points)
      ggplot2::geom_point(size = 2,
                          position = ggplot2::position_dodge(width = 0.5)) } +
    { if (label_points)
      ggplot2::geom_text(ggplot2::aes(label = signif(mean, 3)),
                         vjust = -0.5,
                         position = ggplot2::position_dodge(width = 0.5)) } +
    ggplot2::labs(
      x = "Treatment-start interval s",
      y = expression(Delta(s, K+1)),
      title = expression("Posterior causal contrast " ~ Delta(s, K+1))
    ) +
    theme_fn() +
    ggplot2::theme(legend.position = "bottom")
}
