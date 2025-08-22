#' Posterior Causal Contrast delta(s, K + 1) - Interactive Wrapper
#'
#' @description
#' Convenience wrapper: builds the static plot via
#' \code{\link{plot_posterior_causal_contrast_static}}, optionally saves it, and
#' optionally converts it to an interactive Plotly object.
#'
#' @param contrast_list List of posterior contrasts from
#' @param s_vec Optional vector of treatment levels to plot.
#' @param theme_fn Function to apply a ggplot2 theme (default is
#'  \code{theme_minimal}).
#' @param point_size Numeric; size of points in the plot (default 3).
#' @param error_width Numeric; width of error bars (default 0.15).
#' @param ref_line Optional reference line value (default NULL, no line).
#' @param interactive Logical; if TRUE return a Plotly object.
#' @param save_file Optional filename to save static plot (PNG/PDF, etc.).
#' @param width,height,dpi Device settings when saving.
#' @param ... Additional arguments passed to
#'
#' @return ggplot or plotly object.
#'
#' @examples
#' \dontrun{
#' # interactive plot in RStudio viewer
#' plot_posterior_causal_contrast_interactive(gcomp_out, interactive = TRUE)
#'
#' # save static PNG
#' plot_posterior_causal_contrast_interactive(
#'   gcomp_out, save_file = "delta_static.png", s_vec = 1:K
#' )
#' }
#'
#' @importFrom ggplot2 ggsave
#' @importFrom plotly ggplotly
#' @export

plot_posterior_causal_contrast_interactive <- function(contrast_list,
                                                       s_vec        = NULL,
                                                       theme_fn     = ggplot2::theme_minimal,
                                                       point_size   = 3,
                                                       error_width  = .15,
                                                       ref_line     = NULL,
                                                       interactive  = FALSE,
                                                       save_file    = NULL,
                                                       width        = 8,
                                                       height       = 5,
                                                       dpi          = 300,
                                                       ...) {

  p <- plot_posterior_causal_contrast_static(
    contrast_list = contrast_list,
    s_vec         = s_vec,
    theme_fn      = theme_fn,
    point_size    = point_size,
    error_width   = error_width,
    ref_line      = ref_line,
    ...
  )

  if (!is.null(save_file))
    ggplot2::ggsave(save_file, plot = p,
                    width = width, height = height, dpi = dpi)

  if (interactive) {
    p <- p + ggplot2::labs(
      y = "Delta(s, K+1)",
      title = "Posterior causal contrast delta(s, K+1)"
    )

    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Package 'plotly' required for interactive output. Please install it.")
    return(plotly::ggplotly(p))
  }

  p
}


