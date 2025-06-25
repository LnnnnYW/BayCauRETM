#' Posterior Causal Contrast delta(s, K + 1) - Interactive Wrapper
#'
#' @description
#' Convenience wrapper: builds the static plot via
#' \code{\link{plot_posterior_causal_contrast_static}}, optionally saves it, and
#' optionally converts it to an interactive Plotly object.
#'
#' @inheritParams plot_posterior_causal_contrast_static
#' @param interactive Logical; if TRUE return a Plotly object.
#' @param save_file Optional filename to save static plot (PNG/PDF, etc.).
#' @param line_size Numeric; width of the mean line (default 1).
#' @param width,height,dpi Device settings when saving.
#' @param ribbon_alpha Numeric in 0 to 1; transparency of credible-interval
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
                                                       line_size    = 1,
                                                       ribbon_alpha = 0.3,
                                                       show_points  = TRUE,
                                                       label_points = FALSE,
                                                       ref_line     = NULL,
                                                       interactive  = FALSE,
                                                       save_file    = NULL,
                                                       width        = 8,
                                                       height       = 5,
                                                       dpi          = 300) {

  p <- plot_posterior_causal_contrast_static(
    contrast_list = contrast_list,
    s_vec         = s_vec,
    theme_fn      = theme_fn,
    line_size     = line_size,
    ribbon_alpha  = ribbon_alpha,
    show_points   = show_points,
    label_points  = label_points,
    ref_line      = ref_line
  )

  if (!is.null(save_file))
    ggplot2::ggsave(save_file, plot = p, width = width, height = height, dpi = dpi)

  if (interactive) {
    p <- p + ggplot2::labs(
      y = "Delta(s, K+1)",
      title = "Posterior causal contrast delta(s, K+1)"
    )

    if (!requireNamespace("plotly", quietly = TRUE))
      stop("Package 'plotly' required for interactive output.")
    return(plotly::ggplotly(p))
  }

  p
}

