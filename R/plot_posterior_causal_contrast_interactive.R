#' Posterior causal contrast delta(s, K+1) - interactive wrapper
#'
#' @description
#' Build the static plot via [plot_posterior_causal_contrast_static()], optionally
#' save the static ggplot to disk, and optionally return an interactive Plotly object.
#'
#' @inheritParams plot_posterior_causal_contrast_static
#' @param interactive Logical. If TRUE, convert the ggplot to a Plotly object with
#'   `plotly::ggplotly()` and return it.
#' @param save_file Optional file path to save the static ggplot (PNG/PDF, etc.).
#' @param width,height,dpi Device settings used by `ggplot2::ggsave()` when saving.
#'
#' @return A ggplot object when `interactive = FALSE`; otherwise a Plotly object.
#'
#' @examples
#' \dontrun{
#' # Interactive plot in the viewer
#' plot_posterior_causal_contrast_interactive(gcomp_out, interactive = TRUE)
#'
#' # Save a static PNG (the returned value is still a ggplot)
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


