% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/result_summary_table.R
\name{result_summary_table}
\alias{result_summary_table}
\alias{print.result_summary_table}
\alias{summary.result_summary_table}
\title{Summarise posterior parameters & g-computation contrasts}
\usage{
result_summary_table(
  fit_out,
  gcomp_out,
  pars_to_report = c("beta1", "theta1", "theta_lag"),
  s_vec = NULL,
  filter_pars = NULL,
  sort_by = "Mean",
  sort_desc = TRUE,
  format = "data.frame",
  export_file = NULL
)

\method{print}{result_summary_table}(x, ...)

\method{summary}{result_summary_table}(object, ...)
}
\arguments{
\item{fit_out}{Output list of \code{\link[=fit_causal_recur]{fit_causal_recur()}}.}

\item{gcomp_out}{Output list of \code{\link[=g_computation]{g_computation()}}.}

\item{pars_to_report}{Character; Stan parameter names to keep.}

\item{s_vec}{Integer vector of start intervals to keep.
\code{NULL} = all available.}

\item{filter_pars}{Optional tidy‐expression to filter the parameter table
(e.g. \code{Mean > 0 & Rhat < 1.1}).}

\item{sort_by}{Column to sort the parameter table by.}

\item{sort_desc}{Logical; descending sort?}

\item{format}{\code{"data.frame"}, \code{"kable"}, or \code{"gt"}.}

\item{export_file}{Optional path; if given, writes CSV / XLSX export.}

\item{x}{an object used to select a method.}

\item{...}{Ignored; present for S3 consistency.}

\item{object}{A \code{result_summary_table} object (output of
\code{result_summary_table()}).}
}
\value{
A \strong{result_summary_table} object (see \code{print()} method).

\code{object}, invisibly.
}
\description{
Extracts posterior summaries from a Stan fit (\strong{fit_out}) and
causal contrasts Δ(s, K+1) from a \strong{gcomp_out}, then merges them into
a single object that can be printed, exported, or rendered with
\strong{knitr} / \strong{gt}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(result_summary_table)}: Formatted console/knitr/gt output of the two summary tables.

\item \code{summary(result_summary_table)}: Alias for print.

}}
\examples{
\dontrun{
res <- result_summary_table(fit_out, gcomp_out,
                            pars_to_report = c("beta_Y","beta_A"),
                            s_vec = 1:5, format = "kable")
print(res)   # pretty table
}

}
\seealso{
\code{\link[=print.result_summary_table]{print.result_summary_table()}}, \code{\link[=g_computation]{g_computation()}},
\code{\link[=plot_posterior_causal_contrast_static]{plot_posterior_causal_contrast_static()}}.
}
