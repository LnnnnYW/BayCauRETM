% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_causal_recur.R
\docType{class}
\name{causal_recur_fit}
\alias{causal_recur_fit}
\alias{fit_causal_recur}
\alias{print.causal_recur_fit}
\alias{summary.causal_recur_fit}
\alias{plot.causal_recur_fit}
\title{Fit Bayesian Causal Recurrent - and Terminal - Event Model}
\usage{
fit_causal_recur(
  data,
  K,
  x_cols = NULL,
  formula_T,
  formula_Y,
  prior,
  num_chains = 4,
  iter = 2000,
  stan_model_file = NULL,
  control = list(adapt_delta = 0.95, max_treedepth = 15),
  verbose = TRUE
)

\method{print}{causal_recur_fit}(x, ...)

\method{summary}{causal_recur_fit}(
  object,
  pars_to_report = c("beta_Y", "beta_A", "gamma_Y", "gamma_A"),
  ...
)

\method{plot}{causal_recur_fit}(x, ...)
}
\arguments{
\item{data}{A \strong{long-format} \code{data.frame} that already contains the standard
analysis columns: \code{pat_id}, \code{k_idx}, \code{Y_obs}, \code{T_obs}, \code{A}, and \emph{optionally}
extra covariates.  The function will run a light version of
\code{preprocess_data()} to add \code{Y_prev} and \code{T_prev} (and to fill in any missing
intervals) but \strong{will not} rename columns.}

\item{K}{Integer. Total number of discrete intervals in the study.}

\item{x_cols}{Character vector of additional (static or time - varying)
covariate names to keep; \code{NULL} if none.}

\item{formula_T}{A formula for the terminal-event (death) sub‑model, e.g.
\code{T_obs ~ Y_prev + A + k_idx}.}

\item{formula_Y}{A formula for the recurrent-count sub‑model, e.g.
\code{Y_obs ~ Y_prev + A + k_idx}.}

\item{prior}{Named list of gAR(1) hyperparameters with elements
\code{eta_beta}, \code{sigma_beta}, \code{rho_beta}, \code{eta_gamma}, \code{sigma_gamma},
\code{rho_gamma}.}

\item{num_chains}{Integer. Number of MCMC chains (default \code{4}).}

\item{iter}{Integer. Total iterations \emph{per} chain including warm-up (default
\code{2000}).}

\item{stan_model_file}{Optional path to a pre-compiled Stan model
(\code{.stan} -> \verb{*.rds}).  If \code{NULL}, the package-internal model is used.}

\item{control}{List passed to Stan sampling (see \strong{cmdstanr} docs).}

\item{verbose}{Logical. Print progress messages (default \code{TRUE}).}

\item{x}{an object used to select a method.}

\item{...}{further arguments passed to or from other methods.}

\item{object}{Object returned by this function (used by summary method).}

\item{pars_to_report}{Character vector of parameter names to include in the summary/plot S3 methods.}
}
\value{
An object of class \code{causal_recur_fit} (list) with elements
\code{stan_fit}, \code{data_preprocessed}, \code{n_pat}, \code{K}, \code{design_info}, \code{prior}, and
\code{stan_data_list}.
}
\description{
Fit a discrete - time Bayesian model for recurrent event counts and a terminal
event (death), using gAR(1) smoothing priors on the time - varying intercepts.
}
\details{
Internally, the function first calls \code{\link[=preprocess_data]{preprocess_data()}} (which assumes the
input already uses the \verb{_obs} / \verb{_prev} naming convention), then builds the
required design matrices, compiles or loads the Stan model, and finally runs
MCMC sampling via \strong{rstan}.
}
\section{Methods (by generic)}{
\itemize{
\item \code{print(causal_recur_fit)}: Print a brief summary

\item \code{summary(causal_recur_fit)}: Summarise posterior estimates

\item \code{plot(causal_recur_fit)}: Prompt user to run diagnostics

}}
\examples{
df <- data.frame(
  pat_id = rep(1:2, each = 2),
  k_idx  = rep(1:2, 2),
  Y_obs  = rpois(4, 1),
  T_obs  = rbinom(4, 1, 0.2),
  A      = rbinom(4, 1, 0.5)
)
prior <- list(
  eta_beta    = 0, sigma_beta  = 1, rho_beta   = 0.5,
  eta_gamma   = 0, sigma_gamma = 1, rho_gamma  = 0.5
)
\dontrun{
fit <- fit_causal_recur(
  data       = df,
  K          = 2,
  formula_T  = T_obs ~ Y_prev + A + k_idx,
  formula_Y  = Y_obs ~ Y_prev + A + k_idx,
  prior      = prior,
  num_chains = 1,
  iter       = 100
)
print(fit)
summary(fit)
plot(fit)
}

}
