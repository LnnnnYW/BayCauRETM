---
title: 'BayCauRETM: R package for Bayesian Causal Inference for Recurrent Event Outcomes'
tags:
- R
- Bayesian inference
- Causal inference
- Recurrent events
- timing misalignment
- survival analysis
date: "11 September 2025"
output: pdf_document
authors:
- name: Yuqin Wang
  orcid: "0009-0003-8345-9318"
  affiliation: '1'
- name: Keming Zhang
  orcid: "0009-0001-5495-0058"
  affiliation: '1'
- name: Arman Oganisian
  orcid: "0000-0002-0437-4611"
  affiliation: '1'
  corresponding: true
bibliography: paper.bib
affiliations:
- name: Department of Biostatistics, Brown University, Providence, RI, United States
  index: 1
---

# Summary

Observational studies are often conducted to estimate the effect of different medical treatments on the rate of a recurrent event outcome within a specified follow-up window. Recurrent events are outcomes that can occur multiple times (e.g. multiple hospitalization events, relapses, etc.) the follow-up window. Causal analysis of recurrent event processes is complicated for several reasons: 1) they are typically observed jointly with a terminal event process such as death. This is strictly dependent on the recurrent event process as, for example, the occurrence of death precludes subsequent recurrence of the outcome of interest. 2) Both the event counts and terminal event process are unobserved if a subject drops out of a study before the end of follow-up - a phenomenon known as censoring. 3) Patients may initiate a treatment at different times during the follow-up, meaning that for each treatment there are as many possible strategies as there are possible initiation times. Finally, 4) since treatments are not randomized in observational studies, formal causal methods such as g-computation are needed to adjust for observed confounders recorded in the data.

This paper presents `BayCauRETM`, an `R` package for estimating causal effects of different treatment initiation strategies on a recurrent event outcome in the presence of death and censoring. The user specifies a given treatment initiation time and supplies a `data.frame` with confounders to adjust for, as well as columns for censoring, death, and recurrent event counts. The user then specifies a model for the terminal and recurrent event process using standard `R` regression syntax. With these user inputs, `BayCauRETM` can run a causal adjustment procedure and output adjust expected event rates within a follow-up under the specified treatment initiation time. `BayCauRETM` also provides functions for model diagnostics and visualization.

Intended users include statisticians, epidemiologists, and health-services researchers analyzing observational data.

# Statement of need

Standard software implementing methods for time-to-event and recurrent event data remain valuable for descriptive purposes, but generally do not target causal estimands while dealing with complexities 1-4 described in the Summary section above [@ghoshlin2002; @schaubel2010; @janvin2024]. @Oganisian2024 developed Bayesian statistical methods that accommodate these complexities and conducted a thorough simulation-based validation of these methods. However, due to the focus on methodological development and validation, only proof-of-concept replication code was provided along with the paper. There is need for user-friendly, off-the-shelf software with readable help files that can implement the methods developed in @Oganisian2024.`BayCauRETM` fills this methodological and practical gap by operationalizing the Bayesian approach of @Oganisian2024 in `R`. `BayCauRETM` is designed to have a syntax familiar to base `R` users and which mirror more standard regression functions like `lm()` and `glm`. It comes complete with extensive help files accessed via the `?` command in `R`. Thus, `BayCauRETM` provides the first user-friendly software for analyzing complex recurrent event data in the face of complexities 1)-4) described in the Summary section above.

# Data structure, model, and outputs

In this section, we provide an overview of the expected input data structure, models that are run under-the-hood, and expected outputs. We refer readers to @Oganisian2024 for methodological details.

### Data structure and preprocessing

The package expects longitudinal data that is arranged in an `R` `data.frame` object in long, person-interval format. That is, for some follow-up time, $\tau$, the follow-up window $[0,\tau)$ is partitioned into $K$ equal-length intervals $I_k=[\tau_{k-1},\tau_k)$ for $k=1,\dots,K$ with $\tau_0=0$ and $\tau_K=\tau$. Each row is a patient-interval of time and a subject has as many rows as they have intervals for which they are at-risk (i.e. uncensored and alive).

The required `data.frame` variables are: a subject identifier, an interval index representing $k$, the treatment indicator (0 until the interval in which treatment is initiated, then 1 thereafter), the interval-specific count of recurrent events, a terminal-event indicator (0 up to death and 1 from the first interval after death). Optional fields include baseline covariates, lagged history (e.g., a one-interval lag of the event count). This structure aligns event attribution and at-risk time with treatment history, and it mirrors typical data-collection schemes used in observational data sets like medical insurance claims and or electronic health records (EHR) sources.

### Model specification 

At each row, the `data.frame` should have a monotone binary indicator of death at the start of interval $k$, which we denote with $T_k$. It should also include a binary monotone indicator of whether treatment has been initiated by the end of interval $k$, $A_k$. Finally, it should contain, at each row, the count of the number of events in that interval, $Y_k$, as well as a copy of a set of baseline (i.e. time-constant) covariates, which we denote by $L\in\mathcal{L}$.

Let $a(s)=(\underbrace{0,\dots,0}_{s-1},1,\dots,1)$ define a hypothetical strategy in which we intervene to initiate treatment for everyone in the target population at interval $s\in\{1,2,\dots, K+1\}$. We define potential outcomes $T_k^{a(s)}$ and $Y_k^{a(s)}$ which represent whether a patient would have been dead in interval $k$ and the number of events they would have experienced in interval $k$, respectively, had they - possibly counter to the fact - initiated treatment at interval $s$. This package provides inference for the the difference in average potential incidence rates over the follow-up window under two different initiation strategies $$
\Delta(s,s') =
\mathbb{E}\!\left[\frac{\sum_{k=1}^K Y_k^{a(s)}}{K-\sum_{k=1}^K T_k^{a(s)}}\right]
-
\mathbb{E}\!\left[\frac{\sum_{k=1}^K Y_k^{a^{(s')}}}{K-\sum_{k=1}^K T_k^{a^{(s')}}}\right],
$$

The package runs a pair of discrete-time models conditional on shared treatment and covariate terms:\
1. A discrete-time hazard model for the terminal event (e.g. death) that models death at a given interval conditional on survival up to that interval:\
$$\lambda_k(a_k,\bar y_{k-1},l)=\Pr\!\big(T_k=1\mid T_{k-1}=0,\,a_k,\bar y_{k-1}, l \big)$$ 2. A distribution for the number of event occurrences in a given interval conditional on survival through that interval:\
$$ f(y_k\mid a_k,\bar y_{k-1},l)=\Pr\!\big(Y_k=y_k\mid T_k=0,\,a_k,\bar y_{k-1},l\big) $$ Here, $f(y_k\mid a_k,\bar y_{k-1},\ell)$ represents the Poisson probability mass function with conditional mean/intensity of the event-count $\mu_k(a_k, \bar y_{k-1}, l) = E[Y_k\mid A_k,\bar Y_{k-1},L]$. Together, these two models multiply to form a joint model for the terminal and recurrent event occurrence at a given interval.

The functions in `BayCauRETM` implement the following models for the hazard and intensity, respectively, $$
\text{logit}\,\lambda_k(a_k,\bar y_{k-1},l)=\beta_{0k}+l^\top\beta_L+y_{k-1}\beta_Y+\beta_A a_k,\qquad
\log \mu_k(a_k, \bar y_{k-1}, l)=\theta_{0k}+l^\top\theta_L+y_{k-1}\theta_Y+\theta_A a_k,
$$ The time-varying intercepts $\{\beta_{0k}\}$ and $\{\theta_{0k}\}$ parameterize the baseline hazard and event intensity, respectively. They are assigned a first-order autoregressive (AR1) smoothing prior to improve stability as at-risk counts diminish in later intervals. See @Oganisian2024 for more details.

### Posterior inference and g-computation

`BayCauRETM` conducts full posterior inference for the models described in the previous section. Since the posterior is not available in closed form, it back-ends to Stan [@Stan2017] via the `rstan` package to obtain draws of the model parameters from their joint posterior. Stan is a probabilistic programming language (PPL) that implements cutting edge Hamiltonian Monte Carlo methods to obtain these draws.

Using a given draw of the model parameters obtained via Stan , `BayCauRETM` computes a posterior draw of $\Delta(s, s')$. It does this by approximating both expectations in $\Delta(s, s')$ via Monte Carlo simulations from the model. That is, we simulate the joint death-recurrent event process many times under $a(s)$ and compute the average incidence rate across these simulations. We do the same under $a(s')$ and take the difference to obtain a posterior draw of $\Delta(s, s')$. Doing this across many posterior parameter draws yields a set of posterior draws for $\Delta(s, s')$ as described by @Oganisian2024. The mean of these draws is used as a posterior point estimate while a 95% credible interval is formed by taking the 2.5th and 97.5th percentiles as endpoints.

### Software workflow

Standard workflow using the `BayCauRETM` involves the following functions:

1.  **Model fitting** via `fit_causal_recur()`, which parses formulas for the death hazard and count mean, selects the count family, and configures gAR(1) priors and MCMC settings.

2.  **Posterior g-computation** with `g_computation()` to evaluate a vector of switch times $s$.

3.  **Reporting and diagnostics** using `result_summary_table()` and `plot()` to produce tables/figures, `mcmc_diagnosis()` for $\hat R$, effective sample sizes and traceplots, and `propensity_score_diagnostics()` / `switching_probability_summary()` to help assess validity of the required positivity assumptions.

Sensible defaults are provided, while expert users can override these defaults. Outputs are returned as tidy data frames, facilitating integration into literate analysis workflows.

Detailed usage and example results are available on [GitHub](https://github.com/LnnnnYW/BayCauRETM) (see the [demo PDF](https://github.com/LnnnnYW/BayCauRETM/blob/master/inst/demo_code/demo.pdf)).

# Acknowledgements

This work was partially funded by the Patient Centered Outcomes Research Institute (PCORI) Contract ME-2023C1-31348.

# References
