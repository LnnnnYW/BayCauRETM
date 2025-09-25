---
title: 'BayCauRETM: Bayesian Causal Inference for Recurrent Events with Timing Misalignment
  in R'
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

Many studies track recurrent events (e.g., hospitalizations) alongside a terminal event(e.g., death). In practice, patients often start treatment after they first become clinically eligible. This timing misalignment makes naive treated/untreated or before/after comparisons misleading because person-time and events accrued before initiation may be wrongly attributed to the treated group.

`BayCauRETM` estimate the causal contrast in cumulative recurrent events under a start-by-$s$ intervention versus a control group (never-treated), accounting for a terminal event via a discrete-time g-formula and posterior-predictive g-computation. The package implements a joint, discrete-time Bayesian framework that models the terminal event and recurrent counts with shared covariates/history, uses light temporal smoothing to stabilize late follow-up, and applies posterior-predictive g-computation to estimate counterfactual cumulative-event measures under start-by-$s$ strategies. Typical users include biostatisticians, epidemiologists, and health-services researchers analyzing observational cohorts or pragmatic trials.

# Statement of need

Real-world observational cohorts and pragmatic trials often exhibit a lag between clinical eligibility and actual treatment initiation, while a terminal event (e.g., death) truncates subsequent event accrual. Under these features, naive comparisons, such as treated vs. untreated or pre- vs. post-initiation, can misattribute person-time and events generated before initiation to the treated group, inducing timing-related (e.g., immortal-time) biases. Consequently, many standard tools for recurrent events and survival, frailty or multi-state frameworks, or marginal approaches, remain valuable for description and prediction but generally do not target a decision-relevant causal estimand for start-time policies when a competing terminal event is present and timing is misaligned [@ghoshlin2002; @schaubel2010; @janvin2024]. Analysts therefore need software that (i) encodes a clearly decision-focused estimand (e.g.,start-by-$s$ vs. never-treat), (ii) aligns event attribution and risk sets with the observed treatment history across discrete intervals, thereby avoiding misclassification of exposure and person-time, and (iii) returns posterior summaries that transparently quantify uncertainty.

`BayCauRETM` fills this methodological and practical gap by operationalizing [@Oganisian2024] in R. The package formalizes start-by-$s$ estimands on a discrete-time grid that mirrors common data-collection schedules; establishes identification via a discrete-time g-formula under standard causal assumptions; and delivers decision-oriented contrasts through posterior-predictive g-computation. In practice, users fit a joint Bayesian model for the terminal-event hazard and recurrent-event counts with shared covariates and selected history, with light temporal smoothing to stabilize sparse late follow-up. Counterfactual cumulative-event trajectories are then simulated under user-specified start-time strategies (e.g., start-by-$s$ vs. never-treat), and summarized as tidy tables and publication-ready plots with coherent uncertainty intervals, providing an end-to-end, reproducible workflow for causal start-time questions in the presence of a competing terminal event.

# Data structure, model, and outputs

### Data structure and preprocessing

We work with long-format panels where each row represents one subject and one discrete interval after eligibility. Required fields are: a subject identifier, an interval index, the treatment indicator (0 until the interval in which treatment starts, then 1 thereafter), the interval-specific count of recurrent events, a terminal-event indicator (0 up to death and 1 from the first interval after death). Optional fields include baseline covariateslagged history (e.g., a one-interval lag of the count) and an exposure-time offset if intervals differ marginally in length after data cleaning. This structure aligns event attribution and risk time with treatment history, and it mirrors typical data-collection schemes used in administrative or EHR sources.

### Model specification

We adopt a joint, discrete-time formulation with two components that share treatment, covariates, and selected history terms:\
(i) a hazard model for the terminal event, and\
(ii) a count model for recurrent events observed while at risk.\
Baseline-by-interval effects are parameterized piecewise and lightly smoothed across time to stabilize estimation in sparse late follow-up. Bayesian inference proceeds via Hamiltonian Monte Carlo in Stan, interfaced through `rstan`.

Partition $[0,\tau)$ into $K$ equal-length intervals $I_k=[\tau_{k-1},\tau_k)$ for $k=1,\dots,K$ with $\tau_0=0$. Let - $T_k=\mathbb{I}(\tilde U\le \tau_{k-1},\,\delta_U=1)$ denote the death indicator at the start of interval $k$ (monotone: zeros until death, then ones); - $A_k=\mathbb{I}(\tilde W\le \tau_k,\,A=1)$ be the treatment indicator by the end of interval $k$ (zeros until switching, then ones); - $Y_k=\sum_{j=1}^J \mathbb{I}(V_j\in I_k)$ be the count of recurrent events in interval $k$.

For a switch-time strategy $a^{(s)}=(\underbrace{0,\dots,0}_{s-1},1,\dots,1)$, define potential outcomes $T_k^{a^{(s)}}$ and $Y_k^{a^{(s)}}$. Our primary estimand is the difference in average potential incidence rates $$
\Delta(s,s') =
\mathbb{E}\!\left[\frac{\sum_{k=1}^K Y_k^{a^{(s)}}}{K-\sum_{k=1}^K T_k^{a^{(s)}}}\right]
-
\mathbb{E}\!\left[\frac{\sum_{k=1}^K Y_k^{a^{(s')}}}{K-\sum_{k=1}^K T_k^{a^{(s')}}}\right],
$$ which compares treatment strategies that start by interval $s$ versus $s'$ (including “never treat” as a reference).

Under sequential ignorability and positivity, the joint law of $(\bar Y^{a^{(s)}},\bar T^{a^{(s)}})$ is identified via a g-formula constructed from\
(i) the discrete-time death hazard $$
\lambda_k(a_k,\bar y_{k-1},\ell)=\Pr\!\big(T_k=1\mid T_{k-1}=C_{k-1}=0,\,a_k,\bar y_{k-1},\ell\big),
$$ and\
(ii) a model for interval counts $$
f(y_k\mid a_k,\bar y_{k-1},\ell)=\Pr\!\big(Y_k=y_k\mid T_k=C_k=0,\,a_k,\bar y_{k-1},\ell\big),
$$ with baseline covariates $L\in\mathcal{L}$. We use semiparametric mean models $$
\text{logit}\,\lambda_k=\beta_{0k}+L^\top\beta_L+Y_{k-1}\beta_Y+\beta_A a_k,\qquad
\log \mu_k=\theta_{0k}+L^\top\theta_L+Y_{k-1}\theta_Y+\theta_A a_k,
$$ where $\mu_k=\mathbb{E}(Y_k\mid\cdot)$ follows either a Poisson or a negative-binomial family. Time-varying intercepts $\{\beta_{0k}\}$ and $\{\theta_{0k}\}$ are assigned gAR(1) priors to smooth hazards and intensities over $k$, improving stability when information is limited late in follow-up. This setup respects the composite nature of rate denominators, since death terminates at-risk time.

### Estimation and posterior g-computation

Posterior g-computation simulates the intervention distributions implied by each switch-time strategy $a^{(s)}$. Concretely, we draw parameters from the joint posterior, generate potential $\{T_k^{a^{(s)}},Y_k^{a^{(s)}}\}_{k=1}^K$ forward under the specified $a^{(s)}$ with observed covariates $L$, and evaluate $\Delta(s,s')$ as a functional of the resulting posterior predictive draws. This delivers full posterior summaries (means, intervals) for contrasts of direct decision relevance.

### Software workflow

The package exposes a streamlined path:

1\. **Preprocessing** with `preprocess_data()` to create long-format panels, compute lags such as $Y_{k-1}$, and index intervals consistently.

2\. **Model fitting** via `fit_causal_recur()`, which parses formulas for the death hazard and count mean, selects the count family, and configures gAR(1) priors and MCMC settings.

3\. **Posterior g-computation** with `g_computation()` to evaluate a vector of switch times $s$.

4\. **Reporting and diagnostics** using `result_summary_table()` and `plot()` to produce tables/figures, `mcmc_diagnosis()` for $\hat R$, effective sample sizes and traceplots, and `propensity_score_diagnostics()` / `switching_probability_summary()` to visualize key identification assumptions.

Sensible defaults are provided, while expert users retain control over priors, number of chains, warmup, and iteration counts. Outputs are returned as tidy data frames, facilitating integration into literate analysis workflows.

Detailed usage and example results are available on [GitHub](https://github.com/LnnnnYW/BayCauRETM) (see the [demo PDF](https://github.com/LnnnnYW/BayCauRETM/blob/master/inst/demo_code/demo.pdf)).

# Acknowledgements

This work was partially funded by the Patient Centered Outcomes Research Institute (PCORI) Contract ME-2023C1-31348. Thanks to all the supports.

# References
