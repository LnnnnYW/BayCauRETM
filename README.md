BayCauRETM
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/LnnnnYW/BayCauRETM/actions/workflows/r.yml/badge.svg)](https://github.com/LnnnnYW/BayCauRETM/actions/workflows/r.yml)
[![Coverage
Status](https://coveralls.io/repos/github/LnnnnYW/BayCauRETM/badge.svg)](https://coveralls.io/github/LnnnnYW/BayCauRETM)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)
<!-- badges: end -->

**BayCauRETM** (Bayesian Causal Recurrent-Event and Terminal-Event
Modeling)  
provides a *fully Bayesian* workflow for **discrete-time causal
inference**  
with recurrent-event counts jointly modeled with a terminal event.

- Canonical long-format preprocessing with automated lags  
- Stan-based joint modeling (`fit_causal_recur()`) with gAR(1) priors  
- Posterior *g*-computation for alternative treatment-start strategies  
- Diagnostics: MCMC convergence, propensity-score, switching-hazard  
- Publication-ready tables & plots

## Installation

``` r
# Install development version from GitHub
# install.packages("pak")
pak::pak("LnnnnYW/BayCauRETM")

# Or with devtools:
# devtools::install_github("LnnnnYW/BayCauRETM")
```

## Dependencies

The following R packages are required for `BayCauRETM`:

- [rstan](https://cran.r-project.org/package=rstan)
- [dplyr](https://cran.r-project.org/package=dplyr)
- [ggplot2](https://cran.r-project.org/package=ggplot2)
- [bayesplot](https://cran.r-project.org/package=bayesplot)
- [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/00Index.html)
  *(base package)*
- [magrittr](https://cran.r-project.org/package=magrittr)
- [rlang](https://cran.r-project.org/package=rlang)
- [gt](https://cran.r-project.org/package=gt)
- [knitr](https://cran.r-project.org/package=knitr)
- [plotly](https://cran.r-project.org/package=plotly)
- [tidyr](https://cran.r-project.org/package=tidyr)
- [writexl](https://cran.r-project.org/package=writexl)

These are listed in the `Imports` field of the `DESCRIPTION` file and
will be installed automatically with the package.

## Documentation and Example

The [paper](https://academic.oup.com/biometrics/article/80/4/ujae145/7914699) associated with this package contains the statistical details of the model as well as a detailed walk-through demonstration. 

The code for demostration in the paper is available in the folder [inst/demo_code](https://github.com/LnnnnYW/BayCauRETM/tree/master/inst/demo_code).


## Reporting Issues

If you encounter bugs or have feature suggestions, please  
open an [issue on
GitHub](https://github.com/LnnnnYW/BayCauRETM/issues).  
Pull requests are always welcome!

## Citation

If you use **BayCauRETM**, please cite:

## Contact

The corresponding package author are Yuqin Wang (email:
<yuqin_wang@brown.edu>) and Arman Oganisian (email:
<arman_oganisian@brown.edu>).

> Tip: For reproducibility, knit this README from a `README.Rmd`
> source  
> using `devtools::build_readme()`.
