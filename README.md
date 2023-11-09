
<!-- README.md is generated from README.Rmd. Please edit that file -->

# A data-generation perspective on the use of multiple Fine–Gray models

<!-- badges: start -->
<!-- See hema review example-->
<!-- badges: end -->

**Authors**: Edouard F. Bonneville, Liesbeth C. de Wreede, and Hein
Putter

## Abstract

Studies considering competing risks will often aim to estimate the
cumulative incidence functions conditional on an individual’s baseline
characteristics. While the Fine–Gray subdistribution hazard model is
tailor-made for analysing only one of the competing events, it may still
be used in settings where multiple competing events are of scientific
interest, where it is specified for each cause in turn. In this work, we
provide an overview of data-generating mechanisms where proportional
subdistribution hazards hold for at least one cause. We use these to
motivate why the use of multiple Fine–Gray models should be reconsidered
in favour of better alternatives.

## Usage

Two main files are of interest:

- [`manuscript-figures.R`](./manuscript-figures.R) - R script to
  reproducte the figures in the manuscript.
- [`helpers.R`](./helpers.R) - helper functions that
  `manuscript-figures.R` relies on to compute true cumulative incidence
  function, cause-specific hazards, and subdistribution hazards, for
  both events.

The R environment can be reproduced using
[`{renv}`](https://github.com/rstudio/renv/), by calling

``` r
renv::restore()
```

## Ideas for further extension

- Generalize functions to accommodate user-defined hazard or cumulative
  distribution functions, and vectorize covariate input.

- Create a shiny app to automate plot creation for different
  parametrizations, and show when the total failure probability exceeds
  one/cause-specific hazards become negative.

- Make functions allowing user to actually simulate data (instead of
  just checking the true values).
