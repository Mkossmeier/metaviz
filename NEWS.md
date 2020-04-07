
## Package **metaviz** version 0.1.0

  - February 6, 2017 first release on CRAN:
    <https://CRAN.R-project.org/package=metaviz>

## Package **metaviz** version 0.1.1

  - March 14, 2017: Added a new function for visual funnel plot
    inference: `funnelinf`. Also available as shiny app:
    <https://metaviz.shinyapps.io/funnelinf_app/>
  - June 29, 2017: Version 0.1.1 submitted to CRAN

## Package **metaviz** version 0.2

  - March 16, 2018: Version 0.2 submitted to CRAN
  - Added viz\_funnel, a new function to create funnel plot variants
    (e.g., additional evidence contour enhanced funnel plots)
  - Added viz\_forest, a new function to create forest plot variants
    (classic, thick forest, and rainforest plots)
  - Now includes an additional example dataset of a meta-analysis with
    dichotomous outcomes.

## Package **metaviz** version 0.3

  - January 14, 2019: Version 0.3 submitted to CRAN
  - Added `viz_sunset`, a dedicated function for sunset (power-enhanced)
    funnel plots. Also available as shiny app:
    <https://metaviz.shinyapps.io/sunset/>
  - Now includes an additional example dataset of a meta-analysis with
    standardized mean differences.

## Package **metaviz** version 0.3.1

  - April 7th, 2020: Version 0.3.1 submitted to CRAN
  - Fixed problem of up-side-down funnel plots (`viz_funnel`,
    `viz_sunset`, `funnelinf`) with y\_axis = “se”, which seemed to
    occur with newer versions of **ggplot2** installed.
  - Function `viz_funnel` now allows to customize funnel plot contours
    for random effects models via the `contours_type` argument.
  - Changed behavior of `viz_forest`, `viz_funnel` and `viz_sunset`,
    such that when output of function `rma.uni` from package **metafor**
    is used as input, then the `method` argument is now extracted from
    the `rma.uni` object.
