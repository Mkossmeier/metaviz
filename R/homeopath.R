#' Example data for a meta-analysis of standardized mean differences: Homeopathic treatment vs. Placebo
#'
#' A dataset consisting of 54 randomized control trials (RCTs) from a published meta-analysis
#' (Mathie et al. 2017) comparing homeopathic treatment with placebo. See reference for further details.
#'
#' @format A data frame with 54 rows and 4 variables:
#' \describe{
#'   \item{name}{short name of each RCT}
#'   \item{year}{publication year of each RCT}
#'   \item{d}{observed effect size (Standardized mean differences). Negative effect sizes indicate superiority of homeopathic treatment.}
#'   \item{se}{standard error of observed effect size}
#'   }
#' @references Mathie, R. T., Ramparsad, N., Legg, L. A., Clausen, J., Moss, S.,
#'   Davidson, J. R., ... McConnachie, A. (2017). Randomised, double-blind, placebo-controlled
#'   trials of non-individualised homeopathic treatment: Systematic review and meta-analysis.
#'   \emph{Systematic Reviews}, \emph{6}, 63.
"homeopath"
