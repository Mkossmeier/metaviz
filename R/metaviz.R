#' metaviz: Forest Plots, Funnel Plots, and Visual Funnel Plot Inference for Meta-Analysis
#'
#' The package metaviz is a collection of functions to create visually
#' appealing and information-rich plots of meta-analytic data using ggplot2.
#' Functions to create several variants of forest plots (\code{viz_forest}),
#' funnel plots (\code{viz_funnel}, \code{viz_sunset}), and to conduct visual
#' inference with funnel plots (\code{funnelinf}) are provided.
#'
#' @section Forest plots (\code{viz_forest}):
#' Several different types and variants of forest plots can be created. This includes
#' classic forest plots, subgroup forest plots, cumulative summary forest plots, and leave-one-out sensitivity
#' forest plots. In addition, the function allows to individually label and color studies and to
#' align tables with furhter study-level and summary-level information.
#'
#' In addition to traditional forest plots, rainforest plots as well as thick forest plots can be used.
#' Rainforest and thick forest plots are two recently proposed variants and enhancements of
#' the classic forest plot. Both variants visually emphasize large studies
#' (with short confidence intervals and more weight in the meta-analysis), while small studies
#' (with wide confidence intervals and less weight in the meta-analysis) are visually less dominant.
#' For further details see \code{help(viz_forest)}, \code{help(viz_rainforest)}, and \code{help(viz_thickforest)}.
#'
#' @section Funnel plots (\code{viz_funnel}, \code{viz_sunset}):
#' Numerous different funnel plot variants can be created. Options for several graphical augmentations
#' (e.g., confidence, significance, and additional evidence contours; choice of the ordinate; showing
#' study subgroups), and different statistical information displayed are provided (Egger's regression line,
#' and imputed studies by, as well as the adjusted summary effect from, the trim-and-fill method).
#' Further details and references can be found in the corresponding help file (\code{help(viz_funnel)}).
#'
#' Moreover, a novel variant of the funnel plot is introduced which displays the power of studies to detect an effect
#' of interest (e.g., the meta-analytic summary effect) using a two-sided Wald test. This sunset (power-enhanced)
#' funnel plot uses color-coded regions and a second y axis to visualize study-level power and can help to
#' critically examine the evidentiality and credibility of a set of studies. For further details see \code{help(viz_sunset)}.
#'
#' @section Visual inference with funnel plots (\code{funnelinf}):
#' Funnel plots are widely used in meta-analysis to assess small study effects as potential indicator of publication bias.
#' Visual inference can help to improve the objectivity and validity of conclusions based on funnel plot examinations by
#' guarding the meta-analyst from interpreting patterns in the funnel plot that might be perfectly plausible by chance.
#' Only if the funnel plot showing the real data is distinguishable from simultaneously presented
#' null funnel plots showing data simulated under the null hypothesis, conclusion based on visually inspecting
#' the real-data funnel plot might be warranted. The function \code{funnelinf} provides numerous tailored
#' options to conduct visual inference with the funnel plot graph in the context of meta-analysis. See \code{help(funnelinf)} for further details and relevant references.
#'
#' @section Meta-analytic example datasets:
#' Four different example datasets from published meta-analyses are distributed
#' with the package: \itemize{
#' \item Two datasets for meta-analysis with standardized mean differences (\code{mozart}, \code{homeopath})
#' \item One dataset for meta-analysis with correlation coefficients (\code{brainvol})
#' \item One dataset for meta-analysis with dichotomous outcome data (\code{exrehab}).
#' }
#' More details and corresponding references can be found in the respective help files
#' (\code{help(mozart)}, \code{help(homeopath)}, \code{help(brainvol)}, \code{help(exrehab)}).
#'
#' @docType package
#' @name metaviz-package
#' @aliases metaviz
NULL
