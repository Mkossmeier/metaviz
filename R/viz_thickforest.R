#'Thick forest plots for meta-analyses
#'
#'Creates a thick forest plot, an enhanced variant of the forest plot.
#'
#'The thick forest plot was proposed by Schild and Voracek (2015) as a variant and
#'enhancement of classic forest plots. Thick forest plots use rectengular error bars
#'instead of traditional lines to display confidence intervals (width of the error bar), as well as the relative
#'meta-analytic weight (height of the error bar) of each study. In addtion, study and summary level
#'point estimates are depicted clearly by a specific symbol.
#'
#'Thick forest plots have the following advantages, as compared to classic forest plots:
#'
#'\enumerate{
#'\item Using the height of bars proportional to the (relative) meta-analytic weight
#'causes small studies (with wide confidence intervals and less weight in the meta-analysis) to
#'be visually less dominant.
#'
#'\item In classic forest plots, it is often hard to depict the magnitude of
#'point estimates to a reasonable degree of accuracy, especially for studies
#'with large meta-analytic weights and correspondingly large plotting symbols
#'(commonly squares). Specific symbols within the thick forest plot improve the
#'visualization of study point estimates. }
#'
#'Note that for subgroup analysis the height of each error bar is scaled by the weight of each study within the subgroup divided by
#'the sum of the weights of all studies irrespective of subgroup. Therefore, with subgroups present, the overall impression of error
#'bar heights within a given subgroup compared to other subgroups conveys information about the relative precision of the meta-analyitc
#'estimate within the subgroup.
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}.
#'@param group factor indicating the subgroup of each study to plot a subgroup thick forest plot. Has to be in the same order than \code{x}.
#'@param summary_symbol logical scalar. Should the meta-analytic summary result(s) be computed using \code{method} and shown?
#'@param method Which method should be used to compute the study weights and summary effect(s)?
#'  Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the
#'  DerSimonian-Laird method to estimate  \eqn{\tau^2}{tau squared}).
#'@param study_labels a vector with names/identifiers to annotate each study in the thick forest plot.
#'  Has to be in the same order than \code{x}.
#'@param summary_label a name to annotate the summary effect. If a subgroup
#'  analysis is plotted, \code{summary_label} should be a vector with a name for each
#'  subgroup summary effect, arranged in the order of the levels of \code{group}.
#'@param study_table a data.frame with addtional study-level variables which should be shown in an aligned table.
#'  Has to be in the same order than \code{x}. See vignette('metaviz').
#'@param summary_table a data.frame with addtional summary-level information shown in an aligned table.
#'  If \code{group} is supplied, \code{summary_table} must have a row for each subgroup
#'  summary effect, arranged in the order of the levels of \code{group}. See vignette('metaviz').
#'@param table_headers character vector. Headers for each column of the aligned table if \code{study_table} and/or \code{summary_table} is supplied.
#'  By default the names of \code{study_table}.
#'@param annotate_CI logical scalar. Should the effect size and confidence interval values be annotated?
#'@param confidence_level the confidence level for the plotted confidence bars.
#'@param col character specifying the color used for the error bars and summary diamond.
#'@param tick_col character specifying the color used for the ticks indicating the point estimates.
#'@param text_size numeric value. Values larger than 1 lead to larger text size,
#'  values smaller than 1 to smaller text size than the default.
#'@param xlab character label of the x axis. Also used for the header of the aligned table if \code{annotate_CI} is TRUE.
#'@param x_limit numeric vector of length 2 with the limits (minimum, maximum) of the x axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficents using \code{tanh}. See vignette('metaviz').
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@references Schild, A. H., & Voracek, M. (2015). Finding your way out of the
#'  forest without a trail of bread crumbs: Development and evaluation of two
#'  novel displays of forest plots. \emph{Research Synthesis Methods}, \emph{6},
#'  74-86.
#'@return A thick forest plot is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Plotting a thick forest plot using the mozart data (for details, see help(mozart)):
#' viz_thickforest(x = mozart[, c("d", "se")],
#' study_labels = mozart[, "study_name"], xlab = "Cohen d")
#'
#' # Visualizing a subgroup analysis of published and unpublished studies
#' viz_thickforest(x = mozart[, c("d", "se")], group = mozart[, "unpublished"],
#' study_labels = mozart[, "study_name"],
#' summary_label = c("Summary (published)", "Summary (unpublished)"),
#' xlab = "Cohen d")
#'
#' # Showing additional information in aligned tables. Log risk ratios are labeled
#' # in their original metric (risk ratios) on the x axis.
#' viz_thickforest(x = exrehab[, c("logrr", "logrr_se")],
#' study_labels = exrehab[, "study_name"],
#' annotate_CI = TRUE, xlab = "RR", x_trans_function = exp,
#' study_table = data.frame(
#' eventsT = paste(exrehab$ai, "/", exrehab$ai + exrehab$bi, sep = ""),
#' eventsC = paste(exrehab$ci, "/", exrehab$ci + exrehab$di, sep = "")),
#' summary_table = data.frame(
#' eventsT = paste(sum(exrehab$ai), "/", sum(exrehab$ai + exrehab$bi), sep = ""),
#' eventsC = paste(sum(exrehab$ci), "/", sum(exrehab$ci + exrehab$di), sep = "")))
#'@export
viz_thickforest <- function(x, group = NULL, summary_symbol = TRUE, method = "FE",
                           study_labels = NULL, summary_label = NULL,
                           study_table = NULL, summary_table = NULL, table_headers = NULL, annotate_CI = FALSE,
                           confidence_level = 0.95, col = "Blues", tick_col = "firebrick",
                           text_size = 3, xlab = "Effect", x_limit = NULL,
                           x_trans_function = NULL, x_breaks = NULL) {

  viz_forest(x, group = group, summary_symbol = summary_symbol, method = method,
             study_labels = study_labels, summary_label = summary_label,
             study_table = study_table, summary_table = summary_table,
             table_headers = table_headers, annotate_CI = annotate_CI,
             confidence_level = confidence_level, col = col,  tick_col = tick_col,
             text_size = text_size, xlab = xlab, x_limit = x_limit,
             x_trans_function = x_trans_function, x_breaks = NULL, type = "thick")
}
