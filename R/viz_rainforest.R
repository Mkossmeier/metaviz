#'Rainforest plots for meta-analyses
#'
#'Creates a rainforest plot, an enhanced variant of the forest plot.
#'
#'Rainforest plots were proposed by Schild and Voracek (2015) as a variant and
#'enhancement of classic forest plots. Rainforest plots use (log-)likelihood drops
#'to depict study level results (for details, see Barrowman &
#'Myers, 2003). \code{viz_rainforest} assumes normality of effect sizes to
#'construct these (log-)likelihood drops. The width of each such raindrop is identical to the
#'(Wald-type) confidence interval. For a given (log-)likelihood raindrop,
#'the height can be interpreted as the plausibility (i.e., (log-)likelihood) of different true
#'values given the observed estimate. Moreover, the height of each raindrop is scaled with respect,
#'to its relative meta-analytic weight considering all studies. Therefore, visually comparing heights of
#'different raindrops highlights the relative importance within the meta-analytis.
#'In addtion, color shading is utilized to further visualize statistical uncertainty, as suggested by Jackson (2008).
#'Finally, study and summary level point estimates are depicted clearly by a specific symbol.
#'
#'Rainforest plots have the following advantages, as compared to classic forest plots:
#'
#'\enumerate{ \item The width of the likelihood raindrops corresponds to the
#'confidence intervals, as also shown in the classic forest plot. In addition,
#'for each likelihood drop the height (and color shading) visualizes the
#'plausibility of true values given the observed estimate.
#'
#'\item Low likelihood drops and light color shading causes small studies (with
#'wide confidence intervals and less weight in the meta-analysis) to be visually less
#'dominant.
#'
#'\item In classic forest plots, it is often hard to depict the magnitude of
#'point estimates to a reasonable degree of accuracy, especially for studies
#'with large meta-analytic weights and correspondingly large plotting symbols
#'(commonly squares). Specific symbols within the raindrops improve the
#'visualization of study point estimates. }
#'
#'Note that for subgroup analysis the height of each raindroop is scaled by the weight of each study within the subgroup divided by
#'the total weight sum of all studies irrespective of subgroup. Therefore, with subgroups present, the overall impression of raindrop
#'heights and color shading within a given subgroup compared to other subgroups conveys information about the relative precision of the meta-analyitc
#'subgroup estimate.
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}.
#'@param group factor indicating the subgroup of each study to plot a subgroup rainforest plot. Has to be in the same order than \code{x}.
#'@param summary_symbol logical scalar. Should the meta-analytic summary result(s) be computed using \code{method} and shown?
#'@param method Which method should be used to compute the study weights and summary effect(s)?
#'  Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the
#'  DerSimonian-Laird method to estimate  \eqn{\tau^2}{tau squared}).
#'@param study_labels a vector with names/identifiers to annotate each study in the
#'  rainforest plot. Has to be in the same order than \code{x}
#'@param summary_label a name to annotate the summary effect. If a subgroup
#'  analysis is plotted, \code{summary_label} should be a vector with a name for each
#'  subgroup summary effect, arranged in the order of the levels in \code{group}.
#'@param study_table a data.frame with addtional study-level variables shown in an aligned table on
#'  the left hand side of the rainforest plot. Should be in the same order than \code{x}. See vignette('metaviz').
#'@param summary_table  a data.frame with addtional summary-level information shown in an aligned table on
#'  the left hand side of the rainforest plot. If \code{group} is supplied, \code{summary_table} should have a row for each subgroup
#'  summary effect, arranged in the order of the levels in \code{group}. See vignette('metaviz').
#'@param table_headers character vector. Headers for each column of the aligned table if \code{study_table} and/or \code{summary_table} is supplied.
#'  By default the names of \code{study_table}.
#'@param annotate_CI logical scalar. Should the effect size and confidence interval values be annotated?
#'@param confidence_level the confidence level for the plotted confidence intervals and likelihood raindrops.
#'@param detail_level numeric value. Values larger than 1 lead to a higher
#'  plotting detail (i.e., smoother likelihood raindrop polygons and more fluent
#'  color shading), values smaller than 1 to less plotting detail compared to
#'  the default plot.
#'@param col character specifying the color palette from package \pkg{RColorBrewer} used.
#'  Can be any of "Blues", "Greys", "Oranges", "Greens", "Reds", and "Purples".
#'@param text_size numeric value. Values larger than 1 lead to larger text size,
#'  values smaller than 1 to smaller text size than the default.
#'@param xlab character label of the x axis. Also used for the header of the aligned table if \code{annotate_CI} is TRUE.
#'@param x_limit numeric vector of length 2 with the limits (minimum, maximum) of the x axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficents using \code{tanh}. See vignette('metaviz').
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@param ... deprecated argument names from earlier versions can still be passed to \code{viz_rainforest}.
#'@references Barrowman, N. J., & Myers, R. A. (2003). Raindrop plots: A new way
#'  to display collections of likelihoods and distributions. \emph{American
#'  Statistician}, \emph{57}, 268-274.
#'@references Jackson, C. H. (2008). Displaying uncertainty with shading.
#'  \emph{American Statistician}, \emph{62}, 340-347.
#'@references Schild, A. H., & Voracek, M. (2015). Finding your way out of the
#'  forest without a trail of bread crumbs: Development and evaluation of two
#'  novel displays of forest plots. \emph{Research Synthesis Methods}, \emph{6},
#'  74-86.
#'@return A Rainforest plot is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Plotting a rainforest plot using the mozart data (for details, see help(mozart)):
#' viz_rainforest(x = mozart[, c("d", "se")],
#' study_labels = mozart[, "study_name"], xlab = "Cohen d")
#'
#' # Visualizing a subgroup analysis of published and unpublished studies
#' viz_rainforest(x = mozart[, c("d", "se")], group = mozart[, "unpublished"],
#' study_labels = mozart[, "study_name"],
#' summary_label = c("Summary (published)", "Summary (unpublished)"),
#' xlab = "Cohen d")
#'
#' # Showing additional information in aligned tables. Log risk ratios are labeled
#' # in their original metric (risk ratios) on the x axis.
#' viz_rainforest(x = exrehab[, c("logrr", "logrr_se")],
#' study_labels = exrehab[, "study_name"],
#' annotate_CI = TRUE, xlab = "RR", x_trans_function = exp,
#' study_table = data.frame(
#' eventsT = paste(exrehab$ai, "/", exrehab$ai + exrehab$bi, sep = ""),
#' eventsC = paste(exrehab$ci, "/", exrehab$ci + exrehab$di, sep = "")),
#' summary_table = data.frame(
#' eventsT = paste(sum(exrehab$ai), "/", sum(exrehab$ai + exrehab$bi), sep = ""),
#' eventsC = paste(sum(exrehab$ci), "/", sum(exrehab$ci + exrehab$di), sep = "")))
#'@export
viz_rainforest <- function(x, group = NULL, summary_symbol = TRUE, method = "FE",
                       study_labels = NULL, summary_label = NULL,
                       study_table = NULL, summary_table = NULL, table_headers = NULL, annotate_CI = FALSE,
                       confidence_level = 0.95, detail_level = 1, col = "Blues",
                       text_size = 3, xlab = "Effect", x_limit = NULL,
                       x_trans_function = NULL, x_breaks = NULL, ...) {

  # check if deprecated arguments were supplied
  add_arg <- list(...)
  if(!is.null(add_arg)) {
    if(is.null(study_labels) & "names" %in% names(add_arg)) {
      study_labels <- add_arg$names
    }
    if(is.null(summary_label) & "summary_name" %in% names(add_arg)) {
      summary_label <- add_arg$summary_name
    }
    if(!is.null(names(add_arg)) & any(!(names(add_arg) %in% c("names", "summary_name")))) {
      warning(paste("unknown arguments supplied to viz_rainforest:", paste(names(add_arg)[which(!(names(add_arg) %in% c("names", "summary_name")))], collapse = ", ")))
    }
  }
  if(summary_symbol == "none") summary_symbol <- FALSE

  viz_forest(x, group = group, summary_symbol = summary_symbol, method = method,
             study_labels = study_labels, summary_label = summary_label,
             study_table = study_table, summary_table = summary_table,
             table_headers = table_headers, annotate_CI = annotate_CI,
             confidence_level = confidence_level, detail_level = detail_level, col = col,
             text_size = text_size, xlab = xlab, x_limit = x_limit,
             x_trans_function = x_trans_function, x_breaks = NULL, type = "rain")
}

#' @export
#' @rdname viz_rainforest
rainforest <- viz_rainforest

