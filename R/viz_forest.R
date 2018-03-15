#'Forest plot variants for meta-analyses
#'
#'Creates a rainforest, thick forest, or classic forest plot.
#'
#'The forest plot is the most widely used display to visualize meta-analytic results.
#'The function \code{viz_forest} creates visually appealing and informative-rich forest plots using ggplot2. Many options
#'to flexibly customize the visual appearance and statistical information displayed are provided. In addition,
#'rainforest plots as well as the thick forest plots can be created, two variants and enhancements of the
#'classical forest plot recently proposed by Schild and Voracek (2015). For further details see the documentation of
#'\code{\link[metaviz]{viz_rainforest}}, and \code{\link[metaviz]{viz_thickforest}}.
#'
#'\bold{Available forest plot types}
#'
#'Different aspects of meta-analytic data can be shown in forest plots. Five different types are available in \code{viz_forest} via the \code{type} parameter.
#'Argument \code{"standard"} (default) shows study results as well as summary results in the forest plot. \code{"study_only"} allows to only show study results without the meta-analytic summary estimate.
#'\code{"summary_only"} can be used to only show meta-analytic summary estimate(s), which is primarily useful to visualize several subgroup results (using \code{group}).
#'\code{"cumulative"} shows a cumulative meta-analysis, that is, meta-analytic summary effects are computed sequentially by adding each study one-by-one.
#'Studies are added in the same order than they were supplied in \code{x}. Finally, \code{"sensitivity"} shows for each study the meta-analytic summary
#'effect if that particular study is not considered in the computation of the summary effect (leave-one-out analysis).
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}; then effect sizes and standard errors are extracted from \code{x}.
#'@param group factor indicating the subgroup of each study to plot a subgroup forest plot. Has to be in the same order than \code{x}.
#'@param type character string indicating the type of forest plot to be plotted. Can be "standard" (default), "study_only",
#'  "summary_only", "cumulative", or "sensitivity". See 'Details'.
#'@param variant character string indicating the forest plot variant that should be plotted. Can be "rain" for rainforest plot,
#'  "thick" for a thick forest plot, or "classic" (default) for a traditional forest plot.
#'@param method character string indicating which method should be used to compute the study weights and summary effect(s).
#'  Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the
#'  DerSimonian-Laird method to estimate \eqn{\tau^2}{tau squared}).
#'@param study_labels a character vector with names/identifiers to annotate each study in the forest plot.
#'  Has to be in the same order than \code{x}. Ignored if \code{study_table} and/or \code{summary_table} is supplied.
#'@param summary_label a character string specifying the name to annotate the summary effect. If a subgroup
#'  analysis is plotted, \code{summary_label} should be a character vector with a name for each
#'  subgroup summary effect, arranged in the order of the levels of \code{group}. Ignored if \code{study_table} and/or
#'  \code{summary_table} is supplied.
#'@param confidence_level numeric value. The confidence level for the plotted confidence intervals.
#'@param col character string specifying the main color for plotting. For \code{variant = "rain"} must be one of the following palettes from package
#' \pkg{RColorBrewer}: "Blues", "Greys", "Oranges", "Greens", "Reds", or "Purples".
#'@param text_size numeric value. Size of text in the forest plot. Default is 3.
#'@param xlab character string specifying the label of the x axis. Also used for the header of the aligned table if \code{annotate_CI} is \code{TRUE}.
#'@param x_limit numeric vector of length 2 with the limits (minimum, maximum) of the x axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficients using \code{tanh}. See vignette('metaviz').
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@param annotate_CI logical scalar. Should the effect size and confidence interval values be shown as text in an aligned table on the right-hand side of the forest plot?
#'@param study_table a data.frame with additional study-level variables which should be shown in an aligned table.
#'  Has to be in the same order than \code{x}. See vignette('metaviz').
#'@param summary_table a data.frame with additional summary-level information shown in an aligned table.
#'  If \code{group} is supplied, \code{summary_table} must have a row for each subgroup
#'  summary effect, arranged in the order of the levels of \code{group}. See vignette('metaviz').
#'@param table_headers character vector. Headers for each column of the aligned table if \code{study_table} and/or \code{summary_table} is supplied.
#'  By default the column names of \code{study_table} are used.
#'@param table_layout numeric layout matrix passed to \code{layout_matrx} of \code{\link[gridExtra]{arrangeGrob}}. Can be used to overwrite the default spacing
#'  of the forest plot and aligned tables via \code{study_table}, \code{summary_table}, and \code{annotate_CI}.
#'@param ... further arguments passed to \code{\link[metaviz]{viz_rainforest}} for \code{variant = "rain"}, or
#'  \code{\link[metaviz]{viz_thickforest}} for \code{variant = "thick"}.
#'@references Schild, A. H., & Voracek, M. (2015). Finding your way out of the
#'  forest without a trail of bread crumbs: Development and evaluation of two
#'  novel displays of forest plots. \emph{Research Synthesis Methods}, \emph{6},
#'  74-86.
#'@return A forest plot is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Plotting the mozart data using a classic forest plot
#' viz_forest(x = mozart[, c("d", "se")],
#' study_labels = mozart[, "study_name"], xlab = "Cohen d")
#'
#' # Subgroup analysis of published and unpublished studies shown in a rainforest plot
#' viz_forest(x = mozart[, c("d", "se")], study_labels = mozart[, "study_name"], method = "REML",
#' variant = "rain", summary_label = c("Summary (rr_lab = no)", "Summary (rr_lab = yes)"),
#' group = mozart[, "rr_lab"], xlab = "Cohen d")
#'
#' # Thick forest plot with additional information in aligned tables. Log risk
#' # ratios are labeled in their original metric (risk ratios) on the x axis.
#' viz_forest(x = exrehab[, c("logrr", "logrr_se")], variant = "thick",
#' xlab = "RR", x_trans_function = exp, annotate_CI = TRUE,
#' study_table = data.frame(
#' Name = exrehab[, "study_name"],
#' eventsT = paste(exrehab$ai, "/", exrehab$ai + exrehab$bi, sep = ""),
#' eventsC = paste(exrehab$ci, "/", exrehab$ci + exrehab$di, sep = "")),
#' summary_table = data.frame(
#' Name = "Summary",
#' eventsT = paste(sum(exrehab$ai), "/", sum(exrehab$ai + exrehab$bi), sep = ""),
#' eventsC = paste(sum(exrehab$ci), "/", sum(exrehab$ci + exrehab$di), sep = "")),
#' table_layout = matrix(c(1, 1, 2, 2, 3), nrow = 1))
#'@export
viz_forest <- function(x, group = NULL, type = "standard", variant = "classic", method = "FE",
                            study_labels = NULL, summary_label = NULL,
                            confidence_level = 0.95, col = "Blues",
                            text_size = 3, xlab = "Effect", x_limit = NULL,
                            x_trans_function = NULL, x_breaks = NULL,
                            annotate_CI = FALSE, study_table = NULL, summary_table = NULL,
                            table_headers = NULL, table_layout = NULL, ...) {
  #'@import ggplot2
  #'@import dplyr

# Handle input object -----------------------------------------------------
  # input is output of rma (metafor)
  if("rma" %in% class(x)) {
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    n <- length(es)

    # check if group argument has the right length
    if(!is.null(group) & (length(group) != length(es))) {
      warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored")
      group <- NULL
    }
    if(method != x$method) {
      warning("Note: method argument used differs from input object of class rma.uni (metafor)")
    }
    # If No group is supplied try to extract group from input object of class rma.uni (metafor)
    if(is.null(group) && ncol(x$X) > 1) {
      #check if only categorical moderators were used
      if(!all(x$X == 1 || x$X == 0) || any(apply(as.matrix(x$X[, -1]), 1, sum) > 1))  {
        stop("Can not deal with metafor output object with continuous and/or more than one categorical moderator variable(s).")
      }
      # extract group vector from the design matrix of the metafor object
      no.levels <- ncol(x$X) - 1
      group <- factor(apply(as.matrix(x$X[, -1])*rep(1:no.levels, each = n), 1, sum))
    }
  } else {
    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        study_labels <- study_labels[stats::complete.cases(x[, c(1, 2)])]
        if(!is.null(group)) {
          group <- group[stats::complete.cases(x)]
        }
        x <- x[stats::complete.cases(x), ]
      }
      # check if input is numeric
      if(!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
        stop("Input argument has to be numeric; see help(viz_forest) for details.")
      }
      # check if there are any negative standard errors
      if(!all(x[, 2] >= 0)) {
        stop("Negative standard errors supplied")
      }
      # extract effects and standard errors
      es <- x[, 1]
      se <- x[, 2]
      n <- length(es)
    } else {
      stop("Unknown input argument. See help ('metaviz').")
    }
}

# Preprocess data ---------------------------------------------------------
  # check if group is a factor
  if(!is.null(group) && !is.factor(group)) {
    group <- as.factor(group)
  }
  # check if group vector has the right length
  if(!is.null(group) && (length(group) != length(es))) {
    warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored")
    group <- NULL
  }

  # if no group argument is supplied, use all cases
  if(is.null(group)) {
    group <- factor(rep(1, times = n))
  }

  # drop unused levels of group factor
  group <- droplevels(group)
  k <- length(levels(group))

  # main data
  x <- data.frame(es, se, group)

  if(n <= 1 && type == "sensitivity") {
    stop('For type = "sensitvitiy" there has to be more than 1 study.')
  }

  # Compute meta-analytic summary effect estimates
  if(type %in% c("standard", "summary_only", "sensitivity", "cumulative")) {
    M <- NULL # To avoid "no visible binding for global variable" warning for non-standard evaluation
    # compute meta-analytic summary effect for each group
    M <- x %>%
      group_by(group) %>%
      summarise(M = metafor::rma.uni(yi = es, sei = se, method = method)$b[[1]]) %>%
      select(M)
      summary_es <-  unlist(M)

    # compute standard error of the meta-analytic summary effect for each group
    M <- x %>%
      group_by(group) %>%
      summarise(M = metafor::rma.uni(yi = es, sei = se, method = method)$se[[1]]) %>%
      select(M)
      summary_se <- unlist(M)
    if(type == "sensitivity") {
      loo_es <- function(es, se) {
        res <- numeric(length(es))
        for(i in 1:length(es)) {
          res[i] <- metafor::rma.uni(yi = es[-i], sei = se[-i], method = method)$b[[1]]
        }
        res
      }
      loo_se <- function(es, se) {
        res <- numeric(length(es))
        for(i in 1:length(es)) {
          res[i] <- metafor::rma.uni(yi = es[-i], sei = se[-i], method = method)$se[[1]]
        }
        res
      }
      sens_data <- x %>%
        group_by(group) %>%
        mutate(summary_es = loo_es(es, se),
               summary_se = loo_se(es, se))
    }
    if(type == "cumulative") {
      rollingma_es <- function(es, se) {
        res <- numeric(length(es))
        for(i in 1:length(es)) {
          res[i] <- metafor::rma.uni(yi = es[1:i], sei = se[1:i], method = method)$b[[1]]
        }
        res
      }
      rollingma_se <- function(es, se) {
        res <- numeric(length(es))
        for(i in 1:length(es)) {
          res[i] <- metafor::rma.uni(yi = es[1:i], sei = se[1:i], method = method)$se[[1]]
        }
        res
      }
      cum_data <- x %>%
        group_by(group) %>%
        mutate(summary_es = rollingma_es(es, se),
               summary_se = rollingma_se(es, se))
    }
  } else {
    if(type != "study_only") {
     stop('Argument of type must be one of "standard", "study_only", "summary_only", "cumulative", or "sensitivity".')
    }
  }

  # Compute tau^2 estimate
  if(type %in% c("standard", "study_only")) {
    if(method != "FE") {
    # compute tau squared for each group
    M <- x %>%
      group_by(group) %>%
      summarise(M = metafor::rma.uni(yi = es, sei = se, method = method)$tau2[[1]]) %>%
      select(M)
     summary_tau2 <- unlist(M)
    } else {
      summary_tau2 <- rep(0, times = k)
    }
  }

  if(type %in% c("study_only")) {
    if(!is.null(summary_table) || !is.null(summary_label)) {
      warning('For type "study_only" supplied summary_table and summary_label are ignored.')
    }
    summary_table <- NULL
    summary_label <- NULL
  }
  if(type == "summary_only") {
    if(!is.null(study_table) || !is.null(study_labels)) {
      warning('For type "summary_only" supplied study_table and study_labels are ignored.')
    }
    study_table <- NULL
    study_labels <- NULL
  }


  # if not exactly one name for every study is supplied the default is used (numbers 1 to the number of studies)
  if(is.null(study_labels) || length(study_labels) != n) {
    if(!is.null(study_labels) && length(study_labels) != n) {
      warning("Argument study_labels has wrong length and is ignored.")
    }
    study_labels <- 1:n
  }

  # if not exactly one name for every subgroup is suppied the default is used
  if(is.null(summary_label) || length(summary_label) != k) {
    if(!is.null(summary_label) && length(summary_label) != k) {
      warning("Argument summary_label has wrong length and is ignored.")
    }
    if(k != 1) {
      summary_label <- paste("Subgroup: ", levels(group), sep = "")
    } else {
      summary_label <- "Summary"
    }
  }

  if(confidence_level <= 0 || confidence_level >= 1) {
    stop("Argument confidence_level must be larger than 0 and smaller than 1.")
  }

  if(!is.null(x_trans_function) && !is.function(x_trans_function)) {
    warning("Argument x_trans_function must be a function; input ignored.")
    x_trans_function <- NULL
  }

  if(!is.null(table_layout) && !is.matrix(table_layout)) {
    warning("Agument of table_layout is not a matrix and is ignored.")
    table_layout <- NULL
  }

  # Determine IDs for studies and summary effects which correspond to plotting y coordinates
  ids <- function(group, n) {
    k <- length(levels(group))
    ki_start <- cumsum(c(1, as.numeric(table(group))[-k] + 3))
    ki_end <- ki_start + as.numeric(table(group)) - 1
    study_IDs <- numeric(n)
    for(i in 1:k) {
      study_IDs[group==levels(group)[i]] <- ki_start[i]:ki_end[i]
    }
    summary_IDs <- ki_end + 2
    data.frame("ID" = (n + 3*k-2) - c(study_IDs, summary_IDs),
               "type" = factor(c(rep("study", times = length(study_IDs)),
                                 rep ("summary", times = length(summary_IDs)))))
  }
  if(type != "summary_only") {
  ID <- ids(group, n = n)
  } else {
    ID <- ids(unique(group), n = k)
  }

  if(type %in% c("standard", "study_only")) {
    plotdata <- data.frame("x" = es, "se" = se,
                           "ID" = ID$ID[ID$type == "study"],
                           "labels" = study_labels,
                           "group"= group,
                           "x_min" = es - stats::qnorm(1 - (1 - confidence_level)/2)*se,
                           "x_max" = es + stats::qnorm(1 - (1 - confidence_level)/2)*se)
    if(type == "standard") {
      madata <- data.frame("summary_es" = summary_es,
                           "summary_se" = summary_se,
                           "summary_tau2" = summary_tau2,
                           "ID" = ID$ID[ID$type == "summary"])
    }
    if(type == "study_only") {
      madata <- data.frame("summary_tau2" = summary_tau2)
    }
  } else {
    if(type == "summary_only") {
      plotdata <- data.frame("x" = summary_es, "se" = summary_se,
                             "ID" = ID$ID[ID$type == "summary"],
                             "labels" = summary_label,
                             "group"= levels(group),
                             "x_min" = summary_es - stats::qnorm(1 - (1 - confidence_level)/2)*summary_se,
                             "x_max" = summary_es + stats::qnorm(1 - (1 - confidence_level)/2)*summary_se)
      madata <- NULL
    } else {
      if(type == "cumulative") {
        plotdata <- data.frame("x" = cum_data$summary_es, "se" = cum_data$summary_se,
                               "ID" = ID$ID[ID$type == "study"],
                               "labels" = study_labels,
                               "group"= group,
                               "x_min" = cum_data$summary_es - stats::qnorm(1 - (1 - confidence_level)/2)*cum_data$summary_se,
                               "x_max" = cum_data$summary_es + stats::qnorm(1 - (1 - confidence_level)/2)*cum_data$summary_se)
        madata <- data.frame("summary_es" = summary_es,
                             "summary_se" = summary_se,
                             "ID" = ID$ID[ID$type == "summary"])
      } else {
        if(type == "sensitivity") {
          plotdata <- data.frame("x" = sens_data$summary_es, "se" = sens_data$summary_se,
                                 "ID" = ID$ID[ID$type == "study"],
                                 "labels" = study_labels,
                                 "group"= group,
                                 "x_min" = sens_data$summary_es - stats::qnorm(1 - (1 - confidence_level)/2)*sens_data$summary_se,
                                 "x_max" = sens_data$summary_es + stats::qnorm(1 - (1 - confidence_level)/2)*sens_data$summary_se)
          madata <- data.frame("summary_es" = summary_es,
                               "summary_se" = summary_se,
                               "ID" = ID$ID[ID$type == "summary"])
        }
      }
    }
  }

# Create forest plot variant ------------------------------------------------------
  args <- c(list(plotdata = plotdata, madata = madata,
            type = type,
            study_labels = study_labels, summary_label = summary_label,
            study_table = study_table, summary_table = summary_table,
            annotate_CI = annotate_CI, confidence_level = confidence_level, col = col,
            text_size = text_size, xlab = xlab, x_limit = x_limit,
            x_trans_function = x_trans_function, x_breaks = x_breaks), list(...))

  if(variant == "rain") {
    p <- do.call(viz_rainforest_internal, args)
  } else {
    if(variant == "thick") {
      p <- do.call(viz_thickforest_internal, args)
    } else {
      if(variant == "classic") {
        p <- do.call(viz_classicforest_internal, args)
      } else {
        stop("The argument of variant must be one of rain, thick or classic.")
      }
    }
  }

# Construct tableplots with study and summary information --------
  if(annotate_CI == TRUE || !is.null(study_table) || !is.null(summary_table)) {

    # set limits for the y axis of the table plots
    if(type %in% c("standard", "sensitivity", "cumulative")) {
      y_limit <- c(min(plotdata$ID) - 3, max(plotdata$ID) + 1.5)
    } else {
      y_limit <- c(min(plotdata$ID) - 1, max(plotdata$ID) + 1.5)
    }

    # Function to create table plots
    table_plot <- function(tbl, ID, r = 5.5, l = 5.5, tbl_titles = table_headers) {
      # all columns and column names are stacked to a vector
      df_to_vector <- function(df) {
        v <- vector("character", 0)
        for(i in 1:ncol(df)) v <- c(v, as.vector(df[, i]))
        v
      }
      if(!is.data.frame(tbl)) tbl <- data.frame(tbl)
      tbl <- data.frame(lapply(tbl, as.character), stringsAsFactors = FALSE)
      if(is.null(tbl_titles)) {
        tbl_titles <- names(tbl)
      }
      v <- df_to_vector(tbl)

      area_per_column <- cumsum(c(1, apply(rbind(tbl_titles, tbl), 2, function(x) max(round(max(nchar(x, keepNA = FALSE))/100, 2),  0.03))))
      x_values <- area_per_column[1:ncol(tbl)]
      x_limit <- range(area_per_column)

      lab <- data.frame(y = rep(ID, ncol(tbl)),
                        x = rep(x_values,
                                each = length(ID)),
                        value = v, stringsAsFactors = FALSE)

      lab_title <- data.frame(y = rep(max(plotdata$ID) + 1, times = length(tbl_titles)),
                              x = x_values,
                              value = tbl_titles)

      # To avoid "no visible binding for global variable" warning for non-standard evaluation
      y <- NULL
      value <- NULL
      ggplot(lab, aes(x = x, y = y)) +
        geom_text(aes(label = value), size = text_size, hjust = 0, vjust = 0.5) +
        geom_text(data = lab_title, aes(x = x, y = y, label = value), size = text_size, hjust = 0, vjust = 0.5) +
        coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F) +
        geom_hline(yintercept = max(plotdata$ID) + 0.5) +
        theme_bw() +
        theme(text = element_text(size = 1/0.352777778*text_size),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "none",
              panel.border = element_blank(),
              axis.text.x = element_text(colour="white"),
              axis.text.y = element_blank(),
              axis.ticks.x = element_line(colour="white"),
              axis.ticks.y = element_blank(),
              axis.line.x = element_line(colour="white"),
              axis.line.y = element_blank(),
              plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt")) +
        labs(x = "", y = "")
    }

    # Study and/or summary table left
    if(!is.null(study_table) || !is.null(summary_table)) {
      # Case study table and summary table are both supplied (type standard, cumulative, or sensitivity)
      if(!is.null(study_table) && !is.null(summary_table)) {
        if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
        if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
        study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
        summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
        if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
        if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
        if(ncol(summary_table) < ncol(study_table)) {
          n_fillcol <- ncol(study_table) - ncol(summary_table)
          summary_table <- data.frame(summary_table, matrix(rep("", times = nrow(summary_table) * n_fillcol), ncol = n_fillcol))
          summary_table<- stats::setNames(summary_table, names(study_table))
        } else {
          if(ncol(summary_table) > ncol(study_table)) {
            n_fillcol <- ncol(summary_table) - ncol(study_table)
            study_table <- data.frame(study_table, matrix(rep("", times = nrow(study_table) * n_fillcol), ncol = n_fillcol))
            study_table <- stats::setNames(study_table, names(summary_table))
          }
        }
        if(any(names(study_table) != names(summary_table))) summary_table <- stats::setNames(summary_table, names(study_table))
      } else {
        # Case only study table is supplied
        if(is.null(summary_table)) {
          if(type %in% c("standard", "sensitivity", "cumulative", "study_only")) {
            if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
            study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
            if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
            summary_table <- as.data.frame(matrix(rep("", times = ncol(study_table) * k), ncol = ncol(study_table)), stringsAsFactors = FALSE)
            summary_table <- stats::setNames(summary_table, names(study_table))
          }
        }
        # Case only summary table is supplied
        if(is.null(study_table)) {
          if(type %in% c("standard", "sensitivity", "cumulative")) {
            if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
            summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
            if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
            study_table <- as.data.frame(matrix(rep("", times = ncol(summary_table) * n), ncol = ncol(summary_table)), stringsAsFactors = FALSE)
            study_table <- stats::setNames(study_table, names(summary_table))
          } else {
            if(type %in% c("summary_only")) {
              if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
              summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
              if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
              study_table <- as.data.frame(matrix(rep("", times = ncol(summary_table) * k), ncol = ncol(summary_table)), stringsAsFactors = FALSE)
              study_table <- stats::setNames(study_table, names(summary_table))
            }
          }
        }
      }

    table_left <- data.frame(rbind(study_table, summary_table))

    if(!is.null(table_headers) && length(table_headers) != ncol(table_left)) {
      warning("Argument table_headers has not the right length and is ignored.")
      table_headers <- NULL
    }
  table_left_plot <- table_plot(table_left, ID = ID$ID, r = 0, tbl_titles = table_headers)
  } else {
    table_left <- NULL
  }

  # Textual CI and effect size values right
  if(annotate_CI == TRUE) {
    if(type %in% c("standard", "sensitivity", "cumulative")) {
      x_hat <- c(plotdata$x, madata$summary_es)
      lb <- c(c(plotdata$x, madata$summary_es) - stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, madata$summary_se))
      ub <-  c(c(plotdata$x, madata$summary_es) + stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, madata$summary_se))

      if(!is.null(x_trans_function)) {
        x_hat <- x_trans_function(x_hat)
        lb <- x_trans_function(lb)
        ub <- x_trans_function(ub)
      }

      lb <- format(round(lb, 2), nsmall = 2)
      ub <- format(round(ub, 2), nsmall = 2)
      x_hat <- format(round(x_hat, 2), nsmall = 2)

      CI <- paste(x_hat, " [", lb, ", ", ub, "]", sep = "")
      CI_label <- data.frame(CI = CI, stringsAsFactors = FALSE)

      table_CI <- table_plot(CI_label, ID = c(plotdata$ID, madata$ID), l = 0, r = 11,  tbl_titles =
                               paste(xlab, " [", confidence_level*100, "% CI]", sep = ""))
    } else {
      if(type %in% c("study_only", "summary_only")) {
        x_hat <- plotdata$x
        lb <- plotdata$x - stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*plotdata$se
        ub <-  plotdata$x + stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*plotdata$se

        if(!is.null(x_trans_function)) {
          x_hat <- x_trans_function(x_hat)
          lb <- x_trans_function(lb)
          ub <- x_trans_function(ub)
        }

        lb <- format(round(lb, 2), nsmall = 2)
        ub <- format(round(ub, 2), nsmall = 2)
        x_hat <- format(round(x_hat, 2), nsmall = 2)
        CI <- paste(x_hat, " [", lb, ", ", ub, "]", sep = "")
        CI_label <- data.frame(CI = CI, stringsAsFactors = FALSE)
        table_CI <- table_plot(CI_label, ID = plotdata$ID, l = 0, r = 11, tbl_titles =
                                 paste(xlab, " [", confidence_level*100, "% CI]", sep = ""))
      }
    }
  } else {
    table_CI <- NULL
  }
# Align forest plot and table(s) -----------------------------------
    if(!is.null(table_CI) && !is.null(table_left)) {
      if(is.null(table_layout)) {
        layout_matrix <- matrix(c(rep(1, times = ncol(table_left)), rep(2, times = 3), 3), nrow = 1)
      } else {
        layout_matrix <- table_layout
      }
      p <- gridExtra::arrangeGrob(table_left_plot, p, table_CI, layout_matrix = layout_matrix)
      ggpubr::as_ggplot(p)
    } else {
      if(!is.null(table_CI) && is.null(table_left)) {
        if(is.null(table_layout)) {
          layout_matrix <- matrix(c(1, 1, 1, 1, 2), nrow = 1)
        } else {
          layout_matrix <- table_layout
        }
        p <- gridExtra::arrangeGrob(p, table_CI, layout_matrix = layout_matrix)
        ggpubr::as_ggplot(p)
      } else {
        if(is.null(table_CI) && !is.null(table_left)) {
          if(is.null(table_layout)) {
            layout_matrix <- matrix(c(rep(1, times = 1 + ncol(table_left)), 2, 2, 2, 2, 2), nrow = 1)
          } else {
            layout_matrix <- table_layout
          }
          p <- gridExtra::arrangeGrob(table_left_plot, p, layout_matrix = layout_matrix)
          ggpubr::as_ggplot(p)
        }
      }
    }
  } else {
    p
  }
}



