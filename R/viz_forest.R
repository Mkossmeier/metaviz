#'Forest plot variants for meta-analyses
#'
#'Creates a rainforest, thick forest, or classic forest plot.
#'
#'The forest plot is probably the most widley used display to visualize meta-analytc results.
#'The function \code{viz_forest} creates visually appealing forest plots using ggplot2. Many options
#'to flexibly customize the visual appearence and statistical information displayed are provided. In addtion,
#'rainforest plots as well as the thick forest plots can be created, two variants and enhancements of the
#'classical forest plot recently proposed by Schild and Voracek (2015). For furhter details see the documentation of
#'\code{\link[metaviz]{viz_rainforest}}, and \code{\link[metaviz]{viz_thickforest}}.
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}.
#'@param group factor indicating the subgroup of each study to plot a subgroup forest plot. Has to be in the same order than \code{x}.
#'@param type character indicating the forest plot variant that should be plotted. Can be "rain" for rainforest plot,
#'"thick" for a thick forest plot, or "classic" for a traditional forest plot. Default is "classic".
#'@param summary_symbol logical scalar. Should the meta-analytic summary result(s) be computed (using \code{method}) and shown?
#'@param method Which method should be used to compute the study weights and summary effect(s)?
#'  Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the
#'  DerSimonian-Laird method to estimate \eqn{\tau^2}{tau squared}).
#'@param study_labels a vector with names/identifiers to annotate each study in the forest plot.
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
#'@param col character specifying the main color for plotting. For type = "rain" must be one of the follwoing palettes from package
#' \pkg{RColorBrewer}: "Blues", "Greys", "Oranges", "Greens", "Reds", or "Purples".
#'@param text_size numeric value. Values larger than 1 lead to larger text size,
#'  values smaller than 1 to smaller text size than the default.
#'@param xlab character label of the x axis. Also used for the header of the aligned table if \code{annotate_CI} is TRUE.
#'@param x_limit numeric vector of length 2 with the limits (minimum, maximum) of the x axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficents using \code{tanh}. See vignette('metaviz').
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@param ... furhter arguments passed to \code{\link[metaviz]{viz_rainforest}} for type = "rain", or
#'\code{\link[metaviz]{viz_thickforest}} for type = "thick".
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
#' # Plotting a rainforest plot using the mozart data (for details, see help(mozart)):
#' viz_forest(x = mozart[, c("d", "se")],
#' study_labels = mozart[, "study_name"], xlab = "Cohen d")
#'
#' # Subgroup analysis of published and unpublished studies with a classic forest plot
#' viz_forest(x = mozart[, c("d", "se")], study_labels = mozart[, "study_name"],
#' type = "thick", summary_label = c("Summary (published)", "Summary (unpublished)"),
#' group = mozart[, "unpublished"], xlab = "Cohen d")
#'
#' # Thick forest plot with additional information in aligned tables. Log risk
#' # ratios are labeled in their original metric (risk ratios) on the x axis.
#' viz_forest(x = exrehab[, c("logrr", "logrr_se")], type = "rain",
#' study_labels = exrehab[, "study_name"],
#' annotate_CI = TRUE, xlab = "RR", x_trans_function = exp,
#' study_table = data.frame(
#' eventsT = paste(exrehab$ai, "/", exrehab$ai + exrehab$bi, sep = ""),
#' eventsC = paste(exrehab$ci, "/", exrehab$ci + exrehab$di, sep = "")),
#' summary_table = data.frame(
#' eventsT = paste(sum(exrehab$ai), "/", sum(exrehab$ai + exrehab$bi), sep = ""),
#' eventsC = paste(sum(exrehab$ci), "/", sum(exrehab$ci + exrehab$di), sep = "")))
#'@export
viz_forest <- function(x, group = NULL, type = c("classic", "rain", "thick"), summary_symbol = TRUE, method = "FE",
                            study_labels = NULL, summary_label = NULL,
                            study_table = NULL, summary_table = NULL, table_headers = NULL, annotate_CI = FALSE,
                            confidence_level = 0.95, col = "Blues",
                            text_size = 3, xlab = "Effect", x_limit = NULL,
                            x_trans_function = NULL, x_breaks = NULL, ...) {
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
    if(is.null(group) & ncol(x$X) > 1) {
      #check if only categorical moderators were used
      if(!all(x$X == 1 | x$X == 0) | any(apply(as.matrix(x$X[, -1]), 1, sum) > 1))  {
        stop("Can not deal with metafor output object with continuous and/or more than one categorical moderator variable(s).")
      }
      # extract group vector from the design matrix of the metafor object
      no.levels <- ncol(x$X) - 1
      group <- factor(apply(as.matrix(x$X[, -1])*rep(1:no.levels, each = n), 1, sum))
    }
  } else {
    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) | is.matrix(x)) & ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 | sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing
                values, only complete cases are used.")
        study_labels <- study_labels[stats::complete.cases(x[, c(1, 2)])]
        if(!is.null(group)) {
          group <- group[stats::complete.cases(x)]
        }
        x <- x[stats::complete.cases(x), ]
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
  if(!is.null(group) & !is.factor(group)) {
    group <- as.factor(group)
  }
  # check if group vector has the right length
  if(!is.null(group) & (length(group) != length(es))) {
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
  x <- data.frame(es, se, group)
  if(summary_symbol == TRUE) {
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
  }
  if(method != "FE") {
  # compute tau squared for each group
  M <- x %>%
    group_by(group) %>%
    summarise(M = metafor::rma.uni(yi = es, sei = se, method = method)$tau2[[1]]) %>%
    select(M)
   summary_tau2 <- unlist(M)
  } else {
    summary_tau2 <- 0
  }

  # if not exactly one name for every study is supplied the default is used (numbers 1 to the number of studies)
  if(is.null(study_labels) | length(study_labels) != n) {
    if(!is.null(study_labels) & length(study_labels) != n) {
      warning("Argument study_labels has wrong length and is ignored.")
    }
    study_labels <- 1:n
  }

  # if not exactly one name for every subgroup is suppied the default is used
  if(is.null(summary_label) | length(summary_label) != k) {
    if(!is.null(summary_label) & length(summary_label) != k) {
      warning("Argument summary_label has wrong length and is ignored.")
    }
    if(k != 1) {
      summary_label <- paste("Subgroup: ", levels(group), sep = "")
    } else {
      summary_label <- "Summary"
    }
  }

  if(confidence_level <= 0 | confidence_level >= 1) {
    stop("Argument confidence_level must be larger than 0 and smaller than 1.")
  }

  # Determine IDs for studies and summary effects which correspond to plotting y coordinates
  ids <- function(group) {
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
  ID <- ids(group)

  plotdata <- data.frame("x" = es, "se" = se,
                         "ID" = ID$ID[ID$type == "study"],
                         "labels" = study_labels,
                         "group"= group,
                         "x_min" = es - stats::qnorm(1 - (1 - confidence_level)/2)*se,
                         "x_max" = es + stats::qnorm(1 - (1 - confidence_level)/2)*se)

  if(summary_symbol == TRUE) {
  madata <- data.frame("summary_es" = summary_es,
                       "summary_se" = summary_se,
                       "summary_tau2" = summary_tau2,
                        "ID" = ID$ID[ID$type == "summary"])
  } else {
    madata <- data.frame("summary_tau2" = summary_tau2)
  }

# Create forest plot variant ------------------------------------------------------
  args <- c(list(plotdata = plotdata, madata = madata, ID = ID,
            summary_symbol = summary_symbol,
            study_labels = study_labels, summary_label = summary_label,
            study_table = study_table, summary_table = summary_table,
            annotate_CI = annotate_CI, confidence_level = confidence_level, col = col,
            text_size = text_size, xlab = xlab, x_limit = x_limit,
            x_trans_function = x_trans_function, x_breaks = x_breaks), list(...))
  if(type[1] == "rain") {
    p <- do.call(viz_rainforest_internal, args)
  } else {
    if(type[1] == "thick") {
      p <- do.call(viz_thickforest_internal, args)
    } else {
      if(type[1] == "classic") {
        p <- do.call(viz_classicforest_internal, args)
      } else {
        stop("The type argument must be one of rain, thick or classic.")
      }
    }
  }

# Construct tableplots with study and summary information --------
  if(annotate_CI == TRUE | !is.null(study_table) | !is.null(summary_table)) {
    # set limits for the y axis of the table plots
    if(summary_symbol == TRUE) {
      y_limit <- c(-2, n + 3 * k - 2 + 0.5)
    } else {
      y_limit <- c(0, n + 3 * k - 2 + 0.5)
    }

    # Function to create table plots
    table_plot <- function(tbl, r = 5.5, l = 5.5, tbl_titles = table_headers) {
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

      area_per_column <- cumsum(c(1, apply(rbind(tbl_titles, tbl), 2, function(x) max(floor(max(nchar(x, keepNA = FALSE))/20), 0.4))))
      x_values <- area_per_column[1:ncol(tbl)]
      x_limit <- range(area_per_column)

      lab <- data.frame(y = rep(ID$ID, ncol(tbl)),
                        x = rep(x_values,
                                each = length(ID$ID)),
                        value = v, stringsAsFactors = FALSE)

      lab_title <- data.frame(y = rep(n + 3 * k - 2, times = length(tbl_titles)),
                              x = x_values,
                              value = tbl_titles)
      # To avoid "no visible binding for global variable" warning for non-standard evaluation
      y <- NULL
      value <- NULL
      ggplot(lab, aes(x = x, y = y, label = value)) +
        geom_text(size = text_size, hjust = 0, vjust = 0.5) +
        geom_text(data = lab_title, aes(x = x, y = y, label = value), size = text_size, hjust = 0, vjust = 0.5) +
        coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F) +
        geom_hline(yintercept = n + 3 * k - 2 - 0.5) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
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

    if(annotate_CI == TRUE) {
      x_hat <- c(plotdata$x, summary_es)
      lb <- c(c(plotdata$x, summary_es) - stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, summary_se))
      ub <-  c(c(plotdata$x, summary_es) + stats::qnorm(1 - (1 - confidence_level)/2, 0, 1)*c(plotdata$se, summary_se))

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
      table_CI <- table_plot(CI_label, l = 0, tbl_titles =
                               paste(xlab, " [", confidence_level*100, "% CI]", sep = ""))
    } else {
      table_CI <- NULL
    }

    if(summary_symbol != TRUE) {
      summary_table <- NULL
    }
    if(!is.null(study_table) | !is.null(summary_table)) {
      if(!is.null(study_table) & !is.null(summary_table)) {
        if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
        if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
        study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
        summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
        if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
        if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
        if(ncol(summary_table) < ncol(study_table)) {
          n_fillcol <- ncol(study_table) - ncol(summary_table)
          summary_table <- data.frame(summary_table, matrix(rep("", times = nrow(summary_table) * n_fillcol), ncol = n_fillcol))
          stats::setNames(summary_table, names(study_table))
        }
        if(any(names(study_table) != names(summary_table))) summary_table <- stats::setNames(summary_table, names(study_table))
      } else {
        if(is.null(summary_table) & summary_symbol) {
          if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
          study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
          if(nrow(study_table) != n) stop('study_table must be a data.frame with one row for each study.')
          summary_table <- as.data.frame(matrix(rep("", times = ncol(study_table) * k), ncol = ncol(study_table)), stringsAsFactors = FALSE)
          summary_table <- stats::setNames(summary_table, names(study_table))
        }
        if(is.null(study_table)) {
          if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
          summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
          if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
          study_table <- as.data.frame(matrix(rep("", times = ncol(summary_table) * n), ncol = ncol(summary_table)), stringsAsFactors = FALSE)
          study_table <- stats::setNames(study_table, names(summary_table))
          if(is.null(table_headers)) {
            table_headers <- c("Name", rep("", times = ncol(summary_table)))
          }
        }
      }
      if(summary_symbol == TRUE) {
        table_left <- data.frame(Name = c(as.character(study_labels), as.character(summary_label)), rbind(study_table, summary_table))
      } else {
        table_left <- data.frame(Name = c(as.character(study_labels), as.character(summary_label)), study_table)
      }
      if(!is.null(table_headers) & length(table_headers) != ncol(table_left)) {
        if(length(table_headers) - ncol(table_left) == -1) {
          table_headers <- c("Name", table_headers)
        } else {
          warning("Argument table_headers has not the right length and is ignored.")
          table_headers <- NULL
        }
      }
      table_left_plot <- table_plot(table_left, r = 0)
    } else {
      table_left <- NULL
    }
# Align forest plot and table(s) -----------------------------------
    if(!is.null(table_CI) & !is.null(table_left)) {
      layout_matrix <- matrix(c(1, rep(1, times = ncol(table_left) - 1), 2, 2, 2, 2, 3), nrow = 1)
      p <- gridExtra::arrangeGrob(table_left_plot, p, table_CI, layout_matrix = layout_matrix)
      ggpubr::as_ggplot(p)
    } else {
      if(!is.null(table_CI) & is.null(table_left)) {
        layout_matrix <- matrix(c(1, 1, 1, 1, 2), nrow = 1)
        p <- gridExtra::arrangeGrob(p, table_CI, layout_matrix = layout_matrix)
        ggpubr::as_ggplot(p)
      } else {
        if(is.null(table_CI) & !is.null(table_left)) {
          layout_matrix <- matrix(c(1, rep(1, times = ncol(table_left) - 1), 2, 2, 2, 2), nrow = 1)
          p <- gridExtra::arrangeGrob(table_left_plot, p, layout_matrix = layout_matrix)
          ggpubr::as_ggplot(p)
        }
      }
    }
  } else {
    p
  }
}



