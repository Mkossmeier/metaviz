#'Rainforest plots for meta-analyses
#'
#'Creates a rainforest plot.
#'
#'Rainforest plots were proposed by Schild and Voracek (2015) as a variant and
#'enhancement of classic forest plots. Rainforest plots use likelihood drops
#'to depict study and summary level results (for details, see Barrowman &
#'Myers, 2003). \code{rainforest} assumes normality of effect sizes to
#'construct these likelihood drops. Therefore, the height of a likelihood drop
#'for a hypothetical true value is proportional to the likelihood of that true
#'value given the observed estimate, while the width is identical to the
#'confidence interval when normality of effect sizes is assumed (which, for
#'example, also is the (default) case in \link[metafor]{metafor}, Revman, and CMA). Additionaly,
#'color shading is utilized to further visualize statistical uncertainty, as
#'suggested by Jackson (2008). Finally, study and summary level point estimates
#'are depicted clearly by a specific symbol. Rainforest plots have the
#'following advantages, as compared to classic forest plots:
#'
#'\enumerate{ \item The width of the likelihood raindops corresponds to the
#'confidence intervals, as also shown in the classic forest plot. In addition,
#'the height of the likelihood drops and color shading adequateley visualizes the
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
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \link[metafor]{metafor}.
#'@param names a vector with names/identifiers to annotate each study in the
#'  rainforest plot.
#'@param summary_name a name to annotate the summary effect. If a subgroup
#'  analysis is plotted, \code{summary_name} should be a vector with names for the
#'  summary effect of each subgroup.
#'@param group a factor indicating the subgroup of each study. The group
#'  argument is not necessary and ignored if \code{x} is an output object of package
#'  \link[metafor]{metafor}; in this case, subgroup results are used directly from the \link[metafor]{metafor}
#'  output object (if one categorical moderator was used).
#'@param method Which method should be used to compute the summary effect?
#'  Supported methods are "FEM" for the fixed-effect model and "REM" for the
#'  random-effects model (using the DerSimonian-Laird method to estimate the
#'  between-study variance \eqn{\tau^2}{tau squared}). This argument is ignored if \code{x} is the output
#'  of the function \code{\link[metafor]{rma.uni}} of package \link[metafor]{metafor}; in this
#'  case, meta-analytic results are used directly from the metafor output
#'  object.
#'@param confidence_level the confidence level for the plotted confidence
#'  intervals and likelihood raindrops.
#'@param summary_symbol which symbol should be used to depict the
#'  meta-analytic summary effect and its confidence interval? Available
#'  options are:
#'  \itemize{
#'  \item rain: A likelihood drop is used.
#'  \item diamond:
#'  A classic summary diamond is used.
#'  \item none: No summary effect is displayed.
#'  }
#'@param col character specifying the color palette (from package
#'  \code{\link[RColorBrewer]{RColorBrewer}}) used for shading. Default value is
#'  "Blues". Other options are "Greys", "Oranges", "Greens", "Reds", and
#'  "Purples".
#'@param shading logical scalar indicating if the likelihood drops should be
#'  color shaded. Default is \code{TRUE}.
#'@param xlab character label of the x axis.
#'@param text_size numeric value. Values larger than 1 lead to larger text size,
#'  values smaller than 1 to smaller text size than the default.
#'@param detail_level numeric value. Values larger than 1 lead to a higher
#'  plotting detail (i.e., smoother likelihood raindrop polygons and more fluent
#'  color shading), values smaller than 1 to less plotting detail compared to
#'  the default plot.
#'@references Barrowman, N. J., & Myers, R. A. (2003). Raindrop plots: A new way
#'  to display collections of likelihoods and distributions. \emph{American
#'  Statistician}, \emph{57}, 268-274.
#'@references Jackson, C. H. (2008). Displaying uncertainty with shading.
#'  \emph{American Statistician}, \emph{62}, 340-347.
#'@references Schild, A. H., & Voracek, M. (2015). Finding your way out of the
#'  forest without a trail of bread crumbs: Development and evaluation of two
#'  novel displays of forest plots. \emph{Research Synthesis Methods}, \emph{6},
#'  74-86.
#'@return A Rainforest plot is created using ggplot2. The resulting plot is a
#'  ggplot2 object and is printed by default.
#'@examples
#' library(metaviz)
#' # Plotting a rainforest plot using the mozart data (for details, see help(mozart)):
#' rainforest(x = mozart[, c("d", "se")],
#' names = mozart[, "study_name"], xlab = "Cohen's d")
#'
#' # Visualizing a subgroup analysis of published and unpublished studies
#' rainforest(x = mozart[, c("d", "se")], names = mozart[, "study_name"],
#' summary_name = c("Summary (published)", "Summary (unpublished)"),
#' group = mozart[, "unpublished"], xlab = "Cohen's d")
#'@export
rainforest <- function(x,  names = NULL , summary_name = "Summary", group = NULL, method = "FEM",
                       confidence_level = 0.95, summary_symbol = "rain", col = "Blues", shading = TRUE,
                       xlab = "Effect Size", text_size = 1, detail_level = 1) {
#'@import ggplot2
#'@import stats
  if("rma"%in%class(x)) {  # check if input is output of rma (metafor)
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    n <- length(es)

    if(ncol(x$X) > 1) { # check if subgroups are given
      if(!all(x$X == 1 | x$X == 0)) { #check if only categorical moderators were used
        stop("Not supported input object supplied: Metafor output object with continuous moderator variable(s)")
      }
      # extract group vector from the design matrix of the metafor object
      no.levels <- ncol(x$X) - 1
      group <- factor(apply(as.matrix(x$X[, -1])*rep(1:no.levels, each = n), 1, sum))
      summary_es <- as.numeric(c(x$b[1], x$b[-1] +  x$b[1]))
      summary_se <- as.numeric(c(x$se[1], sqrt(x$se[-1]^2 + x$se[1]^2)))
      k <- length(levels(group))
    } else {
      summary_es <- as.numeric(x$b)
      summary_se <- as.numeric(x$se)
      group <- factor(rep(1, times = n))
      k <- 1
    }

  } else {
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) == 2) { # check if a data.frame or matrix with (exactly) two columns is supplied
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) { # check if there are missing values
        warning("The effect sizes or standard errors contain missing
                values, only complete cases are used")
        x <- x[complete.cases(x), ]
        names = names[complete.cases(x)]
        if(!is.null(group)) {
          group <- group[complete.cases(group), ]
        }
      }

      if(!is.numeric(x[, 1]) || !is.numeric(x[,2])) { # check if effect size and standard error columns are numeric
        x[, 1] <- as.numeric(x[, 1])
        x[, 2] <- as.numeric(x[, 2])
      }

      # check if there are any negative standard errors
      if(!all(x[, 2]>=0)) {
        stop("Negative standard errors supplied")
      }

      es <- x[, 1]
      se <- x[, 2]
      n <- length(es)

      if(!is.null(group) && !is.factor(group)) {
        group <- as.factor(group)
      }

      # check if group vector has the right length
      if(!is.null(group) && (length(group) != length(es)))
      {
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

      # compute summary effect and its standard error
      if(method == "FEM") {
        summary_es <- plyr::ddply(x, plyr::.(group), plyr::summarize, summary_es = sum((1/se^2)*es)/sum(1/se^2))[, 2]
        summary_se <- plyr::ddply(x, plyr::.(group), plyr::summarize, summary_se = sqrt(1/sum(1/se^2)))[, 2]
      } else {
        if(method == "REM") {
          summary_es <- plyr::ddply(x, plyr::.(group), plyr::summarize, rem_es = metaviz::rem_effect(es ,se))[, 2]
          summary_se <- plyr::ddply(x, plyr::.(group), plyr::summarize, rem_se = metaviz::rem_err(es, se))[, 2]
        } else {
          stop("rainforest only supports FEM or REM as method argument")
        }
      }
    } else {
      stop("Unknown input argument. Must be a matrix or data.frame with
           effect sizes and standard errors in two numeric columns or an
           output object of function rma.uni() from the package metafor")
    }
  }

  # if not exactly one name for every study is suppied the default is used (numbers 1 to the number of studies)
  if(is.null(names) || length(es)!=length(names)) {
    names <- 1:n
  }

  # if not exactly one name for every subgroup is suppied the default is used ("Summary.i")
  if(length(summary_name) != k) {
    summary_name <- paste("Summary.", 1:k, sep = "")
  }

  # Determine IDs for studies and summary effects (corresponds to plotting coordinates in the rainforest plot)
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

  # create dataframe for plotting point estimates and confidence intervals
  plotdata <- data.frame("x" = c(es, summary_es), "se" = c(se, summary_se),
                         "ID" = ids(group)$ID, "names" = c(names, summary_name))

  # if no summary_symbol should be plotted, only use study information
  if(summary_symbol == "none") {
    plotdata <- data.frame("x" = es, "se" = se,
                           "ID" = ids(group)$ID[ids(group)$type == "study"],
                           "names" = names)
  }

  # Add confidence intervals
  plotdata <- data.frame(rbind(plotdata, plotdata),
                         "ci_value" = c(plotdata$x - qnorm(1-(1-confidence_level)/2, 0, 1)*plotdata$se,
                                        plotdata$x + qnorm(1-(1-confidence_level)/2, 0, 1)*plotdata$se))

  # function ll constructs a likelihood raindrop of a study. If shading is true, the raindrop is built out of
  # several distinct segments (to adequately color shade the raindrop)
  ll <- function(x, max.range) {
    # width of the region over which the raindop is built
    se.factor <- ceiling(qnorm(1-(1-confidence_level)/2, 0, 1))
    width <- abs((x[1] - se.factor * x[2]) - (x[1] + se.factor * x[2]))

    if(shading) {
      # max.range is internally determined as the width of the broadest raindop.
      # The number of points to construct the raindrop is chosen proportional to
      # the ratio of the width of the raindrop and the max.range,
      # because slim/dense raindrops do not need as many support points as very broad ones.
      # Minimum is 200 points (for the case that width/max.range gets very small)
      length.out <- max(c(floor(1000 * width / max.range), 200))
    } else {
      # if there is no shading, 200 is chosen as (default) number of support points per raindrop
      length.out <- 200
    }
    # Create sequence of points to construct the raindrop. The number of points is chosen by length.out
    # and can be changed by the user with the parameter detail_level.
    # At least 50 support points (for detail level << 1) per raindrop are chosen
    support <- seq(x[1] - se.factor * x[2], x[1] + se.factor * x[2], length.out = max(c(length.out * detail_level, 50)))

    # The values for the likelihood drop are determined: The likelihood for different hypothetical true values
    # minus likelihood for the observed value (i.e. the maximum likelihood) plus the confidence.level quantile of the chi square
    # distribution with one degree of freedeom divided by two.
    # Explanation: -2*(log(L(observed)/L(hypothetical))) is an LRT test and approx. chi^2 with 1 df and significance threshold
    # qchisq(confidence.level, df = 1).
    # That means by adding the confidence.level quantile of the chi square
    # distribution with one degree of freedom (the significance threshold) divided by two,
    # values with l_mu < 0 differ significantly from the observed value.
    threshold <- qchisq(confidence_level, df = 1)/2
    l_mu <- log(dnorm(x[1], mean = support, sd = x[2])) - log(dnorm(x[1], mean = x[1], sd = x[2])) + threshold

    # mirror values for raindrop
    l_mu_mirror <- -l_mu

    # select only likelihood values that are equal or larger than zero,
    # i.e. values that also lie in the confidence interval (using normality assumption)
    sel <- which(l_mu >= 0)

    # Compute area of the likelihood raindrop and scale all values, such that the area is 1 (for every raindrop)
    # That way the height of each raindrop is also porportional to the (unscaled) area of each likelihood raindrop.
    # Otherwise every raindrop would have height equal to confidence.level quantile of the chi square
    # distribution with one degree of freedom divided by two irrespective of its width.
    const <- integrate(function(k) {log(dnorm(x[1], mean = k, sd = x[2])) - log(dnorm(x[1], mean = x[1], sd = x[2])) + threshold},
                       lower = min(support[sel]), upper = max(support[sel]))$value
    l_mu <- l_mu / const
    l_mu_mirror <- l_mu_mirror / const

    # Construct data.frame
    d <- data.frame("support" = c(support[sel], rev(support[sel])), "log_density" = c(l_mu[sel], rev(l_mu_mirror[sel])))

    if(shading) {
      # The number of segments for shading is chosen as follows: 40 segements times the detail_level per drop
      # as default. The minimum count of segments is 20 (If detail_level is << 1), with the exception that
      # if there are too few points for 20 segments then nrow(d)/4 is used (i.e. at least 4 points per segment)
      data.frame(d, "segment" = cut(d$support, max(c(40*detail_level), min(c(20, nrow(d)/4)))))
    } else {
      d
    }
  }

  # compute the max range of all likelihood drops for function ll.
  max.range <- max(abs((es + qnorm(1-(1-confidence_level)/2) * se) - (es - (qnorm(1-(1-confidence_level)/2) * se))))

  # computes all likelihood values and segments (if shading == T). The output is a list, where every element
  # constitutes one study raindop
  res <- apply(cbind(es,se), 1,  FUN = function(x) {ll(x, max.range = max.range)})

  # name every list entry, i.e. raindrop
  names(res) <- ids(group)$ID[ids(group)$type == "study"]

  # The prep.data function prepares the list of raindrops in three ways for plotting (shading of segments):
  # 1) the values are sorted by segments, such that the same segments of each raindrop are joined together
  # 2) segments are renamed with integer values from 1 to the number of segments per raindrop
  # 3) to draw smooth raindops the values at the right hand boundary of each segment have to be the first
  # values at the left hand boundary of the next segment on the right.
  prep.data <- function(res) {
    res <- lapply(res, FUN = function(x) {x <- x[order(x$segment), ]})
    res <- lapply(res, FUN = function(x) {x$segment <- factor(x$segment, labels = 1:length(unique(x$segment))); x})
    res <- lapply(res, FUN = function(x) {
        seg_n <- length(unique(x$segment))
        first <- sapply(2:seg_n, FUN = function(n) {min(which(as.numeric(x$segment)==n))})
        last <-  sapply(2:seg_n, FUN = function(n) {max(which(as.numeric(x$segment)==n))})
        neighbor.top <-   x[c(aggregate(support~segment, FUN = which.max, data=x)$support[1],
                              cumsum(aggregate(support~segment, FUN = length, data=x)$support)[-c(seg_n-1, seg_n)] +
                                aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]), c("support", "log_density")]
        neighbor.bottom <-   x[c(aggregate(support~segment, FUN = which.max, data=x)$support[1],
                                 cumsum(aggregate(support~segment, FUN = length, data=x)$support[-c(seg_n-1, seg_n)])+
                                   aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]) + 1, c("support", "log_density")]
        x[first, c("support", "log_density")] <- neighbor.top
        x[last, c("support", "log_density")] <- neighbor.bottom
        x
      }
    )
  res
  }

  if(shading) {
    res <- prep.data(res)
  }

  # merge the list of raindops in one dataframe for plotting
  res <- plyr::ldply(res)

  # scale the height of each raindrop by the maximum height, such that they
  # fit in their respective plotting region with length 1
  res$log_density <- res$log_density/(2*abs(max(res$log_density)))

  if(summary_symbol == "rain") {
    # prepare likelihood raindrop of the meta-analytical summary effect
    summary.drop <- apply(cbind(summary_es, summary_se), 1,  FUN = function(x) {ll(x, max.range = max.range)})

    # name every list entry, i.e. raindrop
    names(summary.drop) <- ids(group)$ID[ids(group)$type == "summary"]

    if(shading) {
      summary.drop <- prep.data(summary.drop)
    }

    # scale every summary raindrop such that they have maximum size 1
    summary.drop <- lapply(summary.drop, FUN = function (x) {x$log_density <- x$log_density/(2*abs(max(x$log_density))); x})

    # merge the list of raindops in one dataframe for plotting
    summary.drop <- plyr::ldply(summary.drop)

    # alternative scaling of summary drops such that they are scaled relative to all other summary drops (for subgroup analysis),
    # instead of being set to max size (Note: currently not used).
    # summary.drop$log_density <- summary.drop$log_density/(2*abs(max(summary.drop$log_density)))

    # Combine study raindrops and summary raindrop(s) into one data.frame
    res <- rbind(res, summary.drop)
    y_limit <- c(-2, n + 3 * k - 2)
    y_tick_names <- c(as.vector(names), as.vector(summary_name))[order(ids(group)$ID, decreasing = T)]
    y_breaks <- sort(ids(group)$ID, decreasing = T)
  } else {
    if(summary_symbol == "diamond") {
      y_limit <- c(-2, n + 3 * k - 2)
      y_tick_names <- c(as.vector(names), as.vector(summary_name))[order(ids(group)$ID,decreasing = T)]
      y_breaks <- sort(ids(group)$ID, decreasing = T)
      summarydata <- data.frame("x.diamond" = c(summary_es - qnorm(1 - (1 - confidence_level) / 2, 0, 1) * summary_se,
                                      summary_es,
                                      summary_es + qnorm(1 - (1 - confidence_level) / 2, 0, 1) * summary_se,
                                      summary_es),
                                "y.diamond" = c(ids(group)$ID[ids(group)$type=="summary"],
                                                ids(group)$ID[ids(group)$type=="summary"] + 0.3,
                                                ids(group)$ID[ids(group)$type=="summary"],
                                                ids(group)$ID[ids(group)$type=="summary"] - 0.3),
                                "diamond_group" = rep(1:k, times = 4)
                                )
    } else {
      if(summary_symbol == "none") {
        y_limit <- c(0, n + 3 * k - 2)
        y_tick_names <- c(as.vector(names))[order(ids(group)$ID[ids(group)$type=="study"], decreasing =T)]
        y_breaks <- sort(ids(group)$ID[ids(group)$type=="study"], decreasing =T)
      } else {
        warning("Unknown summary_symbol specified. Using default value (rain)")
      }
    }
  }

  # To shade all segments of each raindop symmetrically the min abs(log_density) per raindrop is used
  # as aesthetic to fill the segments. This is necessary because otherwise the first log_density value per
  # segment would be used leading to asymmetrical shading
  if(shading) {
    min.ld <- aggregate(log_density ~ segment + .id, FUN  = function(x) {min(abs(x))}, data = res)
    names(min.ld) <- c("segment", ".id", "min_log_density")
    res <- merge(res, min.ld, sort = F)
  }

  # Set Color palette for shading
  if(!(col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) { # Was a supported color specifed by the user?
    warning("Supported arguments for col are Blues, Greys, Oranges, Greens, Reds, and Purples. Default Blues is used")
    col <- "Blues"
  }
  col <- RColorBrewer::brewer.pal(n=9, name=col)

  # Workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  support <- NULL
  segment <- NULL
  min_log_density <- NULL
  log_density <- NULL
  ci_value <- NULL
  ID <- NULL
  .id <- NULL
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL

  # Create Rainforest plot
  if(shading) {
    p <-
      ggplot(data = res, aes(y = .id, x = support)) +
        geom_polygon(data = res, aes(x = support, y = as.numeric(.id)+log_density,
                                     color = min_log_density, fill = min_log_density,
                                     group = paste(.id, segment)), size = 0.3) +
        geom_point(data = plotdata, shape = "I", aes(x = x, y = ID), col = "snow2", size = 2) +
        geom_line(data = plotdata, col = "snow2", aes(x = ci_value, y = ID, group = ID), size = 0.1)
        if(summary_symbol == "diamond") {
          p <- p + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), color="black", fill = col[9], size = 0.1)
        }
    p <- p +
        scale_fill_gradient(high = col[9], low = col[2], guide = FALSE) +
        scale_color_gradient(high = col[9], low = col[2], guide = FALSE) +
        geom_vline(xintercept = 0, linetype = 2) +
        scale_x_continuous(name = xlab,
                           limits = c(-max(abs(plotdata$ci_value)), max(abs(plotdata$ci_value)))) +
        scale_y_continuous(name = "", limits = y_limit,
                           breaks = y_breaks,
                           labels = y_tick_names,
                           expand = c(0,0)) +
        theme(axis.line = element_line(colour = "black", linetype = "solid"),
              axis.title.x = element_text(size = rel(0.8*text_size), colour = "black"),
              axis.text = element_text(size = rel(0.8*text_size), colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black"),
              panel.background = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line("grey"),
              panel.grid.minor.x = element_line("grey"))
      p
  } else {
    p <-
      ggplot(data = res, aes(y = .id, x = support)) +
        geom_polygon(data = res, aes(x = support, y = as.numeric(.id) + log_density,
                                     group = .id), size = 0.1, col = "black", fill = col [8]) +
        geom_point(data = plotdata, shape = "I", aes(x = x, y = ID), col = "white", size = 3) +
        geom_line(data = plotdata, col = "white", aes(x = ci_value, y = ID, group = ID), size = 0.1)
        if(summary_symbol == "diamond") {
          p <- p + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), color="black", fill=col[9])
        }
      p <- p +
        geom_vline(xintercept = 0, linetype = 2) +
        scale_x_continuous(name = xlab,
                           limits = c(-max(abs(plotdata$ci_value)), max(abs(plotdata$ci_value)))) +
        scale_y_continuous(name = "", limits = y_limit,
                           breaks = y_breaks,
                           labels = y_tick_names,
                           expand = c(0,0)) +
        theme(axis.line = element_line(colour = "black", linetype = "solid"),
              axis.title.x = element_text(size = rel(0.8*text_size), colour = "black"),
              axis.text = element_text(size = rel(0.8*text_size), colour = "black"),
              panel.border = element_rect(fill = NA, colour = "black"),
              panel.background = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.x = element_line("grey"),
              panel.grid.minor.x = element_line("grey"))
      p
  }
}




