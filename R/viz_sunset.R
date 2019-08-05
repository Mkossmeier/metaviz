#'Sunset (power-enhanced) funnel plot
#'
#'Creates a funnel plot with power regions and computes power-related statistics.
#'
#'The funnel plot is the most widely used diagnostic plot in meta-analysis, primarily to assess
#'small-study effects. The sunset (power-enhanced) funnel plot incorporates study-level power information
#'in the funnel display. This directly allows to examine the power studies had to detect an effect of interest
#'(e.g., the observed meta-analytic summary effect), whether funnel plot asymmetry is driven by underpowered but significant studies, and to visually
#'assess if there is an excess of low-powered significant effects in the meta-analysis
#'(conceptually related to the test of excess significance, Ioannidis & Trikalinos, 2007).
#'For effect sizes assumed to be normally distributed (e.g., Cohen d, log OR),
#'the power corresponding to a given standard error is computed by using a two-sided Wald test and (by default) the meta-analytic
#'summary effect as assumed true effect. Colored regions of different power levels and a second axis with study level power are shown in the
#'funnel plot. In addition, power-related statistics are shown: a) The median power of all studies, b) the true effect size necessary such that
#'the median power of the studies would have been 33\% or 66\%, c) results of a test of excess significance (Ioannidis & Trikalinos, 2007), and d)
#'the R-Index for expected replicability (Schimmack, 2016).
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}; then effect sizes and standard errors are extracted from \code{x}.
#'@param y_axis character string indicating which y axis should be used in the funnel plot. Available options are "se" (default) for
#'  standard error and "precision" for the reciprocal of the standard error.
#'@param true_effect numeric scalar. Which true effect should be assumed for power calculations? The default is \code{NULL},
#'  for which the meta-analytic summary effect (fixed-effect model) is used.
#'@param sig_level logical scalar. For which significance level alpha should the study power be computed?
#'@param power_stats logical scalar. Should power-related statistics be computed and printed in the caption of the plot? (see details)
#'@param power_contours character string specifying how different power regions are plotted. Can be either "continuous"
#'  or "discrete" (default).
#'@param contours logical scalar indicating if classic funnel plot confidence contours and the summary effect
#'  should be displayed.
#'@param sig_contours logical scalar. Should significance contours be drawn (at the 0.05 or 0.01 level using a two-sided Wald test)?
#'@param text_size numeric value. Size of text in the funnel plot. Default is 3.
#'@param point_size numeric value. Size of the study points in the funnel plot. Default is 2.
#'@param xlab character string specifying the label of the x axis.
#'@param ylab character string specifying the label of the y axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficients using \code{tanh}.
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@param y_breaks numeric vector of values for the breaks on the y-axis.
#'@param x_limit numeric vector of length two with user specified x axis limits.
#'@param y_limit numeric vector of length two with user specified y axis limits.
#'@references Ioannidis, J. P., & Trikalinos, T. A. (2007). An exploratory test for an excess of significant
#'  findings. \emph{Clinical Trials}, \emph{4}, 245-253.
#'@references Schimmack, U. (2016). The replicability-index: Quantifying statistical research integrity.
#'  Retrieved from https://replicationindex.wordpress.com/2016/01/31/a-revised-introduction-to-the-r-index/
#'@return A power enhanced ("sunset") funnel plot is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Create a power-enhanced ("sunset") funnel plot using confidence and significance contours
#' viz_sunset(x = homeopath[, c("d", "se")], contours = TRUE)
#'@export
viz_sunset <- function(x, y_axis = "se", true_effect = NULL,
                         sig_level = 0.05, power_stats = TRUE,
                         power_contours = "discrete", contours = FALSE, sig_contours = TRUE,
                         text_size = 3, point_size = 2, xlab = "Effect", ylab = NULL,
                         x_trans_function = NULL, x_breaks = NULL, y_breaks = NULL,
                         x_limit = NULL, y_limit = NULL) {
  #'@import ggplot2

  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }
  # input is output of rma (metafor)
  if("rma" %in% class(x)) {
    # extract effect size and standard error
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
  } else {
    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        x <- x[stats::complete.cases(x), ]
      }
      # check if input is numeric
      if(!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
        stop("Input argument has to be numeric; see help(viz_powerfun) for details.")
      }
      # check if there are any negative standard errors
      if(!all(x[, 2] > 0)) {
        stop("Non-positive standard errors supplied")
      }
      # extract effect sizes and standard errors
      es <- x[, 1]
      se <- x[, 2]

    } else {
      stop("Unknown input argument; see help(viz_powerfun) for details.")
    }
  }

  # Compute meta-analytic summary using rma.uni (metafor) and supplied method
  k <- length(es)
  if(is.null(true_effect) | contours == TRUE) {
    summary_es <- metafor::rma.uni(yi = es, sei = se, method = "FE")$b[[1]]
    summary_se <- sqrt(metafor::rma.uni(yi = es, sei = se, method = "FE")$vb[[1]])
  }

  # main data for plotting
  plotdata <- data.frame(es, se)

  # Set Color palette for contours
  contours_col <- "Greys"
  col <- RColorBrewer::brewer.pal(n = 9, name = contours_col)

  # set intial min and max values for x axis limits
  min_x <- min(plotdata$es)
  max_x <- max(plotdata$es)

  if(!is.numeric(sig_level) || length(sig_level) != 1 || sig_level <= 0 || sig_level >= 1) stop("sig_level must be a numeric value greater than 0 and smaller than 1")

  # standard error on the y axis
  if(y_axis =="se") {
    plotdata$y <- se
    if(is.null(y_limit)) {
      max_se <- max(se) + ifelse(diff(range(se)) != 0, diff(range(se))*0.1, max(se)*0.1)
      y_limit <- c(0, max_se)
    } else {
      max_se <- max(y_limit)
    }

    if(is.null(ylab)) {
      ylab <- "Standard Error"
    }

    # determine significance contours
    if(sig_contours == TRUE) {
      sig_funneldata <- data.frame(x = c(-stats::qnorm(0.975) * max_se, 0,
                                         stats::qnorm(0.975) * max_se,
                                         stats::qnorm(0.995) * max_se, 0,
                                         -stats::qnorm(0.995) * max_se),
                                   y = c(max_se, 0, max_se, max_se, 0, max_se))

      min_x <- min(c(min_x, min(sig_funneldata$x)))
      max_x <- max(c(max_x, max(sig_funneldata$x)))
    }
    # determine classic funnel contours
    if(contours == TRUE) {
      funneldata <- data.frame(x = c(summary_es - stats::qnorm(0.975) * sqrt(max_se^2),
                                     summary_es,
                                     summary_es,
                                     summary_es + stats::qnorm(0.975) * sqrt(max_se^2)),
                               y = c(max_se, 0, 0, max_se))
      # update limit values for x axis
      min_x <- min(c(min_x, min(funneldata$x)))
      max_x <- max(c(max_x, max(funneldata$x)))
    }
  } else {
    if(y_axis == "precision") {
      plotdata$y <- 1/se

      if(is.null(y_limit)) {
        # inital value for upper y axis limit
        max_y <- max(1/se) + ifelse(diff(range(se)) != 0, diff(range(1/se))*0.05, 1/se*0.05)
        min_y <- min(1/se) - ifelse(diff(range(se)) != 0, diff(range(1/se))*0.05, 1/se*0.05)
      } else {
        max_y <- max(y_limit)
        min_y <- min(y_limit)
      }
      if(is.null(ylab)) {
        ylab <- "Precision (1/SE)"
      }

      # determine significance contours
      if(sig_contours == TRUE) {
        n_support <- 200

        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec_0.05 <- stats::qnorm(0.975)*(1/prec)
        x_prec_0.01 <- stats::qnorm(0.995)*(1/prec)

        sig_funneldata <- data.frame(x = c(-x_prec_0.01, rev(x_prec_0.01),
                                           x_prec_0.05, rev(-x_prec_0.05)),
                                     y = c(prec, rev(prec), prec, rev(prec)))

        if(is.null(x_limit)) {
          min_x <- min(c(min_x, min(sig_funneldata$x)))
          max_x <- max(c(max_x, max(sig_funneldata$x)))
        } else {
          min_x <- min(x_limit)
          max_x <- max(x_limit)
        }

      }
      # determine classic funnel contours
      if(contours == TRUE) {
        n_support <- 200

        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec <- stats::qnorm(0.975)*sqrt((1/prec)^2)
        funneldata <- data.frame(x = rep(summary_es, times = n_support*2) + c(- x_prec, rev(x_prec)),
                                 y = c(prec, rev(prec)))
        # update x axis limit
        if(is.null(x_limit)) {
          min_x <- min(c(min_x, min(funneldata$x)))
          max_x <- max(c(max_x, max(funneldata$x)))
        } else {
          min_x <- min(x_limit)
          max_x <- max(x_limit)
        }

      }
      if(is.null(y_limit)) y_limit <- c(min_y, max_y)
    } else {
        stop("y_axis argument must be either se or precision")
    }
  }
  if(is.null(x_limit)) {
    x_limit <- c(min_x - diff(c(min_x, max_x))*0.05, max_x + diff(c(min_x, max_x))*0.05)
  }
  # Compute power contours
    if(is.null(true_effect)) {
      true_effect <- summary_es
    }

    yseq <- seq(from = y_limit[1], to = y_limit[2], length.out = 1000)

    if(y_axis == "se") {
      power <- (1 - stats::pnorm(stats::qnorm(1 - sig_level/2) * yseq, abs(true_effect), yseq)) + stats::pnorm(stats::qnorm(sig_level/2) * yseq, abs(true_effect), yseq)
    } else {
      power <- 1 - stats::pnorm(stats::qnorm(1 - sig_level/2) * 1/yseq, abs(true_effect), 1/yseq) + stats::pnorm(stats::qnorm(sig_level/2) * 1/yseq, abs(true_effect), 1/yseq)
    }

    if(power_contours == "discrete") {
      power_y <- numeric(10)
      steps <- c(1:4, 6:10)
      for(i in 1:9) {
        power_y[i+1] <- yseq[which.min(abs(power - steps[i]/10))]
      }

      if(y_axis == "se") {
        power_y[1] <- max(y_limit)
      } else {
        power_y[10] <- max(y_limit)
      }

      if(true_effect == 0) {
        power_y[1] <- y_limit[1]
        power_y[2:10] <- y_limit[2]
      }

      power_recs <- data.frame(xstart = x_limit[1], xend = x_limit[2],
                               ystart = power_y[1:9], yend = power_y[2:10],
                               fill = factor(paste("Power", c(0, steps[-length(steps)])*10, "-", steps*10)))

    } else {
      if(power_contours == "continuous") {
        power_grid <- data.frame(x = rep(x_limit, each = 1000),
                                 y = rep(yseq, times = 2), fill = rep(power, times = 2))
      } else {
        stop('Argument for power_contours must be either "discrete" or "continuous".')
      }
    }

    power_col <- RColorBrewer::brewer.pal(n = 9, name = "RdYlGn")

    # function to find dxx%
    dpower <- function(se, sig_level = 0.05, target_power = 0.33) {
      if(sig_level >= target_power) {
        "n.a."
      } else {
        d <- stats::uniroot(function(x) stats::median((1 - stats::pnorm(stats::qnorm(1 - sig_level/2) * se, abs(x), se)) + stats::pnorm(stats::qnorm(sig_level/2) * se, abs(x), se)) - target_power,
                lower = 0, upper = 10, extendInt = "upX", tol = 0.001, maxiter = 2000)$root
        d <- round(d, 2)
        ifelse(d >= 10, "> 10", ifelse(d <= -10, "< -10", d))
      }
    }

    if(power_stats == TRUE) {
       study_power <- (1 - stats::pnorm(stats::qnorm(1 - sig_level/2) * se, abs(true_effect), se)) + stats::pnorm(stats::qnorm(sig_level/2) * se, abs(true_effect), se)
       med_power <- paste(round(stats::median(study_power)*100, 1), "%", sep = "")
       expected <- sum(study_power)
       observed <- sum(2*(1 - stats::pnorm(abs(es/se))) <= sig_level)
       c2 <- (observed - expected)^2/expected + (observed - expected)^2/(length(study_power) - expected)
       p_tes <- round(1 - stats::pchisq(c2, df = 1), 3)
       R <- 2*stats::median(study_power) - observed/length(se)
       R <- ifelse(R < 0 , 0, ifelse(R > 1, 1, R))
       R <- paste0(round(R*100, 1), "%")
       d33 <- dpower(se, sig_level = sig_level, target_power = 0.33)
       d66 <- dpower(se, sig_level = sig_level, target_power = 0.66)
    }

  if(!is.null(x_trans_function) && !is.function(x_trans_function)) {
    warning("Argument x_trans_function must be a function; input ignored.")
    x_trans_function <- NULL
  }

  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  y <- NULL
  fill <- NULL
  xstart <- NULL
  xend <- NULL
  ystart <- NULL
  yend <- NULL
  x_min <- NULL
  x_max <- NULL

  # Create plot
  p <- ggplot(data = plotdata, aes(x = es, y = y))
  if(power_contours == "continuous") {
    p <- p +
      geom_raster(data = power_grid, aes(x = x, y = y, fill = fill), alpha = 1) +
      scale_fill_gradientn(name = "Power", colours = power_col, limits = c(0.0499, 1), breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1),
                           guide = guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE, barwidth = 10))
  } else {
    p <- p +
      geom_rect(inherit.aes = FALSE, data = power_recs, aes(xmin = xstart, xmax = xend, ymin = ystart,
                                                            ymax = yend, fill = fill), alpha = 1) +
      scale_fill_manual(name = "", values = power_col)
  }

  if(sig_contours == TRUE && y_axis == "se") {
    p <- p +
      geom_polygon(data = sig_funneldata, aes(x = x, y = y), fill = col[9], alpha = 0.6) +
      geom_path(data = sig_funneldata, aes(x = x, y = y))
  } else {
    if(sig_contours == TRUE && y_axis == "precision") {
      p <- p +
        geom_polygon(data = sig_funneldata, aes(x = x, y = y), fill = col[9], alpha = 0.6) +
        geom_path(data = sig_funneldata, aes(x = x, y = y))
    }
  }
  if(contours == TRUE) {
    p <- p +
      geom_path(data = funneldata, aes(x = x, y = y)) +
      geom_vline(xintercept = summary_es)
  }
  if(y_axis == "se") {
    if(is.null(y_breaks)) {
      p <-
        p + scale_y_reverse(name = ylab,
                            sec.axis =
                              dup_axis(~.,
                                       name = "Power",
                                       labels = function(x) {
                                         paste(round((1-stats::pnorm((stats::qnorm(1-sig_level/2)*x - true_effect)/x) + stats::pnorm((-stats::qnorm(1-sig_level/2)*x - true_effect)/x)) * 100, 1), "%", sep = "")}))
      }
    else {
      p <-
        p + scale_y_reverse(name = ylab, breaks = y_breaks,
                            sec.axis =
                              dup_axis(~.,
                                       name = "Power",
                                       labels = function(x) {
                                         paste(round((1-stats::pnorm((stats::qnorm(1-sig_level/2)*x - true_effect)/x) + stats::pnorm((-stats::qnorm(1-sig_level/2)*x - true_effect)/x)) * 100, 1), "%", sep = "")}))

    }
  } else {
    if(y_axis == "precision") {
      if(is.null(y_breaks)) {
        p <-
          p + scale_y_continuous(name = ylab,
                                 sec.axis =
                                   dup_axis(~.,
                                            name = "Power",
                                            labels = function(x) {
                                              paste(round((1-stats::pnorm((stats::qnorm(1-sig_level/2)*1/x - true_effect)/(1/x)) + stats::pnorm((-stats::qnorm(1-sig_level/2)*1/x - true_effect)/(1/x))) * 100, 1), "%", sep = "")}))
      } else {
        p <-
          p + scale_y_continuous(name = ylab, breaks = y_breaks,
                                 sec.axis =
                                   dup_axis(~.,
                                            name = "Power",
                                            labels = function(x) {
                                              paste(round((1-stats::pnorm((stats::qnorm(1-sig_level/2)*1/x - true_effect)/(1/x)) + stats::pnorm((-stats::qnorm(1-sig_level/2)*1/x - true_effect)/(1/x))) * 100, 1), "%", sep = "")}))

      }
    }
  }
  p <- p + geom_point(size = point_size , fill = "white", shape = 21, col = "black", alpha = 1)
  if(!is.null(x_trans_function)) {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)})
    } else {
      p <- p +
        scale_x_continuous(name = xlab,
                           labels = function(x) {round(x_trans_function(x), 3)},
                           breaks = x_breaks)
    }
  } else {
    if(is.null(x_breaks)) {
      p <- p +
        scale_x_continuous(name = xlab)
    } else {
      p <- p +
        scale_x_continuous(breaks = x_breaks,
                           name = xlab)
    }
  }
  p <- p +
    coord_cartesian(xlim = x_limit,
                    ylim = y_limit, expand = F)
  if(power_stats == TRUE) {
    if(is.null(x_trans_function)) {
      p <- p +
        labs(caption = bquote(
               paste(alpha, " = ", .(sig_level), ", ",
                     delta , " = ", .(round(true_effect, 2)), " | ",
                     med[power], " = ", .(med_power), ", ",
                     d[33*'%'], " = ", .(d33), ", ",
                     d[66*'%'], " = ", .(d66), " | ",
                     "E = ",  .(round(expected, 2)), ", ",
                     "O = ", .(observed), ", ",
                     p[TES], .(ifelse(p_tes == 0, " < ", " = ")), .(ifelse(p_tes == 0, "0.001", p_tes)), ", ",
                     "R-Index = ", .(R), sep = "")))
    } else {
      p <- p +
        labs(caption = bquote(
          paste(alpha, " = ", .(sig_level), ", ",
                delta , " = ", .(round(x_trans_function(true_effect), 2)), " | ",
                med[power], " = ", .(med_power), ", ",
                d[33*'%'], " = ", .(round(x_trans_function(d33), 2)), ", ",
                d[66*'%'], " = ", .(round(x_trans_function(d66), 2)), " | ",
                "E = ",  .(round(expected, 2)), ", ",
                "O = ", .(observed), ", ",
                p[TES], .(ifelse(p_tes == 0, " < ", " = ")), .(ifelse(p_tes == 0, "0.001", p_tes)), ", ",
                "R-Index = ", .(R), sep = "")))
    }
  }
  p <- p +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  p
}

