#'Power enhanced funnel plot
#'
#'Creates a funnel plot with power regions.
#'
#'The funnel plot is a widely used diagnostic plot in meta-analysis to assess small study effects
#'and in particular publication bias. Small study effects can be defined as association of effect size and standard error,
#'such that small studies (with large standard errors) tend to report larger effect sizes on average than larger studies
#'with preciser estimates. This is widely seen as indicator for potential publication bias, because small studies need to observe
#'larger effect sizes to obtain significant results. The power of a study is a function of its standard error (for a fixed true effect and significance level).
#'Therefore it can be useful to assess if indeed low powered studies tend to report larger and significant results. For this purpose,
#'the power corresponding to a given standard error is computed by using a two-sided Wald test and (by default) the meta-analytic summary effect as true effect.
#'Colored regions of different power levels and a second axis with study level power is shown within the funnel plot.
#'
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}; then effect sizes and standard errors are extracted from \code{x}.
#'@param method character string indicating the method used to compute the meta-analytic summary effect and, for a random effects model,
#'  the between-study variance \eqn{\tau^2}{tau squared}. Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the DerSimonian-Laird method to estimate \eqn{\tau^2}{tau squared}).
#'  Used for \code{contours}, and \code{addev_contours}.
#'@param y_axis character string indicating which y axis should be used in the funnel plot. Available options are "se" (default) for
#'  standard error and "precision" for the reciprocal of the standard error.
#'@param true_effect numeric scalar. Which true effect should be assumed for power caluclations? By default the meta-analytic summary effect is used.
#'@param power_contours character string specifying how different power regions are plotted. Can be eiter "continuous"
#'  (default) or "discrete".
#'@param contours logical scalar indicating if classic funnel plot confidence contours and the summary effect
#'  should be displayed.
#'@param sig_contours logical scalar. Should significance contours be drawn (at the 0.05 or 0.01 level using a Wald test)?
#'@param detail_level numeric scalar. Allows to increase or decrease the plotting detail of contours. Values can be chosen between 0.1 and 10. Default is 1.
#'@param text_size numeric value. Size of text in the funnel plot. Default is 3.
#'@param point_size numeric value. Size of the study points in the funnel plot. Default is 2.
#'@param xlab character string specifying the label of the x axis.
#'@param ylab character string specifying the label of the y axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficients using \code{tanh}.
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@return A funnel plot with power contours is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Create a funnel plot using confidence and significance contours
#' viz_powerfun(x = mozart[1:10, c("d", "se")])
#'@export
viz_powerfun <- function(x, y_axis = "se", method = "FE", true_effect = NULL, power_contours = "continuous",
                       contours = TRUE, sig_contours = FALSE, detail_level = 1,
                       text_size = 3, point_size = 2,
                       xlab = "Effect", ylab = NULL,
                       x_trans_function = NULL, x_breaks = NULL) {
  #'@import ggplot2

  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }
  # input is output of rma (metafor)
  if("rma" %in% class(x)) {
    # extract effect size and standard error
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    # method <- x$method
    if(method != x$method) {
      message("Note: method argument used differs from input object of class rma.uni (metafor)")
    }
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
  summary_es <- metafor::rma.uni(yi = es, sei = se, method = method)$b[[1]]
  summary_se <- sqrt(metafor::rma.uni(yi = es, sei = se, method = method)$vb[[1]])
  summary_tau2 <- metafor::rma.uni(yi = es, sei = se, method = method)$tau2

  # main data for plotting
  plotdata <- data.frame(es, se)

  # Set Color palette for contours
  contours_col <- "Greys"
  col <- RColorBrewer::brewer.pal(n = 9, name = contours_col)

  # detail_level must be between 0.1 and 10
  if(detail_level < 0.1) {
    detail_level <- 0.1
    warning("Argument detail_level too low. Set to minimum value (0.1)")
  }
  if(detail_level > 10) {
    detail_level <- 10
    warning("Argument detail_level too high. Set to minimum value (10)")
  }

  # set intial min and max values for x axis limits
  min_x <- min(plotdata$es)
  max_x <- max(plotdata$es)

  # standard error on the y axis
  if(y_axis =="se") {
    plotdata$y <- se
    max_se <- max(se) + ifelse(length(se) > 1, diff(range(se))*0.1, max(se)*0.1)
    y_limit <- c(0, max_se)

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
      funneldata <- data.frame(x = c(summary_es - stats::qnorm(0.975) * sqrt(max_se^2 + summary_tau2),
                                     summary_es - stats::qnorm(0.975) * sqrt(summary_tau2),
                                     summary_es + stats::qnorm(0.975) * sqrt(summary_tau2),
                                     summary_es + stats::qnorm(0.975) * sqrt(max_se^2 + summary_tau2)),
                               y = c(max_se, 0, 0, max_se))
      # update limit values for x axis
      min_x <- min(c(min_x, min(funneldata$x)))
      max_x <- max(c(max_x, max(funneldata$x)))
    }
  } else {
    if(y_axis == "precision") {
      plotdata$y <- 1/se

      # inital value for upper y axis limit
      max_y <- max(1/se) + ifelse(length(se) > 1, diff(range(1/se))*0.05, 1/se*0.05)
      min_y <- min(1/se) - ifelse(length(se) > 1, diff(range(1/se))*0.05, 1/se*0.05)

      if(is.null(ylab)) {
        ylab <- "Precision (1/SE)"
      }

      # determine significance contours
      if(sig_contours == TRUE) {
        n_support <- 200 * detail_level

        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec_0.05 <- stats::qnorm(0.975)*(1/prec)
        x_prec_0.01 <- stats::qnorm(0.995)*(1/prec)

        sig_funneldata <- data.frame(x = c(-x_prec_0.01, rev(x_prec_0.01),
                                           x_prec_0.05, rev(-x_prec_0.05)),
                                     y = c(prec, rev(prec), prec, rev(prec)))


        min_x <- min(c(min_x, min(sig_funneldata$x)))
        max_x <- max(c(max_x, max(sig_funneldata$x)))


      }
      # determine classic funnel contours
      if(contours == TRUE) {
        n_support <- 200*detail_level

        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec <- stats::qnorm(0.975)*sqrt((1/prec)^2 + summary_tau2)
        funneldata <- data.frame(x = rep(summary_es, times = n_support*2) + c(- x_prec, rev(x_prec)),
                                 y = c(prec, rev(prec)))
        # update x axis limit
        min_x <- min(c(min_x, min(funneldata$x)))
        max_x <- max(c(max_x, max(funneldata$x)))
      }
      y_limit <- c(min_y, max_y)
    } else {
      stop("y_axis argument must be either se or precision")
    }
  }

  x_limit <- c(min_x - diff(c(min_x, max_x))*0.05, max_x + diff(c(min_x, max_x))*0.05)

  # Compute power contours
    if(is.null(true_effect)) {
      true_effect <- summary_es
    }


    yseq <- seq(from = y_limit[1], to = y_limit[2], length.out = 1000)
    if(y_axis == "se") {
      power <- (1 - stats::pnorm(stats::qnorm(0.975) * yseq, abs(true_effect), yseq)) + stats::pnorm(stats::qnorm(0.025) * yseq, abs(true_effect), yseq)
    } else {
      power <- 1 - stats::pnorm(1.96 * 1/yseq, abs(true_effect), 1/yseq) + stats::pnorm(stats::qnorm(0.025) * 1/yseq, abs(true_effect), 1/yseq)
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
        stop('argument for power_contours must be either "discrete" or "continuous".')
      }
    }

    power_col <- RColorBrewer::brewer.pal(n = 9, name = "RdYlGn")

    # compute statistics
    study_power <- (1 - stats::pnorm(stats::qnorm(0.975) * se, abs(true_effect), se)) + stats::pnorm(stats::qnorm(0.025) * se, abs(true_effect), se)
    med_power <- paste(round(stats::median(study_power)*100, 1), "%", sep = "")
    expected <- sum(study_power)
    observed <- sum(2*(1 - stats::pnorm(abs(es/se))) <= 0.05)
    c2 <- (observed - expected)^2/expected + (observed - expected)^2/(length(study_power) - expected)
    p_tes <- round(1 - stats::pchisq(c2, df = 1), 3)

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
  # Construct plot
  p <- ggplot(data = plotdata, aes(x = es, y = y))
  if(power_contours == "continuous") {
    p <- p +
      geom_raster(data = power_grid, aes(x = x, y = y, fill = fill), alpha = 1) +
      scale_fill_gradientn(name = "Power", colours = power_col, limits = c(0.0499, 1), breaks = c(0.05, 0.2, 0.4, 0.6, 0.8, 1),
                           guide = guide_colorbar(draw.ulim = FALSE, draw.llim = FALSE, barwidth = 10))
  } else {
    p <- p +
      geom_rect(inherit.aes = FALSE, data = power_recs, aes(xmin = xstart, xmax = xend, ymin=ystart,
                                                            ymax = yend, fill = fill), alpha = 1) +
      scale_fill_manual(name ="", values = power_col)
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
    p <-
      p + scale_y_reverse(name = ylab,
                          sec.axis =
                            dup_axis(~.,
                                     name = "Power",
                                     labels = function(x) {
                                       paste(round((1-stats::pnorm((stats::qnorm(0.975)*x - true_effect)/x) + stats::pnorm((-stats::qnorm(0.975)*x - true_effect)/x)) * 100, 1), "%", sep = "")}))
  } else {
    if(y_axis == "precision") {
      p <-
        p + scale_y_continuous(name = ylab,
                               sec.axis =
                                 dup_axis(~.,
                                          name = "Power",
                                          labels = function(x) {
                                            paste(round((1-stats::pnorm((stats::qnorm(0.975)*1/x - true_effect)/(1/x)) + stats::pnorm((-stats::qnorm(0.975)*1/x - true_effect)/(1/x))) * 100, 1), "%", sep = "")}))
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
  p <- p +
    labs(caption = bquote(
           paste(delta, " = ", .(round(true_effect, 3)), ", ",
                 tilde(x)[power], " = ", .(med_power), ", ",
                 "Expected = ",  .(round(expected, 2)), ", ",
                 "Observed = ", .(observed), ", ",
                 chi^2, "(1) = ", .(round(c2, 2)),  ", ",
                 p[TES], " = ", .(p_tes), sep = ""))
    ) +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  p
}

