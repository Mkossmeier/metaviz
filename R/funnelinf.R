#' Visual funnel plot inference for meta-analysis
#'
#' Creates a funnel plot showing the real data alongside funnel plots showing simulated data
#' under the null hypothesis to conduct visual funnel plot inference.
#'
#' Funnel plots are widely used in meta-analysis to detect small study effects and in particular publication bias.
#' However, interpretations of funnel plots often lead to false conclusions (e.g., Terrin, Schmid, and Lau, 2005). Visual inference
#' (Buja et al. 2009; Majumder, Hofmann, and Cook 2013) can help to improve the validity of such funnel plot based conclusions.
#' If the alternative hypothesis is true (e.g., small study effects led to an assymtetric funnel plot), the funnel plot showing the real data
#' should be identifiable when presented alongside funnel plots of simulated data under the null hypothesis. Only if this is possible,
#' conclusion based on the funnel plot might be warranted.
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column.
#'@param group a factor indicating the subgroup of each study.
#'@param group_permut logical scalar indicating if subgroup membership should be permutated
#'  in the null plots. Ignored if no group is supplied.
#'@param n integer specifying the absolute number of plots (null plots plus plot of real data)
#'@param y_axis Which y axis should be used in the funnel plot? Available options are se (default) for
#'  standard error and precision for the reciprocal of the standard error.
#'@param null_model Which meta-analytic model should be used to simulate the effect sizes for the null plots?
#'  Available options are FEM for the fixed effect model and REM for the random-effects model (using the
#'  DerSimonian-Laird method to estimate the between-study variance \eqn{\tau^2}{tau squared})
#'@param contours logical scalar indicating if classic funnel plot contours should be displayed (i.e., summary effect +/-
#'  rnorm(0.975) * SE).
#'@param sig_contours logical scalar. Should significance contours should be drawn? Significance contours show which combination of
#' effect size and standard error lead to p-values smaller than 0.05 or 0.01 (using a Wald test).
#'@param trim_and_fill logical scalar Should imputed studies by the trim and fill method be displayed?
#'@param trim_and_fill_side On which side should studies be imputed by the trim and fill method? Available options are
#'right or left.
#'@param egger logical scalar. Should Egger's regression line be drawn?
#'@param show_solution logical scalar. Should the plot with the real data be highlighted? Default is FALSE.
#'@references Buja, A., Cook, D., Hofmann, H., Lawrence, M., Lee, E. K., Swayne, D. F., & Wickham, H. (2009).
#'  Statistical inference for exploratory data analysis and model diagnostics.
#'  \emph{Philosophical Transactions of the Royal Society of London A: Mathematical, Physical and Engineering Sciences},
#'  \emph{367}, 4361-4383.
#'@references Majumder, M., Hofmann, H., & Cook, D. (2013). Validation of visual statistical inference, applied to linear models.
#'  \emph{Journal of the American Statistical Association}, \emph{108}, 942-956.
#'@references Terrin, N., Schmid, C. H., & Lau, J. (2005). In an empirical evaluation of the funnel plot, researchers could not
#'  visually identify publication bias. \emph{Journal of clinical epidemiology}, \emph{58}, 894-901.
#'@return A lineup of n (20 by default) funnel plots; one showing the real data and n-1 showing
#'  simulated data under the null hypothesis
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Plotting a funnel plot lineup with the mozart data to conduct visual funnel plot inference considering subgroups (for details, see help(mozart)):
#' funnelinf(x = mozart[, c("d", "se")],
#' group = mozart[, "rr_lab"],
#' group_permut = TRUE, null_model = "REM")
#'
#' # Plotting a funnel plot lineup with the brain volume data to conduct visual funnel plot inference considering heterogeneity (for details, see help(brainvol)):
#' funnelinf(x = brainvol[, c("z", "z_se")], null_model = "FEM")
#' @export
funnelinf <- function(x, group = NULL, group_permut = FALSE, n = 20, y_axis = "se", null_model = "FEM",
                      contours = TRUE, sig_contours = TRUE,
                      trim_and_fill = FALSE, trim_and_fill_side = "left",
                      egger = FALSE, show_solution = FALSE) {

  #'@import ggplot2
  #'@import stats
  #'@import nullabor

  es <- x[, 1]
  se <- x[, 2]
  k <- length(es)
  summary_es <- metafor::rma.uni(yi = es, sei=se, method ="FE")$b[1] # for funnel plot contours

  if(null_model == "FEM") {
    summary_es_simul <- summary_es
    se_simul <- se
  } else {
    if(null_model == "REM") {
      summary_es_simul <- metafor::rma.uni(yi = es, sei=se, method ="DL")$b[1]
      se_simul <- sqrt(se^2 + metafor::rma.uni(yi = es, sei=se, method ="DL")$tau2[1])
    } else {
      stop("Supported method arguments are either FEM or REM")
    }
  }

  # Sample effect sizes from a normal distribution mu = summary effect and sd = se_simul
  if(is.null(group)) {
    data <- data.frame(es, se)
  } else {
    data <- data.frame(es, se, group)
  }
  if(n > 1) {
    x <- nullabor::lineup(
      nullabor::null_dist("es", dist="normal", params = list(mean = summary_es_simul, sd = se_simul)),
      data, n = n)
    solution <- attributes(x)$pos
  } else {
    if(n == 1) {
      x <- data.frame(data, .sample = 1)
      solution <- 1
    } else {
      stop("n has to be positive")
    }
  }

  if(!is.null(group) & group_permut) {
    plotdata <- plyr::ddply(x, plyr::.(.sample), plyr::here(plyr::mutate), group_permut = sample(group, size = k, replace = F))
    plotdata[plotdata$.sample == solution, "group_permut"] <- group
    plotdata$group <- plotdata$group_permut
  } else {
    plotdata <- x
  }

  # Add column that shows real data
  plotdata <- data.frame(plotdata, real_data = plotdata$.sample == solution)

  # for lower y axis limit
  max_se <- max(se) + 0.2

  col <- RColorBrewer::brewer.pal(n = 9, name="Blues")

  # set intial min and max values for x axis limits
  min_x <- min(plotdata$es)
  max_x <- max(plotdata$es)

  if(trim_and_fill) {
    trimnfill <- function(es, se, group = NULL, side = "left") {
      if(side == "right") {
        es <- -es
      }
      if(side != "right" & side != "left") {
        stop("trim_and_fill_side argument must be either left or right")
      }

      mean_func <- function(es, se) {
        metafor::rma.uni(yi = es, sei=se, method ="FE")$b[1]
      }

      k0_func <- function(es, se, summary_es) {
        n <- length(es)
        Tn <- sum(rank(abs(es - summary_es))[sign(es - summary_es) > 0])
        round(max((4*Tn - n * (n + 1))/(2*n - 1), 0) , 0) # L0_plus estimator from Duval & Tweedie (2000)
      }

      # trim step
      summary_es_init <- mean_func(es, se)
      k0 <- k0_func(es = es, se = se, summary_es = summary_es_init)
      eps <- 1
      iter <- 0
      while(eps > 0.01 | iter < 20) {
        iter <- iter + 1
        es_ord <- es[order(es, decreasing = T)]
        se_ord <- se[order(es, decreasing = T)]

        if(k0 > 0) {
          es_ord <- es_ord[-(1:k0)]
          se_ord <- se_ord[-(1:k0)]
        }

        summary_es_new <- mean_func(es_ord, se_ord)
        k0 <- k0_func(es = es, se = se, summary_es = summary_es_new)
        eps <- abs(summary_es_init - summary_es_new)
        summary_es_init <- summary_es_new
      }

      if(iter == 19) {
        warning("Trim and fill algorithm did not converge after 10 iterations")
      }

      # fill step
      if(k0 > 0) {
        es_ord <- es[order(es, decreasing = T)]
        se_ord <- se[order(es, decreasing = T)]
        if(!is.null(group)) {
          group_ord <- group[order(es, decreasing = T)]
          group_fill <- group_ord[1:k0]
        }
        if(side == "right") {
          es_fill <- -(summary_es_new + (summary_es_new - es_ord[1:k0]))
          summary_es_init <- - summary_es_init # check if fixed effect mean M of values times -1 is really -M
        } else {
          es_fill <- summary_es_new + (summary_es_new - es_ord[1:k0])
        }

        se_fill <- se_ord[1:k0]

        if(is.null(group)) {
          data.frame(es_fill, se_fill, summary_es_init)
        } else {
          data.frame(es_fill, se_fill, group_fill, summary_es_init)
        }
      }
    }

    side <- trim_and_fill_side

    if(is.null(group)) {
      tnfdata <- plyr::dlply(plotdata, plyr::.(.sample), plyr::here(plyr::summarize), fill = trimnfill(es, se, side = side))
    } else {
      tnfdata <- plyr::dlply(plotdata, plyr::.(.sample), plyr::here(plyr::summarize), fill = trimnfill(es, se, group, side = side))
    }

    buf <- data.frame()
    for(i in 1:n) {
      buf  <- rbind(buf, data.frame(
        tnfdata[[i]][[1]], .sample = rep(i, times = dim(tnfdata[[i]])[1]))
      )
    }

    tnfdata <- buf

    if(dim(tnfdata)[1] > 0) {
      if(is.null(group) ) {
        names(tnfdata) <- c("es", "se", "tnf_summary", ".sample")
      } else {
        names(tnfdata) <- c("es", "se", "group", "tnf_summary", ".sample")
      }



      # update limit values for x axis
      min_x <- min(c(min_x, min(tnfdata$es)))
      max_x <- max(c(max_x, max(tnfdata$es)))
    }

  }



  # standard error on the y axis
  if(y_axis =="se") {
    # determine significance contours
    if(sig_contours) {
      x.1  <- c(-qnorm(0.95) * max_se, 0, qnorm(0.95) * max_se)
      x.05 <- c(-qnorm(0.975) * max_se, 0, qnorm(0.975) * max_se)
      x.01  <- c(-qnorm(0.995) * max_se, 0, qnorm(0.995) * max_se)
      y <- c(max_se, 0, max_se)
      sig_funneldata <- data.frame(x.1 = rep(x.1, times = n), x.05 = rep(x.05, times = n),
                                   x.01 = rep(x.01, times = n), y = rep(y, times = n), .sample = rep(1:n, each = 3))
      # update limit values for x axis
      min_x <- min(c(min_x, x.01[1]))
      max_x <- max(c(max_x, x.01[3]))
    }
    # determine classic funnel contours
    if(contours) {
      M <- plyr::ddply(x, plyr::.(.sample), plyr::summarize, M = sum((1/se^2)*es)/(sum(1/se^2)))[, 2]
      x <- c(M - qnorm(0.975) * max_se, M, M + qnorm(0.975) * max_se)
      y <- rep(c(max_se, 0, max_se), each = n)
      funneldata <- data.frame(x = x, y = y, .sample = rep(1:n, times = 3))

      # update limit values for x axis
      min_x <- min(c(min_x, min(x)))
      max_x <- max(c(max_x, max(x)))

      meandata <- data.frame (M = M, .sample = 1:n)

    }
    # determine egger's regression line
    if(egger) {
      plotdata <- data.frame(plotdata, "z" = (plotdata$es)/plotdata$se)
      plotdata <- data.frame(plotdata, "prec" = 1/plotdata$se)
      radial_intercept <- plyr::ddply(plotdata, plyr::.(.sample), plyr::summarize, coef = coef(lm(z ~ prec))[1])[, 2]
      radial_slope <- plyr::ddply(plotdata, plyr::.(.sample), plyr::summarize, coef = coef(lm(z ~ prec))[2])[, 2]
      # note SE = d*(1/intercept) - slope/intercept, i.e. if egger's regression funnel plot has negative slope, this indicates PB
      eggerdata <- data.frame(intercept = radial_slope/radial_intercept,
                              slope = -1/radial_intercept, #minus slope because ggplot2 does not adjust slope of abline for reversed y axis
                              .sample = 1:n)
    }
  } else {
    if(y_axis == "precision") {
      # inital value for upper y axis limit
      max_y <- max(1/se) + 0.2

      # determine significance contours
      if(sig_contours) {
        n_support <- 200
        x_left <- seq(from = -(max(abs(es))+2), to = -0.1, length.out = n_support)
        x_right <- seq(from = 0.1, to = max(abs(es))+2, length.out = n_support)

        y_left0.1 <- -qnorm(0.95)/x_left
        y_right0.1 <- qnorm(0.95)/x_right

        y_left0.05 <- -qnorm(0.975)/x_left
        y_right0.05 <- qnorm(0.975)/x_right

        y_left0.01 <- -qnorm(0.995)/x_left
        y_right0.01 <- qnorm(0.995)/x_right

        sig_funneldata <- data.frame(x = rep(c(x_left, x_right), times = n),
                                     y0.1 =   rep(c(y_left0.1, y_right0.1), times = n),
                                     y0.05 =  rep(c(y_left0.05, y_right0.05), times = n),
                                     y0.01 =  rep(c(y_left0.01, y_right0.01), times = n),
                                     .sample = rep(1:n, each = 2*n_support)
        )
        # update x axis limit
        min_x <- min(c(min_x, min(x_left)))
        max_x <- max(c(max_x, max(x_right)))

        # update y axis limit
        # max_y <- min(c(max(sig_funneldata$y0.1), max(sig_funneldata$y0.05), max(sig_funneldata$y0.01)))

      }
      # determine classic funnel contour
      if(contours) {
        M <- plyr::ddply(x, plyr::.(.sample), plyr::summarize, M = sum((1/se^2)*es)/(sum(1/se^2)))[, 2]

        n_support <- 200
        x_left <- seq(from = -(max(abs(es)) + 2), to = -0.1, length.out = n_support)
        x_right <- seq(from = 0.1, to = max(abs(es)) + 2, length.out = n_support)

        y_left <- -qnorm(0.975)/x_left
        y_right <- qnorm(0.975)/x_right

        funneldata <- data.frame(x = rep(rep(M, each = n_support*2), times = n) + rep(c(x_left, x_right), times = n),
                                 y = rep(c(y_left, y_right), times = n),
                                 .sample = rep(1:n, each = 2*n_support)
        )

        meandata <- data.frame (M = M, .sample = 1:n)

        # update x axis limit
        min_x <- min(c(min_x, min(funneldata$x)))
        max_x <- max(c(max_x, max(funneldata$x)))
        # update y axis limit
        #max_y <- min(c(max(funneldata$y), max_y))
      }

      if(egger) {
        warning("Note: egger = TRUE ignored: Egger's regression line can only be plotted for y_axis = se")
      }

    } else {
      stop("y_axis argument must be either se or precision")
    }
  }


  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  .sample <- NULL
  real_data <- NULL
  tnf_summary <- NULL
  intercept <- NULL
  slope <- NULL
  y0.01 <- NULL
  y0.05 <- NULL

  # Construct plots
  if(y_axis == "se") {
    p <- ggplot(data = plotdata, aes(x = es, y = se))

    if(show_solution) {
      p <-
        p + geom_rect(data = subset(plotdata, real_data == T), xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "lightgreen")
        # print(paste("True data in position", solution, sep = " "))
    }


    if(sig_contours) {
      p <-
        p +
        geom_polygon(data = sig_funneldata, aes(x = x.01, y = y), fill = col[9], alpha = 0.6) +
        geom_polygon(data = sig_funneldata, aes(x = x.05, y = y), fill = "white", alpha = 0.8) +
        geom_path(data = sig_funneldata, aes(x = x.05, y = y)) +
        geom_path(data = sig_funneldata, aes(x = x.01, y = y))
    }

    if(is.null(group)) {
      p <- p + geom_point(size = 2 , fill = "white", shape = 21, col = "black", alpha = 1)
    } else {
      p <- p + geom_point(aes(fill = group, shape = group), col = "black", size = 2, alpha = 1)
    }

    if(trim_and_fill) {
      if(dim(tnfdata)[1] > 0) {
        if(is.null(group)) {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = se), size = 2 , fill = "black", shape = 21, col = "black", alpha = 1)
        } else {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = se, shape = group),
                              size = 2 , fill = "black", col = "black", alpha = 1)
        }
        if(contours) {
          p <- p + geom_vline(data = tnfdata, aes(xintercept = tnf_summary), lty = "dashed")
        }
      }
    }

    if(egger) {
      p <- p + geom_abline(data = eggerdata, aes(intercept = intercept, slope = slope), lty = "twodash", lwd = 1, color = "firebrick")
    }

    if(contours) {
      p <-
        p +
        geom_path(data = funneldata, aes(x = x, y = y)) +
        geom_hline(yintercept = max_se) +
        geom_vline(data = meandata, aes(xintercept = M))
    }

    p <-
      p +
      scale_y_reverse(name = "Standard Error") +
      coord_cartesian(xlim = c(min_x - 0.2, max_x + 0.2),
                      ylim = c(max(se + 0.1), 0), expand = F) +
      scale_shape_manual(values = 21:25) +
      scale_x_continuous(name = "Effect") +
      scale_fill_brewer(palette = "Set1", type = "qual") +
      theme(legend.position = "none",
            axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
            axis.text = element_text(size = rel(0.8), colour = "black"),
            panel.border = element_rect(fill = NA, colour = "black"),
            panel.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    p <- p + facet_wrap( ~ .sample, ncol = 5, scales = "free")
  }
  if(y_axis == "precision") {
    p <- ggplot(data = plotdata, aes(x = es, y = 1/se))

    if(sig_contours) {
      p <-
        p +
        geom_polygon(data = sig_funneldata, aes(x = x, y = y0.01), fill = col[9], alpha = 0.6) +
        geom_polygon(data = sig_funneldata, aes(x = x, y = y0.05), fill = "white", alpha = 0.8) +
        geom_path(data = sig_funneldata, aes(x = x, y = y0.05)) +
        geom_path(data = sig_funneldata, aes(x = x, y = y0.01))
    }

    if(is.null(group)) {
      p <- p + geom_point(size = 2 , fill = "white", shape = 21, col = "black", alpha = 1)
    } else {
      p <- p + geom_point(aes(fill = group, shape = group), col = "black", size = 2, alpha = 1)
    }

    if(trim_and_fill) {
      if(dim(tnfdata)[1] > 0) {
        if(is.null(group)) {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = 1/se), size = 2 , fill = "black", shape = 21, col = "black", alpha = 1)
        } else {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = 1/se, shape = group),
                              size = 2 , fill = "black", col = "black", alpha = 1)
        }
        p <- p + geom_vline(data = tnfdata, aes(xintercept = tnf_summary), lty = "dashed")
      }
    }

    if(contours) {
      p <-
        p +
        geom_path(data = funneldata, aes(x = x, y = y)) +
        geom_vline(data = meandata, aes(xintercept = M))
    }


    p <-
      p +
      scale_y_continuous(name = "Precision (1/SE)")  +
      coord_cartesian(xlim = c(min_x - 0.2, max_x + 0.2),
                      ylim= c(1/max_se - 0.05, max_y + 0.5), expand = F) +
      scale_shape_manual(values = 21:25) +
      scale_x_continuous(name = "Effect") +
      scale_color_brewer(palette = "Set1", type = "qual") +
      scale_fill_brewer(palette = "Set1", type = "qual") +
      theme(legend.position = "none",
            axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
            axis.text = element_text(size = rel(0.8), colour = "black"),
            panel.border = element_rect(fill = NA, colour = "black"),
            panel.background = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    p <- p + facet_wrap( ~ .sample, ncol = 5, scales = "free")
  }
  p
}
