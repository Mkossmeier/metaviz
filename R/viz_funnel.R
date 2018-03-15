#'Funnel plot variants for meta-analysis
#'
#'Creates a funnel plot. Many options regarding the appearance and
#'statistical information displayed are provided (e.g., significance contours, additional evidence contours, and
#'trim-and-fill analysis).
#'
#'The funnel plot is a widely used diagnostic plot in meta-analysis to assess small study effects
#'and in particular publication bias. The function \code{viz_funnel} is capable to create a large set of different funnel plot variants.
#'Options for several graphical augmentations (e.g., confidence, significance, and additional evidence contours; choice of the ordinate; study subgroups), and
#'different statistical information displayed are provided (Egger's regression line, and imputed studies by, as well as the adjusted summary effect from,
#'the trim-and-fill method).
#'
#'\bold{Contours}
#'
#'Three different contours are available in \code{viz_funnel}:
#'\enumerate{
#'\item \bold{confidence contours} (argument \code{contours}) show the region where one expects 95\% of all studies to fall (assuming the meta-analytic model
#'  applied is true and all estimates are identical to the parameters of interest).
#'  Confidence contours can help to assess the plausibility of observations given the meta-analytic model specified (fixed effect or random effects model).
#'\item \bold{significance contours} (argument \code{sig_contours}) show shaded regions of individual study significance at the 5\% and 1\% level
#'  (using the standard errors supplied and a Wald test). Significance contours were proposed to help distinguish publication bias from other sources of
#'  funnel plot asymmetry (Peters, Sutton, Jones, Abrams, & Rushton, 2008).
#'\item \bold{additional evidence contours: Significance of the summary effect} (argument \code{addev_contours}). These contours define regions
#'  where a new study has to fall such that the updated meta-analytic summary effect is significantly different from zero or not
#'  (using a two-sided test and an alpha level of 5\%). Additional evidence contours allow to assess the robustness of the meta-analysis with respect to
#'  the effect of potentially new published evidence on the significance of the meta-analytic summary effect (Langan, Higgins, Gregory, & Sutton, 2012).
#'}
#'
#'\bold{Measure on the y-axis}
#'
#'  Two different options for the y-axis choice are available. First, to plot the standard errors on a reversed axis
#'  (i.e., studies with small standard errors are at the top). Second, precision (i.e., 1 divided by the standard error) can be used.
#'  Standard errors on the y-axis should be preferred in most situations but precision might have advantages if one or few large studies
#'  (with high precision) should be compared to the results of smaller studies condensed at the bottom of the funnel plot (Sterne & Egger, 2001).
#'
#'\bold{Egger's regression line}
#'
#'  Egger's regression line (Egger, Smith, Schneider & Minder, 1997) can be displayed if the standard error is used on the y axis.
#'  Classic Egger's regression can be computed as the OLS estimator of regressing the standardized effect size (effect size divided by its standard error)
#'  on precision (1 divided by the standard error). Showing this line in the funnel plot can further help to visually assess funnel plot asymmetry.
#'
#'\bold{Trim and fill analysis}
#'
#'  Imputed studies by the trim-and fill method, as well as the adjusted summary effect (Duval & Tweedie, 2000) can be displayed.
#'  The trim-and fill algorithm basically estimates the number of (extreme) studies responsible for funnel plot asymmetry.
#'  It then trims this number of (extreme) studies and computes the adjusted summary effect only considering the remaining studies.
#'  Finally, it imputes studies - presumably missing due to publication bias - by mirroring the trimmed (extreme) studies (driving the funnel plot asymmetry)
#'  around the (adjusted) summary effect. The user has to specify on which side of the funnel plot the trim-and fill method should impute missing studies
#'  (i.e., the direction were studies are presumably missing due to publication bias). To estimate the number of (extreme) studies responsible for funnel plot
#'  asymmetry the L estimator defined in Duval and Tweedie (2000) is used.
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}; then effect sizes and standard errors are extracted from \code{x}.
#'@param group factor indicating the subgroup of each study to show in the funnel plot. Has to be in the same order than \code{x}.
#'@param method character string indicating the method used to compute the meta-analytic summary effect and, for a random effects model,
#'  the between-study variance \eqn{\tau^2}{tau squared}. Can be any method argument from \code{\link[metafor]{rma.uni}}
#'  (e.g., "FE" for the fixed effect model, or "DL" for the random effects model using the DerSimonian-Laird method to estimate \eqn{\tau^2}{tau squared}).
#'  Used for \code{contours}, and \code{addev_contours}.
#'@param y_axis character string indicating which y axis should be used in the funnel plot. Available options are "se" (default) for
#'  standard error and "precision" for the reciprocal of the standard error.
#'@param contours logical scalar indicating if classic funnel plot confidence contours and the summary effect
#'  should be displayed.
#'@param sig_contours logical scalar. Should significance contours be drawn (at the 0.05 or 0.01 level using a Wald test)?
#'@param addev_contours logical scalar. Should approximate additional evidence contours be drawn, showing the significance of the updated summary effect? See Details.
#'  Note: For \code{method} other than "FE" or "DL" runtime is increased significantly. Consider reducing \code{detail_level}.
#'@param contours_col character string indicating the color palette used from package \pkg{RColorBrewer} for
#'  \code{sig_contours}, and \code{addev_contours}. Can be any of "Blues", "Greys", "Oranges", "Greens", "Reds", and "Purples".
#'@param detail_level numeric scalar. Allows to increase or decrease the plotting detail of contours. Values can be chosen between 0.1 and 10. Default is 1.
#'@param trim_and_fill logical scalar. Should studies imputed by the trim and fill method be displayed? Also shows the adjusted summary
#'  effect if \code{contours} is \code{TRUE} as well.
#'@param trim_and_fill_side character string indicating on which side of the funnel plot studies should be imputed by the trim and fill method (i.e., on which side are studies presumably missing due to publication bias).
#'  Must be either "right" or "left" (default).
#'@param egger logical scalar. Should Egger's regression line be drawn? Only available if \code{y_axis} is \code{"se"}.
#'@param text_size numeric value. Size of text in the funnel plot. Default is 3.
#'@param point_size numeric value. Size of the study points in the funnel plot. Default is 2.
#'@param xlab character string specifying the label of the x axis.
#'@param ylab character string specifying the label of the y axis.
#'@param group_legend logical scalar. Should there be a legend shown at the bottom of the graph if \code{group} was supplied?
#'@param group_legend_title a character string specifying the title of the legend if \code{group} was supplied and \code{group_legend} is \code{TRUE}.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficients using \code{tanh}.
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
#'@references Duval, S., & Tweedie, R. (2000). Trim and fill: a simple funnel-plot-based method of testing and adjusting for publication bias
#'  in meta-analysis. \emph{Biometrics}, \emph{56}, 455-463.
#'@references Egger, M., Smith, G. D., Schneider, M., & Minder, C. (1997). Bias in meta-analysis detected by a simple, graphical test. \emph{Bmj}, \emph{315}, 629-634.
#'@references Langan, D., Higgins, J. P., Gregory, W., & Sutton, A. J. (2012). Graphical augmentations to the funnel plot assess the impact of
#'  additional evidence on a meta-analysis. \emph{Journal of clinical epidemiology}, \emph{65}, 511-519.
#'@references Peters, J. L., Sutton, A. J., Jones, D. R., Abrams, K. R., & Rushton, L. (2008). Contour-enhanced meta-analysis
#'  funnel plots help distinguish publication bias from other causes of asymmetry. \emph{Journal of clinical epidemiology}, \emph{61},
#'  991-996.
#'@references  Sterne, J. A., & Egger, M. (2001). Funnel plots for detecting bias in meta-analysis: guidelines on choice of axis.
#'  \emph{Journal of clinical epidemiology}, \emph{54}, 1046-1055
#'@return A funnel plot is created using ggplot2.
#'@author Michael Kossmeier* <michael.kossmeier@univie.ac.at>
#'@author Ulrich S. Tran* <ulrich.tran@univie.ac.at>
#'@author Martin Voracek* <martin.voracek@univie.ac.at>
#'@author *Department of Basic Psychological Research and Research Methods, School of Psychology, University of Vienna
#'@examples
#' library(metaviz)
#' # Create a funnel plot using confidence and significance contours
#' viz_funnel(x = mozart[, c("d", "se")])
#'
#' # Show a trim-and-fill analysis and Egger's regression line:
#' viz_funnel(x = mozart[, c("d", "se")], contours = TRUE,
#' trim_and_fill = TRUE, trim_and_fill_side = "left", egger = TRUE)
#'
#' # Plot log-odds-ratios on the original OR scale and show additional evidence contours:
#' viz_funnel(x = exrehab[, c("logor", "logor_se")], sig_contours = FALSE,
#' addev_contours = TRUE, contours_col = "Greys", xlab = "Odds Ratio",
#' x_trans_function = exp, x_breaks = log(c(0.125, 0.25, 0.5, 1, 2, 4, 8)))
#'
#' # Show study subgroups
#' viz_funnel(x = mozart[, c("d", "se")], group = mozart[, "unpublished"],
#' group_legend_title = "unpublished?")
#'@export
viz_funnel <- function(x, group = NULL, y_axis = "se", method = "FE",
                      contours = TRUE, sig_contours = TRUE, addev_contours = FALSE,
                      contours_col = "Blues", detail_level = 1,
                      egger = FALSE, trim_and_fill = FALSE, trim_and_fill_side = "left",
                      text_size = 3, point_size = 2,
                      xlab = "Effect", ylab = NULL,
                      group_legend = FALSE, group_legend_title = "",
                      x_trans_function = NULL, x_breaks = NULL) {
  #'@import ggplot2
  #'@import dplyr

  # input is output of rma (metafor)
  if("rma" %in% class(x)) {
    # extract effect size and standard error
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    # method <- x$method
    if(method != x$method) {
      message("Note: method argument used differs from input object of class rma.uni (metafor)")
    }
    # If No group is supplied try to extract group from input object of class rma.uni (metafor)
    if(is.null(group) & ncol(x$X) > 1) {
      #check if only categorical moderators were used
      if(!all(x$X == 1 || x$X == 0) || any(apply(as.matrix(x$X[, -1]), 1, sum) > 1))  {
        stop("Can not deal with metafor output object with continuous and/or more than one categorical moderator variable(s).")
      }
      # extract group vector from the design matrix of the metafor object
      no.levels <- ncol(x$X) - 1
      group <- factor(apply(as.matrix(x$X[, -1]) * rep(1:no.levels, each = length(es)), 1, sum))
    }
  } else {
    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        if(!is.null(group)) {
          group <- group[stats::complete.cases(x)]
        }
        x <- x[stats::complete.cases(x), ]
      }
      # check if input is numeric
      if(!is.numeric(x[, 1]) || !is.numeric(x[, 2])) {
        stop("Input argument has to be numeric; see help(viz_funnel) for details.")
      }
      # check if there are any negative standard errors
      if(!all(x[, 2] > 0)) {
        stop("Non-positive standard errors supplied")
      }
      # extract effect sizes and standard errors
      es <- x[, 1]
      se <- x[, 2]

    } else {
      stop("Unknown input argument; see help(viz_funnel) for details.")
    }
  }

  # check if group is a factor
  if(!is.null(group) && !is.factor(group)) {
    group <- as.factor(group)
  }
  # check if group vector has the right length
  if(!is.null(group) && (length(group) != length(es)))
  {
    warning("length of supplied group vector does not correspond to the number of studies; group argument is ignored")
    group <- NULL
  }

  # Compute meta-analytic summary using rma.uni (metafor) and supplied method
  k <- length(es)
  summary_es <- metafor::rma.uni(yi = es, sei = se, method = method)$b[[1]]
  summary_se <- sqrt(metafor::rma.uni(yi = es, sei = se, method = method)$vb[[1]])
  summary_tau2 <- metafor::rma.uni(yi = es, sei = se, method = method)$tau2

  # main data for plotting
  if(is.null(group)) {
    plotdata <- data.frame(es, se)
  } else {
    plotdata <- data.frame(es, se, group)
  }

  # Set Color palette for contours
  if(!(contours_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
    warning("Supported arguments for contours_col are Blues, Greys, Oranges, Greens, Reds, and Purples. Blues is used.")
    contours_col <- "Blues"
  }
  # color palette for significance contours
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

  # Compute data for trim and fill analysis
  if(trim_and_fill == TRUE) {
    trimnfill <- function(es, se, group = NULL, side = "left") {
      if(side == "right") {
        es <- -es
      }
      if(side != "right" && side != "left") {
        stop("trim_and_fill_side argument must be either left or right")
      }
      mean_func <- function(es, se) {
        metafor::rma.uni(yi = es, sei=se, method = method)$b[1]
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
      while(eps > 0.01 || iter < 20) {
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
          summary_es_init <- - summary_es_init
        } else {
          es_fill <- summary_es_new + (summary_es_new - es_ord[1:k0])
        }

        se_fill <- se_ord[1:k0]

        if(is.null(group)) {
          data.frame(es_fill, se_fill, summary_es_init)
        } else {
          data.frame(es_fill, se_fill, group_fill, summary_es_init)
        }
      } else {
        if(is.null(group)) {
          data.frame(es_fill = NULL, se_fill = NULL, summary_es_init = NULL)
        } else {
          data.frame(es_fill = NULL, se_fill = NULL, group_fill = NULL, summary_es_init = NULL)
        }
      }
    }

    side <- trim_and_fill_side

    if(is.null(group)) {
      tnfdata <- trimnfill(es, se, side = side)
    } else {
      tnfdata <- trimnfill(es, se, group, side = side)
    }

    if(nrow(tnfdata) > 0) {
      if(is.null(group)) {
        names(tnfdata) <- c("es", "se", "tnf_summary")
      } else {
        names(tnfdata) <- c("es", "se", "group", "tnf_summary")
      }

      # update limit values for x axis
      min_x <- min(c(min_x, min(tnfdata$es)))
      max_x <- max(c(max_x, max(tnfdata$es)))
    } else {
      trim_and_fill <- FALSE
    }
  }


  # helper function to manually estimate DL estimator for additional evidence contours
  if(method == "DL" && addev_contours == TRUE) {
    rem_dl <- function(es, se) {
      summary_es_FEM <- sum((1/se^2)*es)/sum(1/se^2)
      n <- length(es)
      if(n == 1) {
        # Between study variance cannot be estimated and is set to zero")
        t2 <- 0
      } else {
        Q <- sum((1 / se^2) * (es - summary_es_FEM)^2)
        t2 <- max(c(0, (Q - (n - 1)) / (sum(1 / se^2) - sum((1 / se^2)^2) / sum(1/se^2))))
      }
      w <- 1/(se^2 + t2)
      c(sum(w*es)/sum(w), sqrt(1/sum(w)))
    }
  }

  # standard error on the y axis
  if(y_axis =="se") {
    plotdata$y <- se
    max_se <- max(se) + ifelse(length(se) > 1, diff(range(se))*0.1, max(se)*0.1)
    y_limit <- c(0, max_se)

    if(is.null(ylab)) {
      ylab <- "Standard Error"
    }

    if(trim_and_fill == TRUE) {
      tnfdata$y <- tnfdata$se
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

    # determine egger's regression line
    if(egger == TRUE) {
      plotdata <- data.frame(plotdata, "z" = (plotdata$es)/plotdata$se)
      plotdata <- data.frame(plotdata, "prec" = 1/plotdata$se)
      radial_intercept <- stats::coef(stats::lm(z ~ prec, data = plotdata))[1]
      radial_slope <- stats::coef(stats::lm(z ~ prec, data = plotdata))[2]
      # note SE = d*(1/intercept) - slope/intercept, i.e. if egger's regression funnel plot has negative slope, this indicates PB
      eggerdata <- data.frame(intercept = radial_slope/radial_intercept,
                              slope = -1/radial_intercept #minus slope because ggplot2 does not adjust slope of abline for reversed y axis
                              )
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

      if(trim_and_fill == TRUE) {
        tnfdata$y <- 1/tnfdata$se
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

      if(egger == TRUE) {
        warning("Note: egger = TRUE ignored: Egger's regression line can only be plotted for y_axis = se")
      }
      y_limit <- c(min_y, max_y)
    } else {
      stop("y_axis argument must be either se or precision")
    }
    }

  x_limit <- c(min_x - diff(c(min_x, max_x))*0.05, max_x + diff(c(min_x, max_x))*0.05)

  # Compute additional evidence contours
  if(addev_contours == TRUE) {
    if(y_axis == "se") {
      # set search grid für y_axis == "se"
      y_range <- c(0.001, max_se + diff(range(y_limit))*0.2)
      x_range <- c(min_x - diff(range(x_limit))*0.2, max_x + diff(range(x_limit))*0.2)

      step <- abs(summary_es - x_range[1])/ (150 * detail_level - 1)
      x_add <- c(seq(from = x_range[1],  to = summary_es, length.out = 150 * detail_level), seq(from = summary_es + step, to = x_range[2], by = step))
      y_add <- seq(from = y_range[1], to = y_range[2], length.out = length(x_add))

    } else {
      # set search grid für y_axis == "precision"
      y_range <- c(max_y + diff(range(x_limit))*0.2, min_y - diff(range(x_limit))*0.2)
      x_range <- c(min_x - diff(range(x_limit))*0.2, max_x + diff(range(x_limit))*0.2)

      step <- abs(summary_es - x_range[1])/ (150 * detail_level - 1)
      x_add <- c(seq(from = x_range[1],  to = summary_es, length.out = 150 * detail_level), seq(from = summary_es + step, to = x_range[2], by = step))
      y_add <- 1/seq(from = y_range[1], to = y_range[2], length.out = length(x_add))
    }
    # construct search grid
    study_grid <- expand.grid(x_add, y_add)
    names(study_grid) <- c("x_add", "y_add")


    # Determine the summary effect and significance for the search grid of new studies
    addev_data <- apply(study_grid, 1,
                        function(x) {
                          if(method == "FE") {
                            M_new <- sum((1/c(se, x[2])^2)*c(es, x[1]))/sum(1/c(se, x[2])^2)
                            Mse_new <- sqrt(1/sum(1/c(se, x[2])^2))
                            p.val <- stats::pnorm(M_new/Mse_new)
                            c(M_new, p.val)
                          } else {
                            if(method == "DL") {
                              res_dl <- rem_dl(es = c(es, x[1]), se = c(se, x[2]))
                              M_new <- res_dl[1]
                              p.val <- stats::pnorm(res_dl[1]/res_dl[2])
                              c(M_new, p.val)
                            } else {
                              mod <- metafor::rma.uni(yi = c(es, x[1]), sei = c(se, x[2]), method = method, control = list(stepadj = 0.5, maxiter = 1000))
                              p.val <- stats::pnorm(mod$z)
                              M_new <- mod$b[[1]]
                              c(M_new, p.val)
                            }
                          }
                        }
    )
    addev_data <- t(addev_data)
    addev_data <- data.frame(study_grid,
                               M = addev_data[, 1],
                               sig_group = factor(ifelse(addev_data[, 2] < 0.025, "sig.neg. ", ifelse(addev_data[, 2] > 0.975, "sig.pos. ", "not sig. ")),
                                                  levels = c("sig.neg. ", "not sig. ", "sig.pos. ")))

    addev_data <- addev_data[order(addev_data$x_add, decreasing = F), ]
    if(y_axis == "precision") {
      addev_data$y_add <- 1/addev_data$y_add
    }
  }

  if(!is.null(x_trans_function) && !is.function(x_trans_function)) {
    warning("Argument x_trans_function must be a function; input ignored.")
    x_trans_function <- NULL
  }

  # Construct plot
  y <- NULL
  sig_group <- NULL
  x.01 <- NULL
  x.05 <- NULL
  tnf_summary <- NULL
  intercept <- NULL
  slope <- NULL

  p <- ggplot(data = plotdata, aes(x = es, y = y))
    if(addev_contours == TRUE) {
      p <- p +
        geom_raster(data = addev_data, aes(x = x_add, y = y_add, fill = sig_group), alpha = 0.4) +
        scale_fill_manual(name = "", values = c(col[9], col[1], col[4]), drop = FALSE)
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
      p + scale_y_reverse(name = ylab)
    } else {
      if(y_axis == "precision") {
      p <-
        p + scale_y_continuous(name = ylab)
      }
    }
    if(trim_and_fill == TRUE) {
      if(dim(tnfdata)[1] > 0) {
        if(is.null(group)) {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = y), size = point_size , col = "black", alpha = 1)
        } else {
          p <- p + geom_point(data = tnfdata, aes(x = es, y = y, shape = group),
                              size = point_size , col = "black", alpha = 1)
        }
        if(contours == TRUE) {
          p <- p + geom_vline(data = tnfdata, aes(xintercept = tnf_summary), lty = "dashed")
        }
      }
    }
    if(is.null(group)) {
      p <- p + geom_point(size = point_size , fill = "white", shape = 21, col = "black", alpha = 1)
    } else {
      p <- p + geom_point(aes(col = group, shape = group), size = point_size, alpha = 1)
    }
    if(egger == TRUE && y_axis == "se") {
      p <- p + geom_abline(data = eggerdata, aes(intercept = intercept, slope = slope), lty = "dashed", lwd = 1, color = "firebrick")
    }
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
                      ylim = y_limit, expand = F) +
      scale_shape_manual(values = 15:19, name = group_legend_title) +
      scale_color_brewer(name = group_legend_title, palette = "Set1", type = "qual")
    if(group_legend == FALSE) {
      p <- p +
        guides(color = "none", shape = "none")
    }
    # Add black boxes around legend keys if addev_contour_sig is shown
    if(addev_contours == TRUE) {
      legend.key <-  element_rect(color = "black")
    } else {
      legend.key <-  element_rect(color = "white")
    }
    p <- p +
      theme_bw() +
      theme(text = element_text(size = 1/0.352777778*text_size),
            legend.position = "bottom",
            legend.key = legend.key,
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  p
}

