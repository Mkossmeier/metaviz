#' Visual funnel plot inference for meta-analysis
#'
#' Creates a lineup of funnel plots to conduct visual funnel plot inference. The funnel plot showing the actually observed, supplied data is
#' presented alongside null plots showing simulated data under the null hypothesis.
#'
#' Funnel plots are widely used in meta-analysis to detect small study effects and in particular publication bias.
#' However, interpretations of funnel plots often lead to false conclusions (e.g., Terrin, Schmid, and Lau, 2005). Visual inference
#' (Buja et al. 2009; Majumder, Hofmann, and Cook 2013) can help to improve the validity of conclusions based on the visual inspection of a
#' funnel plot, by saving investigators from interpreting funnel-plot patterns which might be perfectly plausible by chance.
#' Only if the real-data funnel plot is identifiable from null-plots, the null hypothesis is formally rejected and
#' conclusions based on the visual inspection of the real-data funnel plot might be warranted.
#'
#' Function \code{funnelinf} utilizes package \pkg{nullabor} for null plot simulation and \pkg{ggplot2} for
#' plotting the lineup. Several tailored features for visual inference with funnel plots are provided which currently include:
#' \enumerate{
#' \item options for null-plot simulation under both FEM and REM meta-analysis (see below).
#' \item subgroup analysis.
#' \item graphical options specific to the funnel plot (significance and confidence contours, and choice of the ordinate).
#' \item additional options to display various statistical information (Egger's regression line, and imputed studies by, as well as the adjusted summary effect from, the trim-and-fill method).
#' }
#'
#' Null plots are simulated assuming normally distributed effect sizes with expected value equal to the observed summary effect and variance
#' either equal to the observed study variances (\code{null_model = "FEM"}) or the sum of the observed study variances and the estimated
#' between study variance \eqn{\tau^2}{tau squared} (\code{null_model = "REM"}).
#'
#'@param x data.frame or matrix with the effect sizes of all studies (e.g.,
#'  correlations, log odds ratios, or Cohen \emph{d}) in the first column and their
#'  respective standard errors in the second column. Alternatively, x can be the
#'  output object of function \code{\link[metafor]{rma.uni}} from package
#'  \pkg{metafor}; then effect sizes and standard errors are extracted from \code{x}.
#'@param group factor indicating the subgroup of each study to show in the funnel plot. Has to be in the same order than \code{x}.
#'@param group_permut logical scalar indicating if subgroup membership should be permutated
#'  in the null plots. Ignored if no group is supplied.
#'@param n integer specifying the absolute number of plots in the lineup.
#'@param y_axis character string indicating which y axis should be used in the funnel plot. Available options are "se" (default) for
#'  standard error and "precision" for the reciprocal of the standard error.
#'@param null_model character string indicating which meta-analytic model should be used to simulate the effect sizes for the null plots.
#'  Available options are "FEM" for the fixed effect model and "REM" (default) for the random-effects model (using the
#'  DerSimonian-Laird method to estimate the between-study variance \eqn{\tau^2}{tau squared}).
#'@param contours logical scalar indicating if classic funnel plot confidence contours and the summary effect should be displayed (i.e., summary effect +/-
#'  qnorm(0.975) * SE).
#'@param sig_contours logical scalar. Should significance contours be drawn? Significance contours show which combination of
#' effect size and standard error lead to study p-values smaller than 0.05 or 0.01 (using a Wald test).
#'@param contours_col character string indicating the color palette used from package \pkg{RColorBrewer} for
#'  \code{sig_contours}. Can be any of "Blues", "Greys", "Oranges", "Greens", "Reds", and "Purples".
#'@param trim_and_fill logical scalar. Should studies imputed by the trim and fill method be displayed? Also shows the adjusted summary
#'  effect if \code{contours} is \code{TRUE} as well.
#'@param trim_and_fill_side character string indicating on which side of the funnel plot studies should be imputed by the trim and fill method (i.e., on which side are studies presumably missing due to publication bias).
#'  Must be either "right" or "left" (default).
#'@param egger logical scalar. Should Egger's regression line be drawn? Only available if \code{y_axis} is \code{"se"}.
#'@param show_solution logical scalar. Should the real-data plot be highlighted?
#'@param rorschach logical scalar. Should the lineup only consist of null plots?
#'@param text_size numeric value. Size of text in the lineup.
#'@param point_size numeric value. Size of the study points in the funnel plots.
#'@param xlab character string specifying the label of the x axis.
#'@param ylab character string specifying the label of the y axis.
#'@param x_trans_function function to transform the labels of the x axis. Common uses are to transform
#'  log-odds-ratios or log-risk-ratios with \code{exp} to their original scale (odds ratios and risk ratios), or Fisher's z values
#'  back to correlation coefficients using \code{tanh}.
#'@param x_breaks numeric vector of values for the breaks on the x-axis. When used in tandem with \code{x_trans_function}
#'  the supplied values should be not yet transformed.
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
#' \dontrun{
#' # Plotting a funnel plot lineup with the exrehab data to conduct visual funnel plot inference
#' funnelinf(x = exrehab[, c("logrr", "logrr_se")])
#'
#' # Plotting a funnel plot lineup with the mozart data to conduct visual funnel plot inference
#' # considering subgroups
#' funnelinf(x = mozart[, c("d", "se")],
#' group = mozart[, "rr_lab"],
#' group_permut = TRUE, null_model = "REM")
#'
#' # Plotting a funnel plot lineup with the brainvolume data to conduct visual funnel plot inference
#' # considering heterogeneity by using the fixed effect model for null plot simulation
#' funnelinf(x = brainvol[, c("z", "z_se")],
#' null_model = "FEM")
#' }
#' @export
funnelinf <- function(x, group = NULL, group_permut = FALSE, n = 20, null_model = "REM", y_axis = "se",
                      contours = TRUE, sig_contours = TRUE, contours_col = "Blues",
                      trim_and_fill = FALSE, trim_and_fill_side = "left",
                      egger = FALSE, show_solution = FALSE, rorschach = FALSE,
                      point_size = 1.5, text_size = 3, xlab = "Effect", ylab = NULL,
                      x_trans_function = NULL, x_breaks = NULL) {

  #'@import ggplot2
  #'@import dplyr

  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }
  # input is output of rma (metafor)
  if("rma" %in% class(x)) {
    # extract effect size and standard error
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))

    # If No group is supplied try to extract group from input object of class rma.uni (metafor)
    if(is.null(group) && ncol(x$X) > 1) {
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
        stop("Input argument has to be numeric; see help(funnelinf) for details.")
      }
      # check if there are any negative standard errors
      if(!all(x[, 2] > 0)) {
        stop("Non-positive standard errors supplied")
      }
      # extract effect sizes and standard errors
      es <- x[, 1]
      se <- x[, 2]

    } else {
      stop("Unknown input argument; see help(funnelinf) for details.")
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

  k <- length(es)
  summary_es <- metafor::rma.uni(yi = es, sei = se, method = "FE")$b[1] # for funnel plot contours

  if(null_model == "FEM") {
    summary_es_simul <- summary_es
    se_simul <- se
  } else {
    if(null_model == "REM") {
      summary_es_simul <- metafor::rma.uni(yi = es, sei = se, method = "DL")$b[1]
      se_simul <- sqrt(se^2 + metafor::rma.uni(yi = es, sei = se, method = "DL")$tau2[1])
    } else {
      stop("Supported null_model arguments are either FEM or REM")
    }
  }

  if(is.null(group)) {
    data <- data.frame(es, se)
  } else {
    data <- data.frame(es, se, group)
  }

  if(n < 1 || n > 100 || (as.integer(n) != n)) {
    stop("n must be an integer between 1 and 100.")
  }
  # Sample effect sizes from a normal distribution mu = summary effect and sd = se_simul
  if(rorschach != TRUE) {
    if(n > 1) {
      out_message <- NULL
      suppressMessages(
        x <- withCallingHandlers(
          nullabor::lineup(nullabor::null_dist("es", dist="normal", params = list(mean = summary_es_simul, sd = se_simul)),
                           true = data, n = n),
          message = function(w) {
            out_message <<- paste("To see the solution run nullabor::", w$message, sep = "")
          }
        )
      )
      message(out_message, appendLF = FALSE)
      solution <- attributes(x)$pos
    } else {
      if(n == 1) {
        x <- data.frame(data, .sample = 1)
        solution <- 1
      } else {
        stop("n has to be positive")
      }
    }
  } else {
    x <- nullabor::rorschach(nullabor::null_dist("es", dist = "normal", params = list(mean = summary_es_simul, sd = se_simul)),
                       true = data, n = n, p = 0)
    names(x)[1] <- ".sample"
    if(show_solution == TRUE) {
      warning("If rorschach = TRUE the lineup only consists of null plots. Argument show_solution = TRUE ignored.")
      show_solution <- FALSE
    }
  }

  if(!is.null(group) && group_permut == TRUE) {
    plotdata <- x %>%
      group_by(.sample) %>%
      mutate(group_permut = sample(group, size = k, replace = F)) %>%
      ungroup()
    if(rorschach != TRUE) {
      plotdata[plotdata$.sample == solution, "group_permut"] <- group
    }
    plotdata$group <- plotdata$group_permut
  } else {
    plotdata <- x
  }

  # Add column that shows real data
  if(rorschach != TRUE) {
  plotdata <- data.frame(plotdata, real_data = plotdata$.sample == solution)
  }

  # Set Color palette for contours
  if(!(contours_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
    warning("Supported arguments for contours_col are Blues, Greys, Oranges, Greens, Reds, and Purples. Blues is used.")
    contours_col <- "Blues"
  }
  col <- RColorBrewer::brewer.pal(n = 9, name = contours_col)

  # set intial min and max values for x axis limits
  min_x <- min(plotdata$es)
  max_x <- max(plotdata$es)

  if(trim_and_fill == TRUE) {
    trimnfill <- function(es, se, group = NULL, side = "left") {
      if(side == "right") {
        es <- -es
      }
      if(side != "right" && side != "left") {
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
      . <- NULL # avoid no visible binding for global variable warning (non standard evaluation)
      tnfdata <- plotdata %>%
        group_by(.sample) %>%
        do(trimnfill(.$es, .$se, side = side))
    } else {
      . <- NULL # avoid no visible binding for global variable warning (non standard evaluation)
      tnfdata <- plotdata %>%
        group_by(.sample) %>%
        do(trimnfill(.$es, .$se, .$group, side = side))
    }

    if(dim(tnfdata)[1] > 0) {
      # update limit values for x axis
      min_x <- min(c(min_x, min(tnfdata$es_fill)))
      max_x <- max(c(max_x, max(tnfdata$es_fill)))
    } else {
      trim_and_fill <- FALSE
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
      tnfdata$y <- tnfdata$se_fill
    }

    # determine significance contours
    if(sig_contours == TRUE) {
      x.05 <- c(-stats::qnorm(0.975) * max_se, 0, stats::qnorm(0.975) * max_se)
      x.01  <- c(-stats::qnorm(0.995) * max_se, 0, stats::qnorm(0.995) * max_se)
      y <- c(max_se, 0, max_se)
      sig_funneldata <- data.frame(x.05 = rep(x.05, times = n),
                                   x.01 = rep(x.01, times = n),
                                   y = rep(y, times = n),
                                   .sample = rep(1:n, each = 3))
      # update limit values for x axis
      min_x <- min(c(min_x, sig_funneldata$x.01[1]))
      max_x <- max(c(max_x, sig_funneldata$x.01[3]))
    }
    # determine classic funnel contours
    if(contours == TRUE) {
      M <- x %>%
        group_by(.sample) %>%
        summarise(M = sum((1/se^2)*es)/(sum(1/se^2))) %>%
        select(M)
      M <- unlist(M)
      x <- c(M - stats::qnorm(0.975) * max_se, M, M + stats::qnorm(0.975) * max_se)
      y <- rep(c(max_se, 0, max_se), each = n)
      funneldata <- data.frame(x = x, y = y, .sample = rep(1:n, times = 3))

      # update limit values for x axis
      min_x <- min(c(min_x, min(funneldata$x)))
      max_x <- max(c(max_x, max(funneldata$x)))

      meandata <- data.frame (M = M, .sample = 1:n)

    }
    # determine egger's regression line
    if(egger == TRUE) {
      plotdata <- data.frame(plotdata, "z" = (plotdata$es)/plotdata$se)
      plotdata <- data.frame(plotdata, "prec" = 1/plotdata$se)
      radial_coefs <- plotdata %>%
        group_by(.sample) %>%
        summarise(intercept = stats::coef(stats::lm(z ~ prec))[1],
                  slope = stats::coef(stats::lm(z ~ prec))[2])
      # note SE = d*(1/intercept) - slope/intercept, i.e. if egger's regression funnel plot has negative slope, this indicates PB
      eggerdata <- data.frame(intercept = radial_coefs$slope/radial_coefs$intercept,
                              slope = -1/radial_coefs$intercept, #minus slope because ggplot2 does not adjust slope of abline for reversed y axis
                              .sample = 1:n)
    }
  } else {
    if(y_axis == "precision") {
      plotdata$y <- 1/se

      # values for  y limit
      max_y <- max(1/se) + ifelse(length(se) > 1, diff(range(1/se))*0.05, 1/se*0.05)
      min_y <- min(1/se) - ifelse(length(se) > 1, diff(range(1/se))*0.05, 1/se*0.05)
      y_limit <- c(min_y, max_y)

      if(is.null(ylab)) {
        ylab <- "Precision (1/SE)"
      }

      if(trim_and_fill == TRUE) {
        tnfdata$y <- 1/tnfdata$se_fill
      }

      # determine significance contours
      if(sig_contours == TRUE) {
        n_support <- 200

        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec_0.05 <- stats::qnorm(0.975)*(1/prec)
        x_prec_0.01 <- stats::qnorm(0.995)*(1/prec)

        sig_funneldata <- data.frame(x0.05 = rep(c(-x_prec_0.05, rev(x_prec_0.05)), times = n),
                                     x0.01 =  rep(c(-x_prec_0.01, rev(x_prec_0.01)), times = n),
                                     y = rep(c(prec, rev(prec)), times = n),
                                     .sample = rep(1:n, each = 2*n_support)
        )
        # update x axis limit
        min_x <- min(c(min_x, min(sig_funneldata$x0.01)))
        max_x <- max(c(max_x, max(sig_funneldata$x0.01)))
      }
      # determine classic funnel contour
      if(contours == TRUE) {
        M <- x %>%
          group_by(.sample) %>%
          summarise(M = sum((1/se^2)*es)/(sum(1/se^2))) %>%
          select(M)
        M <- unlist(M)
        n_support <- 200
        prec <- seq(from = min_y, to = max_y, length.out = n_support)
        x_prec <- stats::qnorm(0.975)*1/prec
        funneldata <- data.frame(x = rep(rep(M, each = n_support*2), times = n) + rep(c(- x_prec, rev(x_prec)), times = n),
                                 y = rep(c(prec, rev(prec)), times = n),
                                 .sample = rep(1:n, each = 2*n_support)
        )

        meandata <- data.frame (M = M, .sample = 1:n)

        # update x axis limit
        min_x <- min(c(min_x, min(funneldata$x)))
        max_x <- max(c(max_x, max(funneldata$x)))

      }
      if(egger == TRUE) {
        warning("Note: egger = TRUE ignored: Egger's regression line can only be plotted for y_axis = se")
        egger <- FALSE
      }
    } else {
      stop("y_axis argument must be either se or precision")
    }
  }

  if(!is.null(x_trans_function) && !is.function(x_trans_function)) {
    warning("Argument x_trans_function must be a function; input ignored.")
    x_trans_function <- NULL
  }

  x_limit <- c(min_x - diff(c(min_x, max_x))*0.05, max_x + diff(c(min_x, max_x))*0.05)

  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  .sample <- NULL
  real_data <- NULL
  tnf_summary <- NULL
  intercept <- NULL
  slope <- NULL
  x0.01 <- NULL
  x0.05 <- NULL
  es_fill <- NULL
  group_fill <- NULL
  summary_es_init <- NULL

  # Construct plots
  p <- ggplot(data = plotdata, aes(x = es, y = y))
    if(show_solution == TRUE) {
      p <-
        p + geom_rect(data = subset(plotdata, real_data == T), xmin = -Inf, xmax = Inf,
                      ymin = -Inf, ymax = Inf, fill = "lightgreen")
    }
    if(sig_contours == TRUE && y_axis == "se") {
      p <- p +
        geom_polygon(data = sig_funneldata, aes(x = x.01, y = y), fill = col[9], alpha = 0.6) +
        geom_polygon(data = sig_funneldata, aes(x = x.05, y = y), fill = "white", alpha = 0.8) +
        geom_path(data = sig_funneldata, aes(x = x.05, y = y)) +
        geom_path(data = sig_funneldata, aes(x = x.01, y = y))
    } else {
      if(sig_contours == TRUE && y_axis == "precision") {
        p <- p +
          geom_polygon(data = sig_funneldata, aes(x = x0.01, y = y), fill = col[9], alpha = 0.6) +
          geom_polygon(data = sig_funneldata, aes(x = x0.05, y = y), fill = "white", alpha = 0.8) +
          geom_path(data = sig_funneldata, aes(x = x0.01, y = y)) +
          geom_path(data = sig_funneldata, aes(x = x0.05, y = y))
      }
    }
    if(contours == TRUE) {
      p <- p +
        geom_path(data = funneldata, aes(x = x, y = y)) +
        geom_vline(data = meandata, aes(xintercept = M))
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
        if(is.null(group)) {
          p <- p + geom_point(data = tnfdata, aes(x = es_fill, y = y), size = point_size , fill = "black", col = "black", alpha = 1)
        } else {
          p <- p + geom_point(data = tnfdata, aes(x = es_fill, y = y, shape = group_fill),
                              size = point_size , col = "black", fill = "black", alpha = 1)
        }
        if(contours == TRUE) {
          p <- p + geom_vline(data = tnfdata, aes(xintercept = summary_es_init), lty = "dashed")
        }
    }
    if(is.null(group)) {
      p <- p + geom_point(size = point_size , fill = "white", shape = 21, col = "black", alpha = 1)
    } else {
      p <- p + geom_point(aes(fill = group, shape = group), size = point_size, col = "black", alpha = 1)
    }
    if(egger == TRUE) {
      p <- p + geom_abline(data = eggerdata, aes(intercept = intercept, slope = slope), lty = "dashed", lwd = 0.5, color = "firebrick")
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
    p <-
      p +
      coord_cartesian(xlim = x_limit,
                      ylim = y_limit, expand = F) +
      scale_fill_brewer(palette = "Set1", type = "qual") +
      scale_shape_manual(values = 21:25) +
      theme_bw() +
      theme(legend.position = "none",
            text = element_text(size = 1/0.352777778*text_size),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
    p <- p + facet_wrap( ~ .sample, ncol = ceiling(sqrt(n)), scales = "free")
p
}


