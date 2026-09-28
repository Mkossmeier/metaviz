#'@title Lorenz curves for the comparison of study weight concentration of fixed-effect and random-effects models

#'@description Creates two Lorenz curves within the same coordinate system that correspond to the study weight concentration
#'in a fixed-effect and a random-effects model, respectively, based on the same data.

#'@param x metafor rma.uni object conducted with method “FE” or “REML” (the chosen method of the input model makes no difference in the resulting plot)
#'@param type determines which Lorenz curves are shown. The default “FEM_REM” shows Lorenz curves for both the fixed-effect and random-effects models,
#'“FEM” shows only the fixed-effect Lorenz curve and “REM” shows only the random-effects Lorenz curve.
#'@param col boolean argument that determines whether Lorenz curves are colored by model or not.
#'@param tables boolean argument that determines whether a large table with lots of statistical information is plotted below the coordinate system (“TRUE”)
#'or whether a small table with information only on the Gini indices is shown within the coordinate system (“FALSE”)
#' @param seed numeric argument that is used to set a random seed in order to provide reproducible bootstrapped Gini index
#' confidence intervals. If seed == NULL the Gini index confidence intervals will vary slightly between instances of plotting.

#'@details The function wineq_gini creates a plot containing two Lorenz curves (Lorenz, 1905) which represent the concentration of weights among the studies
#'within a fixed-effect model and random-effects model, respectively. An adjoined table provides descriptive statistics about the study weights, as well as
#'Gini indices which correspond to the Lorenz curves (Gini, 1912; Tran et al, 2021). The Gini quotient quantifies the discrepancy between the two Lorenz curves,
#'is directly related to heterogeneity and serves as an effect size for cross-metaanalytic comparisons.
#'@return Two Lorenz curves are plotted in the same coordinate system accompanied by a table of statistical information.
#'@author Verena Pilar <verena.pilar@outlook.com>
#'@references
#'Gini, C. (1912). Variabilità e mutabilità (Variability and Mutability). C. Cuppini, Bologna, 156.
#'
#'Lorenz,M.O. (1905). Methods of measuring the concentration of wealth. \emph{Pub. Am. Stat. Assoc.} 9, 209–219. https://doi.org/10.2307/2276207
#'
#'Tran, U. S., Lallai, T., Gyimesi, M., Baliko, J., Ramazanova, D., & Voracek, M. (2021). Harnessing the fifth element of distributional statistics for psychological science:
#'A practical primer and shiny app for measures of statistical inequality and concentration.\emph{ Frontiers in Psychology}, 12, 716164.
#'@examples
#'library(metafor)
#'
#' # Calculating a random-effects model based on the mozart data
#' mozart_r <- rma(yi = d,
#'                 sei = se,
#'                 data = mozart,
#'                 method = "REML")
#'
#' # using a rma.uni model as input
#' wineq_gini(x = mozart_r, seed = 123)

#'@export

wineq_gini <- function(x, type = "FEM_REM", col = FALSE, tables = TRUE, seed = NULL) {

  #'@import ggplot2
  #'@import dplyr
  #'@import metafor
  #'@import metaviz
  #'@importFrom magrittr %>%
  NULL



  method <- type


  if(missing(x)) {
    stop("argument x is missing, with no default.")
  }

  if("rma" %in% class(x)) {
    es <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    n <- length(es)



  } else {

    # input is matrix or data.frame with effect sizes and standard errors in the first two columns
    if((is.data.frame(x) || is.matrix(x)) && ncol(x) >= 2) { # check if a data.frame or matrix with at least two columns is supplied
      # check if there are missing values
      if(sum(is.na(x[, 1])) != 0 || sum(is.na(x[, 2])) != 0) {
        warning("The effect sizes or standard errors contain missing values, only complete cases are used.")
        #study_labels <- study_labels[stats::complete.cases(x[, c(1, 2)])]

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


  # creating df to use as data for the models
  yi <- NULL # to avoid no visible binding for global variable note

  df <- data.frame (
    yi = es,
    se = se
  )
  # calculating models

  model_fem <- rma(yi = yi,
                   sei = se,
                   data = df,
                   method = "FE")

  model_rem <- rma(yi = yi,
                   sei = se,
                   data = df,
                   method = "REML")



  ### different data processing approach

  weights_fem <- NULL # to avoid no visible binding for global variable note
  weights_rem <- NULL
  share_fem <- NULL
  share_rem <- NULL
  comp_share <- NULL
  share_fem_c <- NULL
  share_rem_c <- NULL


  df_fem <- data.frame(weights_fem = stats::weights(model_fem),
                       ID = seq(from = 1, to = nrow(model_fem$data)))
  df_fem <- df_fem %>%
    mutate(share_fem = weights_fem/sum(weights_fem)*100) %>%
    mutate(across(everything(), ~sort(.))) %>%
    mutate(share_fem_c = cumsum(share_fem)/max(cumsum(share_fem))*100)

  # df for REM
  df_rem <- data.frame(weights_rem = stats::weights(model_rem),
                       ID = seq(from = 1, to = nrow(model_rem$data)))
  df_rem <- df_rem %>%
    mutate(share_rem = weights_rem/sum(weights_rem)*100) %>%
    mutate(across(everything(), ~sort(.))) %>%
    mutate(share_rem_c = cumsum(share_rem)/max(cumsum(share_rem))*100)

  gini_data <- merge(df_fem, df_rem, by = "ID", all = TRUE)
  gini_data <- gini_data %>%
    mutate(comp_share = seq(from = 1, to = nrow(gini_data))/nrow(gini_data)*100)


  if (!is.null(seed)) {
    if (!(is.numeric(seed) &&
          length(seed) == 1 &&
          is.finite(seed))) {
      warning("'seed' must be a single finite numeric value. The argument will be ignored. Bootstrapped Gini index confidence intervals will vary slightly between instances of plotting without a random seed.")
      seed <- NULL
    }
    set.seed(seed) # for reproducible bootstrap CIs
  } else {
    message("Note: Bootstrapped Gini index confidence intervals will vary slightly between instances of plotting. For reproducible bootstrap confidence intervals provide a numeric input for the 'seed' argument.")
  }


  # FEM
  gini_fem_all <- DescTools::Gini(gini_data$weights_fem, conf.level=.95, unbiased = FALSE)
  gini_fem <- gini_fem_all[[1]]
  gini_fem_lb <- gini_fem_all[[2]]
  gini_fem_ub <- gini_fem_all[[3]]
  gini_fem_c_all <- DescTools::Gini(gini_data$weights_fem, conf.level=.95, unbiased = TRUE)
  gini_fem_c <- gini_fem_c_all[[1]]
  gini_fem_c_lb <- gini_fem_c_all[[2]]
  gini_fem_c_ub <- gini_fem_c_all[[3]]

  # REM
  gini_rem_all <- DescTools::Gini(gini_data$weights_rem, conf.level=.95, unbiased = FALSE)
  gini_rem <- gini_rem_all[[1]]
  gini_rem_lb <- gini_rem_all[[2]]
  gini_rem_ub <- gini_rem_all[[3]]
  gini_rem_c_all <- DescTools::Gini(gini_data$weights_rem, conf.level=.95, unbiased = TRUE)
  gini_rem_c <- gini_rem_c_all[[1]]
  gini_rem_c_lb <- gini_rem_c_all[[2]]
  gini_rem_c_ub <- gini_rem_c_all[[3]]

  # DIFF
  gini_dif <- gini_fem - gini_rem
  gini_dif_c <- gini_fem_c - gini_rem_c

  # QUOTIENT
  gini_quot <- gini_rem/gini_fem
  gini_quot_c <- gini_rem_c/gini_fem_c



  gini_data <- rbind(rep(0,length(gini_data)), gini_data)



  if (tables == TRUE) {



    # creating a df for the gini indices
    gini_ind_df <- data.frame(
      Model = c("FEM", "REM", "FEM - REM", "REM/FEM"),
      Gini = c(round(gini_fem,3), round(gini_rem,3), round((gini_fem - gini_rem),3), round(gini_quot,3)),
      `CI` = c(sprintf("[%.3f; %.3f]", gini_fem_lb, gini_fem_ub),
               sprintf("[%.3f; %.3f]", gini_rem_lb, gini_rem_ub),
               "", ""),
      `Gini_corr` = c(round(gini_fem_c,3), round(gini_rem_c,3), round((gini_fem_c - gini_rem_c),3), round(gini_quot_c,3)),
      `CI_corr` = c(sprintf("[%.3f; %.3f]", gini_fem_c_lb, gini_fem_c_ub),
                    sprintf("[%.3f; %.3f]", gini_rem_c_lb, gini_rem_c_ub),
                    "", ""),
      check.names = FALSE     ### to prevent spaces from being turned into dots
    )

    # gathering descriptive data
    FEM <- NULL # to avoid no visible binding for global variable note
    REM <- NULL
    Model <- NULL
    weight <- NULL

    df_weights <- data.frame(
      FEM = stats::weights(model_fem),
      REM = stats::weights(model_rem)
    )

    df_weights_long <- df_weights %>%
      tidyr::pivot_longer(
        cols = c(FEM, REM),
        names_to = "Model",
        values_to = "weight"
      )

    df_weights_sum <- df_weights_long %>%
      group_by(Model) %>%
      reframe(
        Mean = mean(weight),
        Median = stats::median(weight),
        Sd = stats::sd(weight),
        Min = min(weight),
        Max = max(weight),
        IQR = stats::IQR(weight),          # interquartile range
        Skewness = moments::skewness(weight),
        CV = DescTools::CoefVar(weight)     # coefficient of variation
      )





    # table plotting
    table_gini <- gridExtra::tableGrob(
      gini_ind_df %>%
        mutate(across(where(is.numeric), ~ format(round(.x, 3), nsmall = 3))),
      theme = gridExtra::ttheme_minimal(
        core = list(fg_params=list(cex = 0.8)),
        colhead = list(fg_params=list(cex = 0.8))
      ),
      rows = NULL
    )


    table_weights <- gridExtra::tableGrob(
      df_weights_sum %>%
        mutate(across(where(is.numeric), ~ format(round(.x, 3), nsmall = 3))),
      theme = gridExtra::ttheme_minimal(
        core = list(fg_params=list(cex = 0.8)),
        colhead = list(fg_params=list(cex = 0.8))#,
      ),
      rows = NULL
    )

  } else {

    # small version for inside the plot
    gini_mini_df <- data.frame(
      Model = c("FEM", "REM", "FEM - REM", "REM/FEM"),
      `Gini_corr` = c(round(gini_fem_c,3), round(gini_rem_c,3), round((gini_fem_c - gini_rem_c),3), round(gini_quot_c,3)),
      `CI_corr` = c(sprintf("[%.3f; %.3f]", gini_fem_c_lb, gini_fem_c_ub),
                    sprintf("[%.3f; %.3f]", gini_rem_c_lb, gini_rem_c_ub),
                    "", ""),
      check.names = FALSE
    )

    table_mini <- gridExtra::tableGrob(
      gini_mini_df %>%
        mutate(across(where(is.numeric), ~ format(round(.x, 3), nsmall = 3))),
      theme = gridExtra::ttheme_minimal(
        core = list(bg_params = list(fill = "white", col = NA),
                    fg_params = list(cex = 0.8)
        ),
        colhead = list(bg_params = list(fill = "white", col = NA),
                       fg_params=list(cex = 0.8)
        )
      ),
      rows = NULL
    )

    table_mini <- gtable::gtable_add_grob(table_mini,
                                          grobs = grid::rectGrob(gp = grid::gpar(fill = NA, lwd = 2)),
                                          t = 1, b = nrow(table_mini), l = 1, r = ncol(table_mini))
  }

  plot_title <- "Study weight inequality Lorenz curves"

  # creating plot
  p <- ggplot2::ggplot(gini_data, aes(x=comp_share))

  if (method == "FEM_REM") {
    p <- p + ggpattern::geom_ribbon_pattern(aes(ymin = share_fem_c, ymax = share_rem_c),
                                            fill = "grey90",
                                            alpha = 0.5,
                                            pattern="stripe", pattern_angle=90, pattern_size=.4,
                                            pattern_fill="black", pattern_density =.2, pattern_alpha = .1,
                                            pattern_spacing = 0.015, pattern_color = NA)

  }

  if (method %in% c("FEM_REM", "FEM")) {
    p <- p + geom_ribbon(aes(ymin = share_fem_c, ymax = comp_share),
                         fill = "grey90", alpha = 0.4) +
      geom_line(aes(y=share_fem_c), lwd=.5, alpha=.5, color = "black")+
      geom_point(aes(y = share_fem_c, color="FEM", shape = "FEM"))

  } else if (method == "REM") {
    p <- p + geom_ribbon(aes(ymin = share_rem_c, ymax = comp_share),
                         fill = "grey90", alpha = 0.4)
  }


  if (method %in% c("FEM_REM", "REM")) {
    p <- p + geom_line(aes(y=share_rem_c), lwd=.5, alpha=.5, color = "black")+
      geom_point(aes(y = share_rem_c, color="REM", shape = "REM"))
  }


  p <- p + annotate("segment", x = 0, y = 0, xend = 100, yend = 100, lwd = 0.9, alpha = 0.9)

  if (tables == FALSE) {
    p <- p + annotation_custom(grob = table_mini,
                               xmin = 0, xmax = 50, ymin = 80, ymax = 100)
  }

  p <- p + labs(x = "Cumulative component (study) share in %",
                y = "Cumulative unit (weight) share in %")


  if (method == "FEM_REM") {
    p <- p +
      scale_shape_manual(values = c("FEM" = 16, "REM" = 15)) +
      scale_color_manual(values = c("FEM" = ifelse(col == FALSE, "black", "red3"),
                                    "REM" = ifelse(col == FALSE, "black", "steelblue3")))
  } else if (method == "FEM") {
    p <- p +
      scale_shape_manual(values = c("FEM" = 16)) +
      scale_color_manual(values = c("FEM" = ifelse(col == FALSE, "black", "red3")))
  } else if (method == "REM") {
    p <- p +
      scale_shape_manual(values = c("REM" = 15)) +
      scale_color_manual(values = c("REM" = ifelse(col == FALSE, "black", "steelblue3")))
  }




  p <- p + guides(color = guide_legend(title = NULL),
                  shape = guide_legend(title = NULL))

  p <- p + theme_minimal()+
    theme(legend.position = "inside",
          legend.position.inside = c(.85,0.15),
          legend.margin = margin(2,10,2,1, "pt"),
          legend.background = element_rect(fill = "white", linewidth =.4),
          legend.title = element_text(face = "bold", margin = margin(l=3, t=1, unit = "pt")),
          panel.grid.major = element_line(linewidth = .6),
          plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    coord_fixed(ratio = 1)




  if (!(method %in% c("FEM_REM", "FEM", "REM"))) {
    stop("The chosen plot model is not available.")
  }



  if (tables == TRUE) {
    left_vp <- grid::viewport(x = 0, just = "left")
    gridExtra::grid.arrange(p, table_gini, table_weights, nrow = 3,heights = c(10, 3, 2), vp = left_vp)
  } else {
    p
  }

}









