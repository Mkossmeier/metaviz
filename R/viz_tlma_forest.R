
#'@title Forest plot for the visualization of all three levels of a three-level meta-analysis
#'@description Creates a thick or classic forest plot which shows the effects contained in each study,
#'the overall effect of each study and the summary effect of a three-level meta-analysis.
#'
#'@param x metafor rma.mv object
#'@param variant “classic” (default) or “thick” to create a classic or thick TLMA forest plot variant
#'@param annotate_CI adds a right-hand side table to the plot containing the confidence intervals and number of effects of each study
#'@param study_table custom table on the left-hand side of the plot that contains study information. Takes a dataframe as input
#'which has to be of a length equal to the number of studies.
#'@param summary_table custom table on the left-hand side of the plot that contains information about the summary effect.
#'Takes a dataframe as input which contains one row.
#'@param table_headers headers for each column of the left-hand side table. Takes a character vector as input
#'@param median_precision adds grey errorbars that represent the median precision of an effect contained in the respective study.
#'The errorbar’s thickness denotes the number of effects contained in the study.
#'@param median_precision_thick determines whether the thickness of the median precision errorbars represents the number of
#'effect sizes contained within the respective study
#'@param ordered orders the plot by effect size when TRUE
#'@param clouds “TRUE”: shows the effects contained in each study as a cloud around the study effect. “FALSE”: only study effects are shown.
#'@param spread determines how far the single effects spread around the study effect
#'@param col colors single effects by study
#'@param labels y-axis tick labels. Takes a vector the length of the dataset as input.
#'@param xlab x-axis label
#'@param ylab y-axis label
#'@param title plot title
#'@param confidence_level_ci numeric confidence level for the confidence intervals of the study effects.
#'This argument is also used in the calculation of the median precision of an effect included in a study for additional grey error bars.
#'@param prediction_level_pi numeric confidence level for the prediction interval of the summary effect
#'@param show_nr_ES adds columns to the right-hand side table which shows the number of effects contained in the respective study
#'@param x_limit determines the limits of the x-axis. Input is a numeric vector of length 2 (min, max).
#'@param table_layout numeric layout matrix to customize the arrangement of the plot and tables
#'@param line if “TRUE” it creates a line connecting the ordered study effects
#'@param linewidth determines the width of the line connecting the ordered study effects
#'@param text_size determines text size within the plot
#'@param tick_col determines color of study effect markers

#'@details The function viz_tlma_forest creates a forest plot which visualizes all three levels of a three-level meta-analysis.
#'The study effects are most prominently featured and sized according to their weight in the meta-analysis. Each overall effect
#'of a study that contains at least two effects is surrounded by a cloud of the effects that are contained within it. The plot
#'is completed by the overall result with its prediction interval shown at the bottom. It is available in the classic and the thick
#'(Schild & Voracek, 2015) forest plot variant.
#'
#'Note: This function was developed on the basis of the viz_forest code which was created by Michael Kossmeier.
#'Note: This function adapted parts of the forest_plot_3 function code by Fernández-Castilla et al (2020) to generate study-level
#'information of the three-level meta-analysis.
#'@return A forest plot containing point estimates for single effects, study effects and the overall result of a three-level
#'meta-analysis is created using ggplot2.
#'@author Verena Pilar <verena.pilar@outlook.com>
#'@references
#'Fernández-Castilla, B., Declercq, L., Jamshidi, L., Beretvas, N., Onghena, P., & Van den Noortgate, W. (2020). Visual representations
#'of meta-analyses of multiple outcomes: extensions to forest plots, funnel plots, and caterpillar plots. \emph{Methodology}, 16(\emph{4}), 299-315.
#'https://doi.org/10.1002/jrsm.1424
#'
#'Schild, A. H. E., & Voracek, M. (2015). Finding your way out of the forest without a trail of bread crumbs: Development and evaluation
#'of two novel displays of forest plots. \emph{Research Synthesis Methods}, 6(\emph{1}), 74–86. https://doi.org/10.1002/jrsm.1125
#'@examples
#' if (requireNamespace("psymetadata", quietly = TRUE)) {
#'
#'   # Get wibbelink2017 data
#'   testdata <- psymetadata::wibbelink2017
#'
#'
#'   # Calculate the three-level meta-analytic model
#'   testmodel <- metafor::rma.mv(yi,
#'                                vi,
#'                                random = ~ 1 | study_id/es_id,
#'                                tdist = TRUE,
#'                                data = testdata,
#'                                method = "REML")
#'
#'   # Plot the TLMA forest plot
#'   viz_tlma_forest(x = testmodel)
#'   # Plot the thick variant of the TLMA forest plot with a table showing study
#'   # effects plus their confidence intervals
#'   viz_tlma_forest(x = testmodel, variant = "thick", annotate_CI = TRUE)
#' }


#'@export

viz_tlma_forest <- function (x, variant="classic", median_precision = FALSE, median_precision_thick = TRUE,
                             annotate_CI=FALSE, study_table=NULL, summary_table=NULL,
                             table_headers=NULL, ordered=FALSE, clouds=TRUE, spread=0.3,
                             col=FALSE,  labels=NULL, xlab="Effect Size",
                             ylab=NULL, title=NULL, confidence_level_ci = 0.95,
                             prediction_level_pi = 0.95, show_nr_ES = TRUE,
                             x_limit=NULL, table_layout = NULL, line=TRUE,
                             linewidth=0.4, text_size=3, tick_col="firebrick"
) {

  #'@import ggplot2
  #'@import dplyr
  #'@import metaviz
  #'@import metafor
  #'@importFrom magrittr %>%
  NULL
  #'@import ggbeeswarm



  cluster_vars <- all.vars(x$call$random)

  if (!all(cluster_vars %in% names(x$data))) {
    stop("At least one cluster variable (study and/or effect IDs) cannot be matched to the dataset.
         Make sure the cluster variables used for model fitting are named the same as the
         respective columns in the dataset.")
  }

  cluster_ids <- x$data[x$not.na, cluster_vars]
  ID <- cluster_ids[,1]
  ID2 <- cluster_ids[,2]



  n_ID <- max(ID)

  if("rma.mv" %in% class(x)) {
    yi <- as.numeric(x$yi)
    se <- as.numeric(sqrt(x$vi))
    vi <- as.numeric(x$vi)
    n <- length(yi)
    n_ID <- length(unique(ID))

  } else {
    stop("The input has to be a rma.mv model.")
  }

  if(is.null(ID) || is.null(ID2)) {
    stop("Please provide the arguments study_ID and effect_ID")
  } else if (length(ID) != length(yi) || length(ID2) != length(yi)) {
    stop("study_ID and effect_ID have to be of the same length as you model input dataset.")
  } else {
    # ID <- as.factor(ID)
    # ID2 <- as.factor(ID2)
  }

  if (!is.numeric(confidence_level_ci) ||
      length(confidence_level_ci) != 1 ||
      is.na(confidence_level_ci) ||
      confidence_level_ci < 0.01 ||
      confidence_level_ci > 0.99) {
    stop("confidence_level_ci must be a single number between 0.01 and 0.99.")
  }

  if (!is.numeric(prediction_level_pi) ||
      length(prediction_level_pi) != 1 ||
      is.na(prediction_level_pi) ||
      prediction_level_pi < 0.01 ||
      prediction_level_pi > 0.99) {
    stop("prediction_level_pi must be a single number between 0.01 and 0.99.")
  }

  group <- NULL

  if(is.null(group)) {
    group <- as.factor(rep(1, times = length(yi)))
  } else {
    group <- as.factor(group)
  }

  # drop unused levels of group factor
  group <- droplevels(group)
  k <- length(levels(group))

  data <- data.frame(yi, se, vi, ID, ID2, group)

  if (!is.null(labels)) {
    if (length(labels) != length(yi)) {
      warning("Labels must be the same length as the dataset.")
      labels <- NULL
    }
    labels <- as.character(labels)
    data$labels <- labels
  }



  # CI preparation
  alpha <- 1 - confidence_level_ci
  p_upper <- 1 - alpha / 2
  z_crit <- stats::qnorm(p_upper)



  model <- x

  estimate <- round(model$b[1], 2)
  var_bs <- model$sigma2[1]    # between-studies variance
  var_ws <- model$sigma2[2]  # within-study variance



  # PIs for overall effect
  pred <- predict.rma(model, level = prediction_level_pi)
  pi_lb <- pred$pi.lb
  pi_ub <- pred$pi.ub

  # CIs
  data <- data %>%
    mutate(ci_lb = yi - se * stats::qnorm(p_upper)) %>%
    mutate(ci_ub = yi + se * stats::qnorm(p_upper))

  ##############################################################################
  ##############################################################################
  ####        #########  #####################################  ################
  ####        #########  #####################################  ################
  ####        #########                                         ######      ####
  ####        #########  #####################################  ######      ####
  #####################  #####################################  ######      ####
  ##############################################################################
  ##############################################################################

  # creating a separate dataset for study-level information

  yi_ID <- NULL # to avoid no visible binding for global variable note
  se_ID <- NULL
  k <- NULL
  ci_lb_ID <- NULL
  ci_ub_ID <- NULL
  ci_lb_ES <- NULL
  ci_ub_ES <- NULL
  weight_ID <- NULL
  type <- NULL

  studydata = data.frame(
    yi_ID = numeric(n_ID),
    se_ID = numeric(n_ID),
    k = numeric(n_ID),
    ci_lb_ID = numeric(n_ID),
    ci_ub_ID = numeric(n_ID),
    ci_lb_ES = numeric(n_ID),
    ci_ub_ES = numeric(n_ID),
    weight_ID = numeric(n_ID),
    type = numeric(n_ID)
  )

  # add ID column for merging later on
  studydata$ID <- unique(data$ID)

  # appending study_table info to study_length df
  if (!is.null(study_table)) {
    study_table <- study_table %>% rename_with( ~ paste0("studytbl_", .x))
    studytbl_names <- names(study_table)
    studydata <- cbind(studydata, study_table)
  }


  row <- 1


  ###############################


  #for (i in 1:max(data$ID)){
  for (i in unique(data$ID)){
    subdata<-subset(data, ID==i)
    uni=nrow(subdata)

    if (uni==1) {
      studydata$yi_ID[row] <- subdata$yi
      studydata$se_ID[row] <- subdata$se
      studydata$ci_lb_ID[row] <-  subdata$yi - (subdata$se * z_crit)
      studydata$ci_ub_ID[row] <-  subdata$yi + (subdata$se * z_crit)
      studydata$ci_lb_ES[row] <- subdata$yi - (subdata$se * z_crit)
      studydata$ci_ub_ES[row] <- subdata$yi + (subdata$se * z_crit)
      studydata$weight_ID[row] <- 1 / subdata$se^2
      studydata$type[row] <- "singleES"
    }
    else {
      model_ID <- metafor::rma.uni(yi = subdata$yi, sei = subdata$se, method = "REML", data = subdata,
                                   level = confidence_level_ci*100)

      diagonal <- 1/(subdata$vi + var_ws)
      D <- diag(diagonal)
      obs <- nrow(subdata)
      I <- matrix(c(rep(1, (obs^2))), nrow = obs)
      M <- D%*%I%*%D
      inv_sumVar <- sum(1/(subdata$vi + var_ws))
      O <- 1/((1/var_bs) + inv_sumVar)
      V <- D - (O*M)
      T <- as.matrix(subdata$yi)
      X <- matrix(c(rep(1, obs)), ncol=1)
      var_effect <- solve(t(X)%*%V%*%X)


      studydata$yi_ID[row] <- model_ID$b
      studydata$se_ID[row] <- model_ID$se
      studydata$ci_lb_ID[row] <-  model_ID$ci.lb
      studydata$ci_ub_ID[row] <-  model_ID$ci.ub
      studydata$ci_lb_ES[row]<- model_ID$b - z_crit * stats::median(subdata$se)
      studydata$ci_ub_ES[row]<- model_ID$b + z_crit * stats::median(subdata$se)
      studydata$weight_ID[row]<- 1/ var_effect
      studydata$type[row] <- "multiES"
    }

    studydata$k[row]<-nrow(subdata)
    studydata$J[row] <- c(paste("J =",studydata$k[i]))

    row <-  row + 1
  }

  data <- merge(data, studydata, by = "ID", all.x = TRUE)

  # arrange studydata for table plotting later
  if (ordered == TRUE) {
    studydata <- studydata %>%
      arrange(yi_ID)
  }

  # extract study_table info after arranging
  if (!is.null(study_table)) {
    study_table <- select(studydata, all_of(studytbl_names))
    study_table <- study_table %>%
      rename_with( ~ sub("^studytbl_", "", .x))
  }
  #################

  # calculate relative weight of studies
  rel_weight <- NULL # to avoid no visible binding for global variable note

  sum_weight_ID <- sum(unique(data$weight_ID))
  data <- data %>%
    mutate(rel_weight = weight_ID / sum_weight_ID)

  if (ordered == TRUE) {
    data <- data %>%
      arrange(yi_ID)
  }

  # create df for only the study effects
  study_data <- data %>%
    distinct(ID, group)


  #############################################################################################
  ########### ##   #######    ## ## ##  ## ###   ##   ## ###  ###  ## #########################
  ########### ## ## ######  #### ## ## # # ## ###### ### ## ## ## # # #########################
  ########### ##   ####### ######  ### ##  ###   ### ### ###  ### ##  #########################
  #############################################################################################


  group <- as.factor(study_data$group)
  n <- nrow(study_data)


  ids <- function(group, n) {
    k <- length(levels(group))
    ki_start <- cumsum(c(1, as.numeric(table(group))[-k] + 3))
    ki_end <- ki_start + as.numeric(table(group)) - 1

    study_IDs <- numeric(n)
    for (i in 1:k) {
      study_IDs[group == levels(group)[i]] <- ki_start[i]:ki_end[i]
    }

    summary_IDs <- ki_end + 2

    data.frame("y_ID" = (n + 3 * k - 2) - c(study_IDs, summary_IDs),
               "type" = factor(c(
                 rep("study", times = length(study_IDs)), rep ("summary", times = length(summary_IDs))
               )))
  }

  IDs <- ids(group, n = n)
  # called IDs instead of ID to not mix up with ID column in dataframe
  ID_study <-  IDs$y_ID[IDs$type == "study"]
  ID_summary <-  IDs$y_ID[IDs$type == "summary"]

  study_data <- cbind(study_data, ID_study)

  # join data for beeswarm plotting
  data <- data %>%
    inner_join(study_data, by = c("ID", "group"))



  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################
  ##############################################################################


  # for "errorbar" thickness corresponding to study weights
  y_min <- NULL # to avoid no visible binding for global variable note
  y_max <- NULL
  y <- NULL

  data <- data %>%
    mutate (y_max = ID_study + rel_weight / (4 * max(rel_weight)),
            y_min = ID_study - rel_weight / (4 * max(rel_weight)))

  # ticksize and tickdata
  tick_size <- max(data$rel_weight/(6 * max(data$rel_weight)))
  tickdata <- data.frame(x = c(data$yi_ID, data$yi_ID),
                         ID = c(data$ID_study, data$ID_study),
                         y = c(data$ID_study + tick_size,
                               data$ID_study - tick_size))

  # height of the polygon
  est_wd = ifelse(max(as.numeric(data$ID_study) >= 150), 0.8, 0.5)

  # position of the summary effect
  summary_yi <- stats::coef(model)
  summary_se <- model$se

  poly_pos <- ID_summary
  # polygon coordinates
  poly_data <- data.frame(
    y=c(poly_pos,
        poly_pos -est_wd,
        poly_pos,
        poly_pos +est_wd),
    x=c(model$ci.lb,
        summary_yi,
        model$ci.ub,
        summary_yi)
  )

  # set color map
  ran_pltt <- function(n) {
    hues <- stats::runif(n)
    lightness <- stats::runif(n, 0.3, 0.5)
    saturation <- 1

    # Convert HSL to RGB using the grDevices package
    colors <- grDevices::hcl(h = hues * 360, c = saturation * 100, l = lightness * 100)
    return(colors)
  }


  if (col == TRUE) {
    pltt <- ran_pltt(n)
  } else {
    if (variant == "classic") {
      pltt <- rep("grey50", n)
    } else if (variant == "thick") {
      pltt <- rep("black", n)
    }
  }


  color_map <- stats::setNames(pltt, study_data$ID)



  ### y-axis
  summary_label <- "Summary"

  if (is.null(labels)) {
    study_labels <- study_data$ID
  } else {
    study_labels <- data %>%
      distinct(ID, .keep_all = TRUE) %>%
      pull(labels)
  }

  y_breaks <- sort(c(ID_study, ID_summary), decreasing = TRUE)
  y_tick_names <- c(as.vector(study_labels), as.vector(summary_label))[order(c(ID_study, ID_summary), decreasing = T)]
  j_labels <- c(as.vector(studydata$J), "")[order(c(ID_study, ID_summary), decreasing = T)]

  # data for beeswarm
  data_multi <- data %>% filter(type == "multiES")

  # coordinate system axes limits
  y_limit <- c(min(ID_study) - 3, max(ID_study) + 1.5)

  if(is.null(x_limit)) {
    if (median_precision == FALSE) {
      all_x <- c(
        studydata$ci_lb_ID,
        studydata$ci_ub_ID,
        data$yi
      )
    } else {
      all_x <- c(
        studydata$ci_lb_ID,
        studydata$ci_ub_ID,
        studydata$ci_lb_ES,
        studydata$ci_ub_ES,
        data$yi
      )
    }

    x_range <- range(all_x, na.rm = TRUE)

    x_limit <- c(
      x_range[1] - diff(x_range) * 0.05,
      x_range[2] + diff(x_range) * 0.05
    )
  }


  ##############################################################################
  #############    ### #######   ###     #######################################
  ############# ### ## ###### ### #### #########################################
  #############    ### ###### ### #### #########################################
  ############# ###### ###### ### #### #########################################
  ############# ######     ###   ##### #########################################
  ##############################################################################

  p <- ggplot(data = data, aes(y = ID_study))

  if (ordered==TRUE) {
    if (line == TRUE) {
      p <- p + geom_line(aes(x = yi_ID), color = "black", lwd = linewidth)
    }
  }

  if(variant=="classic") {
    if(median_precision == TRUE) {
      # median ES within study
      if (median_precision_thick == TRUE)  {
        p <- p + geom_errorbar(aes(xmin = ci_lb_ES, xmax = ci_ub_ES, lwd = k),
                               color = "grey80", width = 0)
      } else {
        p <- p + geom_errorbar(aes(xmin = ci_lb_ES, xmax = ci_ub_ES),
                               color = "grey80", width = 0, linewidth = 0.6)
      }

    }
    # study ES
    p <- p + geom_errorbar(aes(xmin = ci_lb_ID, xmax = ci_ub_ID), color = "black",
                           width = 0, lwd=0.6)+   # lwd = 0.5
      # study ES
      geom_point(aes(x = yi_ID, size=weight_ID), shape = 15, color = "black")
    # single ES
    if(clouds==TRUE){
      p <- p + ggbeeswarm::geom_quasirandom(data = data_multi, aes(x = yi, y = ID_study, color = as.factor(ID)),
                                            groupOnX = FALSE, orientation = "y", size = ifelse(nrow(data)>50,1,1.8), width = spread)
    }
  } else if (variant=="thick") {
    if(median_precision == TRUE) {
      # median ES within study
      if (median_precision_thick == TRUE)  {
        p <- p + geom_errorbar(aes(xmin = ci_lb_ES, xmax = ci_ub_ES, lwd = k),
                               color = "grey80", width = 0)
      } else {
        p <- p + geom_errorbar(aes(xmin = ci_lb_ES, xmax = ci_ub_ES),
                               color = "grey80", width = 0, linewidth = 0.6)
      }
    }
    p <- p +
      # study ES
      geom_errorbar(aes(xmin = ci_lb_ID, xmax = ci_ub_ID), width = 0)+
      # "errorbar" thickness
      geom_rect(aes(xmin = ci_lb_ID, xmax = ci_ub_ID, ymin = y_min, ymax = y_max,
                    group = ID_study), linewidth = 0.1)
    # single ES
    if(clouds==TRUE){
      p <- p + ggbeeswarm::geom_quasirandom(data = data_multi, aes(x = yi, y = ID_study, color = as.factor(ID)),
                                            groupOnX = FALSE, orientation = "y", size = ifelse(nrow(data)>50,1,1.8), width = spread)
    }
    # ES ticks
    p <- p + geom_line(data = tickdata, aes(x = x, y = y, group = ID), col = tick_col, linewidth = 1.5)
  }



  # PI for summary effect
  p <- p + annotate("errorbar", xmin = pi_lb, xmax = pi_ub, y = poly_pos, width = 0, lwd = 0.6) +
    # summary effects
    geom_polygon(data=poly_data, aes(x=x, y=y), color="black", linewidth=.4, fill="black")+

    geom_vline(xintercept = summary_yi, linetype = "dashed", lwd = .4, alpha = .5) +
    geom_vline(xintercept = 0, lwd = .4, alpha= .5)+
    scale_x_continuous(breaks=c(seq(round(min(data$yi, na.rm=TRUE))-.5,
                                    to = round(max(data$yi, na.rm=TRUE))+.5, by = .5)))+
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand=F)+
    labs(x=xlab, y=ylab, title=title)+
    guides(size="none", linewidth="none", color="none") +
    theme_bw()+
    theme(text = element_text(size = 1/0.352777778*text_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_line("grey"),
          plot.title = element_text(hjust = 0.5)
    ) +
    scale_color_manual(values = color_map)

  if (is.null(study_table)){
    p <- p + scale_y_continuous(breaks = y_breaks,
                                labels = y_tick_names)
  } else {
    p <- p + theme(axis.ticks.y.left = element_blank(),
                   axis.text.y.left  = element_blank())
  }





  ##############################################################################
  ###############################              #################################
  #####################################   ######################################
  #####################################   ######################################
  #####################################   ######################################
  #####################################   ######################################
  #####################################   ######################################
  ##############################################################################


  # Construct tableplots with study and summary information --------
  if(annotate_CI == TRUE || !is.null(study_table) || !is.null(summary_table)) {


    y_limit <- c(min(ID_study) - 3, max(ID_study) + 1.5)

    # Function to create table plots
    table_plot <- function(tbl, ID, r = 5.5, l = 5.5, tbl_titles = NULL) {
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

      # For study labels with newlines in it, the width of the column is now set according to longest line and not the whole label
      nchar2<-function(x){unlist(sapply(strsplit(x,"\n"), function(x) max(nchar(x, keepNA = FALSE))))}
      area_per_column <- cumsum(c(1, apply(rbind(tbl_titles, tbl), 2, function(x) max(round(max(nchar2(x))/100, 2),  0.03))))
      #area_per_column <- cumsum(c(1, apply(rbind(tbl_titles, tbl), 2, function(x) max(round(max(nchar(x, keepNA = FALSE))/100, 2),  0.03))))

      x_values <- area_per_column[1:ncol(tbl)]
      x_limit <- range(area_per_column)

      lab <- data.frame(y = rep(ID, ncol(tbl)),
                        x = rep(x_values,
                                each = length(ID)),
                        value = v, stringsAsFactors = FALSE)


      lab_title <- data.frame(y = rep(max(ID) + 1, times = length(tbl_titles)),
                              x = x_values,
                              value = tbl_titles)

      # To avoid "no visible binding for global variable" warning for non-standard evaluation
      y <- NULL
      value <- NULL
      ggplot(lab, aes(x = x, y = y)) +
        geom_text(aes(label = value), size = text_size, hjust = 0, vjust = 0.5) +
        geom_text(data = lab_title, aes(x = x, y = y, label = value), size = text_size, hjust = 0, vjust = 0.5) +
        coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F) +
        geom_hline(yintercept = max(ID) + 0.5) +
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
        if(nrow(study_table) != n_ID) stop('study_table must be a data.frame with one row for each study.')
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
          if(!is.data.frame(study_table)) study_table <- data.frame(study_table)
          study_table <- data.frame(lapply(study_table, as.character), stringsAsFactors = FALSE)
          if(nrow(study_table) != n_ID) stop('study_table must be a data.frame with one row for each study.')
          summary_table <- as.data.frame(matrix(rep("", times = ncol(study_table) * k), ncol = ncol(study_table)), stringsAsFactors = FALSE)
          summary_table <- stats::setNames(summary_table, names(study_table))
        }
        # Case only summary table is supplied
        if(is.null(study_table)) {
          if(!is.data.frame(summary_table)) summary_table <- data.frame(summary_table)
          summary_table <- data.frame(lapply(summary_table, as.character), stringsAsFactors = FALSE)
          if(nrow(summary_table) != k) stop('summary_table must be a data.frame with one row for each summary effect.')
          study_table <- as.data.frame(matrix(rep("", times = ncol(summary_table) * n), ncol = ncol(summary_table)), stringsAsFactors = FALSE)
          study_table <- stats::setNames(study_table, names(summary_table))
        }
      }

      table_left <- data.frame(rbind(study_table, summary_table))


      # set table headers
      if(!is.null(table_headers)) {
        if(length(table_headers) >= ncol(table_left)) {
          table_headers_left <- table_headers[1:ncol(table_left)]
        } else {
          warning("Argument table_headers has not the right length and is ignored.")
          table_headers_left <- NULL
        }
      } else {
        table_headers_left <- NULL
      }

      table_left_plot <- table_plot(table_left, ID = IDs$y_ID, r = 0, tbl_titles = table_headers_left)
    } else {
      table_left <- NULL
    }

    # Textual CI and effect size values right
    if(annotate_CI == TRUE) {

      # set table headers
      if(!is.null(table_headers)) {
        if(is.null(table_left)) {
          if(length(table_headers) == 1) {
            table_headers_right <- table_headers
          } else {
            warning("Argument table_headers has not the right length and is ignored.")
            table_headers_right <- NULL
          }
        } else {
          if(length(table_headers) == ncol(table_left) + 1) {
            table_headers_right <- table_headers[ncol(table_left) + 1]
          } else {
            table_headers_right <- NULL
          }
        }
      } else {
        table_headers_right <- NULL
      }



      x_hat <- c(studydata$yi_ID, summary_yi)
      lb <- c(c(studydata$yi_ID, summary_yi) - stats::qnorm(1 - (1 - confidence_level_ci)/2, 0, 1)*c(studydata$se_ID, summary_se))
      ub <-  c(c(studydata$yi_ID, summary_yi) + stats::qnorm(1 - (1 - confidence_level_ci)/2, 0, 1)*c(studydata$se_ID, summary_se))


      lb <- format(round(lb, 2), nsmall = 2)
      ub <- format(round(ub, 2), nsmall = 2)
      x_hat <- format(round(x_hat, 2), nsmall = 2)

      # "right-aligning" J for plotting
      J <- as.character(studydata$k)
      max_len_J <- max(nchar(J))
      #J <- sprintf("%3s", J)
      J <- sprintf(paste0("%", max_len_J, "s"), J)
      #J <- stringr::str_pad(J, width = 3, side = "left")
      J <- c(paste(J, "  "), "", "")
      CI <- paste(x_hat, " [", lb, ", ", ub, "]", sep = "")
      PI <- paste(prediction_level_pi*100, "% PI", " [", round(pi_lb,2), ", ", round(pi_ub,2), "]", sep = "")


      # with or without J
      if (show_nr_ES == TRUE) {
        if(is.null(table_headers_right)){
          table_headers_right <- c("k", paste("Study Effect", " [", confidence_level_ci*100, "% CI]", sep = ""))
        }
        CI_label <- data.frame(J = J, CI = c(CI,PI), stringsAsFactors = FALSE)
        table_CI <- table_plot(CI_label, ID = c(ID_study, ID_summary + 0.5, ID_summary - 0.5), l = 0, r = 11,  tbl_titles = table_headers_right)

      } else {
        if(is.null(table_headers_right)){
          table_headers_right <- c(paste("Study Effect", " [", confidence_level_ci*100, "% CI]", sep = ""))
        }
        CI_label <- data.frame(CI = c(CI,PI), stringsAsFactors = FALSE)
        table_CI <- table_plot(CI_label, ID = c(ID_study, ID_summary + 0.5, ID_summary - 0.5), l = 0, r = 11,  tbl_titles = table_headers_right)

      }



    } else {
      table_CI <- NULL
    }
    # Align forest plot and table(s) -----------------------------------
    if(!is.null(table_CI) && !is.null(table_left)) {
      if(is.null(table_layout)) {
        # layout_matrix <- matrix(c(rep(1, times = ncol(table_left)), rep(2, times = 3), 3), nrow = 1)
        layout_matrix <- matrix(c(rep(1, times = ncol(table_left)), rep(2, times = 4), 3, 3), nrow = 1)
      } else {
        layout_matrix <- table_layout
      }
      p <- gridExtra::arrangeGrob(table_left_plot, p, table_CI, layout_matrix = layout_matrix)
      ggpubr::as_ggplot(p)
    } else {
      if(!is.null(table_CI) && is.null(table_left)) {
        if(is.null(table_layout)) {
          #layout_matrix <- matrix(c(1, 1, 1, 1, 2), nrow = 1)
          layout_matrix <- matrix(c(1, 1, 1, 1, 2, 2), nrow = 1)
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



