
#'@title A dataframe containing information on the study effects of a three-level meta-analysis
#'@description Provides a dataframe with information on the study effects of a three-level meta-analysis.
#'
#'@param x metafor rma.mv object
#'@param confidence_level_ci numeric confidence level for the confidence intervals of the study effects
#'@param ordered orders the data by effect size when TRUE
#'
#'@details The function viz_tlma_studyinfo creates a dataframe with information pertaining to the study effects of
#'a three-level meta-analysis. The study effects are the result of random-effects models conducted on the effects
#'contained within each respective study. The dataframe contains the number of effects originating from the study,
#'the study effect with its standard error and confidence intervals, and the weight of the study within the three-level
#'meta-analysis.
#'@return A dataframe containing statistical information about the study effects of a three-level meta-analysis is created.
#'@author Verena Pilar <verena.pilar@outlook.com>
#'@examples
#'
#' if (requireNamespace("psymetadata", quietly = TRUE)) {
#'   # Get wibbelink2017 data
#'   testdata <- psymetadata::wibbelink2017
#'
#'
#'   # Calculate the three-level meta-analytic model
#'   testmodel <- metafor::rma.mv(yi,
#'                               vi,
#'                               random = ~ 1 | study_id/es_id,
#'                               tdist = TRUE,
#'                               data = testdata,
#'                               method = "REML")
#'
#'   # Create a table with study-level information
#'   viz_tlma_studyinfo(testmodel)
#' }

#' @export

viz_tlma_studyinfo <- function (x, confidence_level_ci = 0.95, ordered = FALSE) {


  if (!is.numeric(confidence_level_ci) ||
      length(confidence_level_ci) != 1 ||
      is.na(confidence_level_ci) ||
      confidence_level_ci < 0.01 ||
      confidence_level_ci > 0.99) {
    stop("confidence_level_ci must be a single number between 0.01 and 0.99.")
  }


  group = NULL


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
  }


  if(is.null(group)) {
    group <- as.factor(rep(1, times = length(yi)))
  } else {
    group <- as.factor(group)
  }

  # drop unused levels of group factor
  group <- droplevels(group)
  k <- length(levels(group))

  data <- data.frame(yi, se, vi, ID, ID2, group)


  # CI preparation
  alpha <- 1 - confidence_level_ci
  p_upper <- 1 - alpha / 2
  z_crit <- stats::qnorm(p_upper)


  model <- x

  estimate <- round(model$b[1], 2)
  var_bs <- model$sigma2[1]    # between-studies variance
  var_ws <- model$sigma2[2]  # within-study variance


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


  row <- 1


  ###############################


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


  # arrange studydata
  if (ordered == TRUE) {
    studydata <- studydata %>%
      arrange(yi_ID)
  }

  studylevel_info <- studydata %>%
    select(
      k = k,
      yi_study = yi_ID,
      se_study = se_ID,
      ci_lb_study = ci_lb_ID,
      ci_ub_study = ci_ub_ID,
      weight_study = weight_ID
    )

  return(studylevel_info)


}
