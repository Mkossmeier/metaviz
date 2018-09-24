#'Internal helper function of viz_rainforest and viz_forest to create a rainforest plot.
#'
#'Creates a rainforest plot. Called by viz_rainforest and viz_forest for type = "rain"
#'@keywords internal
viz_rainforest_internal <- function(plotdata, madata,
                                    type = "standard",
                                    study_labels = NULL, summary_label = NULL,
                                    study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                    confidence_level = 0.95, col = "Blues", summary_col = "Blues",
                                    detail_level = 1,
                                    text_size = 3, xlab = "Effect", x_limit = NULL,
                                    x_trans_function = NULL, x_breaks = NULL) {
  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))

  # weight of each study used to scale the height of each raindrop
  if(type %in% c("standard", "study_only")) {
    weight <- 1/(plotdata$se^2 + madata$summary_tau2[as.numeric(plotdata$group)])
  } else {
    weight <- 1/plotdata$se^2
  }
  plotdata$rel_weight <- weight/sum(weight)

  tick_size <- max(plotdata$rel_weight/(6 * max(plotdata$rel_weight)))
  tickdata <- data.frame(x = c(plotdata$x, plotdata$x), ID = c(plotdata$ID, plotdata$ID),
                         y = c(plotdata$ID + tick_size,
                               plotdata$ID - tick_size))

  # function ll constructs a likelihood raindrop of a study. Each raindrop is built out of
  # several distinct segments (to color shade the raindrop)
  ll <- function(x, max.range, max.weight) {
    # width of the region over which the raindop is built
    se.factor <- ceiling(stats::qnorm(1 - (1 - confidence_level)/2))
    width <- abs((x[1] - se.factor * x[2]) - (x[1] + se.factor * x[2]))

    # max.range is internally determined as the width of the broadest raindop.
    # The number of points to construct the raindrop is chosen proportional to
    # the ratio of the width of the raindrop and the max.range,
    # because slim/dense raindrops do not need as many support points as very broad ones.
    # Minimum is 200 points (for the case that width/max.range gets very small)
    length.out <- max(c(floor(1000 * width / max.range), 200))

    # Create sequence of points to construct the raindrop. The number of points is chosen by length.out
    # and can be changed by the user with the parameter detail_level (maximum 100).
    # At least 50 support points (for detail level << 1) per raindrop are chosen
    support <- seq(x[1] - se.factor * x[2], x[1] + se.factor * x[2], length.out = max(c(length.out * min(c(detail_level, 100)), 50)))

    # The values for the likelihood drop are determined: The likelihood for different hypothetical true values
    # minus likelihood for the observed value (i.e. the maximum likelihood) plus the confidence.level quantile of the chi square
    # distribution with one degree of freedeom divided by two.
    # Explanation: -2*(log(L(observed)/L(hypothetical))) is an LRT test and approx. chi^2 with 1 df and significance threshold
    # qchisq(confidence.level, df = 1).
    # That means by adding the confidence.level quantile of the chi square
    # distribution with one degree of freedom (the significance threshold) divided by two,
    # values with l_mu < 0 differ significantly from the observed value.
    threshold <- stats::qchisq(confidence_level, df = 1)/2
    l_mu <- log(stats::dnorm(x[1], mean = support, sd = x[2])) - log(stats::dnorm(x[1], mean = x[1], sd = x[2])) + threshold

    #scale raindrop such that it is proportional to the meta-analytic weight and has height smaller than 0.5
    l_mu <- l_mu/max(l_mu) * x[3]/max.weight * 0.45

    # Force raindrops of studies to have minimum height of 0.05 (i.e. approx. one tenth of the raindrop with maximum height)
    if(max(l_mu) < 0.05) {
      l_mu <- l_mu/max(l_mu) * 0.05
    }

    # mirror values for raindrop
    l_mu_mirror <- -l_mu

    # select only likelihood values that are equal or larger than zero,
    # i.e. values that also lie in the confidence interval (using normality assumption)
    sel <- which(l_mu >= 0)

    # Construct data.frame
    d <- data.frame("support" = c(support[sel], rev(support[sel])), "log_density" = c(l_mu[sel], rev(l_mu_mirror[sel])))

    # The number of segments for shading is chosen as follows: 40 segements times the detail_level per drop
    # as default. The minimum count of segments is 20 (If detail_level is << 1), with the exception that
    # if there are too few points for 20 segments then nrow(d)/4 is used (i.e. at least 4 points per segment)
    data.frame(d, "segment" = cut(d$support, max(c(40*detail_level), min(c(20, nrow(d)/4)))))
  }

  # compute the max range of all likelihood drops for function ll.
  max.range <- max(abs((plotdata$x + stats::qnorm(1 - (1 - confidence_level)/2) * plotdata$se) -
                         (plotdata$x - (stats::qnorm(1 - (1 - confidence_level)/2) * plotdata$se))))

  # computes all likelihood values and segments. The output is a list, where every element
  # constitutes one study raindop
  res <- apply(cbind(plotdata$x, plotdata$se, plotdata$rel_weight), 1,  FUN = function(x) {ll(x, max.range = max.range, max.weight = max(plotdata$rel_weight))})

  # name every list entry, i.e. raindrop, and add id column
  names(res) <- plotdata$ID
  for(i in 1:length(res)) {
    res[[i]] <- data.frame(res[[i]], .id = plotdata$ID[i])
  }

  # The prep.data function prepares the list of raindrops in three ways for plotting (shading of segments):
  # 1) the values are sorted by segments, such that the same segments of each raindrop are joined together
  # 2) segments are renamed with integer values from 1 to the number of segments per raindrop
  # 3) to draw smooth raindrops the values at the right hand boundary of each segment have to be the first
  # values at the left hand boundary of the next segment on the right.
  prep.data <- function(res) {
    res <- lapply(res, FUN = function(x) {x <- x[order(x$segment), ]})
    res <- lapply(res, FUN = function(x) {x$segment <- factor(x$segment, labels = 1:length(unique(x$segment))); x})
    res <- lapply(res, FUN = function(x) {
      seg_n <- length(unique(x$segment))
      first <- sapply(2:seg_n, FUN = function(n) {min(which(as.numeric(x$segment)==n))})
      last <-  sapply(2:seg_n, FUN = function(n) {max(which(as.numeric(x$segment)==n))})
      neighbor.top <-   x[c(stats::aggregate(support~segment, FUN = which.max, data=x)$support[1],
                            cumsum(stats::aggregate(support~segment, FUN = length, data=x)$support)[-c(seg_n-1, seg_n)] +
                              stats::aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]), c("support", "log_density")]
      neighbor.bottom <-   x[c(stats::aggregate(support~segment, FUN = which.max, data=x)$support[1],
                               cumsum(stats::aggregate(support~segment, FUN = length, data=x)$support[-c(seg_n-1, seg_n)])+
                                 stats::aggregate(support~segment, FUN = which.max, data=x)$support[-c(1, seg_n)]) + 1, c("support", "log_density")]
      x[first, c("support", "log_density")] <- neighbor.top
      x[last, c("support", "log_density")] <- neighbor.bottom
      x
    }
    )
    res
  }
  res <- prep.data(res)

  # merge the list of raindops in one dataframe for plotting
  res <- do.call(rbind, res)

  # set limits and breaks for the y axis and construct summary diamond (for type standard and sensitivity)
  if(type %in% c("standard", "sensitivity", "cumulative")) {
    y_limit <- c(min(plotdata$ID) - 3, max(plotdata$ID) + 1.5)
    y_tick_names <- c(as.vector(study_labels), as.vector(summary_label))[order(c(plotdata$ID, madata$ID), decreasing = T)]
    y_breaks <- sort(c(plotdata$ID, madata$ID), decreasing = T)
    summarydata <- data.frame("x.diamond" = c(madata$summary_es - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es,
                                              madata$summary_es + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es),
                              "y.diamond" = c(madata$ID,
                                              madata$ID + 0.3,
                                              madata$ID,
                                              madata$ID - 0.3),
                              "diamond_group" = rep(1:k, times = 4)
    )
  } else {
    y_limit <- c(min(plotdata$ID) - 1, max(plotdata$ID) + 1.5)
    y_tick_names <- plotdata$labels[order(plotdata$ID, decreasing = T)]
    y_breaks <- sort(plotdata$ID, decreasing = T)
  }

  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }

  # To shade all segments of each raindop symmetrically the min abs(log_density) per raindrop is used
  # as aesthetic to fill the segments. This is necessary because otherwise the first log_density value per
  # segment would be used leading to asymmetrical shading
  min.ld <- stats::aggregate(log_density ~ segment + .id, FUN  = function(x) {min(abs(x))}, data = res)
  names(min.ld) <- c("segment", ".id", "min_log_density")
  res <- merge(res, min.ld, sort = F)

  # Set Color palette for shading
  if(type != "summary_only") {
    if(!(col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
      warning("Supported arguments for col for rainforest plots are Blues, Greys, Oranges, Greens, Reds, and Purples. Blues is used.")
      col <- "Blues"
    }
    col <- RColorBrewer::brewer.pal(n = 9, name = col)
    if(summary_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples")) {
      summary_col <- RColorBrewer::brewer.pal(n = 9, name = summary_col)[9]
    }
  } else {
    if(type == "summary_only") {
      if(!(summary_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
        warning("Supported arguments for summary_col for summary-only rainforest plots are Blues, Greys, Oranges, Greens, Reds, and Purples. Blues is used.")
        summary_col <- "Blues"
      }
      summary_col <- RColorBrewer::brewer.pal(n = 9, name = col)
      col <- summary_col
    }
  }

  # Set plot margins. If table is aligned on the left, no y axis breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  support <- NULL
  segment <- NULL
  min_log_density <- NULL
  log_density <- NULL
  .id <-
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL
  ID <- NULL

  # Create Rainforest plot
  p <-
    ggplot(data = res, aes(y = .id, x = support)) +
    geom_errorbarh(data = plotdata, col = col[1], aes(xmin = x_min, xmax = x_max, y = ID, height = 0), inherit.aes = FALSE) +
    geom_polygon(data = res, aes(x = support, y = as.numeric(.id) + log_density,
                                 color = min_log_density, fill = min_log_density,
                                 group = paste(.id, segment)), size = 0.1) +
    geom_line(data = tickdata, aes(x = x, y = y, group = ID), col = "grey", size = 1)
    # geom_errorbarh(data = plotdata, col = "grey", aes(x = x, xmin = x_min, xmax = x_max, y = ID, height = 0))
  if(type %in% c("standard", "sensitivity", "cumulative")) {
    p <- p + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), color="black", fill = summary_col, size = 0.1)
  }
  p <- p +
    scale_fill_gradient(high = col[9], low = col[3], guide = FALSE) +
    scale_color_gradient(high = col[9], low = col[3], guide = FALSE) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
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
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))

  p
}


#'Internal helper function of viz_thickforest and viz_forest to create a thick forest plot
#'
#'Creates a thick forest plot. Called by viz_thickforest and viz_forest for type = "thick"
#'@keywords internal
viz_thickforest_internal <- function(plotdata, madata,
                                     type = "standard",
                                     study_labels = NULL, summary_label = NULL,
                                     study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                     confidence_level = 0.95, col = "Blues", summary_col = "Blues", tick_col = "firebrick",
                                     text_size = 3, xlab = "Effect", x_limit = NULL,
                                     x_trans_function = NULL, x_breaks = NULL) {

  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))

  # weight of each study used to scale the height of each raindrop
  if(type %in% c("standard", "study_only")) {
     weight <- 1/(plotdata$se^2 + madata$summary_tau2[as.numeric(plotdata$group)])
  } else {
    weight <- 1/plotdata$se^2
  }
  rel_weight <- weight/sum(weight)
  plotdata$rel_weight <- rel_weight
  plotdata <- plotdata %>%
    mutate(y_max = ID + rel_weight/(4*max(rel_weight)),
           y_min = ID - rel_weight/(4*max(rel_weight))
    )

  tick_size <- max(plotdata$rel_weight/(6*max(plotdata$rel_weight)))
  tickdata <- data.frame(x = c(plotdata$x, plotdata$x), ID = c(plotdata$ID, plotdata$ID),
                         y = c(plotdata$ID + tick_size,
                               plotdata$ID - tick_size))

  # set limits and breaks for the y axis and construct summary diamond (for type standard and sensitivity)
  if(type %in% c("standard", "sensitivity", "cumulative")) {
    y_limit <- c(min(plotdata$ID) - 3, max(plotdata$ID) + 1.5)
    y_tick_names <- c(as.vector(study_labels), as.vector(summary_label))[order(c(plotdata$ID, madata$ID), decreasing = T)]
    y_breaks <- sort(c(plotdata$ID, madata$ID), decreasing = T)
    summarydata <- data.frame("x.diamond" = c(madata$summary_es - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es,
                                              madata$summary_es + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es),
                              "y.diamond" = c(madata$ID,
                                              madata$ID + 0.3,
                                              madata$ID,
                                              madata$ID - 0.3),
                              "diamond_group" = rep(1:k, times = 4)
    )
  } else {
    y_limit <- c(min(plotdata$ID) - 1, max(plotdata$ID) + 1.5)
    y_tick_names <- plotdata$labels[order(plotdata$ID, decreasing = T)]
    y_breaks <- sort(plotdata$ID, decreasing = T)
  }


  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }

  # Set Color palette for shading
  if(type != "summary_only") {
    if(all(col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
      col <- unlist(lapply(col, function(x) RColorBrewer::brewer.pal(n = 9, name = x)[9]))
    }
  }
  if(type != "study_only") {
    if(all(summary_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
      summary_col <- unlist(lapply(summary_col, function(x) RColorBrewer::brewer.pal(n = 9, name = x)[9]))
    }
    if(type == "summary_only") {
      col <- summary_col
    } else {
      if(length(summary_col) > 1) summary_col <- rep(summary_col, times = 4)
      }
  }



  # Set plot margins. If table is aligned on the left, no y axus breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL
  y_min <- NULL
  y_max <- NULL
  ID <- NULL

  # Create thick forest plot
  p <-
    ggplot(data = plotdata, aes(y = ID, x = x)) +
    geom_errorbarh(data = plotdata, col = col, aes(xmin = x_min, xmax = x_max, y = ID, height = 0)) +
    geom_rect(aes(xmin = x_min, xmax = x_max, ymin = y_min, ymax = y_max,
                  group = ID), fill = col, size = 0.1) +
    geom_line(data = tickdata, aes(x = x, y = y, group = ID), col = tick_col, size = 1.5)
  if(type %in% c("standard", "sensitivity", "cumulative")) {
    p <- p + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col, size = 0.1)
  }
  p <- p +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
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
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))
  p
}



#'Internal helper function of viz_forest to create a classic forest plot
#'
#'Creates a classic forest plot. Called by viz_forest for type = "classic"
#'@keywords internal
viz_classicforest_internal <- function(plotdata, madata,
                                       type = "standard",
                                       study_labels = NULL, summary_label = NULL,
                                       study_table = NULL, summary_table = NULL, annotate_CI = FALSE,
                                       confidence_level = 0.95, col = "Blues", summary_col = "Blues", tick_col = "firebrick",
                                       text_size = 3, xlab = "Effect", x_limit = NULL,
                                       x_trans_function = NULL, x_breaks = NULL) {
  n <- nrow(plotdata)
  k <- length(levels(plotdata$group))

  # weight of each study used to scale the height of each raindrop
  if(type %in% c("standard", "study_only")) {
    weight <- 1/(plotdata$se^2 + madata$summary_tau2[as.numeric(plotdata$group)])
  } else {
    weight <- 1/plotdata$se^2
  }
  plotdata$rel_weight <- weight/sum(weight)

  if(type %in% c("cumulative", "sensitivity")) {
    tick_size <- max(plotdata$rel_weight/(6*max(plotdata$rel_weight)))
    tickdata <- data.frame(x = c(plotdata$x, plotdata$x), ID = c(plotdata$ID, plotdata$ID),
                           y = c(plotdata$ID + tick_size,
                                 plotdata$ID - tick_size))
  }


  # set limits and breaks for the y axis and construct summary diamond (for type standard and sensitivity)
  if(type %in% c("standard", "sensitivity", "cumulative")) {
    y_limit <- c(min(plotdata$ID) - 3, max(plotdata$ID) + 1.5)
    y_tick_names <- c(as.vector(study_labels), as.vector(summary_label))[order(c(plotdata$ID, madata$ID), decreasing = T)]
    y_breaks <- sort(c(plotdata$ID, madata$ID), decreasing = T)
    summarydata <- data.frame("x.diamond" = c(madata$summary_es - stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es,
                                              madata$summary_es + stats::qnorm(1 - (1 - confidence_level) / 2, 0, 1) * madata$summary_se,
                                              madata$summary_es),
                              "y.diamond" = c(madata$ID,
                                              madata$ID + 0.3,
                                              madata$ID,
                                              madata$ID - 0.3),
                              "diamond_group" = rep(1:k, times = 4)
    )
  } else {
    y_limit <- c(min(plotdata$ID) - 1, max(plotdata$ID) + 1.5)
    y_tick_names <- plotdata$labels[order(plotdata$ID, decreasing = T)]
    y_breaks <- sort(plotdata$ID, decreasing = T)
  }

  # set limits for the x axis if none are supplied
  if(is.null(x_limit)) {
    x_limit <- c(range(c(plotdata$x_min, plotdata$x_max))[1] - diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05,
                 range(c(plotdata$x_min, plotdata$x_max))[2] + diff(range(c(plotdata$x_min, plotdata$x_max)))*0.05)
  }

  # Set Color palette for shading
  if(type != "summary_only") {
    if(all(col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
      col <- unlist(lapply(col, function(x) RColorBrewer::brewer.pal(n = 9, name = x)[9]))
    }
  }
  if(type != "study_only") {
    if(all(summary_col %in% c("Blues", "Greys", "Oranges", "Greens", "Reds", "Purples"))) {
      summary_col <- unlist(lapply(summary_col, function(x) RColorBrewer::brewer.pal(n = 9, name = x)[9]))
    }
    if(type == "summary_only") {
      col <- summary_col
    } else {
      if(length(summary_col) > 1) summary_col <- rep(summary_col, times = 4)
    }
  }


  # Set plot margins. If table is aligned on the left, no y axus breaks and ticks are plotted
  l <- 5.5
  r <- 11
  if(annotate_CI == TRUE) {
    r <- 1
  }
  if(!is.null(study_table) || !is.null(summary_table)) {
    l <- 1
    y_tick_names <- NULL
    y_breaks <- NULL
  }
  # workaround for "Undefined global functions or variables" Note in R CMD check while using ggplot2.
  x.diamond <- NULL
  y.diamond <- NULL
  diamond_group <- NULL
  ID <- NULL
  x <- NULL
  y <- NULL
  x_min <- NULL
  x_max <- NULL

  # create classic forest plot
  p <-
    ggplot(data = plotdata, aes(y = ID, x = x)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(data = plotdata, col = "black", aes(xmin = x_min, xmax = x_max, y = ID, height = 0))

  if(type %in% c("cumulative", "sensitivity")) {
    p <- p + geom_line(data = tickdata, aes(x = x, y = y, group = ID), col = col, size = 1)
  } else {
    p <- p + geom_point(aes(size = weight), shape = 22, col = "black", fill = col)
  }

  if(type %in% c("standard", "sensitivity", "cumulative")) {
    p <- p + geom_polygon(data = summarydata, aes(x = x.diamond, y = y.diamond, group = diamond_group), color= "black", fill = summary_col, size = 0.1)
  }
  p <- p +
    scale_y_continuous(name = "",
                       breaks = y_breaks,
                       labels = y_tick_names) +
    coord_cartesian(xlim = x_limit, ylim = y_limit, expand = F)
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
    scale_size_area(max_size = 3) +
    theme_bw() +
    theme(text = element_text(size = 1/0.352777778*text_size),
          legend.position = "none",
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_line("grey"),
          panel.grid.minor.x = element_line("grey"),
          plot.margin = margin(t = 5.5, r = r, b = 5.5, l = l, unit = "pt"))
  p
}
