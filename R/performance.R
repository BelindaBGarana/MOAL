performance <- function(prediction.result, og.result,
                        prediction.name = "Prediction", og.name = "Original",
                        corr.metric = "NES_avg", max.overlap = 20,
                        n.digits = 2) {
  # avoid unclear initialization
  qc_pass <- NULL
  self <- NULL
  percent_overlap.og <- NULL
  pair <- NULL

  ### load themes for plots
  ng.theme.w.legend.title <- ggplot2::theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_line(colour = "black"),
    axis.title.y = element_text(size = 8, colour = "black")
  )

  bg.theme.larger <- ggplot2::theme(
    legend.background = element_rect(), legend.position = "right",
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 14),
    axis.title.x = element_text(size = 40),
    axis.text.x = element_text(size = 20),
    axis.title.y = element_text(size = 40),
    axis.text.y = element_text(size = 20),
    plot.title = element_text(lineheight = .8, face = "bold", size = 36)
  )

  ##### check classes of inputs are as expected
  testthat::expect_is(prediction.result, "data.frame")
  testthat::expect_is(og.result, "data.frame")
  testthat::expect_is(prediction.name, "character")
  testthat::expect_is(og.name, "character")
  testthat::expect_is(corr.metric, "character")
  testthat::expect_is(max.overlap, "numeric")
  testthat::expect_is(n.digits, "numeric")

  ##### check inputs for required column names
  if (!all(c(
    "pair", "moa_i", "moa_j", "self", "qc_pass",
    "percent_overlap", corr.metric
  ) %in%
    colnames(prediction.result))) {
    stop(paste(
      "prediction.result must include column names 'pair', 'moa_i',",
      "'moa_j', 'self', 'qc_pass', 'percent_overlap',",
      "and your corr.metric"
    ))
  } else if (!all(c("pair", "self", "qc_pass", "percent_overlap",
                    corr.metric) %in%
    colnames(og.result))) {
    stop(paste(
      "og.result must include column names 'pair', 'moa_i',",
      "'moa_j', 'self', 'qc_pass', 'percent_overlap',",
      "and your corr.metric"
    ))
  }

  ##### merge data sets
  og.prediction.results <- merge(og.result, prediction.result,
    by = c("pair", "self"),
    suffixes = c(".og", ".prediction")
  )

  ##### prepare data for scatter plot
  ### categorize MOA pairs based on if they pass QC in each data set
  og.prediction.results$qc_pass <- "False"

  if (nrow(
    og.prediction.results[og.prediction.results$qc_pass.prediction, ]
  ) > 0) {
    og.prediction.results[og.prediction.results$qc_pass.prediction,
                          ]$qc_pass <- "True in prediction"
  }

  if (nrow(
    og.prediction.results[og.prediction.results$qc_pass.og, ]
  ) > 0) {
    og.prediction.results[og.prediction.results$qc_pass.og,
                          ]$qc_pass <- "True in og"
  }

  if (nrow(
    og.prediction.results[og.prediction.results$qc_pass.prediction &
      og.prediction.results$qc_pass.og, ]
  ) > 0) {
    og.prediction.results[og.prediction.results$qc_pass.prediction &
      og.prediction.results$qc_pass.og, ]$qc_pass <- "True"
  }

  og.prediction.results$qc_pass <- factor(og.prediction.results$qc_pass,
    levels = c("True", "True in prediction", "True in og", "False")
  )

  #### run Pearson correlation
  Pearson.results <- stats::cor.test(
    og.prediction.results[, c(paste0(corr.metric, ".og"))],
    og.prediction.results[, c(paste0(corr.metric, ".prediction"))],
    method = "pearson"
  )

  #### format correlation parameters for scatter plot
  stats_pearson <- substitute(
    r == est * "," ~ ~"p" ~ "=" ~ p,
    list(
      est = format(as.numeric(Pearson.results$estimate), digits = 3),
      p = format(Pearson.results$p.value, digits = 3)
    )
  )

  #### factor pairs based on if they are self or not
  og.prediction.results$self <- factor(og.prediction.results$self,
    levels = c(TRUE, FALSE)
  )

  ##### generate scatter plot
  og.prediction.plot <- ggplot2::ggplot(
    data = og.prediction.results,
    aes(
      x = paste0(corr.metric, ".og"),
      y = paste0(corr.metric, ".prediction")
    )
  ) +
    ggplot2::geom_point(size = 4, aes(color = qc_pass, shape = self)) +
    ggplot2::scale_color_discrete(
      name = "Pass QC in",
      breaks = c("True", "True in prediction", "True in og", "False"),
      labels = c(
        "Both", paste(prediction.name, "only"),
        paste(og.name, "only"), "Neither"
      )
    ) +
    ggplot2::scale_shape_discrete(
      name = "Self", breaks = c(TRUE, FALSE),
      labels = c("True", "False")
    ) +
    ggplot2::geom_smooth(
      method = "lm", linewidth = 1.5, linetype = "solid",
      color = "blue", se = TRUE, na.rm = TRUE
    ) +
    ggplot2::labs(
      x = paste(og.name, corr.metric),
      y = paste(prediction.name, corr.metric)
    ) +
    ggplot2::geom_text(
      x = min(og.prediction.results[, c(paste0(corr.metric, ".og"))]),
      y = max(og.prediction.results[, c(paste0(corr.metric,
                                               ".prediction"))]),
      vjust = "inward", hjust = "inward",
      colour = "blue", parse = TRUE,
      label = as.character(as.expression(stats_pearson)), size = 14
    ) +
    ng.theme.w.legend.title +
    bg.theme.larger +
    ggrepel::geom_label_repel(
      data = subset(og.prediction.results, qc_pass == "True" & self == "FALSE" &
        percent_overlap.og < max.overlap,
      mapping = aes(label = pair)
      ), box.padding = 8, label.size = 1,
      aes(label = pair)
    )

  ##### generate contingency matrix
  #### prepare data frame for contingency matrix
  og.sig <- c(
    nrow(og.prediction.results[og.prediction.results$qc_pass ==
      "True", ]), # number sig in both
    nrow(og.prediction.results[og.prediction.results$qc_pass ==
      "True in CTRPv2", ])
  ) # number false insig
  og.insig <- c(
    nrow(og.prediction.results[og.prediction.results$qc_pass ==
      "True in CEVIChE", ]), # number false sig
    nrow(og.prediction.results[og.prediction.results$qc_pass ==
      "False", ])
  ) # number not sig in both
  contingency <- as.data.frame(og.sig)
  colnames(contingency)[1] <- "Significant"
  contingency[, c("Insignificant")] <- c(og.insig)
  rownames(contingency) <- colnames(contingency)

  #### prepare table for contingency matrix
  contingencyTable <- gridExtra::tableGrob(contingency)

  ### create plot elements for labels
  colLabel <- grid::textGrob("    Measured",
    x = 0, hjust = 0, gp = grid::gpar(fontsize = 23)
  )
  rowLabel <- grid::textGrob("Predicted",
    x = 1, hjust = 1, gp = grid::gpar(fontsize = 23)
  )

  ### add space on plot for labels
  contingencyTable <-
    gtable::gtable_add_rows(
      contingencyTable,
      unit(2, "line"), 0
    ) # add blank row at top
  contingencyTable <-
    gtable::gtable_add_cols(
      contingencyTable,
      grid::grobWidth(rowLabel), 0
    ) # add blank space on left

  ### add column label
  contingencyTable <-
    gtable::gtable_add_grob(contingencyTable, colLabel,
      t = 1, l = 3, r = ncol(contingency)
    )

  ### add row label
  contingencyTable <-
    gtable::gtable_add_grob(contingencyTable, rowLabel,
      t = 3, l = 1
    )

  ##### evaluate performance of prediction compared to true drug screen
  acc <- 100 * (og.sig[1] + og.insig[2]) / sum(c(og.sig, og.insig))
  sens <- 100 * og.sig[1] / sum(og.sig)
  spec <- 100 * og.insig[2] / sum(og.insig)
  ppv <- 100 * og.sig[1] / sum(c(og.sig[1], og.insig[1]))
  npv <- 100 * og.insig[2] / sum(c(og.sig[2], og.insig[2]))
  fpr <- 100 - spec
  fnr <- 100 - sens
  nPairs <- nrow(og.prediction.results)
  nMOAs <- length(unique(c(og.prediction.results$moa_i,
                           og.prediction.results$moa_j)))

  #### format performance metrics for table
  Metric <- c(
    "Accuracy (%)", "Sensitivity (%)", "Specificity (%)",
    "Positive Predictive Value (%)", "Negative Predictive Value (%)",
    "False Positive Rate (%)", "False Negative Rate (%)",
    "# of MOA Pairs", "# of MOAs"
  )
  perf <- as.data.frame(Metric)
  colnames(perf)[1] <- "Performance Metric"
  perf[, c(prediction.name)] <- c(
    acc, sens, spec, ppv, npv,
    fpr, fnr, nPairs, nMOAs
  )

  #### generate table of performance metrics if n.digits is numeric
  if (is.numeric(n.digits)) {
    perf.df <- perf
    perf.df[, c(prediction.name)] <- format(perf.df[, c(prediction.name)],
      digits = n.digits
    )
    performance.table <- gridExtra::tableGrob(perf.df)
  } else {
    performance.table <- NULL
    warning(paste(
      "n.digits input must be numeric in order to",
      "create performance.table output"
    ))
  }

  return(list(
    scatter.plot = og.prediction.plot,
    scatter.plot.df = og.prediction.results,
    contingency.table = contingencyTable,
    performance.df = perf,
    performance.table = performance.table
  ))
}
