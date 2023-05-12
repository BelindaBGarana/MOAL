performance <- function(prediction.result, og.result,
                       prediction.name = "Prediction", og.name = "Original",
                       venn.name.size = 13, venn.text.size = 8,
                       corr.metric = "NES_avg", max.overlap = 20) {
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
  testthat::expect_is(venn.name.size, "numeric")
  testthat::expect_is(venn.text.size, "numeric")
  testthat::expect_is(corr.metric, "character")
  testthat::expect_is(max.overlap, "numeric")

  ##### check inputs for required column names
  if (!all(c("pair", "self", "qc_pass", "percent_overlap", corr.metric) %in%
        colnames(prediction.result))) {
    stop(paste("prediction.result must include column names 'pair',",
               "'self', 'qc_pass', 'percent_overlap', and your corr.metric"))
  } else if (!all(c("pair", "self", "qc_pass", 'percent_overlap', corr.metric) %in%
               colnames(og.result))) {
    stop(paste("og.result must include column names 'pair',",
               "'self', 'qc_pass', 'percent_overlap', and your corr.metric"))
  }

  ##### merge data sets
  og.prediction.results <- merge(og.result, prediction.result,
                                 by = c("pair", "self"),
                                 suffixes = c(".og", ".prediction"))

  ##### venn diagram
  venn.data <- list(
    "og" = og.prediction.results[og.prediction.results$qc_pass.og, ]$pair,
    "prediction" =
      og.prediction.results[og.prediction.results$qc_pass.prediction, ]$pair
  )
  venn <- ggvenn::ggvenn(venn.data,
    set_name_size = venn.name.size,
    text_size = venn.text.size
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
                            og.prediction.results$qc_pass.og,
                          ]$qc_pass <- "True"
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
                                       levels = c(TRUE, FALSE))

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
      labels = c("Both", paste(prediction.name, "only"),
                 paste(og.name, "only"), "Neither")
    ) +
    ggplot2::scale_shape_discrete(name = "Self", breaks = c(TRUE, FALSE),
                                  labels = c("True", "False")) +
    ggplot2::geom_smooth(method = "lm", linewidth = 1.5, linetype = "solid",
                         color = "blue", se = TRUE, na.rm = TRUE) +
    ggplot2::labs(x = paste(og.name, corr.metric),
                  y = paste(prediction.name, corr.metric)) +
    ggplot2::geom_text(
      x = min(og.prediction.results[, c(paste0(corr.metric, ".og"))]),
      y = max(og.prediction.results[, c(paste0(corr.metric, ".prediction"))]),
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

  ##### evaluate performance of prediction compared to true drug screen
  sig.og <-
    unique(og.prediction.results[og.prediction.results$qc_pass.og, ]$pair)
  sig.prediction <-
    unique(og.prediction.results[og.prediction.results$qc_pass.prediction,
                                 ]$pair)
  all.sig <- unique(c(sig.og, sig.prediction))
  success.rate <- 100 *
    length(sig.og[sig.og %in% sig.prediction]) / length(all.sig)
  false.pos.rate <- 100 *
    length(sig.prediction[!(sig.prediction %in% sig.og)]) / length(all.sig)
  false.neg.rate <- 100 *
    length(sig.og[!(sig.og %in% sig.prediction)]) / length(all.sig)
  performance.prediction <- as.data.frame(prediction.name)
  performance.prediction[, c("success.rate", "false.pos.rate",
                             "false.neg.rate", "N")] <-
    c(success.rate, false.pos.rate, false.neg.rate, length(all.sig))

  #### format performance metrics for table if nonzero
  performance.rounded <- performance.prediction

  if (!is.na(success.rate) & success.rate != 0 & success.rate != 100) {
    performance.rounded$success.rate <-
      format(performance.rounded$success.rate, digits = 3)
  }

  if (!is.na(false.pos.rate) & false.pos.rate != 0 & false.pos.rate != 100) {
    performance.rounded$false.pos.rate <-
      format(performance.rounded$false.pos.rate, digits = 3)
  }

  if (!is.na(false.neg.rate) & false.neg.rate != 0 & false.neg.rate != 100) {
    performance.rounded$false.neg.rate <-
      format(performance.rounded$false.neg.rate, digits = 3)
  }

  performance.table <- gridExtra::tableGrob(performance.rounded,
    cols = c(
      "Algorithm", "Success Rate (%)",
      "False Positive Rate (%)",
      "False Negative Rate (%)",
      "# Shared MOA Pairs Passing QC"
    ),
    rows = NULL
  )

  return(list(
    scatter.plot = og.prediction.plot,
    scatter.plot.df = og.prediction.results,
    venn.diagram = venn,
    performance.df = performance.prediction,
    performance.table = performance.table
  ))
}
