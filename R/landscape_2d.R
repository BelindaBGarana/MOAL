landscape_2d <- function(input, drug.moa = NULL, gmt = NULL, drug = "Drug",
                         moa = "moa", sep = ", ", min.per.set = 5,
                         p = 0.05, FDR = 0.25, rm.qc.fail = FALSE,
                         max.overlap = 20, n.top = 2, mtn.plots = FALSE,
                         suppress.warnings = TRUE) {
  ## suppress warnings if needed
  if (suppress.warnings) {
    # store original warnings option
    og.warnings.option <- getOption("warn")

    # suppress warnings globally
    options(warn = -1)
  }

  # avoid unclear initialization
  NES_avg <- NULL
  Fisher_p <- NULL
  Significance <- NULL
  label <- NULL
  pair <- NULL

  ###### Step 1: Prepare inputs #####
  ##### check classes of inputs are as expected
  testthat::expect_is(input, "data.frame")
  testthat::expect_is(sep, "character")
  testthat::expect_is(min.per.set, "numeric")
  testthat::expect_is(p, "numeric")
  testthat::expect_is(FDR, "numeric")
  testthat::expect_is(rm.qc.fail, "logical")
  testthat::expect_is(max.overlap, "numeric")
  testthat::expect_is(n.top, "numeric")
  testthat::expect_is(mtn.plots, "logical")

  #### make sure required columns are available in input
  if (!all(c(drug, "Drug_set", "NES") %in% colnames(input))) {
    stop(paste("drug, 'Drug_set', and 'NES' are required",
               "columns for input data frame"))
  }

  #### make sure MOA annotations are available
  if (is.null(drug.moa)) {
    if (moa %in% colnames(input)) {
      drug.moa <- dplyr::select(input, all_of(c(drug, moa)))
      drug.moa <- drug.moa[!duplicated(drug.moa), ]
    } else {
      stop(paste("drug.moa is a required input if moa",
                 "is not a column in 'input' parameter"))
    }
  }

  #### remove drugs which did not pass quality control if needed
  if (rm.qc.fail) {
    failures <- unique(input[!input$qc_pass, c(drug)])
    input <- input[!(input[ , c(drug)] %in% failures), ]
  }

  #### generate gmt object if needed
  if (is.null(gmt)) {
    gmt <- DMEA::as_gmt(drug.moa, min.per.set = min.per.set, sep = sep)
  }

  ###### Step 2: Run MOA-level landscape #####
  ### run drugSEA for each MOA
  moas <- unique(gmt$geneset.names)
  moas <- moas[moas %in% input$Drug_set]
  all.moa.results <- list()
  all.mtn.plots <- list()
  for (i in seq_len(length(moas))) {
    message(paste("Evaluating moa", i, "out of", length(moas)))

    moa.drugSEA <- input[input$Drug_set == moas[i], ]
    if (nrow(moa.drugSEA) > 0) {
      if (mtn.plots) {
        drugSEA.moa <- DMEA::drugSEA(data = moa.drugSEA, gmt = gmt,
                                     drug = drug, set.type = moa,
                                     rank.metric = "NES", FDR = FDR,
                                     sep = sep, min.per.set = min.per.set)
        all.mtn.plots[[moas[i]]] <- drugSEA.moa$mtn.plots
      } else {
        drugSEA.moa <- DMEA::drugSEA(data = moa.drugSEA, gmt = gmt,
                                     drug = drug, set.type = moa,
                                     rank.metric = "NES", FDR = 0, sep = sep,
                                     min.per.set = min.per.set)
      }
      all.moa.results[[moas[i]]] <- drugSEA.moa$result
    }
  }
  ### collapse results from list into data frames
  all.moa.drugSEA <- data.table::rbindlist(all.moa.results, use.names = TRUE,
                                           idcol = moa)

  ### also save self-enrichments separately
  self <- all.moa.drugSEA[all.moa.drugSEA[ ,c(moa)] == all.moa.drugSEA$Drug_set, ]

  ###### Step 3: Compile results for each MOA pair #####
  message("Compiling results across moa pairs...")

  ### create list of all possible MOA pairs in alphabetical order
  pairs <- list()
  for (i in seq_len(length(moas))) {
    for (j in seq_len(length(moas))) {
      alpha.pair <- sort(c(moas[i], moas[j]))
      pairs[paste0(alpha.pair[1], " & ", alpha.pair[2])] <-
        paste0(alpha.pair[1], " & ", alpha.pair[2])
    }
  }

  ### reduce list to unique MOA pairs in alphabetical order
  pairs <- unique(pairs)

  ### create data frame to store results for each MOA pair
  pairs.df <- as.data.frame(t(as.data.frame(pairs)))
  colnames(pairs.df)[1] <- "pair"

  ## extract first moa from pair (moa_i)
  pairs.df$moa_i <- sub("\\&.*", "", pairs.df$pair)

  ## remove space at end of moa_i
  pairs.df$moa_i <- substr(pairs.df$moa_i, 1, nchar(pairs.df$moa_i) - 1)

  ## extract second moa from pair (moa_j)
  pairs.df$moa_j <- sub("^[^_]*&", "", pairs.df$pair)

  ## remove space at start of moa_j
  pairs.df$moa_j <- substr(pairs.df$moa_j, 2, nchar(pairs.df$moa_j))

  ### for each pair: calculate percent overlap, store enrichment results,
  ### and determine if MOAs are correctly self-enriched
  pairs.df[, c(
    "ES_ij", "NES_ij", "p_ij", "FDR_ij", "sig_ij", "score_ij",
    "ES_ji", "NES_ji", "p_ji", "FDR_ji", "sig_ji", "score_ji",
    "N_sig", "score_avg", "score_sd", "NES_avg", "NES_sd", "Fisher_p",
    "self", "percent_overlap", "qc1_pass", "qc_pass"
  )] <- NA
  valid.MOAs <- c()
  for (k in seq_len(nrow(pairs.df))) {
    moa.i <- pairs.df$moa_i[k]
    moa.j <- pairs.df$moa_j[k]
    self.ij <- moa.i == moa.j

    ## calculate percent overlap for MOA pair
    if (self.ij) {
      pairs.df$percent_overlap[k] <- 100
    } else {
      drugs.i <- c()
      drugs.j <- c()
      in.both <- c()
      drugs.studied <- drug.moa[drug.moa[ , c(drug)] %in% input[ , c(drug)], ]
      drugs.studied <- drugs.studied[!duplicated(drugs.studied), ]
      for (m in seq_len(nrow(drugs.studied))) {
        n.sets <- 0

        # check if drug is in MOA set i
        if (sjmisc::str_contains(drugs.studied[m, c(moa)], moa.i)) {
          drugs.i <- c(drugs.i, drugs.studied[m, c(drug)])
          n.sets <- n.sets + 1
        }

        # check if drug is in MOA set j
        if (sjmisc::str_contains(drugs.studied[m, c(moa)], moa.j)) {
          drugs.j <- c(drugs.j, drugs.studied[m, c(drug)])
          n.sets <- n.sets + 1
        }

        # count if drug is in both MOA sets
        if (n.sets == 2) {
          in.both <- c(in.both, drugs.studied[m, c(drug)])
        }
      }
      pairs.df$percent_overlap[k] <- 100 * length(unique(in.both)) /
        length(unique(c(drugs.i, drugs.j)))
    }

    ### compile results for each pair
    data.ij <- all.moa.drugSEA[all.moa.drugSEA[ , c(moa)] == moa.i &
      all.moa.drugSEA$Drug_set == moa.j, ]
    data.ji <- all.moa.drugSEA[all.moa.drugSEA[ , c(moa)] == moa.j &
      all.moa.drugSEA$Drug_set == moa.i, ]
    n.sig <- 0
    if (nrow(data.ij) > 0 & nrow(data.ji) > 0) {
      ## store ij data
      pairs.df$ES_ij[k] <- data.ij$ES
      pairs.df$NES_ij[k] <- data.ij$NES
      pairs.df$p_ij[k] <- data.ij$p_value
      pairs.df$FDR_ij[k] <- data.ij$FDR_q_value
      if (pairs.df$p_ij[k] < p & pairs.df$FDR_ij[k] < FDR) {
        pairs.df$sig_ij[k] <- TRUE
        n.sig <- n.sig + 1
      } else {
        pairs.df$sig_ij[k] <- FALSE
      }

      ## store ji data
      pairs.df$ES_ji[k] <- data.ji$ES
      pairs.df$NES_ji[k] <- data.ji$NES
      pairs.df$p_ji[k] <- data.ji$p_value
      pairs.df$FDR_ji[k] <- data.ji$FDR_q_value
      if (pairs.df$p_ji[k] < p & pairs.df$FDR_ji[k] < FDR) {
        pairs.df$sig_ji[k] <- TRUE
        n.sig <- n.sig + 1
      } else {
        pairs.df$sig_ji[k] <- FALSE
      }

      pairs.df$N_sig[k] <- n.sig
      if (self.ij) {
        pairs.df$self[k] <- TRUE
        if (n.sig == 2 & pairs.df$NES_ij[k] > 0 & 
           pairs.df$NES_ji[k] > 0) {
          valid.MOAs <- c(valid.MOAs, moa.i)
        }
      } else {
        pairs.df$self[k] <- FALSE
      }

      ## calculate score_avg
      if (!is.na(pairs.df$ES_ij[k]) & !is.na(pairs.df$ES_ji[k])) {
        # for ij combinations
        # use FDR = 0.00099 for score if FDR = 0
        if (pairs.df$FDR_ij[k] == 0) {
          pairs.df$score_ij[k] <-
            -log(0.00099, base = 10) *
            as.numeric(pairs.df$NES_ij[k])
        } else {
          pairs.df$score_ij[k] <-
            -log(as.numeric(pairs.df$FDR_ij[k]), base = 10) *
            as.numeric(pairs.df$NES_ij[k])
        }

        # for ji combinations
        # use FDR = 0.00099 for score if FDR = 0
        if (pairs.df$FDR_ji[k] == 0) {
          pairs.df$score_ji[k] <-
            -log(0.00099, base = 10) *
            as.numeric(pairs.df$NES_ji[k])
        } else {
          pairs.df$score_ji[k] <-
            -log(as.numeric(pairs.df$FDR_ji[k]), base = 10) *
            as.numeric(pairs.df$NES_ji[k])
        }

        # average across ij & ji combinations
        pairs.df$score_avg[k] <- base::mean(c(pairs.df$score_ij[k],
                                        pairs.df$score_ji[k]))
        pairs.df$score_sd[k] <- stats::sd(c(pairs.df$score_ij[k],
                                     pairs.df$score_ji[k]))

        pairs.df$NES_avg[k] <- base::mean(c(pairs.df$NES_ij[k],
                                      pairs.df$NES_ji[k]))
        pairs.df$NES_sd[k] <- stats::sd(c(pairs.df$NES_ij[k],
                                   pairs.df$NES_ji[k]))

        # calculate Fisher p-value
        p.vals <- c(pairs.df$p_ij[k], pairs.df$p_ji[k])

        # use p = 0.00099 for Fisher p-value calculation if p = 0
        if (p.vals[1] == 0) {
          p.vals[1] <- 0.00099
        }
        if (p.vals[2] == 0) {
          p.vals[2] <- 0.00099
        }
        pairs.df$Fisher_p[k] <- as.numeric(metap::sumlog(p.vals)$p)
      }
    }
  }

  ### annotate pairs where both MOAs are self-enriched (qc1_pass) and pass
  ### quality control (N_sig == 2, sign(NES_ij) == sign(NES_ji), & qc1_pass)
  pairs.df[, c("qc1_pass", "qc_pass")] <- FALSE
  for (k in seq_len(nrow(pairs.df))) {
    if (pairs.df$moa_i[k] %in% valid.MOAs & pairs.df$moa_j[k] %in% valid.MOAs) {
      pairs.df$qc1_pass[k] <- "TRUE"
      if (pairs.df$N_sig[k] %in% 2 & pairs.df$NES_ij[k] /
          pairs.df$NES_ji[k] > 0) {
        pairs.df$qc_pass[k] <- TRUE
      }
    }
  }

  ## plot non-self MOA pairs which pass quality control
  plot.data <- pairs.df[pairs.df$qc_pass & !pairs.df$self, ]
  # define x limits
  limit.x <- ceiling(max(abs(as.numeric(plot.data$NES_avg)), na.rm = TRUE))
  if (nrow(plot.data[plot.data$Fisher_p < p, ]) > 0) {
    # categorize data by significance level if there are significant hits
    plot.data$Significance <- paste0("p > ", p)
    plot.data[plot.data$Fisher_p < p, ]$Significance <- paste0("p < ", p)
    plot.data$Significance <- factor(plot.data$Significance,
                                     levels = c(paste0("p < ", p),
                                                paste0("p > ", p)))

    # add labels for top hits which meet max.overlap threshold
    plot.data$label <- FALSE
    n.top.data <- plot.data[plot.data$percent_overlap <= max.overlap &
      plot.data$Significance == paste0("p < ", p), ] %>%
      dplyr::slice_max(abs(NES_avg), n = n.top)
    plot.data[plot.data$pair %in% n.top.data$pair, ]$label <- TRUE

    volc <- ggplot2::ggplot(data = plot.data, ggplot2::aes(x = NES_avg,
                                                  y = -log(Fisher_p, 10),
                                                  color = Significance)) +
      ggplot2::geom_point(size = 4) +
      ggrepel::geom_text_repel(data = subset(plot.data, label),
                               mapping = aes(label = pair, size = I(4)),
                               nudge_y = 0.25) +
      ggplot2::scale_color_manual(values = c("red", "azure4"),
                                  name = "Significance",
                                  breaks = c(paste0("p < ", p),
                                             paste0("p > ", p))) +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::xlab("Average NES") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          linewidth = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        legend.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        panel.background = element_rect(fill = "white", colour = "white",
                                        linewidth = 0.5, linetype = "solid",
                                        color = "black"),
        text = element_text(size = 10),
        legend.position = "bottom", legend.key = element_blank()
      )
  } else {
    volc <- ggplot2::ggplot(data = plot.data, ggplot2::aes(x = NES_avg,
                                                  y = -log(Fisher_p, 10))) +
      ggplot2::geom_point(size = 4, color = "azure4") +
      ggplot2::xlim(-limit.x, limit.x) +
      ggplot2::xlab("Average NES") +
      ggplot2::ylab("-Log(p-value)") +
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "grey",
                          size = 0.5) +
      ggplot2::theme(
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.line = element_line(colour = "black", linewidth = 0.65),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 20, face = "bold"),
        panel.background = element_rect(fill = "white", colour = "white",
                                        linewidth = 0.5, linetype = "solid",
                                        color = "black"),
        text = element_text(size = 10)
      )
  }

  ## remove suppression of warnings if needed
  if (suppress.warnings) {
    # restore original warnings option
    options(warn = og.warnings.option)
  }

  return(list(results = pairs.df, qc.results = self,
              DMEA.results = all.moa.drugSEA, volcano.plot = volc,
              mtn.plots = all.mtn.plots, gmt = gmt))
}
