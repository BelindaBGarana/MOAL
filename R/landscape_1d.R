landscape_1d <- function(drug.sensitivity, drug.moa, gmt = NULL, drug = "Drug",
                         moa = "moa", sep = ", ", min.per.set = 5,
                         FDR = 0.25, qc.min = 1, rm.self.corr = TRUE,
                         scatter.plots = FALSE, FDR.scatter.plots = 0.05,
                         mtn.plots = FALSE, suppress.warnings = TRUE) {
  ## suppress warnings if needed
  if (suppress.warnings) {
    # store original warnings option
    og.warnings.option <- getOption("warn")

    # suppress warnings globally
    options(warn = -1)
  }

  ###### Step 1: Prepare inputs #####
  ##### check classes of inputs are as expected
  testthat::expect_is(drug.sensitivity, "data.frame")
  testthat::expect_is(drug.moa, "data.frame")
  testthat::expect_is(sep, "character")
  testthat::expect_is(min.per.set, "numeric")
  testthat::expect_is(FDR, "numeric")
  testthat::expect_is(qc.min, "numeric")
  testthat::expect_is(rm.self.corr, "logical")
  testthat::expect_is(scatter.plots, "logical")
  testthat::expect_is(FDR.scatter.plots, "numeric")
  testthat::expect_is(mtn.plots, "logical")

  #### reduce input drug.sensitivity data frame if possible
  ### remove any drugs without MOA annotation
  ## only keep MOA annotations for drugs in input data frame
  drug.moa <- drug.moa[drug.moa[, c(drug)] %in% colnames(drug.sensitivity), ]

  ## identify drugs with MOA annotations
  drugs <- unique(drug.moa[, c(drug)])

  ## reduce input drug.sensitivity data frame
  ## and replace sample.names afterwards
  sample.names <- colnames(drug.sensitivity)[1]
  drug.sensitivity <- cbind(drug.sensitivity[ , 1],
                            drug.sensitivity[ , c(drugs)])
  colnames(drug.sensitivity)[1] <- sample.names

  #### generate gmt object if needed
  if (is.null(gmt)) {
    gmt <- DMEA::as_gmt(drug.moa, element.names = drug,
                        set.names = moa,
                        min.per.set = min.per.set,
                        sep = sep)
  }

  ###### Step 2: Run drug-level landscape #####
  #### for each drug, run correlations with other drugs
  #### followed by DMEA using Pearson estimates
  corr.results <- list()
  corr.plots <- list()
  drugSEA.results <- list()
  all.mtn.plots <- list()
  drugs <- colnames(drug.sensitivity)[2:ncol(drug.sensitivity)]
  for (i in seq_len(length(drugs))) {
    message(paste("Evaluating drug", i, "out of",
                length(drugs)))

    ### extract data for drug of interest
    rank.df <- as.data.frame(cbind(drug.sensitivity[ , 1],
                                   drug.sensitivity[, c(drugs[i])]))
    colnames(rank.df)[1] <- colnames(drug.sensitivity)[1]
    colnames(rank.df)[2] <- drugs[i]
    rank.df <- stats::na.omit(rank.df)
    colnames(rank.df)[2] <- "Rank"

    ### require at least 3 unique data points to run correlations
    if (length(unique(rank.df[, 2])) > 2) {
      ## make sure data for drug of interest is in second column
      drug.df <- merge(rank.df, drug.sensitivity)

      ## run correlations and remove self correlation as appropriate
      corr.AUC <- DMEA::rank_corr(data = drug.df,
                                  plots = scatter.plots,
                                  FDR = FDR.scatter.plots)
      if (rm.self.corr) {
        corr.AUC$result <-
          corr.AUC$result[corr.AUC$result[ , c(drug)] != drugs[i], ]
      }
      corr.results[[drugs[i]]] <- corr.AUC$result
      corr.plots[[drugs[i]]] <- corr.AUC$scatter.plots

      ## annotate drugs with MOAs
      corr.result <- merge(corr.AUC$result, drug.moa)

      ## identify MOA(s) for given drug
      self.MOAs <- drug.moa[drug.moa[ , c(drug)] == drugs[i], c(moa)]
      self.MOAs <- strsplit(self.MOAs, sep)[[1]]

      ## run drugSEA
      if (mtn.plots) {
        drugSEA.AUC <- DMEA::drugSEA(data = corr.result, gmt = gmt, drug = drug,
                                     set.type = moa, FDR = FDR,
                                     sep = sep, min.per.set = min.per.set)
        all.mtn.plots[[drugs[i]]] <- drugSEA.AUC$mtn.plots
      } else {
        drugSEA.AUC <- DMEA::drugSEA(data = corr.result, gmt = gmt, drug = drug,
                                     set.type = moa, FDR = 0,
                                     sep = sep, min.per.set = min.per.set)
      }
      drugSEA.result <- drugSEA.AUC$result

      ## check if expected MOAs are enriched per MOA annotations for given drug
      drugSEA.result$self <- FALSE
      drugSEA.result$qc_pass <- NA
      for (j in seq_len(length(self.MOAs))) {
        # if self MOA was evaluated, identify self-enrichment results
        if (self.MOAs[j] %in% drugSEA.result$Drug_set) {
          drugSEA.result[drugSEA.result$Drug_set == self.MOAs[j], ]$self <- TRUE
          self.enr <- drugSEA.result[drugSEA.result$Drug_set == self.MOAs[j], ]

          # determine if drug meets threshold for quality control
          if (self.enr$NES >= qc.min) {
            drugSEA.result[drugSEA.result$Drug_set == self.MOAs[j],
                           ]$qc_pass <- TRUE
          } else if (self.enr$NES < qc.min) {
            drugSEA.result[drugSEA.result$Drug_set == self.MOAs[j],
                           ]$qc_pass <- FALSE
          }
        }
      }
      drugSEA.results[[drugs[i]]] <- drugSEA.result
    }
  }

  #### collapse results from lists into data frames
  all.corr <- data.table::rbindlist(corr.results, use.names = TRUE,
                                    idcol = "Drug for Rank")
  all.drugSEA <- data.table::rbindlist(drugSEA.results, use.names = TRUE,
                                       idcol = "Drug")

  ###### Step 3: Prepare outputs #####
  #### add MOA annotations
  all.corr <- merge(all.corr, drug.moa, by.x = "Drug for Rank", by.y = "Drug")
  colnames(all.corr)[ncol(all.corr)] <- paste(moa, "for Rank")
  all.drugSEA <- merge(all.drugSEA, drug.moa)

  #### identify drugs whose self-MOA enrichment(s) did not pass quality control
  potential.misannotations <- all.drugSEA[all.drugSEA$self &
                                            !all.drugSEA$qc_pass, ]

  ## remove suppression of warnings if needed
  if (suppress.warnings) {
    # restore original warnings option
    options(warn = og.warnings.option)
  }

  return(list(
    DMEA.results = all.drugSEA, corr.results = all.corr,
    drugs.failing.qc = potential.misannotations,
    scatter.plots = corr.plots, mtn.plots = all.mtn.plots, gmt = gmt
  ))
}
