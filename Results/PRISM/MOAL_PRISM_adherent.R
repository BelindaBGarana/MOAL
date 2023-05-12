# MOAL PRISM adherent
# Author: Belinda B. Garana (BG)
# Created: 2023-03-13; last edit: BG 2023-05-12

rm(list = ls(all = TRUE))
library(dplyr)
library(MOAL)
library(utils)
my.path <- "/Users/belindagarana/Documents/Graham\ Lab/Landscape/Results/PRISM/"

##### Step 1: load inputs #####
# load PRISM info
load.PRISM.info <- function(sample.info = FALSE) {
  library(dplyr)
  library(GSA)
  library(utils)

  # import drug info
  cat(file = stderr(), "About to get drug info", "\n")
  drug.info <- 
    utils::read.csv(
      paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
             "PRISM_secondary-screen-replicate-treatment-info.csv")
      )
  colnames(drug.info)[colnames(drug.info) == "name"] <- "Drug"
  drug.info$Drug <- gsub("[[:punct:][:blank:]]", ".", drug.info$Drug)
  drug.moa <- na.omit(dplyr::distinct(drug.info[, c("Drug", "moa")]))

  # load gmt
  cat(file = stderr(), "About to get gmt", "\n")
  gmt <- 
  GSA::GSA.read.gmt(
    paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
    "Inputs/MOA_gmt_file_n6_no_special_chars.gmt")
    )

  # load cell line info if needed
  if (sample.info) {
    cat(file = stderr(), "About to get sample info", "\n")
    sample.info <- 
      utils::read.csv(
        paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
               "Inputs/CCLE_sample_info.csv")
        )
  }

  return(list(drug.moa = drug.moa, gmt = gmt, sample.info = sample.info))
}

load.PRISM <- function(sample.info = FALSE) {
  library(utils)
  
  # load PRISM info
  PRISM.info <- load.PRISM.info(sample.info)

  # load PRISM drug AUC
  cat(file = stderr(), "About to get PRISM AUC data", "\n")
  PRISM.AUC <- 
    utils::read.csv(
      paste0("https://raw.github.com/BelindaBGarana/DMEA/",
             "shiny-app/Inputs/PRISM_drug_mean_AUC_6-23-21.csv")
      )
  PRISM.AUC$X <- NULL

  return(list(PRISM.AUC = PRISM.AUC, drug.moa = PRISM.info$drug.moa, 
              gmt = PRISM.info$gmt, sample.info = PRISM.info$sample.info))
}

inputs <- load.PRISM(sample.info = TRUE)
sample.info <- inputs$sample.info
drug.moa <- inputs$drug.moa
PRISM.AUC <- inputs$PRISM.AUC

# reduce drug.moa to overlap with PRISM AUC
drug.moa <- drug.moa[drug.moa$Drug %in% 
                       colnames(PRISM.AUC)[2:ncol(PRISM.AUC)], ]

# reduce PRISM data to drugs with moa
PRISM.AUC <- PRISM.AUC %>% select("CCLE_ID", drug.moa$Drug)

# reduce PRISM to only include adherent cell lines
sample.info <- sample.info[sample.info$CCLE_Name %in% PRISM.AUC$CCLE_ID, ]
adh.info <- sample.info[sample.info$culture_type == "Adherent", ]

PRISM.AUC <- PRISM.AUC[PRISM.AUC$CCLE_ID %in% adh.info$CCLE_Name, ]

##### Step 2: perform MOA Landscape analysis for each tissue subset #####
#### with or without drugs failing quality control (qc) ####
subset.adherent <- c(
  "all_adherent", "lung", "skin", "cns", "ovary",
  "pancreas", "esophagus", "breast", "urinary",
  "liver", "ua", "uterus", "gastric", "kidney",
  "crc", "soft", "bone", "pns", "bile_duct", "thyroid"
)
for (i in seq_len(length(subset.adherent))) {
  setwd(my.path)
  dir.create(subset.adherent[i])
  setwd(subset.adherent[i])

  # prepare subset information
  if (subset.adherent[i] == "all_adherent") {
    tissue.info <- adh.info
  } else if (subset.adherent[i] == "cns") {
    tissue.info <- adh.info[adh.info$lineage == "central_nervous_system", ]
  } else if (subset.adherent[i] == "urinary") {
    tissue.info <- adh.info[adh.info$lineage == "urinary_tract", ]
  } else if (subset.adherent[i] == "ua") {
    tissue.info <- adh.info[adh.info$lineage == "upper_aerodigestive", ]
  } else if (subset.adherent[i] == "crc") {
    tissue.info <- adh.info[adh.info$lineage == "colorectal", ]
  } else if (subset.adherent[i] == "soft") {
    tissue.info <- adh.info[adh.info$lineage == "soft_tissue", ]
  } else if (subset.adherent[i] == "pns") {
    tissue.info <- adh.info[adh.info$lineage == "peripheral_nervous_system", ]
  } else {
    tissue.info <- adh.info[adh.info$lineage == subset.adherent[i], ]
  }

  # prepare input drug sensitivity score data frame
  drug.sensitivity <- PRISM.AUC[PRISM.AUC$CCLE_ID %in% tissue.info$CCLE_Name, ]

  # perform drug-level landscape and store results
  landscape.1d.results <- MOAL::landscape_1d(drug.sensitivity, drug.moa, gmt)
  utils::write.csv(landscape.1d.results$DMEA.results,
    "MOAL_drug_DMEA_results.csv",
    row.names = FALSE
  )
  utils::write.csv(landscape.1d.results$drugs.failing.qc,
    "MOAL_drugs_failing_qc.csv",
    row.names = FALSE
  )
  saveRDS(
    landscape.1d.results$corr.results,
    "MOAL_correlations.rds"
  )

  # clear unneeded outputs for memory space
  landscape.1d.results <- landscape.1d.results$DMEA.results

  # allow computer to cool down for 30 min
  Sys.sleep(1800)

  # perform MOA Landscape and store results
  landscape.2d.results <- MOAL::landscape_2d(landscape.1d.results, 
                                             drug.moa, gmt)
  utils::write.csv(landscape.2d.results$results,
    "MOAL_results.csv",
    row.names = FALSE
  )
  utils::write.csv(landscape.2d.results$qc.results,
    "MOAL_qc_results.csv",
    row.names = FALSE
  )
  utils::write.csv(landscape.2d.results$DMEA.results,
    "MOAL_DMEA_results.csv",
    row.names = FALSE
  )
  ggplot2::ggsave(
    "MOAL_volcano_plot.pdf",
    landscape.2d.results$volcano.plot
  )
  saveRDS(
    landscape.2d.results$mtn.plots,
    "MOAL_mountain_plots.rds"
  )

  # clear results for memory space
  landscape.2d.results <- NULL

  # allow computer to cool down for 30 min
  Sys.sleep(1800)

  # also try performing MOA Landscape without drugs failing QC and store results
  landscape.2d.results.wo.qc.fail <- MOAL::landscape_2d(landscape.1d.results,
    drug.moa, gmt, rm.qc.fail = TRUE
  )
  dir.create("wo_drugs_failing_qc")
  setwd("wo_drugs_failing_qc")
  utils::write.csv(landscape.2d.results.wo.qc.fail$results,
    "MOAL_results.csv",
    row.names = FALSE
  )
  utils::write.csv(landscape.2d.results.wo.qc.fail$qc.results,
    "MOAL_qc_results.csv",
    row.names = FALSE
  )
  utils::write.csv(landscape.2d.results.wo.qc.fail$DMEA.results,
    "MOAL_DMEA_results.csv",
    row.names = FALSE
  )
  ggplot2::ggsave(
    "MOAL_volcano_plot.pdf",
    landscape.2d.results.wo.qc.fail$volcano.plot
  )
  saveRDS(
    landscape.2d.results.wo.qc.fail$mtn.plots,
    "MOAL_mountain_plots.rds"
  )

  # clear results for memory space
  landscape.2d.results.wo.qc.fail <- NULL

  # allow computer to cool down for 30 min
  Sys.sleep(1800)
}
