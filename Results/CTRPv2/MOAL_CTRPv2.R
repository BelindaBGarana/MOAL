# MOAL CTRPv2
# Author: Belinda B. Garana (BG)
# Created: 2023-03-13; last edit: BG 2023-05-12

rm(list = ls(all = TRUE))
library(dplyr)
library(GSA)
library(MOAL)
library(na.omit)
library(stats)
library(utils)
my.path <- "/Users/belindagarana/Documents/Graham\ Lab/Landscape/Results/"
CTRPv2.path <- paste0(my.path, "CTRPv2")

##### Step 1: load inputs #####
# load drug sensitivity data frame
CTRPv2.AUC <- utils::read.csv(
  paste0(
    "https://raw.github.com/BelindaBGarana/MOAL/shiny-app/Inputs/",
    "CTRPv2_sensitivity_AUC_accessed_2022-10-13.csv"))
colnames(CTRPv2.AUC)[1] <- "CCLE_ID"

# load drug moa set information
drug.info <- utils::read.csv(
  paste0(
    "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
    "PRISM_secondary-screen-replicate-treatment-info.csv"
  )
)
colnames(drug.info)[colnames(drug.info) == "name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]", ".", drug.info$Drug)
drug.moa <- stats::na.omit(dplyr::distinct(drug.info[, c("Drug", "moa")]))

# load gmt object
gmt <- GSA::GSA.read.gmt(
  paste0(
    "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
    "MOA_gmt_file_n6_no_special_chars.gmt"
  )
)

##### Step 2: perform MOA Landscape analysis #####
#### with or without filtering for adherent culture and ####
#### with or without drugs failing quality control (qc) ####
subset <- c("all_adherent", "wo_filtering_for_adherent_culture")
for (i in seq_len(length(subset))) {
  setwd(CTRPv2.path)
  dir.create(subset[i])
  setwd(subset[i])

  # filter for adherent culture if needed
  if (subset[i] == "all_adherent") {
    # get sample info for adherent cancer cell lines
    sample.info <- 
      utils::read.csv(
        paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/",
               "Inputs/CCLE_sample_info.csv")
      )
    adh.info <- sample.info[sample.info$culture_type == "Adherent", ]
    
    # filter drug sensitivity to only represent adherent cancer cell lines
    drug.sensitivity <- 
      CTRPv2.AUC[CTRPv2.AUC$CCLE_ID %in% adh.info$stripped_cell_line_name, ]
  } else {
    drug.sensitivity <- CTRPv2.AUC
  }
  
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

##### Step 3: compare to PRISM results ##### 
### default (with drugs failing qc) ###
# load results for PRISM (all adherent)
setwd(my.path)
setwd("PRISM/all_adherent")
PRISM.result <- utils::read.csv("MOAL_results.csv")

# load results for CTRPv2 (all adherent)
setwd(my.path)
setwd("CTRPv2/all_adherent")
CTRPv2.result <- utils::read.csv("MOAL_results.csv")

# compare WGV results to original PRISM results
setwd(my.path)
comparison.CTRPv2 <- MOAL::performance(CTRPv2.result,
                                     PRISM.result,
                                     prediction.name = "CTRPv2",
                                     og.name = "PRISM")

# store performance metrics
ggplot2::ggsave("PRISM_vs_CTRPv2_scatter_plot.pdf",
                comparison.CTRPv2$scatter.plot)
utils::write.csv("PRISM_vs_CTRPv2_results.csv",
                 comparison.CTRPv2$scatter.plot.df)
ggplot2::ggsave("PRISM_vs_CTRPv2_venn_diagram.pdf",
                comparison.CTRPv2$venn.diagram)
utils::write.csv("PRISM_vs_CTRPv2_performance.csv",
                 comparison.CTRPv2$performance.df)
ggplot2::ggsave("PRISM_vs_CTRPv2_performance.pdf",
                comparison.CTRPv2$performance.table)

### without drugs failing qc ###
# load results for PRISM (all adherent)
setwd(my.path)
setwd("PRISM/all_adherent/wo_drugs_failing_qc")
PRISM.result <- utils::read.csv("MOAL_results.csv")

# load results for CTRPv2 (all adherent)
setwd(my.path)
setwd("CTRPv2/all_adherent/wo_drugs_failing_qc")
CTRPv2.result <- utils::read.csv("MOAL_results.csv")

# compare WGV results to original PRISM results
setwd(my.path)
comparison.CTRPv2 <- MOAL::performance(CTRPv2.result,
                                       PRISM.result,
                                       prediction.name = "CTRPv2",
                                       og.name = "PRISM")

# store performance metrics
ggplot2::ggsave("PRISM_vs_CTRPv2_scatter_plot.pdf",
                comparison.CTRPv2$scatter.plot)
utils::write.csv("PRISM_vs_CTRPv2_results.csv",
                 comparison.CTRPv2$scatter.plot.df)
ggplot2::ggsave("PRISM_vs_CTRPv2_venn_diagram.pdf",
                comparison.CTRPv2$venn.diagram)
utils::write.csv("PRISM_vs_CTRPv2_performance.csv",
                 comparison.CTRPv2$performance.df)
ggplot2::ggsave("PRISM_vs_CTRPv2_performance.pdf",
                comparison.CTRPv2$performance.table)


