# DMEA landcape: CEVIChE
# Author: Belinda B. Garana (BG)
# Created: 2022-11-30; last edit: BG 2022-12-06

rm(list = ls(all = TRUE))
library(plyr)
library(dplyr)
library(GSA)
library(MOAL)
library(na.omit)
library(readr)
library(reshape2)
library(stats)
library(utils)
my.path <- "/Users/belindagarana/Documents/Graham\ Lab/Landscape/Results/"

##### Step 1: load inputs #####
## load drug sensitivity data frame
# import published data
# downloaded from: https://saezlab.shinyapps.io/ceviche/
# under "Predicted cell viability" tab
# Compound-Ligand selection
# uploaded in compressed form
input.df <- readr::read_csv(
  paste0("https://raw.github.com/BelindaBGarana/MOAL/shiny-app/Inputs/",
         "CEVIChE_selected_viabilities_achilles_ctrp.csv.zip")
)
input.df <- input.df[input.df$pert_type=="trt_cp",]
doses <- unique(input.df$pert_idose)
times <- unique(input.df$pert_itime)

# find dose, time with max nrow
max.treatments <- 0
for(i in 1:length(doses)){
  for(j in 1:length(times)){
    dose.time.data <- input.df[input.df$pert_idose==doses[i] & input.df$pert_itime==times[j],]
    if(nrow(dose.time.data)>max.treatments){
      max.treatments <- nrow(dose.time.data)
      optimal.dose <- doses[i]
      optimal.time <- times[j]
    }
  }
}

# filter data for dose and time with most data
filtered.input <- input.df[input.df$pert_idose==optimal.dose & input.df$pert_itime==optimal.time,]
n.cell.lines <- length(unique(filtered.input$cell_id)) #67

# resolve duplicates
avg.input <- plyr::ddply(filtered.input, .(pert_iname, cell_id), summarize, CTRP_prediction = mean(CTRP_prediction, na.rm=TRUE))

# format drug sensitivity data frame
drug.sensitivity <- reshape2::dcast(avg.input, cell_id ~ pert_iname, value.var='CTRP_prediction')
colnames(drug.sensitivity)[1] <- "CCLE_ID"

## load drug moa set information
drug.info <- utils::read.csv(
  paste0(
    "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
    "PRISM_secondary-screen-replicate-treatment-info.csv"
  )
)
colnames(drug.info)[colnames(drug.info) == "name"] <- "Drug"
drug.info$Drug <- gsub("[[:punct:][:blank:]]", ".", drug.info$Drug)
drug.moa <- stats::na.omit(dplyr::distinct(drug.info[, c("Drug", "moa")]))

## reduce input data sets to just represent drugs with moa annotations
drug.moa <- drug.moa[drug.moa$Drug %in% 
                       colnames(drug.sensitivity[,2:ncol(drug.sensitivity)]), ]
drug.sensitivity <- drug.sensitivity %>% select("CCLE_ID", drug.moa$Drug)

## load gmt object
gmt <- GSA::GSA.read.gmt(
  paste0(
    "https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
    "MOA_gmt_file_n6_no_special_chars.gmt"
  )
)

##### Step 2: perform MOA Landscape analysis #####
subset <- "CEVIChE"
setwd(my.path)
dir.create(subset)
setwd(subset)

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

# clear unnecessary results for memory space
landscape.2d.results <- landscape.2d.results$results

# allow computer to cool down for 30 min
Sys.sleep(1800)

##### Step 3: evaluate performance ##### 
# load results for CTRPv2
setwd(my.path)
setwd("CTRPv2")
og.result <- utils::read.csv("MOAL_results.csv")

# compare WGV results to original PRISM results
performance.CEVIChE <- MOAL::performance(landscape.2d.results,
                                     og.result,
                                     prediction.name = "CEVIChE",
                                     og.name = "CTRPv2")

# store performance metrics
ggplot2::ggsave("CTRPv2_vs_CEVIChE_scatter_plot.pdf",
                performance.CEVIChE$scatter.plot)
utils::write.csv("CTRPv2_vs_CEVIChE_results.csv",
                 performance.CEVIChE$scatter.plot.df)
ggplot2::ggsave("CTRPv2_vs_CEVIChE_venn_diagram.pdf",
                performance.CEVIChE$venn.diagram)
utils::write.csv("CTRPv2_vs_CEVIChE_performance.csv",
                 performance.CEVIChE$performance.df)
ggplot2::ggsave("CTRPv2_vs_CEVIChE_performance.pdf",
                performance.CEVIChE$performance.table)