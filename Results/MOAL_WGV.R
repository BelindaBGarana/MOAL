# MOAL WGV: PRISM adherent
# Author: Belinda B. Garana (BG)
# Created: 2023-03-13; last edit: BG 2023-05-12

rm(list = ls(all = TRUE))
library(dplyr)
library(GEOquery)
library(limma)
library(MOAL)
library(utils)
my.path <- "/Users/belindagarana/Documents/Graham\ Lab/Landscape/Results/"

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

# load CCLE RNAseq for weighted gene voting (WGV)
download.file(
  paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
         "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"), 
  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin"
  )
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_1-200.Rbin")
download.file(
  paste0("https://raw.github.com/BelindaBGarana/DMEA/shiny-app/Inputs/",
         "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"), 
  destfile = "Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin"
  )
load("Normalized_adherent_CCLE_RNAseq_19Q4_samples_in_PRISM_201-327.Rbin")
RNA.df <- rbind(RNA.first200, RNA.rest)

##### Step 2: perform WGV for each drug ##### 
#### run WGV for each drug based on genes differentially expressed ####
#### between top and bottom 10% of sensitive adherent cell lines ####
setwd(my.path)
setwd("WGV")
n.genes <- 500
CCLE_ID <- PRISM.AUC.reduced$CCLE_ID
drugs <- colnames(PRISM.AUC.reduced[,2:ncol(PRISM.AUC.reduced)])
WV.results <- data.frame(CCLE_ID)
for (i in seq_len(length(drugs))){
  #identify top, bottom 10% of samples based on drug AUC
  drug.AUC <- na.omit(PRISM.AUC[,c("CCLE_ID",drugs[i])])
  colnames(drug.AUC)[2] <- "Rank"
  top.SENS <- drug.AUC %>% slice_min(Rank, n=0.1*nrow(drug.AUC))
  top.RES <- drug.AUC %>% slice_max(Rank, n=0.1*nrow(drug.AUC))
  top.SENS.RNAseq <- RNA.df[RNA.df$CCLE_ID %in% top.SENS$CCLE_ID,]
  top.RES.RNAseq <- RNA.df[RNA.df$CCLE_ID %in% top.RES$CCLE_ID,]
  train.RNAseq <- rbind(top.SENS.RNAseq, top.RES.RNAseq)
  top.SENS.RNAseq$CCLE_ID <- NULL
  top.RES.RNAseq$CCLE_ID <- NULL
  
  #use eBayes to identify differentially expressed genes
  if(nrow(top.SENS.RNAseq) > 1 & nrow(top.RES.RNAseq) > 1){
    ### assign samples to groups 
    data.SENS <- t(top.SENS.RNAseq)
    data.RES <- t(top.RES.RNAseq)
    SENS.RES.data <- cbind(data.SENS,data.RES)
    gset <- ExpressionSet(assayData = SENS.RES.data)
    seq0 <- seq(0, 0, length.out=nrow(top.SENS.RNAseq))
    seq1 <- seq(1, 1, length.out=nrow(top.RES.RNAseq))
    sml <- c(seq0, seq1)
    
    ### set up design matrix
    gs <- factor(sml)
    groups <- make.names(c("test","test2"))
    levels(gs) <- groups
    gset$group <- gs
    design <- model.matrix(~group + 0, gset)
    colnames(design) <- levels(gs)
    
    fit <- lmFit(gset, design)  # fit linear model
    
    ### set up contrasts of interest and recalculate model coefficients
    cts <- paste(groups[1], groups[2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    
    ### compute statistics and table of top significant genes
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(fit2))
    tT$Gene <- rownames(tT) # 19,111 genes
    tT.noblanks <- tT[tT$Gene!="",] # 19,111 genes
    colnames(tT.noblanks)[1] <- "Log2FC"
    weights <- stats::na.omit(tT.noblanks[tT.noblanks$Gene %in% 
                                     colnames(RNA.df),]) # 19,111 genes
    utils::write.csv(weights, paste0(drugs[i],"_full_gene_signature.csv"), 
                     row.names = FALSE)
    filtered.weights.q0.05 <- weights[weights$adj.P.Val<0.05,]
    if(nrow(filtered.weights.q0.05)>0){
      avg.filtered.weights.q0.05 <- ddply(filtered.weights.q0.05, .(Gene), summarize,
                                          Log2FC = mean(Log2FC, na.rm=TRUE),
                                          sd.Log2FC = sd(Log2FC, na.rm=TRUE)) # prevent duplicate gene names; 31 genes
      if(nrow(avg.filtered.weights.q0.05) < nrow(filtered.weights.q0.05)){
        filtered.weights <- avg.filtered.weights.q0.05 # columns are already in order of Gene, Log2FC, sd.Log2FC
      }else{
        filtered.weights <- filtered.weights.q0.05
        filtered.weights <- filtered.weights[,c("Gene","Log2FC")]
      }
      utils::write.csv(filtered.weights, file=paste0(drugs[i],"_filtered_gene_signature_no_duplicates.csv"), row.names = FALSE)
      
      #select gene weights using top n.genes based on abs(Log2FC)
      top.genes <- filtered.weights %>% slice_max(abs(Log2FC), n=n.genes)
      
      #perform WV with test set
      test.RNAseq <- RNA.df[!(RNA.df$CCLE_ID %in% train.RNAseq$CCLE_ID) & (RNA.df$CCLE_ID %in% type.list$CCLE_Name), c("CCLE_ID",top.genes$Gene)]
      WV.test <- WV(test.RNAseq, top.genes) #make so these save into 1 df and then can run landscape
      colnames(WV.test)[2] <- drugs[i]
      WV.results <- merge(WV.results, WV.test, all=TRUE)
    }
  }
}
utils::write.csv(WV.results, paste0("WV_matrix_CCLE_RNAseq_PRISM_drug_AUC_",Sys.Date(),".csv"), row.names = FALSE)

##### Step 3: perform drug-level landscape ##### 
drugs <- colnames(WV.results[ , 2:ncol(WV.results)])
corr.results <- list()
drugSEA.results <- list()
qc.min <- 1
for (i in seq_len(length(drugs))){
  # prep correlation input
  rank.list <- WV.results[ , c("CCLE_ID", drugs[i])]
  colnames(rank.list)[2] <- "Rank"
  input.df <- merge(rank.list, WV.results, by="CCLE_ID")
  
  # run correlations
  corr.AUC <- rank.corr(data = input.df, plots = F)
  corr.results[[drugs[i]]] <- corr.AUC$result
  drugSEA.input <- corr.AUC$result
  
  # remove self-correlations
  drugSEA.input <- drugSEA.input[drugSEA.input$Drug != 
                                   drugSEA.input$'Drug for Rank', ]
  
  # run drugSEA and set FDR = 0 to prevent production of mountain plots
  drugSEA.AUC <- DMEA::drugSEA(drugSEA.input, gmt, FDR = 0)
  drugSEA.result <- drugSEA.AUC$result
  
  ## check if expected MOAs are enriched per MOA annotations for given drug
  # identify MOA(s) for given drug
  self.MOAs <- drug.moa[drug.moa[ , c("Drug")] == drugs[i], c("moa")]
  self.MOAs <- strsplit(self.MOAs, sep)[[1]]
  
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
all.corr <- rbindlist(corr.results, use.names = TRUE, idcol = "Drug for Rank")
utils::write.csv(all.corr, "MOAL_correlations.csv")
all.drugSEA <- rbindlist(drugSEA.results, use.names = TRUE, idcol = "Drug")
all.drugSEA <- merge(all.drugSEA, drug.moa, by = "Drug")
utils::write.csv(all.drugSEA, "MOAL_drug_DMEA_results.csv")

# label drugs failing qc
drugs.failing.qc <- all.drugSEA[all.drugSEA$self & !all.drugSEA$qc_pass, ]
utils::write.csv(drugs.failing.qc, "MOAL_drugs_failing_qc.csv")

##### Step 3: perform MOA Landscape ##### 
# clear unneeded outputs for memory space
all.corr <- NULL
drugs.failing.qc <- NULL

# allow computer to cool down for 30 min
Sys.sleep(1800)

# perform MOA Landscape and store results
landscape.2d.results <- MOAL::landscape_2d(all.drugSEA, 
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

##### Step 4: evaluate performance ##### 
# load results for PRISM (all adherent)
setwd(my.path)
setwd("PRISM/all_adherent")
og.result <- utils::read.csv("MOAL_results.csv")

# compare WGV results to original PRISM results
performance.WGV <- MOAL::performance(landscape.2d.results,
                                     og.result,
                                     prediction.name = "WGV",
                                     og.name = "PRISM")

# store performance metrics
ggplot2::ggsave("PRISM_vs_WGV_scatter_plot.pdf",
                performance.WGV$scatter.plot)
utils::write.csv("PRISM_vs_WGV_results.csv",
                 performance.WGV$scatter.plot.df)
ggplot2::ggsave("PRISM_vs_WGV_venn_diagram.pdf",
                performance.WGV$venn.diagram)
utils::write.csv("PRISM_vs_WGV_performance.csv",
                 performance.WGV$performance.df)
ggplot2::ggsave("PRISM_vs_WGV_performance.pdf",
                performance.WGV$performance.table)

