# MOAL simulation study
# Author: Belinda B. Garana (BG)
# Created: 2023-02-23; last edit: BG 2023-03-28
library(dplyr);library(data.table);library(DMEA);library(MOAL);library(ggplot2);
setwd("/Users/belindagarana/Documents/Graham\ Lab/Landscape/Simulation/")

##### Step 1: generate simulated drug AUC data #####
# create names for synthetic cell lines (25) and drugs (60)
drug.set.size <- 6
num.drug.sets <- 10
num.drugs <- num.drug.sets*drug.set.size
synthetic.cell.names <- paste0("Cell_No_", seq(from = 1, to = 25, by = 1))
synthetic.drug.names <- paste0("Drug_No_", seq(from = 1, to = num.drugs, by = 1))

# set drug AUC parameters based on PRISM distribution
frac.low.AUC <- 0.72
num.drugs.low <- round(frac.low.AUC*num.drugs,digits=0) # 72
num.drugs.high <- round((1-frac.low.AUC)*num.drugs,digits=0) # 28
num.synth.drugs.low <- round(frac.low.AUC*drug.set.size, digits=0) # 4
num.synth.drugs.high <- round((1-frac.low.AUC)*drug.set.size, digits=0) # 2
low.mean.AUC <- 0.83
low.sd.AUC <- 0.1657471
high.mean.AUC <- 1.31
high.sd.AUC <- 0.1580992
half.low.sd.AUC <- 0.08
half.high.sd.AUC <- 0.08

##### Step 2: perturb 2 synthetic drug sets #####
## set perturbation values
pert <- 0.2 # max perturbation value

# want drug set 1 to be shifted -0.4, -0.3, -0.2, -0.1, 0, +0.1, +0.2, +0.3, +0.4
# and for each of those cases, drug set 2 to be shifted -0.4, -0.3, -0.2, -0.1, 0, +0.1, +0.2, +0.3, +0.4
values.to.vary1 <- pert*c(rep(-1, times=9),
                          rep(-0.75, times=9),
                          rep(-0.5, times=9),
                          rep(-0.25, times=9),
                          rep(0, times=9),
                          rep(0.25, times=9),
                          rep(0.5, times=9),
                          rep(0.75, times=9),
                          rep(1, times=9))
values.to.vary2 <- pert*rep(c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), times=9)

# create synthetic drug set annotations
Drug <- synthetic.drug.names
drug.moa <- as.data.frame(Drug)
drug.moa$moa <- NA
for(i in 1:num.drug.sets){
  drug.moa[(drug.set.size*(i-1)+1):(drug.set.size*i),]$moa <- paste0("Synthetic drug set ",i)
}
synthetic.drug.set1 <- drug.moa[drug.moa$moa=="Synthetic drug set 1",]$Drug
synthetic.drug.set2 <- drug.moa[drug.moa$moa=="Synthetic drug set 2",]$Drug
gmt <- DMEA::as_gmt(drug.moa, min.per.set=(drug.set.size-1))

# use parallel computing if possible
library(parallel);library(snow);library(doSNOW);

cores <- parallel::detectCores() # number of cores available
if(cores[1] > 1){
  cl <- snow::makeCluster(cores[1]-1) # cluster using all but 1 core
  doSNOW::registerDoSNOW(cl) # register cluster
}

## run 25 replicates
# create list to store results
synth.result.list <- list()
n.reps <- 25
for(i in 1:n.reps){
  # create list to store results for each replicate
  all.synth.results <- list()

  # run analysis with various perturbation amounts for each drug set
  for(j in 1:length(values.to.vary1)){
    varied.value1 <- values.to.vary1[j]
    varied.value2 <- values.to.vary2[j]
    varied.value.name <- paste0("values_added_",varied.value1,"_",varied.value2)

    # determine range of perturbation across cell lines
    if(varied.value1 == 0){
      range.to.vary1 <- rep(0, times=length(synthetic.cell.names))
    }else{
      range.to.vary1 <- seq(from = -varied.value1, to = varied.value1, by = 2*varied.value1 / (length(synthetic.cell.names)-1)) #from0 BG 20210913
    }

    if(varied.value2 == 0){
      range.to.vary2 <- rep(0, times=length(synthetic.cell.names))
    }else{
      range.to.vary2 <- seq(from = -varied.value2, to = varied.value2, by = 2*varied.value2 / (length(synthetic.cell.names)-1)) #from0 BG 20210913
    }

    # create empty AUC data frame with drugs as row names and synthetic cell line names as column names
    AUC.df <- as.data.frame(Drug)
    AUC.df[,synthetic.cell.names] <- NA

    for (k in 1:length(synthetic.cell.names)){
      # fill AUC data frame with distribution matching PRISM for each cell line
      normal.low <- rnorm(num.drugs.low, mean = low.mean.AUC, sd = half.low.sd.AUC)
      normal.high<- rnorm(num.drugs.high, mean = high.mean.AUC, sd = half.high.sd.AUC)
      AUC.df[,synthetic.cell.names[k]] <- c(normal.low, normal.high)

      # perturb synthetic drug set 1
      normal.low1 <- rnorm(num.synth.drugs.low, mean = low.mean.AUC + range.to.vary1[k], sd = half.low.sd.AUC)
      normal.high1<- rnorm(num.synth.drugs.high, mean = high.mean.AUC + range.to.vary1[k], sd = half.high.sd.AUC)
      AUC.df[AUC.df$Drug %in% synthetic.drug.set1, synthetic.cell.names[k]] <- c(normal.low1, normal.high1)

      # perturb synthetic drug set 2
      normal.low2 <- rnorm(num.synth.drugs.low, mean = low.mean.AUC + range.to.vary2[k], sd = half.low.sd.AUC)
      normal.high2<- rnorm(num.synth.drugs.high, mean = high.mean.AUC + range.to.vary2[k], sd = half.high.sd.AUC)
      AUC.df[AUC.df$Drug %in% synthetic.drug.set2, synthetic.cell.names[k]] <- c(normal.low2, normal.high2)
    }

    # transpose AUC data frame so that cell line names are in column 1
    # and drug names are the rest of the column names
    AUC.df <- as.data.frame(t(AUC.df))
    colnames(AUC.df) <- AUC.df[1,]
    AUC.df$Drug <- rownames(AUC.df)
    AUC.df <- AUC.df[2:nrow(AUC.df), c("Drug",Drug)]
    colnames(AUC.df)[1] <- "CCLE_ID"
    AUC.df[, 2:ncol(AUC.df)] <- lapply(AUC.df[, 2:ncol(AUC.df)], as.numeric) # make sure data is numeric

    ##### Step 3: run MOA Landscape #####
    analysis <- MOAL::landscape(AUC.df, drug.moa = drug.moa, gmt = gmt, 
                                min.per.set = (drug.set.size - 1),
                                rm.qc.fail = FALSE, mtn.plots.2d = FALSE)
    synth.result <- analysis$results[analysis$results$pair == "Synthetic drug set 1 & Synthetic drug set 2",]
    synth.result[,c("varied_value1", "varied_value2")] <- c(varied.value1, varied.value2)
    all.synth.results[[varied.value.name]] <- synth.result
  }
  synth.result.list[[i]] <- data.table::rbindlist(all.synth.results, use.names = TRUE, idcol = "varied_value_name")
  Sys.sleep(1800) # let the computer rest for 30 min
}
synth.results <- data.table::rbindlist(synth.result.list, use.names = TRUE, idcol = "Replicate")
write.csv(synth.results, 
          file = paste0(n.reps,"_replicates_synthetic_drug_set_enrichments_",Sys.Date(),".csv"), 
          row.names = FALSE)

if(cores[1] > 1){
  snow::stopCluster(cl) # stop cluster
  rm(cl)
}

##### Step 4: compile results for tile plots #####
# compile NES_avg, percent_qc_pass across all replicates for each varied value combination
plot.df <- as.data.frame(values.to.vary1)
plot.df$values.to.vary2 <- values.to.vary2
plot.df$varied_value_name <- paste0("values_added_",plot.df$values.to.vary1,"_",plot.df$values.to.vary2)
plot.df[,c("NES_avg","score_avg","percent_qc_pass")] <- NA
for(i in 1:nrow(plot.df)){
  varied.values <- plot.df$varied_value_name[i]
  temp.result <- synth.results[synth.results$varied_value_name==varied.values,]
  plot.df$NES_avg[i] <- mean(temp.result$NES_avg, na.rm=TRUE)
  plot.df$score_avg[i] <- mean(temp.result$score_avg, na.rm=TRUE)
  plot.df$percent_qc_pass[i] <- 100*nrow(temp.result[temp.result$qc_pass,])/nrow(temp.result)
}
write.csv(plot.df, 
          file = paste0(n.reps,"_replicates_synthetic_drug_set_enrichment_plot_data_",Sys.Date(),".csv"), 
          row.names = FALSE)

## produce tile plots
# load themes for plots
ng.theme <- ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                           panel.border = element_rect(fill=NA), panel.background = element_blank(),
                           axis.line = element_line(colour = "black"), axis.text.x = element_text(colour = "black"),
                           axis.text.y = element_text(colour = "black"), axis.ticks.x = element_line(colour="black"),
                           axis.ticks.y = element_line(colour="black"), legend.title = element_blank(),
                           axis.title.y = element_text(size=8, colour="black"))

bg.theme <- ggplot2::theme(legend.background = element_rect(), legend.position="top", 
                           legend.text = element_text(size=14), legend.key = element_blank(),
                           axis.title.x = element_text(size=20), axis.text.x  = element_text(size=16),
                           axis.title.y = element_text(size=20), axis.text.y  = element_text(size=16),
                           plot.title = element_text(lineheight=.8, face="bold", size=36))

# plot NES_avg
NES.tile.plot <- ggplot(plot.df, aes(x = values.to.vary1, y = values.to.vary2, fill=NES_avg)) +
  geom_tile() + ng.theme + theme_light() + scale_fill_gradient2() +
  labs(x="Value Added to Synthetic Drug Set 1", y = "Value Added to Synthetic Drug Set 2", fill = "Mean NES")
ggsave(NES.tile.plot, file = paste0("MOAL_simulation_NES_",n.reps,"_replicates_",Sys.Date(),".pdf"))

# plot percent_qc_pass
qc.tile.plot <- ggplot(plot.df, aes(x = values.to.vary1, y = values.to.vary2, fill=percent_qc_pass)) +
  geom_tile() + ng.theme + theme_light() + scale_fill_gradient2() +
  labs(x="Value Added to Synthetic Drug Set 1", y = "Value Added to Synthetic Drug Set 2", fill = "% Passing QC")
ggsave(qc.tile.plot, file = paste0("MOAL_simulation_qc_",n.reps,"_replicates_",Sys.Date(),".pdf"))
