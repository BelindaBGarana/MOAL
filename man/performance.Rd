\name{performance}
\alias{performance}
\title{Network-level Performance}
\usage{
performance(prediction.result, og.result,
            prediction.name = "Prediction", og.name = "Original",
            corr.metric = "NES_avg", max.overlap = 20,
            n.digits = 2)
}
\description{
Evaluates network-level performance of drug sensitivity prediction algorithms
based on Mechanism of Action Landscape results.}

\arguments{
  \item{prediction.result}{Dataframe containing enrichment results for
  mechanism of action (moa) pairs (i.e., results output by landscape or
  landscape_2d functions) based on predicted drug sensitivity scores.
  Required columns: 'pair', 'self', 'qc_pass', 'percent_overlap', and
  corr.metric.}

  \item{og.result}{Dataframe containing enrichment results for
  mechanism of action (moa) pairs (i.e., results output by landscape or
  landscape_2d functions) based on original drug sensitivity scores (i.e.,
  from a real drug screen). Required columns: 'pair', 'self', 'qc_pass',
  'percent_overlap', and corr.metric.}

  \item{prediction.name}{Character string to use as a name for the
  prediction.result. Default: "Prediction".}

  \item{og.name}{Character string to use as a name for the og.result.
  Default: "Original".}

  \item{corr.metric}{Name of column containing metric for correlation.
  Default: "NES_avg".}

  \item{max.overlap}{Maximum percent overlap between moa pairs labeled
  on scatter plot. Default: 20.}

  \item{n.digits}{Numeric value to control the number of digits to use
  when rounding performance metrics for performance.table output. Default: 2.}
}

\value{
  \item{scatter.plot}{ggplot object with a scatter plot and Pearson
  correlation estimate of prediction.result with og.result.}

  \item{scatter.plot.df}{Data frame of overlap of moa pairs evaluated in both
  the prediction.result and the og.result with std used for scatter.plot.}

  \item{contingency.table}{tableGrob object containing the number of
  moa pairs found to be significant or insignificant based on measured drug
  sensitivity values (i.e., in the og.result) compared to the predicted values
  (e.g., the prediciton.result).}

  \item{performance.df}{Data frame containing success rate, false positive
  rate, and false negative rate of prediction.result compared to og.result.}

  \item{performance.table}{tableGrob object containing performance metrics
  in performance.df output rounded to 3 significant digits.}
}

\references{
https://belindabgarana.github.io/MOAL
}

\section{Warning}{The input parameter n.digits must be numeric in order to
get a performance.table output.}

\note{contingency.table and performance.table outputs can each be saved using
the ggsave function from the R package ggplot2.
}

\author{Belinda B. Garana
}

\seealso{
  \link[MOAL]{landscape_1d}
  \link[MOAL]{landscape_2d}
  \link[MOAL]{landscape}
}

\examples{
### Step 1: prepare drug.moa data frame with moa annotations for each drug
# create list of drug names
Drug <- paste0("Drug_", seq(from = 1, to = 12))

# Drug is our default name for the column containing drug names
drug.moa <- as.data.frame(Drug)

# moa is our default name for the column containing drug moa set annotations
drug.moa$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 2)]), 6)

### Step 2: perform moa landscape for original data set
## Step 2a: prepare drug sensitivity data frame to represent original data
# create list of sample names
Sample_ID <- seq(from = 1, to = 21)

og.drug.sensitivity <- as.data.frame(Sample_ID)

# give each drug values representative of AUC sensitivity scores
for(i in seq_len(length(Drug))){
  og.drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 2b: perform moa landscape and store results
og.landscape <- MOAL::landscape(og.drug.sensitivity, drug.moa)

### Step 3: perform moa landscape for prediction data set
## Step 3a: prepare drug sensitivity data frame to represent prediction data
prediction.drug.sensitivity <- as.data.frame(Sample_ID)

# give each drug values representative of AUC sensitivity scores
for(i in seq_len(length(Drug))){
  prediction.drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 3b: perform moa landscape and store results
prediction.landscape <- MOAL::landscape(prediction.drug.sensitivity, drug.moa)

### Step 4: compare moa landscape results from prediction to original
performance.test <- MOAL::performance(prediction.landscape$results,
                                og.landscape$results)
}
