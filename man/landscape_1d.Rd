\name{landscape_1d}
\alias{landscape_1d}
\title{1D Landscape}
\usage{
landscape_1d(drug.sensitivity, drug.moa, gmt = NULL, drug = "Drug",
            moa = "moa", sep = ", ", min.per.set = 5,
            FDR = 0.25, qc.min = 1, rm.self.corr = TRUE,
            scatter.plots = FALSE, FDR.scatter.plots = 0.05,
            mtn.plots = FALSE, suppress.warnings = TRUE)
}
\description{
Evaluates correlations between drugs and enrichment of mechanisms-of-action
(moa) for each input drug.
}

\arguments{
  \item{drug.sensitivity}{Data frame containing drug sensitivity scores (e.g.,
  AUC) with drug names as column names, except for the first column which
  contain sample names (e.g., CCLE IDs).}

  \item{drug.moa}{Data frame with moa annotations for each drug.}

  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from
  GSA.read.gmt function from GSA package). If NULL, moa annotations will be
  gathered from drug.moa. Default: NULL.}

  \item{drug}{Column name in drug.moa containing the drug names. Default:
  "Drug".}

  \item{moa}{Column name in drug.moa containing the drug moa annotations.
  Default: "moa".}

  \item{sep}{Character(s) to use for splitting drug moa set names. Default:
  ", ".}

  \item{min.per.set}{The minimum number of drugs required per set. Default:
  5.}

  \item{FDR}{False discovery rate cut-off value to produce mountain plots
  for moa enrichments. Default: 0.25.}

  \item{qc.min}{Minimum NES value to meet quality control (qc) expectation
  for self moa enrichment. Drugs with NES < qc.min for their own moa(s) will
  be labeled as failing qc. Default: 1.}

  \item{rm.self.corr}{Logical operator describing whether or not to remove
  self correlations. Default: TRUE.}

  \item{scatter.plots}{Logical operator describing whether or not to produce
  scatter plots for correlations passing FDR.scatter.plots threshold.
  Default: FALSE.}

  \item{FDR.scatter.plots}{False discovery rate cut-off value to produce
  scatter plots for correlations. Default: 0.05.}

  \item{mtn.plots}{Logical operator describing whether or not to produce
  mountain plots for enrichments passing FDR threshold. Default: FALSE.}

  \item{suppress.warnings}{Logical operator describing whether or not to
  suppress warnings. Default: TRUE.}
}

\value{
  \item{DMEA.results}{Data frame containing results for enrichment of each drug
  moa set for each drug.}

  \item{corr.results}{Data frame with Pearson & Spearman correlation estimates,
  p-values, FDR q-values (Benjamini-Hochberg method) and linear regression
  parameters for each drug-drug pair.}

  \item{drugs.failing.qc}{Data frame with self moa enrichment results which did
  not meet the minimum NES value (qc.min) to pass quality control.}

  \item{scatter.plots}{ggplot object with a scatter plot for each significant
  correlation.}

  \item{mtn.plots}{ggplot object list of mountain plots for drug moa sets with
  enrichment which meet the FDR threshold.}

  \item{gmt}{gmt object containing list of elements (genesets), set names
  (geneset.names), and set descriptions (geneset.descriptions).}
}

\references{
https://belindabgarana.github.io/MOAL
}

\author{Belinda B. Garana
}

\note{
To minimize computation time, do not change default input parameters
scatter.plots = FALSE and mtn.plots = FALSE.

If your input drug.sensitivity data frame has ~1k drugs, it is
recommended to clear your R environment after storing results from
landscape_1d to avoid running out of memory.

If rm.self.corr == TRUE, your moa sets of interest must contain 1
more drug than the min.per.set (e.g., 6 drugs if min.per.set = 5).
}

\examples{
## Step 1: prepare drug.sensitivity data frame
# create list of sample names
Sample_ID <- seq(from = 1, to = 21)

drug.sensitivity <- as.data.frame(Sample_ID)

# create list of drug names
Drug <- paste0("Drug_", seq(from = 1, to = 12))

# give each drug values representative of AUC sensitivity scores
for(i in seq_len(length(Drug))){
  drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 2: prepare drug.moa data frame with moa annotations for each drug
# Drug is our default name for the column containing drug names
drug.moa <- as.data.frame(Drug)

# moa is our default name for the column containing drug moa set annotations
drug.moa$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 2)]), 6)

## Step 3: perform 1D landscape and store results
landscape.1d.test <- MOAL::landscape_1d(drug.sensitivity, drug.moa)
}
