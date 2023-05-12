\name{landscape}
\alias{landscape}
\title{Mechanism-of-action Landscape}
\usage{
landscape(drug.sensitivity, drug.moa, gmt = NULL, drug = "Drug",
          moa = "moa", sep = ", ", min.per.set = 5, p = 0.05,
          FDR = 0.25, qc.min = 1, rm.self.corr = TRUE,
          rm.qc.fail = FALSE, max.overlap = 20, n.top = 2,
          scatter.plots = FALSE, FDR.scatter.plots = 0.05,
          mtn.plots.1d = FALSE, mtn.plots.2d = TRUE,
          suppress.warnings = TRUE)
}
\description{
Uses an input drug sensitivity score data frame to evaluate enrichment
of drug mechanisms-of-action (moa) for each input moa.
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

  \item{p}{Fisher p-value cut-off for enrichments to be labeled as
  significant in volcano plot. Default: 0.05.}

  \item{FDR}{False discovery rate cut-off value to produce mountain plots
  for significant enrichments. Default: 0.25.}

  \item{qc.min}{Minimum NES value to meet quality control (qc) expectation
  for self moa enrichment. Drugs with NES < qc.min for their own moa(s) will
  be labeled as failing qc. Default: 1.}

  \item{rm.self.corr}{Logical operator describing whether or not to remove
  self correlations. Default: TRUE.}

  \item{rm.qc.fail}{Logical operator describing whether or not to remove
  input drugs which did not pass quality control based on their enrichment
  for their own moa(s). Default: FALSE.}

  \item{max.overlap}{Maximum percent overlap between moa pairs labeled on
  volcano plot. Default: 20.}

  \item{n.top}{Number of top enrichments to label on volcano plot.
  Default: 2.}

  \item{scatter.plots}{Logical operator describing whether or not to produce
  scatter plots for correlations passing FDR.scatter.plots threshold.
  Default: FALSE.}

  \item{FDR.scatter.plots}{False discovery rate cut-off value to produce
  scatter plots for correlations. Default: 0.05.}

  \item{mtn.plots.1d}{Logical operator describing whether or not to produce
  mountain plots for drug-level moa enrichments passing FDR threshold.
  Default: FALSE.}

  \item{mtn.plots.2d}{Logical operator describing whether or not to produce
  mountain plots for moa-level moa enrichments passing FDR threshold.
  Default: FALSE.}

  \item{suppress.warnings}{Logical operator describing whether or not to
  suppress warnings. Default: TRUE.}
}

\value{
  \item{results}{Data frame containing enrichment results for moa pairs.}

  \item{qc.results}{Data frame containing self-enrichment results for each moa.}

  \item{DMEA.results}{Data frame containing enrichment results for each moa.}

  \item{volcano.plot }{ggplot object for volcano plot summarizing enrichment
  of all moa evaluated.}

  \item{mtn.plots}{ggplot object list of mountain plots for drug sets with
  enrichment which meet the FDR threshold.}

  \item{drug.results}{Data frame containing enrichment results for each drug.}

  \item{corr.result}{Data frame with Pearson & Spearman correlation estimates,
  p-values, FDR q-values (Benjamini-Hochberg method) and linear regression
  parameters for each drug-drug pair.}

  \item{drugs.failing.qc}{Data frame with drug-level self moa enrichment
  results which did not meet the minimum NES value (qc.min) to pass quality
  control (qc).}

  \item{scatter.plots}{ggplot object with a scatter plot for each
  significant correlation on each page.}

  \item{drug.mtn.plots}{ggplot object list of mountain plots for each drug
  with enrichment which meets the FDR threshold.}

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
scatter.plots = FALSE, mtn.plots.1d = FALSE, and mtn.plots.2d - FALSE.

If your input drug.sensitivity data frame has ~1k drugs, it is
recommended to run the landscape_1d function followed by the
landscape_2d function instead of the landscape function and clear
your R environment after storing results from landscape_1d to avoid
running out of memory.
}

\seealso{
  \link[MOAL]{landscape_1d}
  \link[MOAL]{landscape_2d}
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

## Step 3: perform moa landscape and store results
landscape.test <- MOAL::landscape(drug.sensitivity, drug.moa)
}
