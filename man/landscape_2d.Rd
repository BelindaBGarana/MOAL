\name{landscape_2d}
\alias{landscape_2d}
\title{2D Landscape}
\usage{
landscape_2d(input, drug.moa = NULL, gmt = NULL, drug = "Drug",
            moa = "moa", sep = ", ", min.per.set = 5,
            p = 0.05, FDR = 0.25, rm.qc.fail = FALSE,
            max.overlap = 20, n.top = 2, mtn.plots = FALSE,
            suppress.warnings = TRUE)
}
\description{
Evaluates enrichment of drug mechanisms-of-action (moa) for each input moa.
}

\arguments{
  \item{input}{Data frame containing normalized enrichment scores (i.e., NES)
  of moa for various drugs (i.e., DMEA.results output from landscape_1d
  function. Required columns: drug, "Drug_set", and "NES".}

  \item{drug.moa}{Data frame with moa annotations for each drug. If NULL, the
  input data frame must contain a moa column. Default: NULL.}

  \item{gmt}{Optional: gmt object containing drug sets (e.g., output from
  GSA.read.gmt function from GSA package). If NULL, moa annotations will be
  gathered from drug.moa or the input data frame. Default: NULL.}

    \item{drug}{Column name in drug.moa containing the drug names. Default:
  "Drug".}

  \item{moa}{Column name in drug.moa containing the drug moa annotations.
  Default: "moa".}

  \item{sep}{Character(s) to use for splitting drug set names. Default:
  ", ".}

  \item{min.per.set}{The minimum number of drugs required per set. Default:
  5.}

  \item{p}{Fisher p-value cut-off for enrichments to be labeled as
  significant in volcano plot. Default: 0.05.}

  \item{FDR}{False discovery rate cut-off value to produce mountain plots
  for significant enrichments. Default: 0.25.}

  \item{rm.qc.fail}{Logical operator describing whether or not to remove
  input drugs which did not pass quality control based on their enrichment
  for their own moa(s). Default: FALSE.}

  \item{max.overlap}{Maximum percent overlap between moa pairs labeled on
  volcano plot. Default: 20.}

  \item{n.top}{Number of top enrichments to label on volcano plot.
  Default: 2.}

  \item{mtn.plots}{Logical operator describing whether or not to produce
  mountain plots for enrichments passing FDR threshold. Default: FALSE.}

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

  \item{gmt}{gmt object containing list of elements (genesets), set names
  (geneset.names), and set descriptions (geneset.descriptions).}
}

\references{
https://belindabgarana.github.io/MOAL
}

\author{Belinda B. Garana
}

\note{
To minimize computation time, do not change default input parameter
mtn.plots = FALSE.
}

\seealso{
  \link[MOAL]{landscape_1d}
}

\examples{
## Step 1: prepare drug.sensitivity data frame for 1D landscape
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

## Step 3: perform drug-level landscape and store results
landscape.1d.test <- MOAL::landscape_1d(drug.sensitivity, drug.moa)

## Step 4: perform moa landscape and store results
landscape.2d.test <- MOAL::landscape_2d(landscape.1d.test$DMEA.results)
}
