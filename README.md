# Mechanism-of-action Landscape


## Introduction

Resistance is a major barrier to cancer remission. However, it is not
practical to test all the millions of potential drug combinations due to
cost limitations. To help prioritize drug pairings for further
evaluation, we studied the similarity of drugs grouped by
mechanisms-of-action (MOA) in published drug screens. With this process,
termed the Mechanism-of-action Landscape (MOAL), 131 drug MOA pairings
were prioritized from more than 1.8 million potential drug pairings.
Though this method does not evaluate synergy, it has been used to study
more than a thousand drugs whereas other methods have only studied on
the order of ten to one hundred drugs. In addition to recovering
established drug combinations (e.g., MEK and RAF inhibitors in skin
cancer), this method has also detected potential drug pairings which
have not been well studied and demonstrated that results can vary
depending on tissue type. Furthermore, we show that the MOAL method is
beneficial over a simple UMAP visualization because it uses statistics
to eliminate more than 96% of potential drug MOA pairings from
consideration and identifies significant MOA pairings which may not be
obvious from their distance on a UMAP visualization (e.g., BET and HDAC
inhibitors). This method can also be used to evaluate the global
performance of drug sensitivity prediction algorithms. For more
information, please visit our website:
<https://belindabgarana.github.io/MOAL>

## Installation

To install this package, start R and enter:

``` r
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("MOAL")
```

Alternatively, to install this package directly from GitHub, start R and
enter:

``` r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("BelindaBGarana/MOAL")
```

If you are using Windows OS, you may need to change the code above to:

``` r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

devtools::install_github("BelindaBGarana/MOAL", build = FALSE)
```

## Input Option 1: Drug Screen Data Frame (\<1k Drugs)

With the landscape function, you can rank drug moa’s based on similarity
to each other. For example, if your input data frame has sample names
(e.g., CCLE IDs) in the first column and drug sensitivity scores (e.g.,
AUC) for various drugs in the rest of the columns, then moa pairs with
NES_avg \> 0 may be toxic to the same samples and pairs with NES_avg \<
0 may be toxic to opposite samples. In other words, if the moa pair “moa
X & moa Y” passes quality control (qc) and NES_avg \< 0, then drugs in
moa X may be toxic to samples resistant to Y and vice-versa. The only
requirements are that: \* At least 2 drug sets are represented in your
drug sensitivity data frame \* At least 3 samples (e.g., cell lines) are
in the drug sensitivity data frame

Required inputs: \* drug.sensitivity: Data frame containing drug
sensitivity scores (e.g., AUC) with drug names as column names, except
for the first column which contain sample names (e.g., CCLE IDs) \*
drug.moa: Dataframe with moa annotations for each drug.

Example: 1. Prepare drug.sensitivity data frame

``` r
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
```

2.  Prepare drug.moa data frame with moa annotations for each drug

``` r
# Drug is our default name for the column containing drug names
drug.moa <- as.data.frame(Drug)

# moa is our default name for the column containing drug moa set annotations
drug.moa$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 2)]), 6)
```

3.  Perform moa landscape and store results

``` r
landscape.test <- MOAL::landscape(drug.sensitivity, drug.moa)
```

    ## Generating gmt object for enrichment analysis...

    ## Loading required namespace: parallel

    ## Loading required namespace: snow

    ## Loading required namespace: doSNOW

    ## Evaluating drug 1 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 2 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 3 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 4 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 5 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 6 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 7 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 8 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 9 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 10 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 11 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 12 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Generating gmt object for enrichment analysis...

    ## Evaluating moa 1 out of 2

    ## Running enrichment analysis...

    ## Evaluating moa 2 out of 2

    ## Running enrichment analysis...

    ## Compiling results across moa pairs...

## Input Option 2: Drug Screen Data Frame (\>=1k Drugs)

If you are analyzing a large drug screen (\>= 1k drugs), you may run out
of memory if you perform the full landscape analysis all at once. To
prevent this issue, it is recommended to perform each of the 2 parts of
the landscape analysis separately and clear your R environment after
running the first part of the analysis (i.e., landscape_1d). Again, the
only requirements are that: \* At least 2 drug sets are represented in
your drug sensitivity data frame \* At least 3 samples (e.g., cell
lines) are in the drug sensitivity data frame

Required inputs: \* drug.sensitivity: Data frame containing drug
sensitivity scores (e.g., AUC) with drug names as column names, except
for the first column which contain sample names (e.g., CCLE IDs) \*
drug.moa: Dataframe with moa annotations for each drug.

Example: 1. Prepare drug.sensitivity data frame

``` r
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
```

2.  Prepare drug.moa data frame with moa annotations for each drug

``` r
# Drug is our default name for the column containing drug names
drug.moa <- as.data.frame(Drug)

# moa is our default name for the column containing drug moa set annotations
drug.moa$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 2)]), 6)
```

3.  Perform drug-level landscape, store results, and clear unnecessary
    data before continuing

``` r
landscape.1d.test <- MOAL::landscape_1d(drug.sensitivity, drug.moa)
```

    ## Generating gmt object for enrichment analysis...

    ## Evaluating drug 1 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 2 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 3 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 4 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 5 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 6 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 7 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 8 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 9 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 10 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 11 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 12 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

``` r
# after storing results, clear unnecessary data to free memory
landscape.1d.test <- landscape.1d.test$DMEA.results
```

4.  Perform final step of MOA Landscape and store results

``` r
landscape.2d.test <- MOAL::landscape_2d(landscape.1d.test)
```

    ## Generating gmt object for enrichment analysis...

    ## Evaluating moa 1 out of 2

    ## Running enrichment analysis...

    ## Evaluating moa 2 out of 2

    ## Running enrichment analysis...

    ## Compiling results across moa pairs...

## Input Option 3: Prediction Results

The MOA Landscape can also be used to evaluate the network-level
performance of drug sensitivity prediction algorithms. If you use one of
these prediction algorithms to generate a predicted drug.sensitivity
data frame, you can use that data frame as an input to the MOA Landscape
(i.e., using one of the input options detailed above) and compare the
landscape results to that of the original drug screen used to generate
the prediction. Again, the only requirements are that: \* At least 2
drug sets are represented in your drug sensitivity data frame \* At
least 3 samples (e.g., cell lines) are in the drug sensitivity data
frame

Required inputs: \* prediction.result: Dataframe containing enrichment
results for mechanism-of-action (moa) pairs (i.e., results output by
landscape or landscape_2d functions) based on predicted drug sensitivity
scores. \* og.result: Dataframe containing enrichment results for
mechanism-of-action (moa) pairs (i.e., results output by landscape or
landscape_2d functions) based on original drug sensitivity scores (i.e.,
from a real drug screen).

Example: 1. Prepare drug.moa data frame with moa annotations for each
drug

``` r
# create list of drug names
Drug <- paste0("Drug_", seq(from = 1, to = 12))

# Drug is our default name for the column containing drug names
drug.moa <- as.data.frame(Drug)

# moa is our default name for the column containing drug moa set annotations
drug.moa$moa <- rep(paste("Set", LETTERS[seq(from = 1, to = 2)]), 6)
```

2.  Perform MOA landscape for original data set

``` r
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
```

    ## Generating gmt object for enrichment analysis...

    ## Evaluating drug 1 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 2 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 3 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 4 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 5 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 6 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 7 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 8 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 9 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 10 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 11 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 12 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Generating gmt object for enrichment analysis...

    ## Evaluating moa 1 out of 2

    ## Running enrichment analysis...

    ## Evaluating moa 2 out of 2

    ## Running enrichment analysis...

    ## Compiling results across moa pairs...

3.  Perform MOA landscape for prediction data set

``` r
## Step 3a: prepare drug sensitivity data frame to represent prediction data
prediction.drug.sensitivity <- as.data.frame(Sample_ID)

# give each drug values representative of AUC sensitivity scores
for(i in seq_len(length(Drug))){
  prediction.drug.sensitivity[,c(Drug[i])] <- rnorm(length(Sample_ID),
                                        mean = 0.83, sd = 0.166)
}

## Step 3b: perform moa landscape and store results
prediction.landscape <- MOAL::landscape(prediction.drug.sensitivity, drug.moa)
```

    ## Generating gmt object for enrichment analysis...

    ## Evaluating drug 1 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 2 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 3 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 4 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 5 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 6 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 7 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 8 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 9 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 10 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 11 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Evaluating drug 12 out of 12

    ## Running correlations and regressions...

    ## Running enrichment analysis...

    ## Generating gmt object for enrichment analysis...

    ## Evaluating moa 1 out of 2

    ## Running enrichment analysis...

    ## Evaluating moa 2 out of 2

    ## Running enrichment analysis...

    ## Compiling results across moa pairs...

4.  Compare MOA landscape results from prediction to original

``` r
performance.test <- MOAL::performance(prediction.landscape$results,
                                og.landscape$results)
```

## Session Info

``` r
sessionInfo()
```

    ## R version 4.2.1 (2022-06-23)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Big Sur ... 10.16
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Biobase_2.56.0      pkgload_1.3.2       splines_4.2.1      
    ##  [4] foreach_1.5.2       qqconf_1.3.0        brio_1.1.3         
    ##  [7] sn_2.1.0            Rdpack_2.4          stats4_4.2.1       
    ## [10] yulab.utils_0.0.6   TFisher_0.2.0       yaml_2.3.7         
    ## [13] ggrepel_0.9.3       numDeriv_2016.8-1.1 pillar_1.8.1       
    ## [16] lattice_0.20-45     MOAL_0.99.0         glue_1.6.2         
    ## [19] doSNOW_1.0.20       digest_0.6.31       rbibutils_2.2.10   
    ## [22] qvalue_2.28.0       colorspace_2.1-0    sandwich_3.0-2     
    ## [25] ggfun_0.0.9         cowplot_1.1.1       htmltools_0.5.4    
    ## [28] Matrix_1.5-3        plyr_1.8.8          pkgconfig_2.0.3    
    ## [31] mvtnorm_1.1-3       snow_0.4-4          patchwork_1.1.2    
    ## [34] scales_1.2.1        ggplotify_0.1.0     DMEA_0.99.0        
    ## [37] tibble_3.2.1        generics_0.1.3      ggplot2_3.4.1      
    ## [40] sjlabelled_1.2.0    withr_2.5.0         TH.data_1.1-1      
    ## [43] BiocGenerics_0.42.0 cli_3.6.0           mnormt_2.1.1       
    ## [46] survival_3.4-0      magrittr_2.0.3      evaluate_0.20      
    ## [49] fansi_1.0.4         MASS_7.3-58.1       tools_4.2.1        
    ## [52] ggvenn_0.1.9        data.table_1.14.6   lifecycle_1.0.3    
    ## [55] multcomp_1.4-20     mutoss_0.1-12       stringr_1.5.0      
    ## [58] aplot_0.1.10        munsell_0.5.0       plotrix_3.8-2      
    ## [61] compiler_4.2.1      gridGraphics_0.5-1  rlang_1.1.0        
    ## [64] grid_4.2.1          iterators_1.0.14    rstudioapi_0.14    
    ## [67] labeling_0.4.2      rmarkdown_2.20      testthat_3.1.7     
    ## [70] gtable_0.3.2        codetools_0.2-18    multtest_2.52.0    
    ## [73] reshape2_1.4.4      sjmisc_2.8.9        R6_2.5.1           
    ## [76] gridExtra_2.3       zoo_1.8-11          knitr_1.42         
    ## [79] dplyr_1.1.0         fastmap_1.1.0       utf8_1.2.3         
    ## [82] rprojroot_2.0.3     mathjaxr_1.6-0      insight_0.19.1     
    ## [85] desc_1.4.2          metap_1.8           stringi_1.7.12     
    ## [88] parallel_4.2.1      Rcpp_1.0.10         vctrs_0.6.0        
    ## [91] tidyselect_1.2.0    xfun_0.38
