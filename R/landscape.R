landscape <- function(drug.sensitivity, drug.moa, gmt = NULL, drug = "Drug",
                      moa = "moa", sep = ", ", min.per.set = 5, p = 0.05,
                      FDR = 0.25, qc.min = 1, rm.self.corr = TRUE,
                      rm.qc.fail = FALSE, max.overlap = 20, n.top = 2,
                      scatter.plots = FALSE, FDR.scatter.plots = 0.05,
                      mtn.plots.1d = FALSE, mtn.plots.2d = TRUE,
                      suppress.warnings = TRUE) {
  # perform drug (1d) landscape
  landscape_1d.results <- landscape_1d(
    drug.sensitivity, drug.moa, gmt, drug, moa, sep, min.per.set, FDR,
    qc.min, rm.self.corr, scatter.plots, FDR.scatter.plots, mtn.plots.1d,
    suppress.warnings
  )

  # perform MOA (2d) landscape
  landscape_2d.results <- landscape_2d(
    landscape_1d.results$DMEA.results, drug.moa, gmt, drug, moa, sep,
    min.per.set, p, FDR, rm.qc.fail, max.overlap, n.top, mtn.plots.2d,
    suppress.warnings
  )

  return(list(
    results = landscape_2d.results$results,
    qc.results = landscape_2d.results$qc.results,
    DMEA.results = landscape_2d.results$DMEA.results,
    volcano.plot = landscape_2d.results$volcano.plot,
    mtn.plots = landscape_2d.results$mtn.plots,
    drug.results = landscape_1d.results$DMEA.results,
    corr.results = landscape_1d.results$corr.results,
    drugs.failing.qc = landscape_1d.results$drugs.failing.qc,
    scatter.plots = landscape_1d.results$scatter.plots,
    drug.mtn.plots = landscape_1d.results$mtn.plots,
    gmt = landscape_1d.results$gmt
  ))
}
