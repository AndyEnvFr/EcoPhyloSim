# PhyloSim 0.3.1

## New features

* New pipeline for large simulation experiments: `runPhylosimAnalysis()` runs many simulations with fixed and prior-drawn parameters and saves them in chunks; `chunkwiseCalculateSummaries()`, `chunkwiseCleanup()`, `combineCunks()`, `prepareForAnalysis()` and `extractParams()` process the results (#85).
* New functions `rejectionSample()` and `rejectionSamplePhylosimAnalysis()` for ABC rejection sampling on simulation outputs (#85).
* New functions `getSpecies()` and `getSpeciesTraits()` extract the extant species and their traits from a simulation (#85).
* `calculateSummaryStatistics()` returns additional statistics (Clark-Evans index, 8-neighbor index, environmental response coefficients) (#85).

## Changes

* `calculateSummaryStatistics()` returns `NA` if any statistic cannot be calculated; set `strict = FALSE` to get the previous behavior (#85).
* `createCompletePar()` and `runSimulation()` stop with an error if `environment` or the values of `airmat` are outside [0, 1], or if `fitnessBaseMortalityRatio < 1` (#85).
* `createCompletePar()` stops with an error if `fitnessActsOn` is `"reproduction"` or `"both"`, because these settings crashed R (#98, #101).
* `runSimulationBatch()`: errors in single simulations are no longer silently turned into `NULL` during parallel runs, and `backup = TRUE` saves each result as `<scenario>.rds` (#106).
* Default plot sizes in `sac()` are log-spaced, as documented (#115).
* kulife, apTreeshape and PhyloMeasures were removed from Imports because they are archived on CRAN, but `calculateSummaryStatistics()` still requires apTreeshape and PhyloMeasures (#89); sf and spatstat were added to Imports, rmarkdown to Suggests.
* License changed from GPL (>= 2) to GPL (>= 3).

## Bug fixes

* Speciation by fission (`fission = 1` or `2`) assigned the wrong individuals to the new species (#85).
* `runSimulationBatch()`: `strip = "summaries"` was ignored without parallelization, and the results were not named after their scenarios (#106).
* Installation from GitHub failed because `src/RcppExports.cpp` was missing (#102).

## Documentation

* The vignette was converted from Sweave (PDF) to R Markdown (HTML).
* Package help converted to roxygen; examples moved to `inst/examples/` (#87).
* The help of `createCompletePar()` explains vectors in `runs` and the valid ranges of `environment`, `airmat` and `fitnessBaseMortalityRatio` (#85).
