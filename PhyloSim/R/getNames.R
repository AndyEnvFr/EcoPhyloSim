#' Generate simulation names based on model parameters
#'
#' @title Generate simulation names for simulation objects
#' @param runs Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Character vector of names (length 1 for PhyloSim, or 1 per element for PhylosimList)
#' @description Constructs standardized names based on model parameters for consistent labeling.
#' 
#' See details for abbreviation rules.
#'
#' @details
#' \itemize{
#'   \item \code{ddX}: Density dependence with competition strength X
#'   \item \code{dispX}: Dispersal type, X is G or numeric
#'   \item \code{srX}: Speciation rate
#'   \item \code{eX}: Environment dependence with strength X
#'   \item \code{fbmrX}: Fitness-based mortality ratio
#'   \item \code{dcX}: Density cutoff
#'   \item \code{faoX}: Fitness acts on (M = mortality, R = reproduction, B = both)
#'   \item \code{fiX}, \code{rqX}, \code{rqsX}, \code{pX}: Optional parameters if > 0
#' }
#'
#' @example /inst/examples/getConNeigh-help.R
#'
#' @seealso \code{\link{PhyloSim}}
#' @export
getNames <- function(runs) {
  UseMethod("getNames")
}

#' @rdname getNames
#' @method getNames PhyloSim
#' @export
getNames.PhyloSim <- function(runs) {
  m <- runs$Model
  name <- paste0(
    if (isTRUE(m$density)) paste0("dd", m$compStrength, "_"),
    "disp", ifelse(m$dispersal == "global", "G_", paste0(m$dispersal, "_")),
    "sr", m$specRate, "_",
    if (isTRUE(m$environment)) paste0("e", m$envStrength, "_"),
    "fbmr", m$fitnessBaseMortalityRatio, "_",
    "dc", m$densityCut, "_",
    "fao", switch(m$fitnessActsOn, mortality = "M", reproduction = "R", both = "B"),
    if (isTRUE(m[["fission"]] > 0)) paste0("_fi", m$fission),
    if (isTRUE(m[["redQueen"]] > 0)) paste0("_rq", m$redQueen),
    if (isTRUE(m[["redQueenStrength"]] > 0)) paste0("_rqs", m$redQueenStrength),
    if (isTRUE(m[["protracted"]] > 0)) paste0("_p", m$protracted)
  )
  return(name)
}

#' @rdname getNames
#' @method getNames PhylosimList
#' @export
getNames.PhylosimList <- function(runs) {
  sapply(runs, getNames)
}


