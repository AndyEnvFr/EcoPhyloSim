#' Generate simulation names for a PhylosimList
#'
#' @title Generate simulation names for a PhylosimList
#' @param runs Object of class \code{PhylosimList}
#' @return A character vector of names, one for each simulation in the list
#' @description Constructs standardized names for a list of simulations based on model parameters.
#' These names are useful for labeling and comparing simulation outputs systematically.
#' @details
#' The following abbreviations and encodings are used in the generated names:
#'
#' \itemize{
#'   \item \code{ddT}: Density dependence is enabled.
#'   \item \code{compStrength}: Competition strength (numeric, appended after ddT).
#'   \item \code{_dispX}: Dispersal type; \code{X} is \code{G} for global or the specific dispersal value.
#'   \item \code{_sr}: Speciation rate.
#'   \item \code{_eT}: Environment dependence is enabled.
#'   \item \code{envStrength}: Environmental effect strength (numeric, appended after _eT).
#'   \item \code{_fbmr}: Fitness-based mortality ratio.
#'   \item \code{_dc}: Density cutoff.
#'   \item \code{_fao}: Fitness acts on (e.g., birth, death).
#'   \item \code{_fi}: Fission probability or intensity (only shown if > 0).
#'   \item \code{_rq}: Red Queen effect parameter (only shown if > 0).
#'   \item \code{_rqs}: Red Queen strength (only shown if > 0).
#'   \item \code{_p}: Protracted speciation delay (only shown if > 0).
#' }
#'
#' Example name:
#' \code{ddT0.5_dispG_sr2_eT0.8_fbmr10_dc4_fao1_fi0.1_rq0.3_rqs0.2_p3}
#'
#' @examples
#' \dontrun{
#'   names(my_PhyloSim_list) <- getNames(my_PhyloSim_list)
#' }
#' @seealso \code{\link{PhyloSim}}
#' @export
getNames <- function(runs) {
  if (!inherits(runs, "PhylosimList")) {
    stop("Input must be of class 'PhylosimList'")
  }
  
  sapply(runs, function(x) {
    m <- x$Model
    paste0(
      if (isTRUE(m$density)) paste0("dd", m$compStrength, "_"), # only return the abbreviations if the parameter is turned on !
      "disp", ifelse(m$dispersal == "global", "G_", paste0(m$dispersal,"_")),
      "sr", paste0(m$specRate,"_"),
      if (isTRUE(m$environment)) paste0("e", m$envStrength, "_"), # only return the abbreviations if the parameter is turned on !
      "fbmr", paste0(m$fitnessBaseMortalityRatio,"_"),
      "dc", paste0(m$densityCut,"_"),
      "fao", ifelse(m$fitnessActsOn == "mortality", "M_", "R_"),
      if (m$fission != 0) paste0("fi", m$fission, "_"), # only return the abbreviations if the parameter is turned on !
      if (m$redQueen != 0) paste0("rq", m$redQueen, "_"), # only return the abbreviations if the parameter is turned on !
      if (m$redQueenStrength != 0) paste0("rqs", m$redQueenStrength, "_"), # only return the abbreviations if the parameter is turned on !
      if (m$protracted != 0) paste0("p", m$protracted, "_") # only return the abbreviations if the parameter is turned on !
    )
  })
}
