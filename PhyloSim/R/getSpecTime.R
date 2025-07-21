#' Generate simulation names based on model parameters
#'
#' @title Generate simulation names for simulation objects
#' @param runs Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return For PhyloSim: the modified object with getName set. For PhylosimList: character vector of names.
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
#'   \item \code{faoX}: Fitness acts on (M = mortality, R = reproduction)
#'   \item \code{fiX}, \code{rqX}, \code{rqsX}, \code{pX}: Optional parameters if > 0
#' }
#'
#' @examples
#' \dontrun{
#' # For list of simulations
#' names(my_PhyloSim_list) <- getNames(my_PhyloSim_list)
#'
#' # For single simulation
#' my_simulation <- getNames(my_simulation)  # Returns updated object
#' my_simulation$Model$getName  # Now contains the name
#' }
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
    ifelse(isTRUE(m$positiveDensity), paste0("pdd", m$pDDStrength, "Var", m$pDDNicheWidth, "Cut", m$pDensityCut, "_"), "pdd0_"),
    ifelse(isTRUE(m$negativeDensity), paste0("ndd", m$nDDStrength, "Var", m$nDDNicheWidth, "Cut", m$nDensityCut, "_"), "ndd0_"),
    "disp", ifelse(m$dispersal == "global", "G_", paste0(m$dispersal, "_")),
    "sr", m$specRate, "_",
    if (isTRUE(m$environment)) paste0("e", m$envStrength, "_"),
    "fbmr", m$fitnessBaseMortalityRatio, "_",
    "fao", ifelse(m$fitnessActsOn == "mortality", "M", "R"),
    if (m$fission > 0) paste0("_fi", m$fission),
    if (m$redQueen > 0) paste0("_rq", m$redQueen),
    if (m$redQueenStrength > 0) paste0("_rqs", m$redQueenStrength),
    if (m$protracted > 0) paste0("_p", m$protracted)
  )
  
  # Set the name in the object and return the modified object
  runs$Model$getName <- name
  return(runs)
}
#' @rdname getNames
#' @method getNames PhylosimList
#' @export
getNames.PhylosimList <- function(runs) {
  names <- sapply(runs, function(x) {
    getNames.PhyloSim(x)$Model$getName
  })
  
  return(names)
}


#' Plot or extract species richness over time
#'
#' @title Get species richness time series
#' @param runs A \code{PhyloSim} or \code{PhylosimList} object
#' @param thinning_factor Optional integer to reduce the number of time points for plotting (smootes the line). Turned off by default.
#' @param ymax Optional numeric; fixed y-axis max for plotting. Usefull for comparing different runs. Turned off by default. 
#' @param plot Logical; whether to generate plots. Default is TRUE
#' @return A \code{data.frame} (for single object) or list of data.frames (for list input)
#' @description Computes species richness per generation using \code{PhyloSim::specRich()}.
#' Can be visualized or just returned. Richness values can optionally be thinned.
#'
#' @export
getSpecTime <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  UseMethod("getSpecTime")
}
#' @rdname getSpecTime
#' @method getSpecTime PhyloSim
#' @export
# TODO: merge getSpecTime and specRich: when specRich is called with "which.result" = "all", getSpecTime can be called.
getSpecTime.PhyloSim <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  sr <- sapply(seq_along(runs$Output), function(i) {
    PhyloSim::specRich(runs, which.result = i) # uses specRich function. This gives only the richness for one generation, not for a time series.
  })
  yr <- runs$Model$runs
  result <- data.frame(year = yr, spec_rich = sr)
  
  if (!is.null(thinning_factor) && thinning_factor > 1) {
    idx <- seq(1, nrow(result), thinning_factor)
    result <- result[idx, ]
  }
  
  if (plot) {
    # Get title - use getName if it exists, otherwise generate it
    title <- if (!is.null(runs$Model$getName)) {
      runs$Model$getName
    } else {
      # Generate name on the fly if not set
      temp_runs <- getNames(runs)
      temp_runs$Model$getName
    }
    
    plot(result$year, result$spec_rich, type = "l", ylab = "richness", xlab = "generation",
         ylim = if (!is.null(ymax)) c(0, ymax) else c(0, max(result$spec_rich)),
         main = title)
  }
  
  return(result)
}
#' @rdname getSpecTime
#' @method getSpecTime PhylosimList
#' @export
getSpecTime.PhylosimList <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  results <- lapply(runs, getSpecTime, thinning_factor = thinning_factor, ymax = ymax, plot = plot)
  names(results) <- names(runs)
  return(results)
}