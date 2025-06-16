#' Convert simulation matrices to long-format table
#'
#' @title Extract Matrices as Long-Format Table
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return A \code{data.frame} with census, individual ID, species ID, mortality status, and conspecific count
#' @description Extracts key spatial matrices (ID, species, mortality, conspecific neighborhood) 
#' from all available generations of a simulation object and converts them to a tabular format
#' suitable for statistical analysis.
#'
#' The output table includes:
#' \itemize{
#'   \item \code{census} - Generation label (current generation)
#'   \item \code{indId} - Individual ID from \code{idMat} (current generation)
#'   \item \code{specId} - Species ID from \code{specMat} (current generation)
#'   \item \code{mortNextGen} - Mortality status from the NEXT generation's \code{mortMat}
#'   \item \code{con} - Number of conspecific neighbors from \code{conNeighMat} (current generation)
#' }
#'
#' @details
#' The input must be preprocessed using \code{\link{getConNeigh}}, which automatically computes and
#' assigns the required matrices (\code{idMat}, \code{mortMat}, \code{conNeighMat}).
#' 
#' Note: All values except mortality status refer to the current generation. Mortality status
#' is pulled from the subsequent generation's mortality matrix, representing whether each
#' individual survived to the next time step.
#'
#' All available generations will be included in the tabular output, regardless of spacing or pattern.
#'
#' @seealso \code{\link{getConNeigh}}, \code{\link{getMortality}}, \code{\link{getID}}, \code{\link{getTorus}}
#' @export
getMatToTab <- function(simu) {
  UseMethod("getMatToTab")
}

#' @rdname getMatToTab
#' @method getMatToTab PhyloSim
#' @export
getMatToTab.PhyloSim <- function(simu) {
  if (!"conNeighMat" %in% names(simu$Output[[1]])) {
    stop("conNeighMat not found. Please preprocess using getConNeigh().")
  } # the number of conspecific nerighbors must be calculated before !
  
  census <- names(simu$Output)
  genN <- length(census)
  dimInner <- dim(simu$Output[[1]]$specMat)
  
  result <- data.frame(
    census = character(),
    indId = integer(),
    specId = integer(),
    mortNextGen = logical(),
    con = integer(),
    stringsAsFactors = FALSE
  ) # get these cols in the returninng df
  
  for (cidx in seq_len(genN - 1)) {
    cen <- census[cidx]
    
    interim <- data.frame(
      census = rep(cen, dimInner[1] * dimInner[2]),
      indId = as.vector(simu$Output[[cidx]]$idMat),
      specId = as.vector(simu$Output[[cidx]]$specMat),
      con = as.vector(simu$Output[[cidx]]$conNeighMat),
      mortNextGen = as.vector(simu$Output[[cidx + 1]]$mortMat)
    ) #intermin results, for one generation is appendend (rbind) consecutively
    
    result <- rbind(result, interim)
  }
  
  return(result)
}

#' @rdname getMatToTab
#' @method getMatToTab PhylosimList
#' @export
getMatToTab.PhylosimList <- function(simu) {
  do.call(rbind, lapply(simu, getMatToTab)) # return a list with one dataframe per parameter setting. One dataframe contains all generations of one parameter setting.
}