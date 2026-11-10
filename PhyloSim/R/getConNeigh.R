#' Calculate conspecific neighborhood matrices
#'
#' @title Get number of conspecific neighbors
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Input object with \code{conNeighMat} matrices added to all generations.
#'   Each \code{conNeighMat} contains integer counts of conspecific neighbors
#'   for every cell in the simulation grid.
#' @description 
#' Calculates the number of conspecific neighbors within a circular neighborhood
#' for all generations in a simulation. The neighborhood radius
#' is defined by the \code{densityCut} parameter from the simulation model.
#' Conspecific neighbors are cells containing individuals of the same species
#' as the focal cell.
#' 
#' The function uses circular neighborhoods as defined by 
#' \code{\link{getCircularOffsets}}, which generates coordinate offsets
#' for all cells within the specified radius distance from each focal cell.
#' 
#' @details
#' The simulation grid is treated as a torus (periodic boundaries), as in
#' the C++ simulation. If \code{mortMat} or \code{idMat} are missing, they are
#' added with \code{\link{getMortality}} and \code{\link{getID}}, as needed by
#' \code{\link{getMatToTab}}. The parameter name from \code{\link{getNames}}
#' is stored in \code{Model$getName}. For a \code{PhylosimList}, existing list
#' names are kept; otherwise the list is named with \code{\link{getNames}}.
#' 
#' @seealso 
#' \code{\link{getCircularOffsets}} for the neighborhood definition used
#' \code{\link{getID}} for extraction of ID
#' \code{\link{getMortality}} for extraction of mortalities
#' 
#' @example /inst/examples/getConNeigh-help.R
#' @export
getConNeigh <- function(simu) {
  UseMethod("getConNeigh")
}

#' @rdname getConNeigh
#' @method getConNeigh PhyloSim
#' @export
getConNeigh.PhyloSim <- function(simu) {
  simu$Model$getName <- getNames(simu)
  if (is.null(simu$Output[[1]]$mortMat)) simu <- getMortality(simu)
  if (is.null(simu$Output[[1]]$idMat)) simu <- getID(simu)
  
  r <- simu$Model$densityCut
  offsets <- getCircularOffsets(r)
  
  for (i in seq_along(simu$Output)) {
    spec <- simu$Output[[i]]$specMat
    big <- getTorusMatrix(spec, r)
    rows <- r + seq_len(nrow(spec))
    cols <- r + seq_len(ncol(spec))
    con <- matrix(0L, nrow(spec), ncol(spec))
    for (k in seq_len(nrow(offsets))) {
      con <- con + (big[rows + offsets$dx[k], cols + offsets$dy[k]] == spec)
    }
    simu$Output[[i]]$conNeighMat <- con
  }
  
  return(simu)
}

#' @rdname getConNeigh
#' @method getConNeigh PhylosimList
#' @export
getConNeigh.PhylosimList <- function(simu){
  processed_list <- lapply(simu, getConNeigh)
  class(processed_list) <- "PhylosimList"
  if (is.null(names(simu))) names(processed_list) <- getNames(processed_list)
  return(processed_list)
}
