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
#' The function operates by:
#' \enumerate{
#'   \item Expanding all matrices to torus topology (periodic boundary conditions)
#'   \item For each generation, counting conspecific neighbors using circular offsets
#'   \item Cropping all matrices back to original dimensions
#'   \item Preparing missing \code{idMat} and \code{mortMat} by calling 
#'         \code{getID} and \code{getMortality} functions if needed
#' }
#' 
#' The torus expansion ensures that edge effects are minimized by treating
#' the simulation grid as if it wraps around at the boundaries. After
#' neighborhood calculations are complete, all matrices are automatically
#' cropped back to their original dimensions.
#' 
#' @seealso 
#' \code{\link{getCircularOffsets}} for the neighborhood definition used
#' \code{\link{getID}} for extraction of ID
#' \code{\link{getMortality}} for extraction of mortalities
#' 
#' @examples
#' \dontrun{
#' # Apply to single simulation
#' simu_with_neighbors <- getConNeigh(my_simulation)
#' my_simulation$Model$getName <- getNames(my_simulation)
#' 
#' # Apply to list of simulations
#' simu_list_with_neighbors <- getConNeigh(my_simulation_list)
#' names(simu_list_with_neighbors) <- getNames(simu_list_with_neighbors)
#' simu_list_with_neighbors[[1]]$Model$getName
#' }
#' 
#' @export
getConNeigh <- function(simu) {
  UseMethod("getConNeigh")
}

#' @rdname getConNeigh
#' @method getConNeigh PhyloSim
#' @export
getConNeigh.PhyloSim <- function(simu) {
  # Add name to Model$getName for downstream tracking
  simu$Model$getName <- getNames(simu)
  
  simu <- getTorus(simu, overwrite = TRUE)
  
  r <- max(simu$Model$pDensityCut,simu$Model$nDensityCut)
  offsets <- getCircularOffsets(r)
  
  lx <- nrow(simu$Output[[1]]$specMat)
  ly <- ncol(simu$Output[[1]]$specMat)
  
  census <- names(simu$Output)
  for (cen in census) {
    sx <- c(r + 1, lx - r)
    sy <- c(r + 1, ly - r)
    con <- matrix(0, lx - 2 * r, ly - 2 * r)
    inner <- simu$Output[[cen]]$specMat[sx[1]:sx[2], sy[1]:sy[2]]
    
    for (xy in seq_len(nrow(offsets))) {
      xshift <- offsets$dx[xy]
      yshift <- offsets$dy[xy]
      X <- sx + xshift
      Y <- sy + yshift
      shifted <- simu$Output[[cen]]$specMat[X[1]:X[2], Y[1]:Y[2]]
      con <- con + ifelse(shifted == inner, 1, 0)
    }
    
    simu$Output[[cen]]$conNeighMat <- con
  }
  
  crop <- function(mat) mat[(r + 1):(lx - r), (r + 1):(ly - r)]
  for (i in seq_along(simu$Output)) {
    for (mat in c("specMat", "traitMat", "envMat", "compMat", "neutMat", "mortMat", "idMat")) {
      simu$Output[[i]][[mat]] <- crop(simu$Output[[i]][[mat]])
    }
  }
  
  return(simu)
}

#' @rdname getConNeigh
#' @method getConNeigh PhylosimList
#' @export
getConNeigh.PhylosimList <- function(simu){
  processed_list <- lapply(simu, getConNeigh)
  class(processed_list) <- "PhylosimList"
  names(processed_list) <- getNames(processed_list)
  return(processed_list)
}
