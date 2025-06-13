#' Calculate conspecific neighborhood matrices
#'
#' @title Get number of conspecific neighbors
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Input object with \code{conNeighMat} matrices added to all generations
#' @description Calculates the number of conspecific neighbors in a circular radius 
#' for all generations based on torus-expanded \code{specMat}. 
#' Automatically calls \code{\link{getTorus}} with \code{overwrite = TRUE}, 
#' which also prepares \code{idMat} and \code{mortMat} if missing.
#' 
#' The torus is automatically undone after computing neighborhood matrices.
#'
#' @export
getConNeigh <- function(simu) {
  UseMethod("getConNeigh")
}

#' @rdname getConNeigh
#' @method getConNeigh PhyloSim
#' @export
getConNeigh.PhyloSim <- function(simu) {
  simu <- getTorus(simu, overwrite = TRUE)
  
  radius <- simu$Model$densityCut
  offsets <- getCircularOffsets(radius)
  
  r <- radius
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
  
  # Undo torus for all matrices
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
getConNeigh.PhylosimList <- function(simu) {
  structure(lapply(simu, getConNeigh), class = "PhylosimList")
}
