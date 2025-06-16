#' Assign individual IDs across generations
#'
#' @title Generate ID matrix for individuals
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Input object with \code{idMat} matrices added to each generation.  Also, adds the years of generations (e.g., \code{simu$Output$test$Output$`98`}.
#' @description Each individual gets a persistent ID. IDs are unique per run.
#' If the trait stays the same between generations, the ID persists. Otherwise, a new one is assigned.
#'
#' @export
getID <- function(simu) {
  UseMethod("getID")
}

#' @rdname getID
#' @method getID PhyloSim
#' @export
getID.PhyloSim <- function(simu) {
  names(simu$Output) <- as.character(simu$Model$runs) # adds generation names in the list
  
  rXc <- simu$Model$x * simu$Model$y # get dimensions
  ids <- lapply(seq_along(simu$Output), function(i) (1:rXc) + i * rXc)
  
  for (i in seq_along(simu$Output)) {
    simu$Output[[i]]$idMat <- matrix(NA, simu$Model$x, simu$Model$y)
  }
  
  simu$Output[[1]]$idMat[] <- 1:rXc
  
  for (i in 2:length(simu$Output)) {
    unchanged <- simu$Output[[i]]$traitMat == simu$Output[[i - 1]]$traitMat # ids are based on traits, like mortalities
    simu$Output[[i]]$idMat[unchanged] <- simu$Output[[i - 1]]$idMat[unchanged]
    simu$Output[[i]]$idMat[!unchanged] <- ids[[i]][!unchanged]
  }
  
  return(simu)
}

#' @rdname getID
#' @method getID PhylosimList
#' @export
getID.PhylosimList <- function(simu) {
  structure(lapply(simu, getID), class = "PhylosimList")
}
