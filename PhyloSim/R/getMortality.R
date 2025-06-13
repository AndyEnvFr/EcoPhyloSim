#' Compute individual mortality across generations
#'
#' @title Get mortality of individuals
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Same object, with \code{mortMat} matrices added to \code{Output[[i]]}
#' @description Compares each generation's \code{traitMat} to the next. If the trait
#' differs, the individual is considered dead. First generation gets NA. Generation
#' names are set using \code{Model$runs}.
#'
#' @export
getMortality <- function(simu) {
  UseMethod("getMortality")
}

#' @rdname getMortality
#' @method getMortality PhyloSim
#' @export
getMortality.PhyloSim <- function(simu) {
  names(simu$Output) <- as.character(simu$Model$runs)
  
  for (i in seq_along(simu$Output)) {
    simu$Output[[i]]$mortMat <- matrix(NA, nrow = simu$Model$x, ncol = simu$Model$y)
  }
  
  for (i in 1:(length(simu$Output) - 1)) {
    current <- simu$Output[[i]]$traitMat
    nex <- simu$Output[[i + 1]]$traitMat
    simu$Output[[i + 1]]$mortMat <- ifelse(current != nex, TRUE, FALSE)
  }
  
  return(simu)
}

#' @rdname getMortality
#' @method getMortality PhylosimList
#' @export
getMortality.PhylosimList <- function(simu) {
  lapply(simu, getMortality)
}
