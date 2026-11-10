#' Compute individual mortality across generations
#'
#' @title Get mortality of individuals
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return Same object, with \code{mortMat} matrices added to \code{Output[[i]]}. Also, adds the years of generations (e.g., \code{simu$Output$test$Output$`98`}.
#' @description Compares each generation's \code{traitMat} to the next. If the trait
#' differs, the individual is considered dead. First generation gets NA. Generation
#' names are set using \code{Model$runs}.
#'
#' Only the recorded generations (\code{runs} in \code{\link{createCompletePar}})
#' are compared, so \code{TRUE} means that the individual in this cell was replaced
#' at least once since the previous recorded generation, not that it died in the
#' last generation.
#'
#' @example /inst/examples/getConNeigh-help.R
#' @export
getMortality <- function(simu) {
  UseMethod("getMortality")
}

#' @rdname getMortality
#' @method getMortality PhyloSim
#' @export
getMortality.PhyloSim <- function(simu) {
  names(simu$Output) <- as.character(simu$Model$runs) # adds generation names in the list
  
  for (i in seq_along(simu$Output)) {
    simu$Output[[i]]$mortMat <- matrix(NA, nrow(simu$Output[[i]]$traitMat), ncol(simu$Output[[i]]$traitMat))
  }
  
  for (i in seq_along(simu$Output)[-1]) { # mortality based on trait values
    simu$Output[[i]]$mortMat <- simu$Output[[i - 1]]$traitMat != simu$Output[[i]]$traitMat
  }
  
  return(simu)
}

#' @rdname getMortality
#' @method getMortality PhylosimList
#' @export
getMortality.PhylosimList <- function(simu) {
  structure(lapply(simu, getMortality), class = "PhylosimList")
}
