#' Convert simulation matrices to long-format table
#'
#' @title Extract matrices as long-format table
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @return A \code{data.frame} with census, individual ID, species ID, mortality, and conspecific count
#' @description Extracts key spatial matrices (ID, species, mortality, conspecific neighborhood) 
#' from all available generations of a simulation object and flattens them into a tabular format
#' suitable for statistical analysis.
#'
#' The output table includes:
#' \itemize{
#'   \item \code{census} – generation label
#'   \item \code{ind_id} – individual ID from \code{idMat}
#'   \item \code{spec_id} – species ID from \code{specMat}
#'   \item \code{mort} – mortality status from the next generation's \code{mortMat}
#'   \item \code{con} – number of conspecific neighbors from \code{conNeighMat}
#' }
#'
#' @details
#' The input must be preprocessed using \code{\link{getConNeigh}}, which automatically computes and
#' assigns the required matrices (\code{idMat}, \code{mortMat}, \code{conNeighMat}).
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
  }
  
  census <- names(simu$Output)
  gen_n <- length(census)
  dim_inner <- dim(simu$Output[[1]]$specMat)
  
  result <- data.frame(
    census = character(),
    ind_id = integer(),
    spec_id = integer(),
    mort = logical(),
    con = integer(),
    stringsAsFactors = FALSE
  )
  
  for (cidx in seq_len(gen_n - 1)) {
    cen <- census[cidx]
    
    interim <- data.frame(
      census = rep(cen, dim_inner[1] * dim_inner[2]),
      ind_id = as.vector(simu$Output[[cidx]]$idMat),
      spec_id = as.vector(simu$Output[[cidx]]$specMat),
      con = as.vector(simu$Output[[cidx]]$conNeighMat),
      mort = as.vector(simu$Output[[cidx + 1]]$mortMat)
    )
    
    result <- rbind(result, interim)
  }
  
  return(result)
}

#' @rdname getMatToTab
#' @method getMatToTab PhylosimList
#' @export
getMatToTab.PhylosimList <- function(simu) {
  do.call(rbind, lapply(simu, getMatToTab))
}
