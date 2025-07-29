#' Convert simulation matrices to long-format table
#'
#' @title Extract Matrices as Long-Format Table
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @param detailedParams Logical. If \code{FALSE} (default), returns single \code{params} column.
#' If \code{TRUE}, returns separate columns for each parameter.
#' @return A \code{data.frame} with census, individual ID, species ID, mortality status, 
#' conspecific count, and simulation parameter label(s)
#' @description Extracts key spatial matrices (ID, species, mortality, conspecific neighborhood) 
#' from all available generations of a simulation object and converts them to a tabular format
#' suitable for statistical analysis.
#'
#' The output table includes:
#' \itemize{
#'   \item \code{census} - Generation label (current generation)
#'   \item \code{indId} - Individual ID from \code{idMat}
#'   \item \code{specId} - Species ID from \code{specMat}
#'   \item \code{mortNextGen} - Mortality status from the NEXT generation's \code{mortMat}
#'   \item \code{con} - Number of conspecific neighbors from \code{conNeighMat}
#'   \item \code{params} - Parameter string label from \code{Model$getName} (if \code{detailedParams = FALSE})
#'   \item Individual parameter columns (if \code{detailedParams = TRUE})
#' }
#'
#' @details
#' The input must be preprocessed using \code{\link{getConNeigh}}, which computes and
#' assigns the required matrices (\code{idMat}, \code{mortMat}, \code{conNeighMat}).
#' 
#' Note: All values except mortality refer to the current generation. Mortality status
#' is pulled from the subsequent generation.
#'
#' When \code{detailedParams = TRUE}, individual parameter columns are extracted from the
#' simulation object and added to the output table.
#'
#' @seealso \code{\link{getConNeigh}}, \code{\link{getMortality}}, \code{\link{getID}}, \code{\link{getTorus}}
#' @export
getMatToTab <- function(simu, detailedParams = FALSE) {
  UseMethod("getMatToTab")
}


#' @rdname getMatToTab
#' @method getMatToTab PhyloSim
#' @export
getMatToTab.PhyloSim <- function(simu, detailedParams = FALSE) {
  if (!"conNeighMat" %in% names(simu$Output[[1]])) {
    stop("conNeighMat not found. Please preprocess using getConNeigh().")
  }
  
  census <- names(simu$Output)
  genN <- length(census)
  dimInner <- dim(simu$Output[[1]]$specMat)
  paramName <- simu$Model$getName
  
  # Base columns
  base_cols <- data.frame(
    census = character(),
    indId = integer(),
    specId = integer(),
    mortNextGen = logical(),
    con = integer(),
    stringsAsFactors = FALSE
  )
  
  # Add parameter columns based on detailedParams flag
  if (detailedParams) {
    # Extract individual parameter values
    param_cols <- data.frame(
      pDD = numeric(),
      nDD = numeric(),
      pDDVar = numeric(),
      nDDVar = numeric(),
      pDC = numeric(),
      nDC = numeric(),
      disp = numeric(),
      fao = character(),
      fbmr = numeric(),
      sr = numeric(),
      stringsAsFactors = FALSE
    )
  } else {
    # Single params column
    param_cols <- data.frame(
      params = character(),
      stringsAsFactors = FALSE
    )
  }
  
  result <- cbind(base_cols, param_cols)
  
  for (cidx in seq_len(genN - 1)) {
    cen <- census[cidx]
    n_rows <- dimInner[1] * dimInner[2]
    
    # Base data
    interim <- data.frame(
      census = as.numeric(rep(cen, n_rows)),
      indId = as.vector(simu$Output[[cidx]]$idMat),
      specId = as.vector(simu$Output[[cidx]]$specMat),
      mortNextGen = as.vector(simu$Output[[cidx + 1]]$mortMat),
      con = as.vector(simu$Output[[cidx]]$conNeighMat),
      stringsAsFactors = FALSE
    )
    
    # Add parameter columns
    if (detailedParams) {
      
      # Compute abundance per species using base R
      spec_vec <- as.vector(simu$Output[[cidx]]$specMat)
      abund_table <- table(spec_vec)
      abund_lookup <- as.numeric(abund_table[as.character(spec_vec)])
      
      param_data <- data.frame(
        abund = abund_lookup,
        pDD = rep(ifelse(simu$Model$positiveDensity, simu$Model$pDDStrength, 0), n_rows),
        nDD = rep(ifelse(simu$Model$negativeDensity, simu$Model$nDDStrength, 0), n_rows),
        pDDVar = rep(simu$Model$pDDNicheWidth, n_rows),
        nDDVar = rep(simu$Model$nDDNicheWidth, n_rows),
        pDC = rep(simu$Model$pDensityCut, n_rows),
        nDC = rep(simu$Model$nDensityCut, n_rows),
        disp = rep(simu$Model$dispersal, n_rows),
        fao = rep(simu$Model$fitnessActsOn, n_rows),
        fbmr = rep(simu$Model$fitnessBaseMortalityRatio, n_rows),
        sr = rep(simu$Model$specRate, n_rows),
        stringsAsFactors = FALSE
      )
    } else {
      param_data <- data.frame(
        params = rep(paramName, n_rows),
        stringsAsFactors = FALSE
      )
    }
    
    interim <- cbind(interim, param_data)
    result <- rbind(result, interim)
  }
  
  return(result)
}

#' @rdname getMatToTab
#' @method getMatToTab PhylosimList
#' @export
getMatToTab.PhylosimList <- function(simu, detailedParams = FALSE) {
  do.call(rbind, lapply(simu, function(x) getMatToTab(x, detailedParams = detailedParams)))
}