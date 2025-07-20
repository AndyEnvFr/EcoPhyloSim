#' Expand all matrices toroidally
#' @keywords internal
#' @title Toroidal expansion of simulation matrices
#' @param simu Object of class \code{PhyloSim} or \code{PhylosimList}
#' @param overwrite Logical; if \code{TRUE}, original matrices will be overwritten
#' @param radius Integer; neighborhood radius. Defaults to \code{densityCut}
#' @return Input object with torus-expanded matrices added or replaced
#' @description Adds torus-padding to all relevant matrices mainly to examine nieghbors also of edge cases. This can be done, because the simulation in c++ also uses toroidal boundaries. By default, generates \code{mortMat} and \code{idMat} if missing.
#' Names generations based on \code{Model$runs}.
#'
getTorus <- function(simu, overwrite = FALSE, radius = NULL) {
  UseMethod("getTorus")
}

#' @rdname getTorus
#' @method getTorus PhyloSim
getTorus.PhyloSim <- function(simu, overwrite = TRUE, radius = NULL) {
  # take max dc of ndd and pdd
  r <- if (is.null(radius)) max(simu$Model$pDensityCut,simu$Model$nDensityCut) else radius # by default take density cut as the radius
  
  if (!is.matrix(simu$Output[[1]]$mortMat)) {
    simu <- getMortality(simu) # adds mortalities and generation names to simu$Output$__generationName_
  }
  if (!is.matrix(simu$Output[[1]]$idMat)) {
    simu <- getID(simu) # adds ID matrix
  }
  
  names(simu$Output) <- as.character(simu$Model$runs)
  
  for (i in seq_along(simu$Output)) {
    for (matname in c("specMat", "traitMat", "envMat", "compMat", "neutMat", "mortMat", "idMat")) {
      original <- simu$Output[[i]][[matname]]
      padded <- getTorusMatrix(original, r)
      if (overwrite) { # be default replaces the old matrix with the new bigger one 
        simu$Output[[i]][[matname]] <- padded
      } else { # consumes more space, keeps small and big matrix
        simu$Output[[i]][[paste0(matname, "Big")]] <- padded
      }
    }
  }
  return(simu)
}

#' @rdname getTorus
#' @method getTorus PhylosimList
getTorus.PhylosimList <- function(simu, overwrite = FALSE, radius = NULL) {
  lapply(simu, getTorus, overwrite = overwrite, radius = radius)
}

#' @title Toroidal matrix expansion (helper)
#' @param mat Input matrix
#' @param r Radius to expand
#' @return Expanded matrix with toroidal wrapping
#' @keywords internal
getTorusMatrix <- function(mat, r) { # internal function is called within getTorus
  mrow <- nrow(mat)
  mcol <- ncol(mat)
  mbig <- matrix(NA, mrow + 2 * r, mcol + 2 * r) # creates big torus Matrix
  
  sr <- r + 1 # based on the radius enlargen the matrix
  er <- mrow + r
  sc <- r + 1
  ec <- mcol + r
  
  mbig[sr:er, sc:ec] <- mat # copies small matrix
  mbig[1:r, sc:ec] <- mat[(mrow - r + 1):mrow, ] # adds alle edges and corners to enlargen the matrix
  mbig[(er + 1):(er + r), sc:ec] <- mat[1:r, ]
  mbig[sr:er, 1:r] <- mat[, (mcol - r + 1):mcol]
  mbig[sr:er, (ec + 1):(ec + r)] <- mat[, 1:r]
  mbig[1:r, 1:r] <- mat[(mrow - r + 1):mrow, (mcol - r + 1):mcol]
  mbig[1:r, (ec + 1):(ec + r)] <- mat[(mrow - r + 1):mrow, 1:r]
  mbig[(er + 1):(er + r), 1:r] <- mat[1:r, (mcol - r + 1):mcol]
  mbig[(er + 1):(er + r), (ec + 1):(ec + r)] <- mat[1:r, 1:r]
  
  return(mbig)
}
