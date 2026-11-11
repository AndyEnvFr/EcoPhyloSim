#' Get mortality data for density dependence analyses
#'
#' @description Returns a table with one row per cell and census interval: the
#'   individual alive at a recorded generation, its conspecific neighbors and
#'   whether it died before the next recorded generation. The table can be used
#'   to analyze conspecific negative density dependence (CNDD), e.g., with a
#'   binomial GLM or GLMM.
#'
#' @param simu simulation result of class \code{PhyloSim} from \code{\link{runSimulation}},
#'   or a list of them, e.g., from \code{\link{runSimulationBatch}}. At least two
#'   generations must be recorded (see \code{runs} in \code{\link{createCompletePar}}).
#' @param radius radius of the neighborhood, a non-negative integer. If \code{NULL},
#'   the \code{densityCut} of the simulation is used.
#'
#' @details The neighborhood is the same as for density dependence in the model:
#'   all cells with \code{dx^2 + dy^2 <= radius^2} around the focal cell, without
#'   the focal cell, on a torus.
#'
#'   An individual counts as dead if its cell was replaced at least once until the
#'   next recorded generation. This is detected by a change of \code{traitMat}, as
#'   the traits change at every birth. On average, a cell is replaced about once
#'   per generation, so the recorded generations should be consecutive, e.g.,
#'   \code{runs = 1000:1010} after a burn-in.
#'
#'   In the model, density dependence acts through the competition trait
#'   (\code{compTrait}), not through species identity: fitness depends on the
#'   absolute difference in this trait to the neighbors. The number of
#'   conspecific neighbors \code{con} is therefore only a proxy. As the grid is
#'   always full, the number of heterospecific neighbors is the number of cells
#'   in the neighborhood minus \code{con}. With \code{environment > 0}, fitness
#'   also depends on the difference between \code{envTrait} and \code{env}.
#'
#'   Surviving individuals appear in several rows; account for this, e.g., with a
#'   random effect for \code{id}.
#'
#' @return A data frame with the columns
#'   \describe{
#'     \item{run}{only for lists: name or index of the simulation}
#'     \item{census}{recorded generation}
#'     \item{interval}{number of generations until the next recorded generation}
#'     \item{id}{individual ID, constant while the individual survives; for lists, unique only together with \code{run}}
#'     \item{x, y}{column and row of the cell in the output matrices}
#'     \item{species}{species ID}
#'     \item{con}{number of conspecific neighbors within \code{radius}}
#'     \item{compTrait, envTrait}{competition and environmental trait of the individual}
#'     \item{env}{environment of the cell}
#'     \item{died}{\code{TRUE} if the individual died before the next recorded generation}
#'   }
#'
#' @example inst/examples/getMortalityData-help.R
#' @export
getMortalityData <- function(simu, radius = NULL) {
  if (!inherits(simu, "PhyloSim")) { # list of simulations, e.g. from runSimulationBatch
    if (!is.list(simu) || !length(simu) || !all(sapply(simu, inherits, "PhyloSim"))) stop("simu must be a PhyloSim object or a list of them")
    tabs <- lapply(simu, getMortalityData, radius = radius)
    run <- if (is.null(names(simu))) seq_along(simu) else names(simu)
    out <- cbind(run = rep(run, sapply(tabs, nrow)), do.call(rbind, tabs))
    rownames(out) <- NULL
    return(out)
  }
  out <- simu$Output
  if (length(out) < 2) stop("at least two recorded generations are needed, e.g. runs = c(1000, 1001)")
  if (is.null(radius)) radius <- simu$Model$densityCut
  ny <- nrow(out[[1]]$specMat); nx <- ncol(out[[1]]$specMat)
  off <- getCircularOffsets(radius)
  off <- off[off$dx %% nx != 0 | off$dy %% ny != 0, , drop = FALSE] # skip the focal cell also if the neighborhood wraps around it
  row <- c(row(out[[1]]$specMat)); col <- c(col(out[[1]]$specMat))
  nb <- vapply(seq_len(nrow(off)), function(k) (row - 1L + off$dy[k]) %% ny + 1L + ((col - 1L + off$dx[k]) %% nx) * ny, integer(length(row)))
  id <- seq_along(row)
  tabs <- vector("list", length(out) - 1)
  for (i in seq_along(tabs)) {
    o <- out[[i]]
    died <- c(out[[i + 1]]$traitMat != o$traitMat)
    tabs[[i]] <- data.frame(census = simu$Model$runs[i], interval = simu$Model$runs[i + 1] - simu$Model$runs[i],
                            id = id, x = col, y = row, species = c(o$specMat),
                            con = as.integer(rowSums(matrix(o$specMat[nb] == c(o$specMat), nrow = length(row)))),
                            compTrait = c(o$compMat), envTrait = c(o$traitMat), env = c(o$envMat), died = died)
    id[died] <- max(id) + seq_len(sum(died)) # new individuals get new IDs
  }
  do.call(rbind, tabs)
}
