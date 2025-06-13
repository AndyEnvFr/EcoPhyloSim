#' Plot or extract species richness over time
#'
#' @title Get species richness time series
#' @param runs A \code{PhyloSim} or \code{PhylosimList} object
#' @param thinning_factor Optional integer to reduce the number of time points (for plotting)
#' @param ymax Optional numeric; fixed y-axis max for plotting
#' @param plot Logical; whether to generate plots
#' @return A \code{data.frame} (for single object) or list of data.frames (for list input)
#' @description Computes species richness per generation using \code{PhyloSim::specRich()}.
#' Can be visualized or just returned. Richness values can optionally be thinned.
#'
#' @export
getSpecTime <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  UseMethod("getSpecTime")
}

#' @rdname getSpecTime
#' @method getSpecTime PhyloSim
#' @export
getSpecTime.PhyloSim <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  sr <- sapply(seq_along(runs$Output), function(i) {
    PhyloSim::specRich(runs, which.result = i)
  })
  yr <- runs$Model$runs
  result <- data.frame(year = yr, spec_rich = sr)
  
  if (!is.null(thinning_factor) && thinning_factor > 1) {
    idx <- seq(1, nrow(result), thinning_factor)
    result <- result[idx, ]
  }
  
  if (plot) {
    plot(result$year, result$spec_rich, type = "b", ylab = "richness", xlab = "generation",
         ylim = if (!is.null(ymax)) c(0, ymax) else c(0, max(result$spec_rich)),
         main = "Species Richness Over Time")
  }
  
  return(result)
}

#' @rdname getSpecTime
#' @method getSpecTime PhylosimList
#' @export
getSpecTime.PhylosimList <- function(runs, thinning_factor = NULL, ymax = NULL, plot = TRUE) {
  results <- lapply(runs, getSpecTime, thinning_factor = thinning_factor, ymax = ymax, plot = plot)
  names(results) <- names(runs)
  return(results)
}
