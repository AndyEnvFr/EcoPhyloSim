#' @title Rank Abundance Curve
#' @description Plots the Rank Abundance Curve for a given community. 
#' @param runs A simulation object of class "PhyloSim" or "PhylosimList".
#' @param which.result Integer or "all". For PhyloSim, selects the result(s) to use. For PhylosimList, applies to all.
#' @param plot_type Type of plot: "line" (default) or "bar".
#' @param title Optional plot title. For PhylosimList, a character vector or NULL.
#' @param ymax Numeric, y-axis max. If NULL, uses automatic scaling.
#' @param xmax Numeric, x-axis max. If NULL, uses automatic scaling.
#' @return Dataframe(s) with rank-abundance values (only if not which.result="all").
#' @export
rac <- function(runs, which.result = NULL, plot_type = "line", title = NULL, ymax = NULL, xmax = NULL) {
  UseMethod("rac")
}

#' @rdname rac
#' @method rac PhyloSim
#' @export
rac.PhyloSim <- function(runs, which.result = NULL, plot_type = "line", title = NULL, ymax = NULL, xmax = NULL) {
  if (is.null(which.result)) which.result <- length(runs$Output)
  
  if (which.result == "all") {
    if (plot_type == "bar") stop("Argument 'bar' not possible for which.result='all'")
    simulations <- seq_along(runs$Output)
  } else {
    simulations <- which.result
  }
  
  colfunc <- colorRampPalette(c("blue", "red"))
  cols <- colfunc(length(simulations))
  
  RAC <- list()
  max_abundance <- 0
  max_rank <- 0
  
  for (i in simulations) {
    matrix <- runs$Output[[i]]$specMat
    Abundances <- as.data.frame(table(matrix))
    sel <- order(Abundances$Freq, decreasing = TRUE)
    RAC[[i]] <- data.frame(
      Rank = seq_len(nrow(Abundances)),
      Abundance = Abundances$Freq[sel],
      Species = Abundances$matrix[sel]
    )
    max_abundance <- max(max_abundance, max(RAC[[i]]$Abundance))
    max_rank <- max(max_rank, max(RAC[[i]]$Rank))
  }
  
  ylim_max <- if (is.null(ymax)) max_abundance else ymax
  xlim_max <- if (is.null(xmax)) max_rank else xmax
  
  plot_title <- if (!is.null(title)) {
    title
  } else if (!is.null(runs$Model$getName)) {
    runs$Model$getName
  } else {
    getNames(runs)$Model$getName
  }
  
  if (plot_type == "bar") {
    barplot(RAC[[i]]$Abundance, log = "y", ylab = "Log Abundance", xlab = "Rank",
            main = plot_title, names.arg = RAC[[i]]$Rank, ylim = c(1, ylim_max))
  }
  
  if (plot_type == "line") {
    if (length(simulations) == 1) {
      plot(RAC[[i]]$Rank, RAC[[i]]$Abundance, type = "l", log = "y", ylab = "Log Abundance",
           xlab = "Rank", main = plot_title, lwd = 2, xlim = c(0, xlim_max), ylim = c(1, ylim_max))
    } else {
      for (i in simulations) {
        if (i == simulations[1]) {
          plot(RAC[[i]]$Rank, RAC[[i]]$Abundance, type = "l", log = "y", ylab = "Log Abundance",
               xlab = "Rank", main = plot_title, col = cols[i], xlim = c(0, xlim_max), ylim = c(1, ylim_max))
        } else {
          lines(RAC[[i]]$Rank, RAC[[i]]$Abundance, col = cols[i])
        }
      }
    }
  }
  
  if (length(simulations) == 1) return(RAC[[simulations]])
}

#' @rdname rac
#' @method rac PhylosimList
#' @export
rac.PhylosimList <- function(runs, which.result = NULL, plot_type = "line", title = NULL, ymax = NULL, xmax = NULL) {
  if (!is.null(title) && length(title) != length(runs)) {
    warning("Length of title vector does not match length of runs. Using automatic titles.")
    title <- NULL
  }
  
  results <- lapply(seq_along(runs), function(i) {
    current_title <- if (!is.null(title)) title[i] else NULL
    rac(runs[[i]], which.result = which.result, plot_type = plot, title = current_title, ymax = ymax, xmax = xmax)
  })
  
  names(results) <- names(runs)
  return(results)
}
