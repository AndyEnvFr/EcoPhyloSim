#' @title Species Area Curve
#' @description Plots the species area curve for a given community. 
#' @param runs A simulation object of class "PhyloSim" or "PhylosimList".
#' @param which.result Integer or "all". For PhyloSim, selects the result(s) to use. For PhylosimList, applies to all.
#' @param size a single value or a vector determining the edge length of the subplots. If not provided, 10 logarithmic plotsizes will be used.
#' @param rep The number of repetitions per size to calculate the mean 
#' @param plot Logical determining whether to plot the SAC or not
#' @param nested Logical, determining whether nested subplots should be used.
#' @param title Optional plot title. For PhylosimList, a character vector or NULL.
#' @details Displays the accumulated species richness as a function of plot size or the amount of equally sized plots. It serves as an indicator for the clustering of a species community. A positively bent curve usually indicates clustering since an increase in plot size or number leads to an increase in species richness while a negatively bent curve indicates a more neutral distribution of species within the community. \cr\cr If which.result = "all" all intermediate results are shown in one plot. The colors of the lines are plotted as a gradient from blue (first results) to red (end result).
#' @return A list containing the mean species richness for each size and the respective standard deviation. If which.result = "all" only the plot will be returned. For PhylosimList, returns a list of results.
#' @example /inst/examples/plotSAC-help.R
#' @export
sac <- function(runs, which.result = NULL, size = NULL, rep = 50, plot = TRUE, nested = FALSE, title = NULL) {
  UseMethod("sac")
}

#' @rdname sac
#' @method sac PhyloSim
#' @export
sac.PhyloSim <- function(runs, which.result = NULL, size = NULL, rep = 50, plot = TRUE, nested = FALSE, title = NULL) {
  
  meanSpeciesRichnessUpperCI <- numeric()
  meanSpeciesRichnessLowerCI <- numeric()
  
  if(is.null(which.result)) which.result = length(runs$Output) 
  
  if(is.null(which.result) == FALSE){
    if(which.result == "all"){
      simulations <- c(1:length(runs$Output))
    } else {
      simulations <- which.result
    } 
  }
  
  colfunc <- colorRampPalette(c("blue", "red"))
  cols <- colfunc(length(simulations))
  
  # Determine plot title
  plot_title <- if (!is.null(title)) {
    title
  } else if (!is.null(runs$Model$getName)) {
    runs$Model$getName
  } else {
    getNames(runs)$Model$getName
  }
  
  for(t in simulations){
    
    simu_t <- runs$Output[[t]]
    matrix <- simu_t$specMat
    
    if(length(simulations)==1) t<-1
    
    landscapeDim = dim(matrix)
    landscapeArea = landscapeDim[1] * landscapeDim[2]
    if (is.null(size)){
      size = as.integer(seq(2, sqrt(landscapeArea),len = 10))
    }
    
    meanSpeciesRichness <- numeric()
    
    subPlots <- localPlots(size=size, n = rep, simu=runs, which.result = t, nested = nested)$subPlots
    
    speciesRichness <- sapply(subPlots, SR)
    speciesRichnesslist <-list()
    
    f <- rep(1:length(size), length(speciesRichness)/length(size))
    
    for(i in 1: length(speciesRichness)){
      if(i <= length(size)){
        speciesRichnesslist[[f[i]]] <- speciesRichness[i]
      }
      else speciesRichnesslist[[f[i]]] <- c(speciesRichnesslist[[f[i]]],speciesRichness[i])
    }
    
    for(i in 1:length(speciesRichnesslist)){
      meanSpeciesRichness[i] <- mean(speciesRichnesslist[[i]])
      
      if(t == 1){
        meanSpeciesRichnessUpperCI[i] <- quantile(speciesRichnesslist[[i]], probs = 0.95, na.rm = T)
        meanSpeciesRichnessLowerCI[i] <- quantile(speciesRichnesslist[[i]], probs = 0.05, na.rm = T)
      } else{
        if(meanSpeciesRichnessUpperCI[i]< quantile(speciesRichnesslist[[i]], probs = 0.95, na.rm = T)){
          meanSpeciesRichnessUpperCI[i] <- quantile(speciesRichnesslist[[i]], probs = 0.95, na.rm = T)
        }
        if(meanSpeciesRichnessLowerCI[i] > quantile(speciesRichnesslist[[i]], probs = 0.05, na.rm = T)){
          meanSpeciesRichnessLowerCI[i] <- quantile(speciesRichnesslist[[i]], probs = 0.05, na.rm = T)
        }
      }
    }
    
    if(plot == TRUE){
      if(t == 1){
        if(length(simulations)==1){
          plot(size^2, meanSpeciesRichness, log="xy", xlab="Area (n cells)",
               ylab="Number of Species", main=plot_title, col=t, lwd=2, pch=4, type="b")
          
          polygon(x=c(size^2,rev(size^2)), y=c(meanSpeciesRichnessUpperCI,
                                               rev(meanSpeciesRichnessLowerCI)), col="#00000030", border=NA)
          
          lines(size^2, (meanSpeciesRichnessUpperCI), col="red", lty=2, lwd=2)
          lines(size^2, (meanSpeciesRichnessLowerCI), col="red", lty=2, lwd=2)
          
        } else {
          plot(size^2, meanSpeciesRichness, type="l", log="xy", xlab="Area (n cells)", 
               ylab="Number of Species", main=plot_title, col=cols[t])
        }
      } else {
        lines(size^2, meanSpeciesRichness, type="l", col=cols[t])
        
        if(t == length(runs$Output)){
          
          polygon(x=c(size^2,rev(size^2)), y=c(meanSpeciesRichnessUpperCI,
                                               rev(meanSpeciesRichnessLowerCI)), col="#00000030", border=NA)
          
          lines(size^2, (meanSpeciesRichnessUpperCI), col="red", lty=2, lwd=2)
          lines(size^2, (meanSpeciesRichnessLowerCI), col="red", lty=2, lwd=2)
        }
      }
    }
  }
  
  if(length(simulations)==1) {
    return(data.frame(size = size^2, sr.Mean = meanSpeciesRichness,
                      sr.UpperCI = meanSpeciesRichnessUpperCI, sr.LowerCI = meanSpeciesRichnessLowerCI))
  }
}

#' @rdname sac
#' @method sac PhylosimList
#' @export
sac.PhylosimList <- function(runs, which.result = NULL, size = NULL, rep = 50, plot = TRUE, nested = FALSE, title = NULL) {
  if (!is.null(title) && length(title) != length(runs)) {
    warning("Length of title vector does not match length of runs. Using automatic titles.")
    title <- NULL
  }
  
  results <- lapply(seq_along(runs), function(i) {
    current_title <- if (!is.null(title)) title[i] else NULL
    sac(runs[[i]], which.result = which.result, size = size, rep = rep, 
        plot = plot, nested = nested, title = current_title)
  })
  
  names(results) <- names(runs)
  return(results)
}

## Internal function called by SAC
SR <- function(matrix){  
  sr <- length(unique(c(matrix)))
  return(speciesRichness = sr)
}