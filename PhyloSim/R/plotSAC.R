
#' @title Species Area Curve
#' @description Plots the species area curve for a given community. 
#' @param simu Simulation output of the class "Phylosim", usually consisting out of several lists. Needs at least the spatial distribution of the species stored in a matrix ($specMat)
#' @param which.result Integer, determines which result should be used. This argument is only usefull if your 'runs' argument in \code{\link{createCompletePar}} contains more than one element. By default (NULL), the last result is used. If you choose "all" all results are shown in one plot (see 'Details).
#' @param size a single value or a vector determining the edge length of the subplots. If not provided, 10 logarithmic plotsizes will be used.
#' @param rep The number of repetitions per size to calculate the mean 
#' @param plot Logical determining whether to plot the SAC or not
#' @param nested Logical, determining whether nested subplots should be used.
#' @param title String, determining the title of the plot
#' @details Displays the accumulated species richness as a function of plot size or the amount of equally sized plots. It serves as an indicator for the clustering of a species community. A positively bent curve usually indicates clustering since an increase in plot size or number leads to an increase in species richness while a negatively bent curve indicates a more neutral distribution of species within the community. \cr\cr If which.result = "all" all intermediate results are shown in one plot. The colors of the lines are plotted as a gradient from blue (first results) to red (end result).
#' @return A list containing the mean species richness for each size and the respective standard deviation. If which.result = "all" only the plot will be returned.
#' @example /inst/examples/plotSAC-help.R

#' @export
#' 
## Modification of the SAC function to plot all recorded time steps
sac <- function(simu, which.result = NULL, size = NULL, rep=50, plot=T, title="SAC" , nested = FALSE){
  
  meanSpeciesRichnessUpperCI <- numeric()
  meanSpeciesRichnessLowerCI <- numeric()
  
  if(is.null(which.result)) which.result = length(simu$Output) 
  
  if(is.null(which.result) == FALSE){
    if(which.result == "all"){
      simulations <- c(1:length(simu$Output))
    
    }else{
      simulations<-which.result
    } 
  }
  
  colfunc <- colorRampPalette(c("blue", "red"))
  cols<-colfunc(length(simulations))
  
  for(t in simulations){
    
    simu_t <- simu$Output[[t]]
    matrix <- simu_t$specMat
    
    if(length(simulations)==1) t<-1
    
    landscapeDim = dim(matrix)
    landscapeArea = landscapeDim[1] * landscapeDim[2]
    if (is.null(size)){
      size = unique(as.integer(round(exp(seq(log(2), log(sqrt(landscapeArea)), length.out = 10)))))
    }
    
    
     meanSpeciesRichness <- numeric()

      subPlots <- localPlots(size=size, n = rep, simu=simu, which.result = t, nested = nested)$subPlots

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
      }}
    
    if(plot == T){
      if(t == 1){
        if(length(simulations)==1){
          plot(size^2, meanSpeciesRichness, log="xy", xlab="Area (n cells)",
               ylab="Number of Species", main=title, col=t, lwd=2, pch=4, type="b")
          
          polygon(x=c(size^2,rev(size^2)), y=c(meanSpeciesRichnessUpperCI,
                  rev(meanSpeciesRichnessLowerCI)), col="#00000030", border=NA)
          
          lines(size^2, (meanSpeciesRichnessUpperCI), col="red", lty=2, lwd=2)
          lines(size^2, (meanSpeciesRichnessLowerCI), col="red", lty=2, lwd=2)
          
        }else  plot(size^2,meanSpeciesRichness,type="l", log="xy", xlab="Area (n cells)", ylab="Number of Species", main=title, col=cols[t])
      } else{
        lines(size^2,meanSpeciesRichness,type="l", col=cols[t])
        
        if(t == length(simu$Output)){
         
          polygon(x=c(size^2,rev(size^2)), y=c(meanSpeciesRichnessUpperCI,
                  rev(meanSpeciesRichnessLowerCI)), col="#00000030", border=NA)
          
          lines(size^2, (meanSpeciesRichnessUpperCI), col="red", lty=2, lwd=2)
          lines(size^2, (meanSpeciesRichnessLowerCI), col="red", lty=2, lwd=2)
        }
      }
    }
  }
  if(length(simulations)==1) return(data.frame(size = size^2, sr.Mean = meanSpeciesRichness,
                                    sr.UpperCI = meanSpeciesRichnessUpperCI, sr.LowerCI = meanSpeciesRichnessLowerCI))
}


## Internal function called by SAC
SR <- function(matrix){  
  sr <- length(unique(c(matrix)))
  return(speciesRichness = sr)
}

