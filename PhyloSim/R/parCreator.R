#' @title  Parameter Generator 
#' @description Function to create a list of parameters for biogeographical simulations with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}
#' @param x  Integer, Dimension of the model landscape in x-direction
#' @param y  Integer, Dimension of the model landscape in y-direction
#' @param dispersal Integer. Type 0 or "global" for global dispersion. For local dispersion all integers >=1 set the dispersal distance.
#' @param runs  Integer or vector of Integers, Number of generations or sequence of generations the model runs over (see Details). 
#' @param specRate Integer, Number of Individuals introduced to the community in each generation
#' @param density Float, determining whether or how strong the density dependence influences the model. By default (density=0) there is no density dependence. The higher the value of the parameter, the stronger is the density dependence.
#' @param environment Float, determining whether or how strong the environment influences the model.  By default (environment=0) there is no influence of the environment. The higher the value of the parameter, the stronger is the influence of the environment. Environment value must be between 0 and 1.
#' @param fitnessActsOn Character, determining how the fitness influences the individuals. Possible inputs are "mortality" (default), "reproduction" or "both"
#' @param fitnessBaseMortalityRatio Integer, determines the fitness based mortality ratio. Must be greater than or equal to 1.
#' @param densityCut Integer, defines the effective range of the competition (ignored if density = FALSE)
#' @param fission Integer, determining which fission type should be used. Options are 0 (none = default), 1 (every second individual becomes part of new species) and 2 (population is geographically split in two parts).
#' @param protracted Integer, determining the time span in generations a new species stays 'incipient' before turning into a 'good' species. Default is 0.
#' @param redQueenStrength Float, determining the strength of the Red Queen effect. A value > 0 mean a new species gets a fitness boost due to its novelty.
#' @param redQueen Float, determining the strength of the fitness decline of an aging species. 
#' @param airmat Matrix, deteriming the environment of the simulation. airmat needs to be a matrix with the same dimensions as the grid. Must be scaled between 0 and 1.
#' @param seed numerical, sets the random seed
#' @param type Character, determining which model should be used. "base" is running the default model. Other possibilities are "Leipzig" and "Rneutral" which will run a neutral model purely in R.
#' @param scenario String, further information you want to add to the parameter set in order to refer to a model run more conveniently. 
#' @param calculateSummaries Logical, determining wheter summary statistics should be calculated
#' @param convertToBinaryTree Logical, determining if the phylogeny should be converted into a binary tree
#' @param prunePhylogeny Logical, determining whether the phylogeny should be prune by the internal pruning function
#' @details If runs is a sequence of generations, intermediate results are saved. E.g. when runs is c(500, 600, 700), the simulation runs 700 generations in total, and the intermediate results at generations 500 and 600 are saved additionally. The intermediate and end results are saved in the output of \code{runSimulation}. \cr\cr If type = "Rneutral" the model will run entirely in R. This model is to be seen only for test and teaching purpose. To be used in practice it is far too slow. Also the output is reduced. Only the species landscape and the parameter settings will be displayed in the output.
#' @return A List with parameters
#' @example /inst/examples/parCreator-help.R
#' @export

createCompletePar <- function(x = 50, y = 50, dispersal = "global", runs = 100, specRate = 1.0, density = 0, environment = 0, fitnessActsOn = "mortality" , fitnessBaseMortalityRatio = 10, densityCut = 1, seed=NULL,  type = "base", fission = 0, redQueen = 0, redQueenStrength = 0, protracted = 0, airmat = 1, scenario = NULL, calculateSummaries = FALSE, convertToBinaryTree = TRUE, prunePhylogeny = TRUE){
    
  soilmat <- 1 # Needs to be defined here as long as it is not
               # implemented in the model
  if (length(runs)>1){
    if (any(runs[-length(runs)] > runs[-1])) stop("When passing a vector for runs, the last element must be the largest.")
  }
  
  if(length(airmat)!= 1){
    if((nrow(airmat) != y) | (ncol(airmat)!= x)) stop("Environment and matrix size do not match")
  }
  if (max(airmat) > 1 || min(airmat) < 0) stop("Values of airmat must be between 0 and 1")
  
  if(length(soilmat) != 1){
    if((nrow(soilmat) != y) | (ncol(soilmat)!= x)) stop("Environment and matrix size do not match")
  }
  
  if (environment > 1 || environment < 0) stop("Parameter environment must be between 0 and 1")
  
  if (fitnessBaseMortalityRatio < 1) stop("Parameter fitnessBaseMortalityRation must be greater than or equal to 1")
  
  airmat <- as.numeric(airmat)
  
  # The soilmat is not implemented in the model so far.
  soilmat = 1
  soilmat <- as.numeric(soilmat)
  
  if (is.null(seed)) seed = sample(1:10000,1)
  
  par = list(x=x,y=y,dispersal = dispersal, runs = runs, specRate = specRate, 
             density = density, environment = environment, fitnessActsOn=fitnessActsOn,
             fitnessBaseMortalityRatio=fitnessBaseMortalityRatio, densityCut = densityCut, 
             seed = seed, type = type, scenario = scenario, fission = fission, 
             redQueen = redQueen, redQueenStrength = redQueenStrength, protracted = protracted,
             airmatR = airmat, soilmatR=soilmat, calculateSummaries=calculateSummaries,
             convertToBinaryTree = convertToBinaryTree, prunePhylogeny = prunePhylogeny)
  
  
  if (par$fitnessActsOn %in% c("both", "reproduction") && par$environment >= 0) {
    stop(paste0("Incompatible settings: fitnessActsOn = 'both' and environment = ", par$environment,
                ". This will result in R crashing. Either set environment = 0 or set fitnessActsOn back to default."))
  }
  
  


  return(par)
  
}

