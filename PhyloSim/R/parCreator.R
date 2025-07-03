#' @title  Parameter Generator
#' @description Function to create a list of parameters for biogeographical simulations with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}
#'
#' @param x Integer, dimension of the model landscape in x-direction
#' @param y Integer, dimension of the model landscape in y-direction
#' @param dispersal Integer. Type 0 or "global" for global dispersion. For local dispersion, all integers >=1 set the dispersal distance.
#' @param runs Integer or vector of Integers, number of generations or sequence of generations the model runs over (see Details).
#' @param specRate Integer, number of individuals introduced to the community in each generation
#'
#' @param negativeDensity,positiveDensity,environment Float, defines the strength of the respective ecological process. By default all are 0 (no effect). Higher values increase process strength. See Details for process-specific information.
#' @param nDDNicheWidth,pDDNicheWidth,envNicheWidth Double, width (σ) of the Gaussian fitness kernel for the respective ecological process. Smaller values imply stronger trait-specific filtering. See Details for defaults and process-specific effects.
#' @param nDensityCut,pDensityCut Integer, define the spatial cutoff radius for negative and positive density dependence (i.e., the local neighborhood over which competition or facilitation is calculated). Only used if the respective density dependence process is active.
#'
#' @param fitnessActsOn Character, determining how the fitness influences the individuals. Possible inputs are "mortality" (default), "reproduction" or "both"
#' @param fitnessBaseMortalityRatio Integer, determines the fitness-based mortality ratio. Must be greater than or equal to 1.
#'
#' @param fission Integer, determining which fission type should be used. Options are 0 (none = default), 1 (every second individual becomes part of new species), and 2 (population is geographically split in two parts).
#' @param protracted Integer, determining the time span in generations a new species stays 'incipient' before turning into a 'good' species. Default is 0.
#' @param redQueenStrength Float, determining the strength of the Red Queen effect. A value > 0 means a new species gets a fitness boost due to its novelty.
#' @param redQueen Float, determining the strength of the fitness decline of an aging species.
#' @param airmat Matrix, determining the environment of the simulation. Must have the same dimensions as the landscape grid and values scaled between 0 and 1.
#' @param seed Numeric, sets the random seed
#' @param type Character, determining which model should be used. "base" is running the default model. Other possibilities are "Leipzig" and "Rneutral" which will run a neutral model purely in R.
#' @param scenario String, additional identifier to label a specific parameter set / simulation run.
#' @param calculateSummaries Logical, determining whether summary statistics should be calculated
#' @param convertToBinaryTree Logical, determining if the phylogeny should be converted into a binary tree
#' @param prunePhylogeny Logical, determining whether the phylogeny should be pruned by the internal pruning function
#'
#' @details
#' If \code{runs} is a sequence of generations, intermediate results are saved. E.g., when \code{runs = c(500, 600, 700)}, the simulation runs 700 generations in total, and results at generations 500 and 600 are stored as checkpoints.
#'
#' The model incorporates three main ecological processes, each controlled by a strength parameter and a corresponding niche width parameter:
#' \itemize{
#'   \item Negative Density Dependence (Competition): Controlled by \code{negativeDensity}, \code{nDDNicheWidth}, and \code{nDensityCut}. Higher \code{negativeDensity} increases competitive pressure. Smaller \code{nDDNicheWidth} results in more trait-specific filtering. \code{nDensityCut} defines the spatial scale of competition.
#'   \item Positive Density Dependence (Facilitation): Controlled by \code{positiveDensity}, \code{pDDNicheWidth}, and \code{pDensityCut}. Higher \code{positiveDensity} increases facilitative interactions. Smaller \code{pDDNicheWidth} increases trait specificity. \code{pDensityCut} defines the spatial scale of facilitation.
#'   \item Environmental Selection: Controlled by \code{environment} and \code{envNicheWidth}. Higher \code{environment} increases environmental filtering. Smaller \code{envNicheWidth} leads to narrower environmental tolerances.
#' }
#'
#' If \code{type = "Rneutral"} the model will run entirely in R. This mode is for testing/teaching only and is significantly slower. The output is reduced and includes only the species landscape and parameter set.
#'
#' @return A list of parameter settings to be used with \code{\link{runSimulation}} or \code{\link{runSimulationBatch}}.
#'
#' @example /inst/examples/parCreator-help.R
#' @export


# changed variables names, when implemented positive density dependence [Andy]
# dens -> negativeDens
# density -> negativeDensity
# compStrength -> nDDStrength
# nicheWidth -> envNicheWidth

createCompletePar <- function(x = 50, y = 50, type = "base", dispersal = "global", runs = 100, specRate = 1.0, negativeDensity = 0, nDensityCut = 1, nDDNicheWidth = 0.1, positiveDensity = 0, pDensityCut = 1, pDDNicheWidth = 0.1, environment = 0, envNicheWidth = 0.03659906, fitnessActsOn = "mortality", fitnessBaseMortalityRatio = 10, seed = NULL, fission = 0, redQueen = 0, redQueenStrength = 0, protracted = 0, airmat = 1, scenario = NULL, calculateSummaries = FALSE, convertToBinaryTree = TRUE, prunePhylogeny = TRUE) {
  soilmat <- as.numeric(1) # not implemented yet
  airmat <- as.numeric(1) # not implemented yet

  if (length(runs) > 1) {
    if (any(runs[-length(runs)] > runs[-1])) stop("When passing a vector for runs, the last element must be the largest.")
  }

  if (length(airmat) != 1) {
    if ((nrow(airmat) != y) | (ncol(airmat) != x)) stop("Environment and matrix size do not match")
  }
  if (max(airmat) > 1 || min(airmat) < 0) stop("Values of airmat must be between 0 and 1")

  if (length(soilmat) != 1) {
    if ((nrow(soilmat) != y) | (ncol(soilmat) != x)) stop("Environment and matrix size do not match")
  }

  if (environment > 1 || environment < 0) stop("Parameter environment must be between 0 and 1")

  if (fitnessBaseMortalityRatio < 1) stop("Parameter fitnessBaseMortalityRation must be greater than or equal to 1")


  if (is.null(seed)) seed <- sample(1:10000, 1)

  par <- list(
    x = x, y = y, dispersal = dispersal, runs = runs, specRate = specRate,
    negativeDensity = negativeDensity, nDDNicheWidth = nDDNicheWidth, nDensityCut = nDensityCut,
    positiveDensity = positiveDensity, pDDNicheWidth = pDDNicheWidth, pDensityCut = pDensityCut,
    environment = environment, envNicheWidth = envNicheWidth, fitnessActsOn = fitnessActsOn,
    fitnessBaseMortalityRatio = fitnessBaseMortalityRatio,
    seed = seed, type = type, scenario = scenario, fission = fission,
    redQueen = redQueen, redQueenStrength = redQueenStrength, protracted = protracted,
    airmat = airmat, soilmat = soilmat, calculateSummaries = calculateSummaries,
    convertToBinaryTree = convertToBinaryTree, prunePhylogeny = prunePhylogeny
  )



  return(par)
}
