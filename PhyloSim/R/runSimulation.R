#' @title  Species Community Simulation
#' @description A model of species community assembly under different assembly mechanisms.  
#' @param par, a list of model parameters created with \link{createCompletePar}
#' @return An object of class "Phylosim". This objet contains the species matrix, the trait matrix, the environmental matrix, the competition matrix and the neutral matrix, as well as the phlogeny and the parameter settings of the simulation. 
#' @details If your parameterset contains more than one runs argument, each interim step is saved in the Phylosim object. \cr\cr For larger simularions consider \link{runSimulationBatch} to make use of parallel computing. \cr\cr If you are using `type="Rneutral"` only one runs argument can be processed.\cr\cr It is possible that more than one new species arises per generation. This leads to a multifurcated phylogeny, yet many packages such as "ape" can only work with bifurcated tree. Setting `converToBinaryTree = TRUE` converts the generated phylogeny to a bifurcate one using `multi2di()` from the "ape" package. 
#' @importFrom adephylo distTips
#' @import Rcpp
#' @useDynLib PhyloSim, .registration = TRUE
#' @export
#' @example /inst/examples/runSimulation-help.R

runSimulation <- function(par)
{
  
  # GENERAL SETTINGS
  
  if (par$type == "base"){
    
    if (par$dispersal == "global" | par$dispersal == 0){
      dispersal = 1
      dispersalCut = 1
    }else if (is.numeric(par$dispersal)){
      if (par$dispersal < 0.5) stop("dispersal parameter too small") 
      dispersalCut = 2*par$dispersal
      dispersal = 3
    }else stop("wrong dispersal parameter") 
    
    if(par$fitnessActsOn == "mortality"){
      reproductiveFitness = F
      mortalityFitness = T    
    }else if(par$fitnessActsOn == "reproduction"){
      reproductiveFitness = T
      mortalityFitness = F    
    }else if (par$fitnessActsOn == "both"){
      reproductiveFitness = T
      mortalityFitness = T    
    }else stop("wrong fitness parameters")
    
    
    if (is.numeric(par$density)) {
      if (par$density == 0) {
        par$density = F
        par$compStrength = 1
      } else{
        par$compStrength = par$density
        par$density = T
      }
    }
    
    
    if (is.numeric(par$positiveDensity)) {
      if (par$positiveDensity == 0) {
        par$positiveDensity = F
        par$pDDStrength = 1
      } else{
        par$pDDStrength = par$positiveDensity
        par$positiveDensity = T
      }
    }
    
    if (is.numeric(par$environment)) {
      if (par$environment > 1 || par$environment < 0) stop("Parameter environment must be between 0 and 1")
      if (par$environment == 0) {
        par$environment = F
        par$envStrength = 1
      } else {
        par$envStrength = par$environment
        par$environment = T
      }
    }

    if (!is.null(par$nicheWidth) && par$nicheWidth <= 0) {
      stop("nicheWidth must be a positive number (used as standard deviation in the environmental kernel)")
    }    

    if (!is.null(par$densityNicheWidth) && par$densityNicheWidth <= 0) {
      stop("densityNicheWidth must be a positive number (used as standard deviation in the density kernel)")
    }    

    if (!is.null(par$positiveDensityDependenceNicheWidth) && par$positiveDensityDependenceNicheWidth <= 0) {
      stop("positiveDensityDependenceNicheWidth must be a positive number (used as standard deviation in the positive density dependence kernel)")
    }    
    
    if(par$density == 0 & par$environment == 0 & par$positiveDensity == 0)
    {
      neutral = TRUE
    }else{
      neutral = F
    }    
    
    if (par$fitnessBaseMortalityRatio < 1) stop("Parameter fitnessBaseMortalityRation must be greater than or equal to 1")
    
    if (max(par$airmat) > 1 || min(par$airmat) < 0) stop("Values of airmat must be between 0 and 1")
    
    ptm <- proc.time()
    
    out <- callModel( x = par$x, 
                      y = par$y, 
                      dispersal = dispersal, 
                      runs = round(par$runs), 
                      specRate = par$specRate, 
                      dens = par$density, 
                      positiveDens = par$positiveDensity, 
                      env = par$environment, 
                      neutral = neutral, 
                      mort = mortalityFitness, 
                      mortStrength = par$fitnessBaseMortalityRatio, 
                      repro = reproductiveFitness, 
                      dispersalCutoff = dispersalCut, 
                      densityCutoff = par$densityCut, 
                      seed = par$seed, 
                      envStrength = par$envStrength, 
                      compStrength = par$compStrength, 
                      pDDStrength = par$pDDStrength, 
                      fission = par$fission, 
                      redQueen = par$redQueen, 
                      redQueenStrength = par$redQueenStrength, 
                      protracted = par$protracted, 
                      airmatR = par$airmat, 
                      soilmatR = par$soilmat,
                      prunePhylogeny = par$prunePhylogeny,
                      nicheWidth = par$nicheWidth,
                      densityNicheWidth = par$densityNicheWidth,
                      pDDNicheWidth = par$positiveDensityDependenceNicheWidth
                      )  
    
    runtime <- as.numeric((proc.time() - ptm)[3])
    
    print (paste("Core simulation finished after",floor(((proc.time() - ptm)[3])/60), "minute(s) and", ((proc.time() - ptm)[3])%%60, "second(s). Converting data"))
    
    #     Rcpp::List listResults = Rcpp::List::create(Rcpp::Named("Species") = specOut,
    #                                                 Rcpp::Named("EnvTrait") = traitOut,
    #                                                 Rcpp::Named("NeutralTrait") = neutralOut,
    #                                                 Rcpp::Named("CompetitionTrait") = compOut,
    #                                                 Rcpp::Named("Environment") = envOut,
    #                                                 Rcpp::Named("Phylogeny") = phyloOut);
    
    # The following code reshuffles the output to the structure 
    # of a phylosim class
    output<-list()
    for (i in 1:length(par$runs)){
      
      specMati = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y)
      if(length(unique(c(specMati)))==1){
        phyloi <- 0
        warning("Cannot build Phylogeny")
      }else{
        # TODO: ape can not read the tree if it is not pruned by prunePhylogenyR:
        # Error in if (tp[3] != "") obj$node.label <- tp[3] : missing value where TRUE/FALSE needed 
        phyloi <- ape::read.tree(text= out[[i]]$Phylogeny)
      }
       
      output$Output[[i]] = list(
        specMat = matrix(out[[i]]$Species,ncol=par$x, nrow=par$y), 
        traitMat= matrix(out[[i]]$EnvTrait,ncol=par$x, nrow=par$y), 
        envMat = matrix(out[[i]]$Environment,ncol=par$x, nrow=par$y), 
        compMat = matrix(out[[i]]$CompetitionTrait,ncol=par$x, nrow=par$y), 
        neutMat = matrix(out[[i]]$NeutralTrait,ncol=par$x, nrow=par$y),
        phylogeny = if (par$convertToBinaryTree && !is.double(phyloi)) ape::multi2di(phyloi) else phyloi, 
        phyloTXT = out[[i]]$Phylogeny,
        summaries = NA)
    }
    cat("done! \n")
    
    par$runtime <- runtime
    
    output$Model <- par
    
    class(output) <- "PhyloSim"
    
    if(par$calculateSummaries) {
      output$Output[[1]]$summaries = calculateSummaryStatistics(output)
    }
    
    return(output)
    
  ##################################################################
  # The following two options are other models that are implemented in Phylosiom
  # Leipzig class C code similar to Chave 2002 - probably currently not functional
  # Rneutral calls R code for a neutralm model similar to Chave 2002 that was implemented for 
  # Testing purposes 
  
  }else if(par$type == "Leipzig"){
    
    # TODO ... NOT FUNCTIONAL YET !!!
    out <- .C(callLeipzig, as.integer(x),as.integer(y), as.integer(dispersal), as.integer(nSpec),as.integer(specRate), as.integer(runs), as.logical(density), as.numeric(densityStrength), specOut = as.integer(outVec), densOut = as.numeric(outVec))[9:10]
    
    specMat = matrix(out[[1]],ncol=x, nrow=y)
    densMat = matrix(out[[2]],ncol=x, nrow=y)
    print (paste("Finished after",floor(((proc.time() - ptm)[3])/60), "minute(s) and", ((proc.time() - ptm)[3])%%60, "second(s)."))
    return(list(specMat = specMat, densMat=densMat))
    
    ###################################################################
  }else if (par$type == "Rneutral"){
    ptm <- proc.time()
    result = NeutralMod(xdim = par$x, ydim=par$y, specRate = par$specRate,  seed = par$seed, runs = par$runs)
    runtime <- as.numeric((proc.time() - ptm)[3])
    par$runtime <- runtime
    out = list(Output = list(list("specMat"=result)), Model = par)
    class(out) <- "PhyloSim"
    return(out)
  }
}


# neutral = T
# x = 50
# y = 50
# dispersal = 0.5
# runs = c(100,200)
# specRate = 1.0
# density = F
# environment = F
# seed = 1
# fitnessBaseMortalityRatio = 5
# fitnessActsOn = "mortality" 
# densityCut = 1
# saveTimes = "last"
# 
# 
# if (dispersal == "global" | dispersal == 0){
#   dispersal = 1
#   dispersalCut = 1
# }else if (is.numeric(dispersal)){
#   if (dispersal < 0.5) stop("dispersal parameter too small") 
#   dispersalCut = 2*dispersal
#   dispersal = 3
# }else stop("wrong dispersal parameter") 
# 
# if(fitnessActsOn == "mortality"){
#   reproductiveFitness = F
#   mortalityFitness = T    
# }else if(fitnessActsOn == "reproduction"){
#   reproductiveFitness = T
#   mortalityFitness = F    
# }else if (fitnessActsOn == "both"){
#   reproductiveFitness = T
#   mortalityFitness = T    
# }else stop("wrong fitness parameters")
# 
# if(density == FALSE & environment == FALSE)
# {
#   neutral = TRUE
# }else{
#   neutral = F
# }
# if (length(runs)>1){
#   if (any(runs[-length(runs)] > runs[-1])) stop( "wrong argument to runs")
# }
# 
# 
# out <- callModel( x,  y,  dispersal,  runs,  specRate, density, environment, neutral, mortalityFitness, fitnessBaseMortalityRatio, reproductiveFitness, dispersalCut, densityCut, seed)  


