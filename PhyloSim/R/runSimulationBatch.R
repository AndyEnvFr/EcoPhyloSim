#' @title Batch runner
#' @description A model of species community assembly under different assembly mechanisms, using parallel computing to make use of multi core CPUs and clusters in order to reduce computation time. The function is an extension to \link{runSimulation} in order to accelerate the simulation of multiple scenarios.
#' @param pars A list of parameter sets as created by \code{\link{createCompletePar}} 
#' @param parallel Integer, determining the number of cores used for parallel computing. If parallel = "auto", n-1 cores are used.
#' @param backup Logical, determining whether the results of the individual scenario runs should be saved as a workspace image (advised if the simulation takes a long time, or if the individual scenarios vary greatly in runtime. Default is FALSE)
#' @param strip NULL or "summaries", determining whether the whole simulation or only the summary statistics should be returned.
#' @return An object of class "PhylosimList".
#' @details The "PhylosimList" object is a list of "Phylosim" objects. They can be accessed by indexing (see Example).\cr\cr This function uses the \code{\link{foreach}} and \code{\link{doParallel}} package to compute the model scenarios parallel on several cores. \cr\cr The phylogeny is passed to R in the newick format and parsed to an object of class "phylo" with the function \code{\link[ape]{read.tree}} from the \code{\link{ape}} package. 
#' @example /inst/examples/runSimulationBatch-help.R
#' @export
#' 
runSimulationBatch <- function(pars, parallel = FALSE, backup = FALSE, strip = NULL){
  #start timing
  ptm <- proc.time() 
  library(foreach)
  
  # TODO getParametersXML(XMLfile)
  
  if (parallel != FALSE){
    cat("running", length(pars), "batch simulations with parallelization\n")
    
    if (parallel == TRUE | parallel == "auto") cores <- parallel::detectCores() - 1
    if (is.numeric(parallel)) cores <- parallel
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    
    out <- foreach(i=1:length(pars), .packages = c("PhyloSim")) %dopar%{
      
      if(backup == TRUE){
        name <- paste(pars[[i]]$scenario, "_params.RData", sep="")
        par <- pars[[i]]
        save(par, file = name)
      }
      
      # TODO: remove this after test
      OUT <- tryCatch(runSimulation(pars[[i]]),
                      error = function(e) NULL)
      
      if(backup == TRUE){
        name <- paste(pars[[i]]$scenario, "_out.RData", sep="")
        save(OUT, file = name)
      }
      if(!is.null(strip) && !is.null(OUT)) {
        if(strip == "summaries") {
          OUT <- OUT$Output[[1]]$summaries
          class(OUT) <- "list"
        }
        else stop("Phylosim::runSimulationBatch unrecognized argument to strip")
      }
      OUT
    }
    parallel::stopCluster(cl)
  } else {
    cat("running", length(pars), "batch simulations without parallelization\n")
    out <- foreach(i=1:length(pars), .packages = c("PhyloSim")) %do%{
      
      cat("running parameter", i, "\n")
      
      OUT <- runSimulation(pars[[i]])
      
      if(backup == TRUE){
        name <- paste(pars[[i]]$scenario, ".RData", sep="")
        save(OUT, file = name)
      }
      if(!is.null(strip)) {
        if(strip == "summaries") {
          OUT <- OUT$Output[[1]]$summaries
          class(OUT) <- "list"
        }
        else stop("Phylosim::runSimulationBatch unrecognized argument to strip")
      }
      OUT
    }
  }
  
  # Fix the naming assignment
  for (i in 1:length(pars)) names(out)[i] <- pars[[i]]$scenario
  
  if(is.null(strip)) {
    class(out) <- "PhylosimList"
  } else {
    class(out) <- "list"
  }
  
  # Fix timing output
  time <- proc.time() - ptm
  cat("Finished after", floor(time[3]/60), "minute(s) and", round(time[3]%%60, 2), "second(s).\n")
  
  return(out)
}
