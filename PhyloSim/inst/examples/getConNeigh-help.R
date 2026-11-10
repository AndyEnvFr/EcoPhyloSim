## Small simulation with four recorded generations
par <- createCompletePar(x = 20, y = 20, dispersal = 1, runs = c(10, 20, 30, 40),
                         density = 1, densityCut = 2, scenario = "example")
simu <- runSimulation(par)

getNames(simu)                  # parameter label of the run
getSpecTime(simu)               # species richness of each recorded generation

simu <- getMortality(simu)      # adds mortMat: replaced since previous recorded generation?
simu <- getID(simu)             # adds idMat: individual IDs across generations
simu <- getConNeigh(simu)       # adds conNeighMat: number of conspecific neighbours
head(getMatToTab(simu))         # long table for statistical analysis
