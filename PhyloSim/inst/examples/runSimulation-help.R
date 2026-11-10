 library(PhyloSim)
# Define a parameter set
par <- createCompletePar(x = 50, y = 50, dispersal = 1 , runs = 1000,
        density = 0, environment = 0, specRate = 1, fission = 0, redQueen=0, redQueenStrength=0,
        protracted=0)

# Run the model
simu <- runSimulation(par)

plot(simu)

# Look at the phylogeny (requires package 'ape')
require(ape)

# Get the phylogeny of the last run. In this example this is after 1000 runs.
phylogeny <- simu$Output[[2]]$phylogeny

# Only extant taxa
extantPhylogeny <- drop.fossil(phylogeny)

# Display the results
plot(extantPhylogeny)

#Look at the species area relation
sac(simu, area = c(1,10,100,1000), rep = 100, plot= TRUE)