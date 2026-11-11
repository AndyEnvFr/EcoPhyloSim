# Local dispersal and density dependence, burn-in of 300 generations,
# then 6 consecutive recorded generations
par <- createCompletePar(x = 30, y = 30, dispersal = 1, density = 1, densityCut = 2,
                         runs = 300:305, seed = 1)
simu <- runSimulation(par)

dat <- getMortalityData(simu)
head(dat)
mean(dat$died) # share of individuals that died within one generation

# Mortality depending on the number of conspecific neighbors
fit <- glm(died ~ con, family = binomial, data = dat)
summary(fit)
