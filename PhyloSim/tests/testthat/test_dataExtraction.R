context("getMortalityData")

set.seed(1)
par <- createCompletePar(x = 20, y = 15, dispersal = 1, runs = c(10, 11, 13), # non-square grid
                         density = 1, densityCut = 2, scenario = "test")
simu <- runSimulation(par)

# conspecific neighbors as in Landscape::densityUpdate (Grid.cpp), cell (x, y) = mat[y + 1, x + 1]
bruteCon <- function(sp, r) {
  X <- ncol(sp); Y <- nrow(sp); con <- matrix(0, Y, X)
  for (x in 0:(X - 1)) for (y in 0:(Y - 1)) for (dx in -r:r) {
    yl <- floor(sqrt(r^2 - dx^2))
    for (dy in -yl:yl) {
      nx <- (x + dx) %% X; ny <- (y + dy) %% Y
      if (!(nx == x && ny == y)) con[y + 1, x + 1] <- con[y + 1, x + 1] + (sp[y + 1, x + 1] == sp[ny + 1, nx + 1])
    }
  }
  c(con)
}

test_that("getCircularOffsets matches the C++ neighborhood", {
  expect_equal(nrow(PhyloSim:::getCircularOffsets(1)), 4)
  expect_equal(nrow(PhyloSim:::getCircularOffsets(3)), 28)
  expect_equal(nrow(PhyloSim:::getCircularOffsets(0)), 0)
  expect_error(PhyloSim:::getCircularOffsets(1.5))
})

test_that("getMortalityData counts conspecific neighbors on the torus", {
  dat <- getMortalityData(simu)
  expect_equal(nrow(dat), 2 * 20 * 15)
  expect_equal(unique(dat$census), c(10, 11))
  expect_equal(unique(dat$interval), c(1, 2))
  for (i in 1:2) {
    expect_equal(dat$con[dat$census == par$runs[i]], bruteCon(simu$Output[[i]]$specMat, 2))
  }
  expect_equal(getMortalityData(simu, radius = 1)$con[1:300], bruteCon(simu$Output[[1]]$specMat, 1))
  expect_true(all(getMortalityData(simu, radius = 0)$con == 0))
  # many species on a 7 x 4 grid; radius 4 wraps onto the focal cell, which is not counted, as in C++
  set.seed(2)
  mat <- function() matrix(sample(1:3, 28, replace = TRUE), nrow = 4, ncol = 7)
  fake <- structure(list(Model = list(runs = c(1, 2), densityCut = 2),
                         Output = lapply(1:2, function(i) list(specMat = mat(), traitMat = mat(), compMat = mat(), envMat = mat()))),
                    class = "PhyloSim")
  for (r in 1:4) expect_equal(getMortalityData(fake, radius = r)$con, bruteCon(fake$Output[[1]]$specMat, r))
})

test_that("getMortalityData gives deaths and individual IDs", {
  dat <- getMortalityData(simu)
  d1 <- dat[dat$census == 10, ]; d2 <- dat[dat$census == 11, ]
  expect_equal(d1$died, c(simu$Output[[2]]$traitMat != simu$Output[[1]]$traitMat))
  expect_equal(d1$envTrait, c(simu$Output[[1]]$traitMat))
  expect_equal(d1$species, c(simu$Output[[1]]$specMat))
  expect_true(any(d1$died) && any(!d1$died))
  expect_equal(d2$id[!d1$died], d1$id[!d1$died]) # survivors keep their ID
  expect_false(any(d2$id[d1$died] %in% d1$id)) # new individuals get new IDs
  expect_false(anyDuplicated(d2$id) > 0)
})

test_that("getMortalityData works with lists and needs two recorded generations", {
  lst <- getMortalityData(list(a = simu, b = simu))
  expect_equal(names(lst)[1], "run")
  expect_equal(lst$run, rep(c("a", "b"), each = 600))
  expect_equal(unique(getMortalityData(list(simu, simu))$run), 1:2)
  one <- runSimulation(createCompletePar(x = 10, y = 10, runs = 10))
  expect_error(getMortalityData(one), "two recorded generations")
  expect_error(getMortalityData(5), "PhyloSim object")
})
