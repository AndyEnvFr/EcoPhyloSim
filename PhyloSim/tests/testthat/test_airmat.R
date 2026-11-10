# fixes #122 (airmat was read with the wrong index when x != y)
test_that("airmat is used cell by cell also for non-square grids", {
  for (xy in list(c(20, 30), c(30, 20))) {
    x <- xy[1]; y <- xy[2]
    air <- matrix(seq(0.01, 0.99, length.out = x * y), nrow = y, ncol = x)
    p <- createCompletePar(x = x, y = y, runs = 100, environment = 0.5, airmat = air,
                             specRate = 1, convertToBinaryTree = FALSE, seed = 1)
    env <- suppressWarnings(runSimulation(p))$Output[[1]]$envMat
    expect_equal(env, air)
  }
})
