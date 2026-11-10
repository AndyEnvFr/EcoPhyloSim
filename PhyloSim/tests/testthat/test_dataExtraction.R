context("data extraction functions")

set.seed(1)
par <- createCompletePar(x = 20, y = 15, dispersal = 1, runs = c(10, 11, 12), # non-square grid
                         density = 1, densityCut = 2, scenario = "test")
simu <- runSimulation(par)

test_that("getCircularOffsets matches the C++ neighbourhood", {
  expect_equal(nrow(getCircularOffsets(1)), 4)
  expect_equal(nrow(getCircularOffsets(3)), 28)
  expect_equal(nrow(getCircularOffsets(0)), 0)
})

test_that("getTorusMatrix wraps around, also for r = 0 and r > size", {
  m <- matrix(1:6, 2, 3)
  expect_identical(PhyloSim:::getTorusMatrix(m, 0), m)
  big <- PhyloSim:::getTorusMatrix(m, 1)
  expect_identical(dim(big), c(4L, 5L))
  expect_identical(big[1, 1], m[2, 3])
  expect_identical(dim(PhyloSim:::getTorusMatrix(m, 4)), c(10L, 11L))
})

test_that("getConNeigh counts conspecific neighbours on the torus", {
  res <- getConNeigh(simu)
  sp <- simu$Output[[1]]$specMat
  off <- getCircularOffsets(2)
  ref <- matrix(0, nrow(sp), ncol(sp))
  for (i in seq_len(nrow(sp))) for (j in seq_len(ncol(sp))) for (k in seq_len(nrow(off))) {
    ref[i, j] <- ref[i, j] + (sp[(i - 1 + off$dx[k]) %% nrow(sp) + 1,
                                 (j - 1 + off$dy[k]) %% ncol(sp) + 1] == sp[i, j])
  }
  expect_equal(res$Output[[1]]$conNeighMat, ref, check.attributes = FALSE)
  expect_identical(res$Output[[1]]$specMat, sp)
  expect_identical(res$Model$getName, getNames(simu))
})

test_that("getMortality, getID and getMatToTab are consistent", {
  res <- getConNeigh(simu)
  expect_identical(names(res$Output), c("10", "11", "12"))
  expect_true(all(is.na(res$Output[[1]]$mortMat)))
  changed <- res$Output[[2]]$traitMat != res$Output[[1]]$traitMat
  expect_identical(res$Output[[2]]$mortMat, changed)
  expect_true(any(changed) && any(!changed))
  expect_equal(res$Output[[2]]$idMat[!changed], res$Output[[1]]$idMat[!changed])
  expect_false(any(res$Output[[2]]$idMat[changed] %in% res$Output[[1]]$idMat))
  tab <- getMatToTab(res)
  expect_equal(nrow(tab), 2 * 300)
  expect_identical(dim(res$Output[[1]]$idMat), dim(res$Output[[1]]$specMat))
  expect_identical(dim(res$Output[[1]]$mortMat), dim(res$Output[[1]]$specMat))
  expect_identical(tab$mortNextGen[1:300], as.vector(res$Output[[2]]$mortMat))
  expect_error(getMatToTab(simu), "getConNeigh")
})

test_that("functions work with one recorded generation and with lists", {
  one <- runSimulation(createCompletePar(x = 10, y = 10, dispersal = 1, runs = 5,
                                         density = 1, scenario = "one"))
  expect_equal(nrow(getMatToTab(getConNeigh(one))), 0)
  lst <- createPhylosimList(list(simu, one))
  expect_length(getNames(lst), 2)
  expect_s3_class(getConNeigh(lst), "PhylosimList")
  expect_length(getSpecTime(lst), 2)
})

test_that("getNames works if optional parameters are missing", {
  s <- simu
  s$Model$fission <- NULL
  s$Model$redQueen <- NULL
  expect_identical(getNames(s), getNames(simu))
})
