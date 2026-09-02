test_that("evoParams works with defaults", {
  params <- evoParams(no_sp = 5)
  expect_s3_class(params, "MizerParams")
  expect_true("mizerEvolution" %in% names(getMetadata(params)$extensions))
  expect_true("lineage" %in% names(species_params(params)))
  expect_true("zeta" %in% names(species_params(params)))
  expect_true("pop" %in% names(species_params(params)))
  expect_true("ext" %in% names(species_params(params)))
  expect_equal(nrow(species_params(params)), 5)
})

test_that("evoParams handles deprecated arguments", {
  expect_warning(
    params <- evoParams(no_sp = 4, min_w_inf = 5),
    "deprecated"
  )
  expect_s3_class(params, "MizerParams")

  expect_warning(
    params <- evoParams(no_sp = 4, max_w_inf = 5000),
    "deprecated"
  )
  expect_s3_class(params, "MizerParams")

  expect_warning(
    params <- evoParams(no_sp = 4, R_factor = 10),
    "deprecated"
  )
  expect_s3_class(params, "MizerParams")
})

test_that("evoParams works with updateParams", {
  params <- evoParams(updateParams = NS_params)
  expect_s3_class(params, "MizerParams")
  expect_true("mizerEvolution" %in% names(getMetadata(params)$extensions))
  expect_true("lineage" %in% names(species_params(params)))
  expect_equal(nrow(species_params(params)), nrow(species_params(NS_params)))
})

test_that("extinctionRDD works", {
  sp <- data.frame(
    species = c("1", "2"),
    lineage = factor(c("1", "2")),
    R_max = c(100, 100)
  )
  rdi <- c(10, 10)
  rdd <- extinctionRDD(rdi, sp)
  expect_equal(length(rdd), 2)
  expect_true(all(rdd > 0))

  rdi_zero <- c(1e-35, 10)
  rdd_zero <- extinctionRDD(rdi_zero, sp)
  expect_equal(unname(rdd_zero[1]), 0)
  expect_true(unname(rdd_zero[2]) > 0)
})

test_that("alien_synthesis and nameGenerator and tempFun work", {
  aliens <- alien_synthesis(n = 3)
  expect_equal(nrow(aliens), 3)
  expect_true(all(c("w_inf", "w_mat", "beta", "sigma") %in% names(aliens)))

  names <- nameGenerator(iteration = 4)
  expect_equal(length(names), 4)
  expect_type(names, "character")

  w <- 10^seq(-3, 3, length.out = 10)
  scalar <- tempFun(w = w, temperature = 15, t_ref = 10, Ea = 0.6, c_a = 0.01)
  expect_true(is.matrix(scalar))
  expect_equal(dim(scalar), c(1, 10))
})
