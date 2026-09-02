test_that("evoProject runs without events", {
  saveFolder <- file.path(tempdir(), "test_sim_none")
  params <- evoParams(no_sp = 3)
  sim <- evoProject(params = params, t_max = 2, mutation = 0, alien = 0,
                    saveFolder = saveFolder)
  expect_s3_class(sim, "MizerSim")
  expect_equal(dim(sim@n)[1], 3)
})

test_that("evoProject runs with mutation and alien events", {
  saveFolder <- file.path(tempdir(), "test_sim_events")
  params <- evoParams(no_sp = 3)

  mut_df <- data.frame(species = "1", time = 2)
  sim <- evoProject(params = params, t_max = 4, mutation = mut_df, alien = 0,
                    saveFolder = saveFolder)
  expect_s3_class(sim, "MizerSim")
  expect_true(dim(sim@n)[1] >= 4)
  expect_true(nrow(species_params(sim@params)) >= 3)
})
