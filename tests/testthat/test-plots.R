test_that("plotting functions produce ggplots", {
  saveFolder <- file.path(tempdir(), "test_sim_plots")
  params <- evoParams(no_sp = 3)
  sim <- evoProject(params = params, t_max = 3, mutation = 0, alien = 0,
                    saveFolder = saveFolder)

  expect_s3_class(plotDynamics(sim), "ggplot")
  expect_s3_class(plotSS(sim), "ggplot")
  expect_s3_class(plotevoFeeding(sim), "ggplot")
  expect_s3_class(plotevoGrowth(sim), "ggplot")
  expect_s3_class(plotStarvation(sim), "ggplot")
  expect_s3_class(plotevoMortality(sim), "ggplot")
  expect_s3_class(plotSpawn(sim), "ggplot")
  expect_s3_class(plotPredRate(sim), "ggplot")
  expect_s3_class(plotPredMort(sim), "ggplot")
  expect_s3_class(plotevoTrait(sim, returnData = FALSE), "ggplot")
})
