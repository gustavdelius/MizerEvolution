test_that("evoParams works", {
  params <- evoParams(no_sp = 2, info_level = 0)
  expect_s4_class(params, "MizerParams")
  expect_equal(unname(reproduction_level(params)), rep(1 / 4, 2))
  expect_true("mizerEvolution" %in% names(params@extensions))
})

test_that("evoParams converts a model without rebuilding its resource", {
  original <- newTraitParams(no_sp = 2, info_level = 0)
  original <- renameSpecies(original, c(`1` = "Alpha", `2` = "Beta"))

  converted <- evoParams(updateParams = original, info_level = 0)

  expect_equal(as.numeric(resource_rate(converted)),
               as.numeric(resource_rate(original)))
  expect_equal(as.numeric(resource_capacity(converted)),
               as.numeric(resource_capacity(original)))
  expect_equal(resource_dynamics(converted), resource_dynamics(original))
  expect_equal(as.character(species_params(converted)$species), c("1", "2"))
  expect_equal(unname(species_params(converted)$name), c("Alpha", "Beta"))
})

test_that("the species addition path preserves the model and requested abundance", {
  params <- evoParams(no_sp = 2, info_level = 0)
  new_species <- species_params(params)[1, , drop = FALSE]
  new_species$species <- "3"
  target_n <- rbind(initialN(params), initialN(params)[1, ] * 0.05)
  rownames(target_n)[3] <- "3"

  updated <- mizerEvolution:::addSpecies(
    params,
    new_species,
    init_n = target_n,
    info_level = 0
  )

  expect_equal(as.numeric(initialN(updated)), as.numeric(target_n))
  expect_equal(getRateFunction(updated, "RDD"), "extinctionRDD")
  expect_true("mizerEvolution" %in% names(updated@extensions))
})

test_that("evoProject can create an initial phenotype pool", {
  params <- evoParams(no_sp = 2, info_level = 0)
  sim <- evoProject(
    params = params,
    t_max = 1,
    mutation = 0,
    initPool = 1,
    trait = "w_mat",
    saveFolder = file.path(tempdir(), "mizer-evolution-init-pool"),
    info_level = 0
  )

  expect_equal(nrow(species_params(sim@params)), 4)
  expect_true(all(species_params(sim@params)$w_mat25 <
                    species_params(sim@params)$w_mat))
})

test_that("evoProject assembles scheduled invasion runs by dimension name", {
  set.seed(21)
  params <- evoParams(no_sp = 3, info_level = 0)
  alien <- species_params(params)[1, , drop = FALSE]
  alien$time <- 1
  alien$lineage <- factor(4, levels = 4)
  alien$species <- "4"
  save_folder <- file.path(tempdir(), "mizer-evolution-scheduled-invasion")
  dir.create(save_folder, showWarnings = FALSE)
  marker <- file.path(save_folder, "keep-me")
  writeLines("user file", marker)

  sim <- evoProject(
    params = params,
    t_max = 2,
    mutation = 0,
    alien = alien,
    saveFolder = save_folder,
    info_level = 0
  )

  expect_equal(dim(N(sim)), c(3, 4, length(sim@params@w)))
  expect_equal(getTimes(sim), 0:2)
  expect_equal(dimnames(getEffort(sim))$gear, "knife_edge_gear")
  expect_true(file.exists(marker))
})

test_that("evoProject supports generated invaders", {
  set.seed(22)
  params <- evoParams(no_sp = 3, info_level = 0)

  sim <- evoProject(
    params = params,
    t_max = 2,
    mutation = 0,
    alien = 101,
    saveFolder = tempdir(),
    info_level = 0
  )

  expect_equal(dim(N(sim))[2], 4)
  expect_equal(getTimes(sim), 0:2)
})

test_that("North Sea mutation aligns abundance after adding species", {
  set.seed(23)
  params <- evoParams(updateParams = NS_params, info_level = 0)

  sim <- evoProject(
    params = params,
    t_max = 2,
    mutation = 101,
    trait = "beta",
    saveFolder = tempdir(),
    info_level = 0
  )

  expect_equal(dim(N(sim))[2], 2 * nrow(species_params(params)))
  expect_equal(getTimes(sim), 0:2)
  expect_equal(dim(getEffort(sim))[2], length(initial_effort(params)))
})

test_that("plotDynamics selections contain the same curves as the full plot", {
  set.seed(24)
  params <- evoParams(no_sp = 3, info_level = 0)
  sim <- evoProject(
    params,
    t_max = 1,
    mutation = 0,
    initPool = 1,
    saveFolder = tempdir(),
    info_level = 0
  )

  full_plot <- plotDynamics(sim)
  full_totals <- full_plot$data
  full_phenotypes <- full_plot$layers[[2]]$data

  species_plot <- plotDynamics(sim, species = 2)
  expect_equal(
    species_plot$data,
    dplyr::filter(full_totals, lineage == "2")
  )
  expect_equal(
    species_plot$layers[[2]]$data,
    dplyr::filter(full_phenotypes, lineage == "2")
  )

  subset_plot <- plotDynamics(sim, SpIdx = c(1, 3))
  expect_equal(
    subset_plot$data,
    dplyr::filter(full_totals, lineage %in% c("1", "3"))
  )
  expect_equal(
    subset_plot$layers[[2]]$data,
    dplyr::filter(full_phenotypes, lineage %in% c("1", "3"))
  )
})

test_that("plotSS aggregates by lineage and removes numerical extinction dust", {
  set.seed(25)
  params <- evoParams(no_sp = 3, info_level = 0)
  sim <- evoProject(
    params,
    t_max = 1,
    mutation = 0,
    initPool = 1,
    saveFolder = tempdir(),
    info_level = 0
  )

  spectrum <- plotSS(sim, biomass = FALSE, returnData = TRUE)
  final_n <- as.matrix(N(sim)[dim(N(sim))[1], , ])
  sp <- species_params(finalParams(sim))
  lineage <- as.character(sp$lineage[
    match(rownames(final_n), as.character(sp$species))
  ])
  expected <- rowsum(final_n, lineage, reorder = FALSE)
  expected_data <- data.frame(
    value = c(expected),
    Species = rep(rownames(expected), times = ncol(expected)),
    w = rep(sim@params@w, each = nrow(expected)),
    row.names = NULL
  )
  expected_data <- expected_data[expected_data$value >= 1e-20, ]

  expect_equal(spectrum, expected_data)
  expect_gte(min(spectrum$value), 1e-20)

  phenotype_spectrum <- plotSS(
    sim,
    biomass = FALSE,
    species = FALSE,
    returnData = TRUE
  )
  lineage_by_phenotype <- setNames(lineage, rownames(final_n))
  expect_equal(
    phenotype_spectrum$lineage,
    unname(lineage_by_phenotype[phenotype_spectrum$Species])
  )
})

test_that("current mortality and starvation plotting APIs are supported", {
  skip_if_not_installed("mizerStarvation", minimum_version = "0.1.1")
  suppressPackageStartupMessages(library(mizerStarvation))
  params <- evoParams(no_sp = 2, info_level = 0)
  params <- mizerStarvation::setStarvation(params)
  sim <- project(params, t_max = 1, progress_bar = FALSE)

  expect_s3_class(
    plotevoMortality(sim, returnData = TRUE),
    "data.frame"
  )
  expect_s3_class(
    plotStarvation(sim, returnData = TRUE),
    "data.frame"
  )
})
