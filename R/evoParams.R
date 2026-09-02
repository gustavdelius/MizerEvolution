#' Set up parameters for an evolutionary model
#'
#' Sets up a `MizerParams` object for an evolutionary or invasion model.
#'
#' @param no_sp The number of species in the model.
#' @param min_w_max The maximum size of the smallest species in the
#'   community. This will be rounded to lie on a grid point.
#' @param max_w_max The maximum size of the largest species in the community.
#'   This will be rounded to lie on a grid point.
#' @param min_w The size of the the egg of the smallest species. This also
#'   defines the start of the community size spectrum.
#' @param max_w The largest size in the model. By default this is set to the
#'   largest maximum size `max_w_max`. Setting it to something larger
#'   only makes sense if you plan to add larger species to the model later.
#' @param eta Ratio between maturity size and maximum size of a species.
#'   Ignored if `min_w_mat` is supplied. Default is 10^(-0.6),
#'   approximately 1/4.
#' @param min_w_mat The maturity size of the smallest species. Default value is
#'   \code{eta * min_w_max}. This will be rounded to lie on a grid point.
#' @param no_w The number of size bins in the community spectrum. These bins
#'   will be equally spaced on a logarithmic scale. Default value is such that
#'   there are 50 bins for each factor of 10 in weight.
#' @param min_w_pp The smallest size of the resource spectrum. By default this
#'   is set to the smallest value at which any of the consumers can feed.
#' @param w_pp_cutoff The largest size of the resource spectrum. Default value
#'   is max_w_max unless \code{perfect_scaling = TRUE} when it is Inf.
#' @param n Scaling exponent of the maximum intake rate.
#' @param p Scaling exponent of the standard metabolic rate. By default this is
#'   equal to the exponent `n`.
#' @param lambda Exponent of the abundance power law.
#' @param r_pp Growth rate parameter for the resource spectrum.
#' @param kappa Coefficient in abundance power law.
#' @param alpha The assimilation efficiency of the community.
#' @param ks Standard metabolism coefficient. If not provided, default will be
#'   calculated from critical feeding level argument `fc`.
#' @param fc Critical feeding level. Used to determine `ks` if it is not given
#'   explicitly.
#' @param h Maximum food intake rate.
#' @param beta Preferred predator prey mass ratio.
#' @param sigma Width of prey size preference.
#' @param f0 Expected average feeding level. Used to set `gamma`, the
#'   coefficient in the search rate. Ignored if `gamma` is given
#'   explicitly.
#' @param gamma Volumetric search rate. If not provided, default is determined
#'   by [mizer::get_gamma_default()] using the value of `f0`.
#' @param zeta Amplitude of lognormal distribution around trait when generating
#'   new phenotypes.
#' @param ext_mort_prop The proportion of the total mortality that comes from
#'   external mortality, i.e., from sources not explicitly modelled. A number in
#'   the interval [0, 1).
#' @param reproduction_level Level of density dependence in reproduction.
#' @param R_factor `r lifecycle::badge("deprecated")` Use `reproduction_level = 1 / R_factor` instead.
#' @param gear_names The names of the fishing gears for each species. A
#'   character vector, the same length as the number of species.
#' @param knife_edge_size The minimum size at which the gear or gears select
#'   fish. A single value for each gear or a vector with one value for each
#'   gear.
#' @param RDD Stock recruitment relationship function name. Default is `"extinctionRDD"`.
#' @param egg_size_scaling If TRUE, the egg size is a constant fraction of the
#'   maximum size of each species. This fraction is \code{min_w / min_w_max}. If
#'   FALSE, all species have the egg size `w_min`.
#' @param resource_scaling If TRUE, the carrying capacity for larger resource
#'   is reduced to compensate for the fact that fish eggs and larvae are
#'   present in the same size range.
#' @param perfect_scaling If TRUE then parameters are set so that the community
#'   abundance, growth before reproduction and death are perfect power laws. In
#'   particular all other scaling corrections are turned on.
#' @param updateParams Optional `MizerParams` object to convert into evolutionary parameters.
#' @param ea_int Temperature parameter Ea for intake.
#' @param ca_int Temperature parameter ca for intake.
#' @param ed_int Temperature parameter Ed for metabolism.
#' @param cd_int Temperature parameter cd for metabolism.
#' @param t_ref Reference temperature.
#' @param t_d Deactivation temperature.
#' @param min_w_inf `r lifecycle::badge("deprecated")` Use `min_w_max` instead.
#' @param max_w_inf `r lifecycle::badge("deprecated")` Use `max_w_max` instead.
#' @param info_level Controls the amount of information messages shown.
#' @export
#' @return An object of type `MizerParams`
#' @family functions for setting up models
#' @examples
#' \dontrun{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 5)
#' plotDynamics(sim)
#' }
evoParams <- function(no_sp = 11,
                      min_w_max = 10,
                      max_w_max = 10 ^ 4,
                      min_w = 10 ^ (-3),
                      max_w = max_w_max,
                      eta = 10^(-0.6),
                      min_w_mat = min_w_max * eta,
                      no_w = log10(max_w_max / min_w) * 20 + 1,
                      min_w_pp = 1e-10,
                      w_pp_cutoff = min_w_mat,
                      n = 2 / 3,
                      p = n,
                      lambda = 2.05,
                      r_pp = 0.1,
                      kappa = 0.005,
                      alpha = 0.4,
                      h = 40,
                      beta = 100,
                      sigma = 1.3,
                      f0 = 0.6,
                      fc = 0.25,
                      ks = NA,
                      gamma = NA,
                      zeta = .2,
                      ext_mort_prop = 0,
                      reproduction_level = 1 / 4,
                      R_factor = deprecated(),
                      gear_names = "knife_edge_gear",
                      knife_edge_size = 1000,
                      RDD = "extinctionRDD",
                      egg_size_scaling = FALSE,
                      resource_scaling = FALSE,
                      perfect_scaling = FALSE,
                      updateParams = NULL,
                      ea_int = 0,
                      ca_int = 0,
                      ed_int = 0,
                      cd_int = 0,
                      t_ref = 0,
                      t_d = 0,
                      min_w_inf = deprecated(),
                      max_w_inf = deprecated(),
                      info_level = default_info_level()
) {
  if (lifecycle::is_present(min_w_inf)) {
    lifecycle::deprecate_warn("3.4.0", "evoParams(min_w_inf)", "evoParams(min_w_max)")
    min_w_max <- min_w_inf
  }
  if (lifecycle::is_present(max_w_inf)) {
    lifecycle::deprecate_warn("3.4.0", "evoParams(max_w_inf)", "evoParams(max_w_max)")
    max_w_max <- max_w_inf
  }
  if (lifecycle::is_present(R_factor)) {
    lifecycle::deprecate_warn("3.4.0", "evoParams(R_factor)", "evoParams(reproduction_level)",
                              details = "Set `reproduction_level = 1 / R_factor`.")
    reproduction_level <- 1 / R_factor
  }

  with_info_level(info_level = info_level, {
    if (is.null(updateParams)) {
      params <- newTraitParams(
        no_sp = no_sp,
        min_w_max = min_w_max,
        max_w_max = max_w_max,
        min_w = min_w,
        max_w = max_w,
        eta = eta,
        min_w_mat = min_w_mat,
        no_w = no_w,
        min_w_pp = min_w_pp,
        w_pp_cutoff = w_pp_cutoff,
        n = n,
        p = p,
        lambda = lambda,
        r_pp = r_pp,
        kappa = kappa,
        alpha = alpha,
        h = h,
        beta = beta,
        sigma = sigma,
        f0 = f0,
        fc = fc,
        ks = ks,
        gamma = gamma,
        ext_mort_prop = ext_mort_prop,
        reproduction_level = reproduction_level,
        gear_names = gear_names,
        knife_edge_size = knife_edge_size,
        egg_size_scaling = egg_size_scaling,
        resource_scaling = resource_scaling,
        perfect_scaling = perfect_scaling,
        info_level = info_level)

      species_params(params)$lineage <- as.factor(1:no_sp)
      species_params(params)$name <- nameGenerator(no_sp)
    } else {
      new_sp_params <- species_params(updateParams)
      speciesName <- new_sp_params$species
      new_sp_params$species <- as.factor(1:nrow(new_sp_params))

      params <- newMultispeciesParams(new_sp_params,
                                      interaction = interaction_matrix(updateParams),
                                      min_w_pp = min(w_full(updateParams)),
                                      info_level = info_level)
      resource_params(params) <- resource_params(updateParams)

      gear_params(params)$species <- as.factor(gear_params(params)$species)
      species_params(params)$lineage <- as.factor(1:nrow(species_params(params)))
      species_params(params)$name <- speciesName
    }

    species_params(params)$zeta <- zeta
    species_params(params)$pop <- 1
    species_params(params)$ext <- 0
    params <- setReproduction(params, RDD = RDD)

    # temperature parameters
    species_params(params)$ea_int <- ea_int
    species_params(params)$ca_int <- ca_int
    species_params(params)$ed_int <- ed_int
    species_params(params)$cd_int <- cd_int
    species_params(params)$t_d <- t_d
    species_params(params)$t_ref <- t_ref

    params <- recordExtension(
      params, "mizerEvolution",
      version = as.character(utils::packageVersion("mizerEvolution")),
      requirement = "sizespectrum/mizerEvolution"
    )

    return(params)
  })
}

neighbourDistance <- function(x) {
  y <- vector("numeric", length = length(x))
  x <- c(0, x)
  for (i in 1:(length(x) - 1)) y[i] <- x[i + 1] - x[i]
  return(y)
}

finalTouch <- function(saveFolder, params, t_max) {
  no_w <- length(w(params))
  no_w_pp <- length(w_full(params))
  no_phen <- nrow(species_params(params))

  biomass <- array(0, dim = c(t_max + 1, no_phen, no_w),
                   dimnames = list("time" = 1:(t_max + 1),
                                   "sp" = as.character(species_params(params)$species),
                                   "w" = names(w(params))))
  biomassPP <- array(0, dim = c(t_max + 1, no_w_pp),
                     dimnames = list("time" = 1:(t_max + 1),
                                     "w" = names(w_full(params))))
  gears <- unique(gear_params(params)$gear)
  effort <- array(0, dim = c(t_max + 1, length(gears)),
                  dimnames = list("time" = 1:(t_max + 1), "gear" = gears))

  sim_start <- 1
  runs <- dir(saveFolder, pattern = "^run[0-9]+\\.rds$")
  runs <- runs[order(as.numeric(gsub("[^0-9]", "", runs)))]

  for (runFile in runs) {
    tempRun <- readRDS(file.path(saveFolder, runFile))
    n_time <- dim(tempRun@n)[1]
    n_sp <- dim(tempRun@n)[2]
    biomass[sim_start:(sim_start - 1 + n_time), 1:n_sp, ] <- tempRun@n
    biomassPP[sim_start:(sim_start - 1 + n_time), ] <- tempRun@n_pp
    effort[sim_start:(sim_start - 1 + n_time), ] <- tempRun@effort
    sim_start <- sim_start - 1 + n_time
  }

  biomass[is.na(biomass)] <- 0
  tempRun@n <- biomass
  tempRun@effort <- effort
  tempRun@n_pp <- biomassPP

  component_names <- names(params@other_dynamics)
  no_components <- length(component_names)
  no_t <- t_max + 1

  list_n_other <- rep(list(NA), no_t * no_components)
  dim(list_n_other) <- c(no_t, no_components)
  dimnames(list_n_other) <- list(time = 1:(t_max + 1),
                                 component = component_names)
  tempRun@n_other <- list_n_other

  for (iSpecies in as.character(tempRun@params@species_params$species)) {
    itime <- dim(tempRun@n)[1]
    count <- 0
    while (sum(tempRun@n[itime, iSpecies, ]) <= 1e-30 && count < dim(tempRun@n)[1]) {
      count <- count + 1
      itime <- itime - 1
      tempRun@params@species_params$ext[as.character(tempRun@params@species_params$species) == iSpecies] <- itime
    }
  }

  return(tempRun)
}

#' Project evolutionary simulation forward in time
#'
#' Runs an evolutionary simulation by alternating projection steps and phenotype
#' mutation / species invasion events.
#'
#' @param params A `MizerParams` object created by [evoParams()].
#' @param initCondition An optional previous `MizerSim` object to continue from.
#' @param t_max Total simulation time in years. Default is 100.
#' @param dt Time step for projection. Default is 0.1.
#' @param mutation Mutation rate per time step per thousand (if numeric), or a data frame specifying mutation events.
#' @param trait The trait subject to mutation. Default is `"w_mat"`.
#' @param initPool Initial number of mutant phenotypes created at start. Default is 0.
#' @param initSpread Spread factor for initial phenotype trait distribution. Default is 5.
#' @param alien Rate of alien invasions per time step per thousand (if numeric), or a data frame of invading species.
#' @param alien_init_n Multiplier for initial alien abundance density.
#' @param trait_range Optional trait range data frame for [alien_synthesis()].
#' @param saveFolder Directory for storing temporary simulation step files.
#' @param effort Fishing effort. Default is 0.
#' @return A `MizerSim` object containing the full evolutionary trajectory.
#' @export
#' @examples
#' \donttest{
#' saveFolder <- file.path(tempdir(), "simTemp")
#' dir.create(saveFolder, showWarnings = FALSE)
#' params <- evoParams(no_sp = 5, RDD = "extinctionRDD")
#' sim <- evoProject(params = params, t_max = 50, mutation = 5,
#'                   saveFolder = saveFolder)
#' plotDynamics(sim)
#' }
evoProject <- function(params = NULL, initCondition = NULL, t_max = 100, dt = 0.1,
                       mutation = 2, trait = "w_mat", initPool = 0, initSpread = 5, alien = 0,
                       alien_init_n = NULL, trait_range = NULL,
                       saveFolder = file.path(tempdir(), "simTemp"), effort = 0)
{
  if (is.null(params) && is.null(initCondition)) {
    stop("Either `params` or `initCondition` must be provided.")
  }

  if (is.null(initCondition)) {
    if (initPool > 0) {
      for (iSpecies in sort(unique(species_params(params)$lineage))) {
        for (iPhenotype in seq(1, initPool)) {
          newSp <- species_params(params)[species_params(params)$species == iSpecies, , drop = FALSE]
          new_id <- max(as.numeric(as.character(species_params(params)$species))) + 1
          newSp$species <- factor(as.character(new_id), levels = new_id)
          zeta <- if (!is.null(newSp$zeta)) newSp$zeta else 0.2

          switch(trait,
                 size = {
                   sd <- as.numeric(zeta * initSpread * newSp$w_inf)
                   newSp$w_inf <- abs(newSp$w_inf + rnorm(1, 0, sd))
                   if (!is.null(newSp$w_mat)) newSp$w_mat <- newSp$w_inf * (newSp$w_mat / (newSp$w_inf))
                 },
                 Beta = ,
                 beta = {
                   sd <- as.numeric(zeta * initSpread * newSp$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd))
                   while (newSp$beta < 10) newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd))
                 },
                 sigma = {
                   sd <- as.numeric(zeta * initSpread * newSp$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd))
                   while (newSp$sigma < 0.5 || newSp$sigma > 5) newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd))
                 },
                 predation = {
                   sd_beta <- as.numeric(zeta * initSpread * newSp$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd_beta))
                   while (newSp$beta < 10) newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd_beta))
                   sd_sigma <- as.numeric(zeta * newSp$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd_sigma))
                   while (newSp$sigma < 0.5 || newSp$sigma > 5) newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd_sigma))
                 },
                 eta = {
                   eta_val <- if (!is.null(newSp$eta)) newSp$eta else (newSp$w_mat / newSp$w_inf)
                   sd <- as.numeric(zeta * initSpread * eta_val)
                   new_eta <- abs(eta_val + rnorm(1, 0, sd))
                   if (new_eta >= 1) new_eta <- 0.95
                   newSp$w_mat <- newSp$w_inf * new_eta
                   if ("eta" %in% names(newSp)) newSp$eta <- new_eta
                 },
                 ed_int = {
                   sd <- as.numeric(zeta * initSpread * newSp$ed_int)
                   newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                   while (newSp$ed_int < 2.5) newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                 },
                 t_d = {
                   sd <- as.numeric(zeta * initSpread * newSp$t_d)
                   newSp$t_d <- abs(newSp$t_d + rnorm(1, 0, sd))
                 },
                 temperature = {
                   sd <- as.numeric(zeta * initSpread * newSp$ed_int)
                   newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                   while (newSp$ed_int < 2.5) newSp$ed_int <- abs(newSp$ed_int + rnorm(1, 0, sd))
                   sd2 <- as.numeric(zeta * initSpread * newSp$t_d)
                   newSp$t_d <- abs(newSp$t_d + rnorm(1, 0, sd2))
                 },
                 all = {
                   sd <- as.numeric(zeta * initSpread * newSp$w_inf)
                   newSp$w_inf <- abs(newSp$w_inf + rnorm(1, 0, sd))
                   if (!is.null(newSp$w_mat)) newSp$w_mat <- newSp$w_inf * (newSp$w_mat / (newSp$w_inf))
                   sd_beta <- as.numeric(zeta * initSpread * newSp$beta)
                   newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd_beta))
                   while (newSp$beta < 10) newSp$beta <- abs(newSp$beta + rnorm(1, 0, sd_beta))
                   sd_sigma <- as.numeric(zeta * initSpread * newSp$sigma)
                   newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd_sigma))
                   while (newSp$sigma < 0.5 || newSp$sigma > 5) newSp$sigma <- abs(newSp$sigma + rnorm(1, 0, sd_sigma))
                 },
                 {
                   sd <- as.numeric(zeta * initSpread * as.numeric(newSp[[trait]]))
                   newSp[[trait]] <- abs(newSp[[trait]] + rnorm(1, 0, sd))
                 })

          n_init <- initialN(params)
          n_newSp <- array(0, dim = c(1, dim(n_init)[2]),
                           dimnames = list("sp" = as.character(newSp$species),
                                           "w" = dimnames(n_init)$w))
          n_init <- rbind(n_init, n_newSp)
          dimnames(n_init)[[1]][nrow(n_init)] <- as.character(newSp$species)

          params <- addSpecies(params = params, species_params = newSp, init_n = n_init)
        }
      }

      original_n <- initialN(params)[unique(as.character(species_params(params)$lineage)), , drop = FALSE]
      init_n <- initialN(params)
      for (iSpecies in unique(as.character(species_params(params)$lineage))) {
        biomRandom <- runif(initPool + 1, 0, 1)
        biomFrac <- biomRandom / sum(biomRandom)
        position <- 0
        for (iPhen in which(as.character(species_params(params)$lineage) == iSpecies)) {
          position <- position + 1
          init_n[iPhen, ] <- original_n[iSpecies, ] * biomFrac[position]
        }
      }
      initialN(params) <- init_n
    }
  } else {
    params <- initCondition@params
    initialN(params) <- initCondition@n[dim(initCondition@n)[1], , ]
    initialNResource(params) <- initCondition@n_pp[dim(initCondition@n_pp)[1], ]
  }

  SpIdx <- unique(as.character(species_params(params)$lineage))

  if (!dir.exists(saveFolder)) dir.create(saveFolder, recursive = TRUE)

  ## Evolution / mutation schedule
  if (is.numeric(mutation)) {
    t_mutation <- matrix(0, nrow = length(SpIdx), ncol = (t_max / 1),
                         dimnames = list("species" = SpIdx, "time" = 1:(t_max / 1)))
    for (iSpecies in SpIdx) {
      for (iTime in 1:(t_max / 1)) {
        if (mutation > sample(seq(0, 100, .1), 1))
          t_mutation[iSpecies, iTime] <- 1
      }
    }
    t_phen <- apply(t_mutation, 2, sum)
    t_phen <- which(t_phen >= 1)
  } else if (is.data.frame(mutation)) {
    t_phen <- mutation$time
    t_mutation <- matrix(0, nrow = length(SpIdx), ncol = (t_max / 1),
                         dimnames = list("species" = SpIdx, "time" = 1:(t_max / 1)))
    for (iRow in 1:nrow(mutation)) {
      sp_name <- as.character(mutation$species[iRow])
      if (sp_name %in% rownames(t_mutation)) {
        t_mutation[sp_name, mutation$time[iRow]] <- 1
      }
    }
  } else {
    stop("The mutation argument has to be numeric or dataframe. Try again.")
  }

  ## Invasion schedule
  if (is.numeric(alien)) {
    t_invasion <- matrix(0, nrow = length(SpIdx), ncol = (t_max / 1),
                         dimnames = list("alien" = SpIdx, "time" = 1:(t_max / 1)))
    for (iTime in 1:(t_max / 1)) {
      if (alien > sample(seq(0, 100, .1), 1))
        t_invasion[, iTime] <- 1
    }
    t_alien <- apply(t_invasion, 2, sum)
    t_alien <- which(t_alien >= 1)
  } else if (is.data.frame(alien)) {
    t_alien <- alien$time
    t_invasion <- matrix(0, nrow = nrow(alien), ncol = (t_max / 1),
                         dimnames = list("alien" = as.character(alien$species), "time" = 1:(t_max / 1)))
    for (iRow in 1:nrow(alien)) t_invasion[as.character(alien$species[iRow]), alien$time[iRow]] <- 1
  } else {
    stop("The alien argument has to be numeric or dataframe. Try again.")
  }

  t_event <- sort(unique(c(t_phen, t_alien)))
  t_max_vec <- neighbourDistance(x = c(t_event, t_max))
  if (t_max_vec[length(t_max_vec)] == 0) t_max_vec <- t_max_vec[-length(t_max_vec)]

  mySim <- project(params, t_max = t_max_vec[1], dt = dt, progress_bar = FALSE, effort = effort)

  if (length(t_max_vec) > 1) {
    saveRDS(mySim, file = file.path(saveFolder, "run1.rds"))

    for (iSim in 2:length(t_max_vec)) {
      lastBiom_updated <- FALSE

      # alien invasion
      if (sum(t_invasion[, t_event[iSim - 1]])) {
        if (is.numeric(alien)) {
          newSp <- alien_synthesis(trait_range = trait_range)
          new_id <- max(as.numeric(as.character(species_params(mySim@params)$species))) + 1
          newSp$species <- factor(as.character(new_id), levels = new_id)
          newSp$name <- nameGenerator()
          if (is.infinite(species_params(mySim@params)$R_max[1])) {
            newSp$R_max <- Inf
          } else {
            newSp$R_max <- resource_params(mySim@params)$kappa * newSp$w_inf^-1
          }
          newSp$ca_int <- 0
          newSp$cd_int <- 0
          newSp$zeta <- .2
          newSp$fc <- .25
          newSp$f0 <- .6
          new_lineage <- max(as.numeric(as.character(species_params(mySim@params)$lineage))) + 1
          newSp$lineage <- factor(as.character(new_lineage), levels = new_lineage)
          newSp$pop <- t_event[iSim - 1]
          newSp$ext <- 0
        } else if (is.data.frame(alien)) {
          alien_row <- alien[iSim - 1, , drop = FALSE]
          newSp <- alien_synthesis()
          new_id <- max(as.numeric(as.character(species_params(mySim@params)$species))) + 1
          newSp$species <- factor(as.character(new_id), levels = new_id)
          newSp$name <- nameGenerator()
          new_lineage <- max(as.numeric(as.character(species_params(mySim@params)$lineage))) + 1
          newSp$lineage <- factor(as.character(new_lineage), levels = new_lineage)
          newSp$pop <- t_event[iSim - 1]
          newSp$ext <- 0
          for (col in names(alien_row)) {
            if (!is.na(alien_row[[col]][1])) {
              newSp[[col]] <- alien_row[[col]][1]
            }
          }
        }

        lastBiom <- array(mySim@n[dim(mySim@n)[1], , ],
                          dim = c(dim(mySim@n)[2], dim(mySim@n)[3]),
                          dimnames = dimnames(mySim@n)[2:3])
        lastBiom_updated <- TRUE

        n0_mult <- alien_init_n
        a <- 0.35
        no_w <- length(w(mySim@params))
        n <- 2 / 3
        q <- 0.8
        if (is.null(n0_mult)) {
          lambda <- resource_params(mySim@params)$lambda
          kappa <- resource_params(mySim@params)$kappa
          n0_mult <- kappa / 1000
        }

        w_inf_val <- if (!is.null(newSp$w_inf) && !is.na(newSp$w_inf)) newSp$w_inf else 100
        initial_n <- unlist(tapply(w(mySim@params), 1:no_w, function(wx, n0_mult, w_inf, a, n, q) {
          n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a)
        }, n0_mult = n0_mult, w_inf = w_inf_val, a = a, n = n, q = q))

        initial_n[w(mySim@params) > w_inf_val] <- 0
        if (!is.null(newSp$w_min) && !is.na(newSp$w_min)) initial_n[w(mySim@params) < newSp$w_min] <- 0
        n_newSp <- t(initial_n)

        init_n <- rbind(lastBiom, n_newSp)
        rownames(init_n)[nrow(init_n)] <- as.character(newSp$species)

        mySim@params <- addSpecies(params = mySim@params, species_params = newSp, init_n = init_n)
      }

      # mutation / new phenotype
      if (sum(t_mutation[, t_event[iSim - 1]])) {
        resident_list <- as.character(which(t_mutation[, t_event[iSim - 1]] > 0))

        for (resident_lineage in resident_list) {
          if (lastBiom_updated) {
            lastBiom <- initialN(mySim@params)
          } else {
            lastBiom <- array(mySim@n[dim(mySim@n)[1], , ],
                              dim = c(dim(mySim@n)[2], dim(mySim@n)[3]),
                              dimnames = dimnames(mySim@n)[2:3])
          }
          lastBiom_updated <- TRUE

          if (is.numeric(mutation)) {
            candidate_residents <- as.character(species_params(mySim@params)$species[which(as.character(species_params(mySim@params)$lineage) == resident_lineage)])
            if (length(candidate_residents) == 0) candidate_residents <- as.character(species_params(mySim@params)$species)
            resident <- if (length(candidate_residents) == 1) candidate_residents else sample(candidate_residents, 1)

            phenCount <- 0
            while (sum(mySim@n[t_max_vec[iSim - 1], resident, ]) < 2e-28 && phenCount < length(candidate_residents) * 2) {
              resident <- if (length(candidate_residents) == 1) candidate_residents else sample(candidate_residents, 1)
              phenCount <- phenCount + 1
            }

            newSp <- species_params(mySim@params)[as.character(species_params(mySim@params)$species) == resident, , drop = FALSE]
            newSp$pop <- t_event[iSim - 1]
            zeta <- if (!is.null(newSp$zeta)) newSp$zeta else 0.2

            switch(trait,
                   size = {
                     sd <- as.numeric(zeta * newSp$w_inf)
                     newSp$w_inf <- newSp$w_inf + rnorm(1, 0, sd)
                   },
                   predation = {
                     sd_beta <- as.numeric(zeta * newSp$beta)
                     newSp$beta <- newSp$beta + rnorm(1, 0, sd_beta)
                     while (newSp$beta < 10) newSp$beta <- newSp$beta + rnorm(1, 0, sd_beta)

                     sd_sigma <- as.numeric(zeta * newSp$sigma)
                     newSp$sigma <- newSp$sigma + rnorm(1, 0, sd_sigma)
                     while (newSp$sigma < 0.5 || newSp$sigma > 5) newSp$sigma <- newSp$sigma + rnorm(1, 0, sd_sigma)
                   },
                   {
                     sd <- as.numeric(zeta * as.numeric(newSp[[trait]]))
                     newSp[[trait]] <- abs(newSp[[trait]] + rnorm(1, 0, sd))
                   })

            n_newSp <- 0.05 * lastBiom[dimnames(lastBiom)[[1]] == resident, , drop = FALSE]
            lastBiom[dimnames(lastBiom)[[1]] == resident, ] <- lastBiom[dimnames(lastBiom)[[1]] == resident, ] - n_newSp
          } else if (is.data.frame(mutation)) {
            resident_id <- as.character(mutation$species[iSim - 1])
            parent_sp <- species_params(mySim@params)[as.character(species_params(mySim@params)$species) == resident_id, , drop = FALSE]
            if (nrow(parent_sp) == 0) parent_sp <- species_params(mySim@params)[1, , drop = FALSE]
            newSp <- parent_sp
            for (col in names(mutation)) {
              if (col %in% names(newSp) && !is.na(mutation[[col]][iSim - 1])) {
                newSp[[col]] <- mutation[[col]][iSim - 1]
              }
            }
            newSp$pop <- t_event[iSim - 1]
            newSp$ext <- 0

            n0_mult <- if ("init_n_multiplier" %in% names(mutation)) mutation$init_n_multiplier[iSim - 1] else NULL
            a <- 0.35
            no_w <- length(w(mySim@params))
            n <- 2 / 3
            q <- 0.8
            if (is.null(n0_mult) || length(n0_mult) == 0 || is.na(n0_mult)) {
              kappa <- resource_params(mySim@params)$kappa
              n0_mult <- kappa / 1000
            }
            w_inf_val <- if (!is.null(newSp$w_inf) && !is.na(newSp$w_inf)) newSp$w_inf else 100
            initial_n <- unlist(tapply(w(mySim@params), 1:no_w, function(wx, n0_mult, w_inf, a, n, q) {
              n0_mult * w_inf^(2 * n - q - 2 + a) * wx^(-n - a)
            }, n0_mult = n0_mult, w_inf = w_inf_val, a = a, n = n, q = q))
            initial_n[w(mySim@params) > w_inf_val] <- 0
            if (!is.null(newSp$w_min) && !is.na(newSp$w_min)) initial_n[w(mySim@params) < newSp$w_min] <- 0
            n_newSp <- t(initial_n)
          }

          new_id <- max(as.numeric(as.character(species_params(mySim@params)$species))) + 1
          newSp$species <- factor(as.character(new_id), levels = new_id)

          init_n <- rbind(lastBiom, n_newSp)
          rownames(init_n)[nrow(init_n)] <- as.character(newSp$species)

          mySim@params <- addSpecies(params = mySim@params, species_params = newSp, init_n = init_n)
        }
      }

      if (t_max_vec[iSim] > 0) {
        mySim <- project(mySim@params, t_max = t_max_vec[iSim], dt = dt, progress_bar = FALSE, effort = effort)
        saveRDS(mySim, file = file.path(saveFolder, paste0("run", iSim, ".rds")))
      }
    }

    sim <- finalTouch(saveFolder = saveFolder, params = mySim@params, t_max = t_max)
    unlink(saveFolder, recursive = TRUE)
    return(sim)
  }

  return(mySim)
}

addSpecies <- function(params, species_params, interaction, defaultInteraction = 1, init_n) {
  if (any(as.character(species_params$species) %in% as.character(params@species_params$species))) {
    stop("You can not add species that are already there.")
  }
  no_old_sp <- nrow(params@species_params)
  old_sp <- 1:no_old_sp
  no_new_sp <- nrow(species_params)
  new_sp <- (no_old_sp + 1):(no_old_sp + no_new_sp)
  no_sp <- no_old_sp + no_new_sp

  if (missing(interaction) || is.null(interaction)) {
    inter <- matrix(defaultInteraction, nrow = no_sp, ncol = no_sp)
    inter[old_sp, old_sp] <- params@interaction
  } else if (all(dim(interaction) == c(no_new_sp, no_new_sp))) {
    inter <- matrix(defaultInteraction, nrow = no_sp, ncol = no_sp)
    inter[old_sp, old_sp] <- params@interaction
    inter[new_sp, new_sp] <- interaction
  } else if (all(dim(interaction) != c(no_sp, no_sp))) {
    stop("interaction matrix has invalid dimensions.")
  } else {
    inter <- interaction
  }

  missing <- setdiff(names(params@species_params), names(species_params))
  species_params[missing] <- NA
  missing <- setdiff(names(species_params), names(params@species_params))
  params@species_params[missing] <- NA
  combi_species_params <- rbind(params@species_params, species_params,
                                stringsAsFactors = FALSE)

  # fishing params
  parent_lineage <- as.character(species_params$lineage)
  parent_gear_idx <- which(as.character(params@gear_params$species) %in% parent_lineage)
  if (length(parent_gear_idx) == 0) parent_gear_idx <- 1
  new_gear_params <- params@gear_params[parent_gear_idx[1], , drop = FALSE]
  new_gear_params$species <- as.character(species_params$species)
  combi_gear_params <- rbind(params@gear_params, new_gear_params, stringsAsFactors = FALSE)

  r_pp_val <- if (!is.null(params@resource_params$r_pp)) params@resource_params$r_pp else 0.1

  p <- newMultispeciesParams(
    combi_species_params,
    interaction = inter,
    min_w = min(params@w),
    max_w = max(params@w),
    min_w_pp = min(params@w_full),
    no_w = length(params@w),
    initial_effort = params@initial_effort,
    gear_params = combi_gear_params,
    RDD = params@rates_funcs$RDD,
    n = params@resource_params$n,
    resource_rate = r_pp_val,
    kappa = params@resource_params$kappa,
    lambda = params@resource_params$lambda,
    w_pp_cutoff = params@resource_params$w_pp_cutoff,
    info_level = 0
  )

  p@initial_n_pp <- params@initial_n_pp
  p@cc_pp <- params@cc_pp
  p@rr_pp <- params@rr_pp
  p@resource_dynamics <- params@resource_dynamics
  p@resource_params <- params@resource_params

  p@other_params <- params@other_params
  p@rates_funcs <- params@rates_funcs

  comment(p) <- comment(params)
  p@initial_n <- init_n
  p@psi[old_sp, ] <- params@psi
  p@mu_b[old_sp, ] <- params@mu_b
  p@mu_b[new_sp, ] <- rep(params@mu_b[1, ], each = no_new_sp)

  p <- recordExtension(
    p, "mizerEvolution",
    version = as.character(utils::packageVersion("mizerEvolution")),
    requirement = "sizespectrum/mizerEvolution"
  )

  return(p)
}

#' Stock recruitment relationship enabling extinction of species
#'
#' When a species reaches an abundance threshold its recruitment gets disabled.
#'
#' @param rdi Vector of density-independent reproduction rates.
#' @param species_params Data frame with species parameters.
#' @param ... Additional arguments.
#' @export
extinctionRDD <- function(rdi, species_params, ...) {
  if (!("R_max" %in% names(species_params))) {
    stop("The R_max column is missing in species_params.")
  }
  rdiNormal <- vector(mode = "numeric", length = length(rdi))
  names(rdi) <- as.character(species_params$lineage)
  for (iSpecies in sort(unique(as.character(species_params$lineage)))) {
    rdiSp <- rdi
    rdiSp[which(names(rdi) != iSpecies)] <- 0

    for (i in seq_along(rdiSp)) {
      if (is.na(rdiSp[i])) rdiSp[i] <- 1e-30
    }

    if (sum(rdiSp) != 0) {
      rdiNormal <- rdiNormal + rdiSp / sum(rdiSp)
    }
  }

  if (sum(is.infinite(species_params$R_max))) {
    rdd <- rdi
  } else {
    r_maxN <- species_params$R_max * rdiNormal
    for (i in seq_along(r_maxN)) {
      if (r_maxN[i] == 0) r_maxN[i] <- 1
    }
    rdd <- rdi / (1 + rdi / r_maxN)
  }
  rdd[rdd <= 1e-30] <- 0
  return(rdd)
}

#' Generate new alien species parameters
#'
#' Uses a range of traits to produce parameters for an invading species.
#'
#' @param trait_range A data frame containing parameter ranges, or `NULL` for defaults.
#' @param n Number of alien species to generate. Default is 1.
#' @return A data frame of species parameters.
#' @export
alien_synthesis <- function(trait_range = NULL, n = 1) {
  if (is.null(trait_range)) {
    trait_range <- data.frame(
      "trait" = c("w_inf", "beta", "sigma", "k_vb", "ks", "eta", "ea_int", "ed_int", "t_ref", "t_d"),
      "distribution" = c("lnorm", "uform", "uform", "uform", "uform", "norm", "uform", "uform", "uform", "uform"),
      "var1" = c(6.7648459, 10, 0.8, 0.1, 2.8, 0.12071530, 0.1, 1, 0, 5),
      "var2" = c(2.2524833, 400000, 3.2, 1, 13, 0.11480856, 5, 15, 25, 35),
      row.names = 1
    )
  }

  result_list <- vector("list", n)
  for (iAlien in 1:n) {
    w_inf <- sigma <- k_vb <- ks <- eta <- beta <- ea_int <- ed_int <- t_ref <- t_d <- -1

    while (w_inf <= 0 || w_inf > 10000) w_inf <- rlnorm(1, trait_range["w_inf", ]$var1, trait_range["w_inf", ]$var2)
    while (beta <= 0) beta <- runif(1, log10(trait_range["beta", ]$var1), log10(trait_range["beta", ]$var2))
    while (sigma <= 0) sigma <- runif(1, trait_range["sigma", ]$var1, trait_range["sigma", ]$var2)
    while (k_vb <= 0) k_vb <- runif(1, trait_range["k_vb", ]$var1, trait_range["k_vb", ]$var2)
    while (ks <= 0) ks <- runif(1, trait_range["ks", ]$var1, trait_range["ks", ]$var2)
    while (eta < 0.027) eta <- rnorm(1, trait_range["eta", ]$var1, trait_range["eta", ]$var2)
    while (ea_int <= 0) ea_int <- runif(1, trait_range["ea_int", ]$var1, trait_range["ea_int", ]$var2)
    while (ed_int < ea_int) ed_int <- runif(1, trait_range["ed_int", ]$var1, trait_range["ed_int", ]$var2)
    while (t_ref <= 0) t_ref <- runif(1, trait_range["t_ref", ]$var1, trait_range["t_ref", ]$var2)
    while (t_d < t_ref) t_d <- runif(1, trait_range["t_d", ]$var1, trait_range["t_d", ]$var2)

    w_mat <- w_inf * eta
    beta <- 10^beta

    result_list[[iAlien]] <- data.frame(
      "w_inf" = w_inf,
      "w_mat" = w_mat,
      "beta" = beta,
      "sigma" = sigma,
      "k_vb" = k_vb,
      "ks" = ks,
      "ea_int" = ea_int,
      "ed_int" = ed_int,
      "t_ref" = t_ref,
      "t_d" = t_d
    )
  }

  return(do.call(rbind, result_list))
}

#' Generate random species names
#'
#' Creates random pronounceable names composed of consonants and vowels.
#'
#' @param iteration Number of names to generate. Default is 1.
#' @return A character vector of names.
#' @export
nameGenerator <- function(iteration = 1) {
  name_vec <- character(iteration)
  vowel <- letters[c(1, 5, 9, 15, 21, 25)]
  consonant <- letters[-c(1, 5, 9, 15, 21, 25)]
  CONSONANT <- LETTERS[-c(1, 5, 9, 15, 21, 25)]

  for (i in 1:iteration) {
    firstLetter <- sample(CONSONANT, 1)
    vowels <- sample(vowel, 2, TRUE)
    consonants <- sample(consonant, 2, TRUE)
    name_vec[i] <- paste0(firstLetter, vowels[1], consonants[1], vowels[2], consonants[2])
  }

  return(name_vec)
}

#' Calculate temperature scalar
#'
#' Calculates a physiological temperature scalar based on the Padfield (2016)
#' temperature performance equation.
#'
#' @param w Vector of body weights.
#' @param temperature Ambient temperature in Celsius.
#' @param t_ref Reference temperature in Celsius.
#' @param Ea Activation energy.
#' @param c_a Temperature scaling parameter.
#' @return A numeric matrix of temperature scalars.
#' @export
tempFun <- function(w, temperature, t_ref, Ea, c_a) {
  k <- 8.617332e-5
  temperature <- temperature + 273
  t_ref <- t_ref + 273

  temperatureScalar <- t(sapply(w, FUN = function(x) { x^(c_a * (temperature - t_ref)) }) *
                           exp((-Ea / k) * ((1 / temperature) - (1 / t_ref))))
  return(temperatureScalar)
}
