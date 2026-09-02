#' Plot biomass variations through time of species and phenotypes
#'
#' Plots the total biomass of each species and optionally its individual phenotypes
#' over the course of the simulation.
#'
#' @param object An object of class [mizer::MizerSim].
#' @param phenotype A boolean value to display the phenotype's
#'   abundance on the plot. Default is TRUE.
#' @param species Plot only a specific species and its phenotypes.
#'   Default is NULL.
#' @param trait Display the trait value as color gradient. Works
#'   only if `species` is not NULL.
#' @param SpIdx Display a subset of species.
#' @param time_range The time range (either a vector of values, a vector of min
#'   and max time, or a single value) to average the abundances over. Default is
#'   the entire simulation time range.
#' @param print_it A boolean value that displays the plot. Default is TRUE.
#' @param returnData A boolean value that controls the output. Either a dataframe
#'   ready for ggplot (TRUE) or the plot itself (FALSE). Default is FALSE.
#' @param save_it A boolean value to save the plot. Default is FALSE.
#' @param nameSave File name for saving the plot. Default is `"Biomass.png"`.
#' @param ylimit To edit the y axis scale on the plot.
#'
#' @return A ggplot object, or a list of data frames if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotDynamics(sim)
#' }
plotDynamics <- function(object, phenotype = TRUE, species = NULL, trait = NULL, SpIdx = NULL,
                         time_range = c(min(as.numeric(dimnames(object@n)$time)),
                                        max(as.numeric(dimnames(object@n)$time))),
                         print_it = TRUE, returnData = FALSE, save_it = FALSE,
                         nameSave = "Biomass.png", ylimit = c(NA, NA)) {

    jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                     "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    min_value <- 1e-30
    biomass <- getBiomass(object)
    time_elements <- get_time_elements(object, time_range)
    biomass <- biomass[time_elements, , drop = FALSE]

    if (is.null(SpIdx)) {
        for (i in unique(as.character(object@params@species_params$lineage))) {
            matching_cols <- which(as.character(object@params@species_params$lineage) == i)
            if (length(matching_cols) > 0 && sum(biomass[, matching_cols, drop = FALSE]) != 0) {
                SpIdx <- c(SpIdx, i)
            }
        }
    }

    biomassSp <- NULL
    biomassTemp <- biomass
    colnames(biomassTemp) <- as.character(object@params@species_params$lineage)
    for (i in SpIdx) {
        biomassPhen <- biomassTemp[, which(colnames(biomassTemp) == i), drop = FALSE]
        if (!is.null(dim(biomassPhen))) biomassPhen <- apply(biomassPhen, 1, sum, na.rm = TRUE)
        biomassSp <- cbind(biomassSp, biomassPhen)
    }
    colnames(biomassSp) <- SpIdx

    spSub <- as.character(object@params@species_params$species[as.character(object@params@species_params$lineage) %in% SpIdx])
    biomass <- biomass[, as.character(dimnames(biomass)$sp) %in% spSub, drop = FALSE]

    plotBiom <- function(x, params) {
        Biom <- reshape2::melt(x)
        colnames(Biom) <- c("time", "phen", "value")
        if (is.null(params@species_params$name)) {
            Biom$sp <- sapply(Biom$phen, function(p) {
                params@species_params$lineage[as.character(params@species_params$species) == as.character(p)][1]
            })
        } else {
            Biom$sp <- sapply(Biom$phen, function(p) {
                params@species_params$name[as.character(params@species_params$species) == as.character(p)][1]
            })
        }
        return(Biom)
    }

    BiomSp <- plotBiom(biomassSp, object@params)
    BiomSp <- BiomSp[BiomSp$value >= min_value, ]

    color_legend <- if (!is.null(object@params@linecolour)) object@params@linecolour[SpIdx] else NULL
    if (!is.null(object@params@species_params$name) && !is.null(color_legend)) {
        names(color_legend) <- object@params@species_params$name[as.numeric(SpIdx)]
    }

    if (phenotype) {
        BiomPhen <- plotBiom(biomass, object@params)
        if (!is.null(trait)) {
            BiomPhen$trait <- rep(trait, each = length(unique(BiomPhen$time)))
        }
        BiomPhen <- BiomPhen[BiomPhen$value >= min_value, ]

        p <- ggplot(BiomSp) +
            geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), linewidth = 1.2) +
            geom_line(data = BiomPhen, aes(x = time, y = value, colour = as.factor(sp), group = phen), alpha = 0.2) +
            scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30, 10, 2)))) +
            scale_x_continuous(name = "Time in years") +
            labs(color = "Species") +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white")) +
            ggtitle("Community biomass")

        if (!is.null(color_legend) && !any(is.na(color_legend))) {
            p <- p + scale_colour_manual(values = color_legend)
        }

        if (!is.null(species)) {
            plotTitle <- paste("Species", species)
            p <- ggplot(dplyr::filter(BiomPhen, sp == species)) +
                geom_line(aes(x = time, y = value, group = phen), alpha = .75) +
                scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30, 4, 2)))) +
                scale_x_continuous(name = "Time in years") +
                theme(panel.background = element_rect(fill = "white", color = "black"),
                      panel.grid.minor = element_line(colour = "grey92"),
                      legend.key = element_rect(fill = "white")) +
                ggtitle(plotTitle)

            if (!is.null(trait)) {
                plotTitle <- paste("Trait of species", species)
                p <- ggplot(dplyr::filter(BiomPhen, sp == species)) +
                    geom_line(aes(x = time, y = value, colour = trait, group = trait), alpha = 1) +
                    scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(seq(-30, 4, 2)))) +
                    scale_x_continuous(name = "Time in years") +
                    scale_colour_gradientn(colours = jet.colors(9), limits = c(NA, NA)) +
                    theme(panel.background = element_rect(fill = "white", color = "black"),
                          panel.grid.minor = element_line(colour = "grey92"),
                          legend.key = element_rect(fill = "white")) +
                    ggtitle(plotTitle)
            }
        }
    } else {
        p <- ggplot(BiomSp) +
            geom_line(aes(x = time, y = value, colour = as.factor(sp), group = sp), linewidth = 1.2) +
            scale_y_log10(name = "Biomass in g.m^-3", limits = ylimit, breaks = c(1 %o% 10^(-30:4))) +
            scale_x_continuous(name = "Time in years") +
            labs(color = "Species") +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white")) +
            ggtitle("Community biomass")

        if (!is.null(color_legend) && !any(is.na(color_legend))) {
            p <- p + scale_colour_manual(values = color_legend)
        }
    }

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(list(BiomSp, BiomPhen)) else if (print_it) return(p)
}


#' Plot size spectrums of species and phenotypes
#'
#' Plots the size spectrum of species, phenotypes, and the resource spectrum.
#'
#' @inheritParams plotDynamics
#' @param min_w Minimum weight to display on size axis.
#' @param ylim Limits for the y axis.
#' @param biomass If TRUE, plots biomass density; otherwise abundance density. Default is TRUE.
#' @param community If TRUE, plots community total spectrum. Default is FALSE.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotSS(sim)
#' }
plotSS <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                   min_w = min(object@params@w) / 100, ylim = c(NA, NA),
                   biomass = TRUE, print_it = TRUE, species = TRUE, community = FALSE,
                   save_it = FALSE, nameSave = "SizeSpectrum.png", returnData = FALSE) {

    cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    time_elements <- get_time_elements(object, time_range)
    spec_n <- apply(object@n[time_elements, , , drop = FALSE], c(2, 3), mean)
    pkt_n <- apply(object@n_pp[time_elements, , drop = FALSE], 2, mean)

    y_axis_name <- "Abundance density in individuals.m^-3"
    if (biomass) {
        spec_n <- sweep(spec_n, 2, object@params@w, "*")
        pkt_n <- pkt_n * object@params@w_full
        y_axis_name <- "Biomass in g.m^-3"
    }

    plot_datSP <- data.frame(value = c(spec_n), Species = rep(dimnames(spec_n)[[1]], times = ncol(spec_n)),
                             w = rep(object@params@w, each = nrow(spec_n)),
                             lineage = rep(object@params@species_params$lineage, times = ncol(spec_n)))
    plot_datPkt <- data.frame(value = c(pkt_n), Species = "Phytoplankton", w = object@params@w_full)

    if (community) {
        plot_datSP <- data.frame(value = apply(spec_n, 2, sum), w = object@params@w)
    } else if (species) {
        SpIdx <- unique(as.character(object@params@species_params$lineage))
        spec_sp <- matrix(data = NA, ncol = dim(spec_n)[2], nrow = length(SpIdx),
                          dimnames = list(SpIdx, dimnames(spec_n)[[2]]))
        names(dimnames(spec_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$lineage) == sp_id)
            if (length(matching_rows) > 0) {
                spec_sp[i, ] <- apply(spec_n[matching_rows, , drop = FALSE], 2, sum)
            }
        }
        plot_datSP <- data.frame(value = c(spec_sp), Species = rep(dimnames(spec_sp)[[1]], times = ncol(spec_sp)),
                                 w = rep(object@params@w, each = length(SpIdx)))
    }

    plot_datSP <- plot_datSP[(plot_datSP$value > 0) & (plot_datSP$w >= min_w), ]
    plot_datPkt <- plot_datPkt[(plot_datPkt$value > 0) & (plot_datPkt$w >= min_w), ]

    if (community) {
        p <- ggplot(plot_datSP) +
            geom_line(aes(x = w, y = value)) +
            geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species), alpha = 0.5, color = "blue", linewidth = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5))) +
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white")) +
            ggtitle("Community spectrum")
    } else if (species) {
        p <- ggplot(plot_datSP) +
            geom_line(aes(x = w, y = value, colour = as.factor(Species), group = Species)) +
            geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species), alpha = 0.5, color = "blue", linewidth = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5))) +
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
            labs(color = "Species") +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white")) +
            ggtitle("Size spectrum")
    } else {
        p <- ggplot(plot_datSP) +
            geom_line(aes(x = w, y = value, colour = as.factor(lineage), group = Species)) +
            geom_line(data = plot_datPkt, aes(x = w, y = value, group = Species), alpha = 0.5, color = "blue", linewidth = 1.5) +
            scale_x_continuous(name = "Size in g", trans = "log10", breaks = c(1 %o% 10^(-6:5))) +
            scale_y_continuous(name = y_axis_name, limits = ylim, trans = "log10") +
            labs(color = "Species") +
            theme(panel.background = element_rect(fill = "white", color = "black"),
                  panel.grid.minor = element_line(colour = "grey92"),
                  legend.key = element_rect(fill = "white")) +
            ggtitle("Size spectrum per phenotypes")
    }

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_datSP) else if (print_it) return(p)
}


#' Plot feeding level
#'
#' Plots the feeding level of species against body size.
#'
#' @inheritParams plotDynamics
#' @param start Initial time index if plotting series.
#' @param every Step interval if plotting series.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotevoFeeding(sim)
#' }
plotevoFeeding <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                           species = TRUE, start = 1000, every = 1000,
                           print_it = TRUE, returnData = FALSE, save_it = FALSE,
                           nameSave = "Feeding.png") {

    feed_time <- getFeedingLevel(object = object, time_range = time_range, drop = FALSE)
    feed <- apply(feed_time, c(2, 3), mean)
    Cfeed <- getCriticalFeedingLevel(object@params)

    SpIdx <- sort(unique(as.character(object@params@species_params$lineage)))

    if (species) {
        feed_sp <- matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx),
                          dimnames = list(SpIdx, dimnames(feed)$w))
        Cfeed_sp <- matrix(data = NA, ncol = dim(feed)[2], nrow = length(SpIdx),
                           dimnames = list(SpIdx, dimnames(feed)$w))

        names(dimnames(feed_sp)) <- list("species", "size")
        names(dimnames(Cfeed_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$lineage) == sp_id)
            if (length(matching_rows) > 0) {
                feed_sp[i, ] <- apply(feed[matching_rows, , drop = FALSE], 2, mean)
                Cfeed_sp[i, ] <- if (is.matrix(Cfeed)) apply(Cfeed[matching_rows, , drop = FALSE], 2, mean) else Cfeed
            }
        }
        feed <- feed_sp
        Cfeed <- Cfeed_sp
    }

    a <- object@params@species_params$w_inf[as.numeric(SpIdx)]
    vlines <- data.frame(xint = a, grp = SpIdx)
    color_legend <- if (!is.null(object@params@linecolour)) object@params@linecolour[SpIdx] else NULL

    plot_dat <- data.frame(value = c(feed), critical = c(Cfeed),
                           species = as.factor(rep(dimnames(feed)[[1]], times = ncol(feed))),
                           size = rep(object@params@w, each = nrow(feed)))

    p <- ggplot(plot_dat) +
        geom_line(aes(x = size, y = value, colour = species)) +
        geom_line(aes(x = size, y = critical), linetype = "dashed", color = "red") +
        geom_vline(data = vlines, aes(xintercept = xint, colour = as.factor(grp)), linetype = "dashed") +
        scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "Feeding Level", limits = c(0, 1)) +
        theme(panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.key = element_rect(fill = "white")) +
        ggtitle(NULL)

    if (!is.null(color_legend) && !any(is.na(color_legend))) {
        p <- p + scale_colour_manual(values = color_legend, name = "Species")
    }

    if (save_it) ggsave(plot = p, filename = nameSave, width = 18, height = 18, units = "cm")

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot growth
#'
#' Plots the instantaneous growth rates of species against body size.
#'
#' @inheritParams plotDynamics
#' @param ylim Limits for the y axis.
#' @param ... Additional arguments passed to methods.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotevoGrowth(sim)
#' }
plotevoGrowth <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                          species = TRUE, print_it = TRUE, returnData = FALSE, save_it = FALSE,
                          ylim = c(NA, NA), nameSave = "Growth.png", ...) {

    time_elements <- get_time_elements(object, time_range)
    growth_time <- plyr::aaply(which(time_elements), 1, function(x) {
        n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
        dimnames(n) <- dimnames(object@n)[2:3]
        growth <- getEGrowth(object@params, n = n, n_pp = object@n_pp[x, ])
        return(growth)
    })

    growth <- growth_time

    if (species) {
        SpIdx <- sort(unique(as.character(object@params@species_params$lineage)))
        growth_sp <- matrix(data = NA, ncol = dim(growth)[2], nrow = length(SpIdx),
                            dimnames = list(SpIdx, dimnames(growth)$w))
        names(dimnames(growth_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$lineage) == sp_id)
            if (length(matching_rows) > 0) {
                growth_sp[i, ] <- apply(growth[matching_rows, , drop = FALSE], 2, mean)
            }
        }
        growth <- growth_sp
    }

    plot_dat <- data.frame(value = c(growth), Species = rep(dimnames(growth)[[1]], times = ncol(growth)),
                           w = rep(object@params@w, each = nrow(growth)))
    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species)) +
        scale_x_continuous(name = "Size", trans = "log10", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "instantaneous growth", trans = "log10", limits = ylim) +
        theme(legend.title = element_blank(),
              legend.justification = c(1, 1),
              legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92")) +
        ggtitle(NULL)

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot starvation level
#'
#' Plots the instantaneous starvation mortality rate against body size.
#'
#' @inheritParams plotDynamics
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotStarvation(sim)
#' }
plotStarvation <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                           species = TRUE, print_it = TRUE, returnData = FALSE, save_it = FALSE,
                           nameSave = "Starvation.png") {

    cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    time_elements <- get_time_elements(object, time_range)
    death_time <- plyr::aaply(which(time_elements), 1, function(x) {
        n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
        dimnames(n) <- dimnames(object@n)[2:3]
        if (exists("getStarvMort", mode = "function")) {
            starv <- get("getStarvMort")(object@params, n = n, n_pp = object@n_pp[x, ])
        } else {
            starv <- matrix(0, nrow = nrow(object@params@species_params), ncol = length(object@params@w),
                            dimnames = list(object@params@species_params$species, names(object@params@w)))
        }
        return(starv)
    })

    if (species) {
        SpIdx <- sort(unique(as.character(object@params@species_params$species)))
        death_sp <- matrix(data = NA, ncol = dim(death_time)[2], nrow = length(SpIdx),
                          dimnames = list(SpIdx, dimnames(death_time)$w))
        names(dimnames(death_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$species) == sp_id)
            if (length(matching_rows) > 0) {
                death_sp[i, ] <- apply(death_time[matching_rows, , drop = FALSE], 2, mean)
            }
        }
        death_time <- death_sp
    }

    plot_dat <- data.frame(value = c(death_time), species = rep(dimnames(death_time)[[1]], times = ncol(death_time)),
                           size = rep(object@params@w, each = nrow(death_time)))

    name <- paste("Starvation at time", time_range, sep = " ")

    p <- ggplot(plot_dat) +
        geom_line(aes(x = size, y = value, colour = as.factor(species))) +
        scale_x_log10(name = "Size", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "Instantaneous starvation mortality") +
        scale_colour_manual(values = cbPalette, name = "Species") +
        theme(panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.key = element_rect(fill = "white")) +
        ggtitle(name)

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot total mortality
#'
#' Plots the total mortality rate of species across sizes.
#'
#' @inheritParams plotDynamics
#' @param comments Logical indicating whether to print diagnostic comments.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotevoMortality(sim)
#' }
plotevoMortality <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                             print_it = TRUE, returnData = FALSE, comments = TRUE) {

    if (is.matrix(object@effort)) {
        effort <- object@effort[time_range, ]
    } else {
        effort <- object@effort[time_range]
    }

    z <- getMort(object@params, n = object@n[time_range, , ],
                 n_pp = object@n_pp[time_range, ], effort = effort)
    rownames(z) <- as.character(object@params@species_params$lineage)
    SpIdx <- as.numeric(as.character(sort(unique(object@params@species_params$lineage))))

    z_sp <- matrix(data = NA, ncol = dim(z)[2], nrow = length(SpIdx),
                   dimnames = list(SpIdx, dimnames(z)[[2]]))
    names(dimnames(z_sp)) <- list("prey", "w_prey")

    for (i in seq_along(SpIdx)) {
        sp_id <- SpIdx[i]
        matching_rows <- which(rownames(z) == as.character(sp_id))
        if (length(matching_rows) > 0) {
            z_sp[i, ] <- apply(z[matching_rows, , drop = FALSE], 2, mean)
        }
    }
    z <- z_sp

    plot_dat <- data.frame(value = c(z), Species = rep(dimnames(z)[[1]], times = ncol(z)),
                           w = rep(object@params@w, each = nrow(z)))

    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species)) +
        scale_x_continuous(name = "Size", trans = "log10", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "Mortality", limits = c(0, max(plot_dat$value))) +
        theme(legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92")) +
        ggtitle(NULL)

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot spawn
#'
#' Plots the energy allocated to reproduction for each species across sizes.
#'
#' @inheritParams plotDynamics
#' @param ... Additional arguments passed to methods.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotSpawn(sim)
#' }
plotSpawn <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                      species = TRUE, print_it = TRUE, returnData = FALSE, save_it = FALSE,
                      nameSave = "Spawn.png", ...) {

    time_elements <- get_time_elements(object, time_range)
    spawn_time <- plyr::aaply(which(time_elements), 1, function(x) {
        n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
        dimnames(n) <- dimnames(object@n)[2:3]
        spawn <- getERepro(object@params, n = n, n_pp = object@n_pp[x, ])
        return(spawn)
    })

    spawn <- spawn_time

    if (species) {
        SpIdx <- sort(unique(as.character(object@params@species_params$species)))
        spawn_sp <- matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx),
                           dimnames = list(SpIdx, dimnames(spawn)$w))
        names(dimnames(spawn_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$species) == sp_id)
            if (length(matching_rows) > 0) {
                spawn_sp[i, ] <- apply(spawn[matching_rows, , drop = FALSE], 2, mean)
            }
        }
        spawn <- spawn_sp
    }

    sp_idx_num <- seq_len(min(9, nrow(object@params@species_params)))
    a <- object@params@species_params$w_inf[sp_idx_num]
    vlines <- data.frame(xint = a, grp = sp_idx_num)

    name <- paste("Spawn level at time", time_range, sep = " ")
    plot_dat <- data.frame(value = c(spawn), Species = rep(dimnames(spawn)[[1]], times = ncol(spawn)),
                           w = rep(object@params@w, each = nrow(spawn)))
    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species)) +
        scale_x_continuous(name = "Size", trans = "log10", breaks = c(1 %o% 10^(-3:5))) +
        scale_y_continuous(name = "Energy allocated to spawning", trans = "log10") +
        geom_vline(data = vlines, aes(xintercept = xint, colour = as.factor(grp)), linetype = "dashed") +
        theme(legend.title = element_blank(),
              legend.justification = c(1, 1),
              legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92")) +
        ggtitle(name)

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot predation rate
#'
#' Plots the predation rate against size.
#'
#' @inheritParams plotDynamics
#' @param ylim Limits for the y axis.
#' @return A ggplot object, or data frame if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotPredRate(sim)
#' }
plotPredRate <- function(object, time_range = max(as.numeric(dimnames(object@n)$time)),
                         species = TRUE, ylim = c(NA, NA),
                         print_it = TRUE, returnData = FALSE, save_it = FALSE,
                         nameSave = "PredRate.png") {

    time_elements <- get_time_elements(object, time_range)
    spawn_time <- plyr::aaply(which(time_elements), 1, function(x) {
        n <- array(object@n[x, , ], dim = dim(object@n)[2:3])
        dimnames(n) <- dimnames(object@n)[2:3]
        spawn <- getPredRate(object@params, n = n, n_pp = object@n_pp[x, ])
        return(spawn)
    })

    spawn <- spawn_time

    if (species) {
        SpIdx <- sort(unique(as.character(object@params@species_params$species)))
        spawn_sp <- matrix(data = NA, ncol = dim(spawn)[2], nrow = length(SpIdx),
                           dimnames = list(SpIdx, dimnames(spawn)$w))
        names(dimnames(spawn_sp)) <- list("species", "size")

        for (i in seq_along(SpIdx)) {
            sp_id <- SpIdx[i]
            matching_rows <- which(as.character(object@params@species_params$species) == sp_id)
            if (length(matching_rows) > 0) {
                spawn_sp[i, ] <- apply(spawn[matching_rows, , drop = FALSE], 2, mean)
            }
        }
        spawn <- spawn_sp
    }

    sp_idx_num <- seq_len(min(9, nrow(object@params@species_params)))
    a <- object@params@species_params$w_inf[sp_idx_num]
    vlines <- data.frame(xint = a, grp = sp_idx_num)

    name <- paste("Predation rate at time", time_range, sep = " ")
    plot_dat <- data.frame(value = c(spawn), Species = rep(dimnames(spawn)[[1]], times = ncol(spawn)),
                           w = rep(object@params@w_full, each = nrow(spawn)))
    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species)) +
        scale_x_continuous(name = "Size", trans = "log10", breaks = c(1 %o% 10^(-7:5))) +
        scale_y_continuous(name = "Potential death rate from predator", trans = "log10", limits = ylim) +
        geom_vline(data = vlines, aes(xintercept = xint, colour = as.factor(grp)), linetype = "dashed") +
        theme(legend.title = element_blank(),
              legend.justification = c(1, 1),
              legend.key = element_rect(fill = "white"),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92")) +
        ggtitle(name)

    if (save_it) ggsave(plot = p, filename = nameSave, scale = 1.5)

    if (returnData) return(plot_dat) else if (print_it) return(p)
}


#' Plot predation mortality rate of each species against size
#'
#' After running a projection, plot the predation mortality rate of each species
#' by size.
#'
#' @inheritParams plotDynamics
#' @param highlight Vector of species names to highlight.
#' @return A ggplot object.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions], [mizer::getPredMort()]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotPredMort(sim)
#' }
plotPredMort <- function(object, species = NULL, time_range, highlight = NULL) {
    if (inherits(object, "MizerSim")) {
        if (missing(time_range)) {
            time_range <- max(as.numeric(dimnames(object@n)$time))
        }
        params <- object@params
    } else {
        params <- validParams(object)
    }
    pred_mort <- getPredMort(object, time_range = time_range, drop = FALSE)
    if (length(dim(pred_mort)) == 3) {
        pred_mort <- apply(pred_mort, c(2, 3), mean)
    }

    if (is.null(species)) {
        species <- as.character(params@species_params$species)
    }
    species_levels <- c(as.character(params@species_params$species),
                        "Background", "Resource", "Total")
    pred_mort <- pred_mort[as.character(dimnames(pred_mort)[[1]]) %in% species, , drop = FALSE]
    plot_dat <- data.frame(value = c(pred_mort),
                           Species = factor(rep(dimnames(pred_mort)[[1]], times = ncol(pred_mort)),
                                            levels = species_levels),
                           w = rep(params@w, each = nrow(pred_mort)))
    p <- ggplot(plot_dat) +
        geom_line(aes(x = w, y = value, colour = Species,
                      linetype = Species))

    p <- p +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(name = "Predation mortality [1/year]",
                           limits = c(0, max(plot_dat$value)))

    if (!is.null(params@linecolour) && !any(is.na(params@linecolour))) {
        p <- p + scale_colour_manual(values = params@linecolour)
    }
    if (!is.null(params@linetype) && !any(is.na(params@linetype))) {
        p <- p + scale_linetype_manual(values = params@linetype)
    }

    return(p)
}


#' Plot trait evolution through time
#'
#' Plots the mean and standard deviation of selected traits through time for each species lineage.
#'
#' @inheritParams plotDynamics
#' @param Normalisation A boolean indicating whether to normalize the mean trait by initial trait value. Default is FALSE.
#' @param traitID The trait column name to plot. Default is `"w_mat"`.
#' @return A ggplot object, or data list if `returnData = TRUE`.
#' @export
#' @family plotting functions
#' @seealso [mizer::plotting_functions]
#' @examples
#' \donttest{
#' params <- evoParams(no_sp = 5)
#' sim <- evoProject(params, t_max = 20)
#' plotevoTrait(sim)
#' }
plotevoTrait <- function(object, SpIdx = NULL, Normalisation = FALSE,
                         returnData = TRUE, traitID = NULL) {
    cbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

    if (is.null(SpIdx)) SpIdx <- unique(as.character(object@params@species_params$lineage))
    TimeMax <- dim(object@n)[1]
    if (is.null(traitID)) traitID <- "w_mat"

    SumPar <- object@params@species_params
    TT <- data.frame("Species" = as.character(SumPar$lineage),
                     "Phenotype" = as.numeric(as.character(SumPar$species)),
                     "Apparition" = SumPar$pop,
                     "Extinction" = SumPar$ext,
                     SumPar[traitID])

    for (i in seq_len(nrow(TT))) {
        if (TT$Extinction[i] == 0) TT$Extinction[i] <- TimeMax
    }

    biomass <- object@n
    for (iTime in seq_len(dim(biomass)[1])) {
        for (iPhen in seq_len(dim(biomass)[2])) {
            w_lim <- SumPar$w_mat[iPhen]
            S <- numeric(length(object@params@w))
            S[which.min(abs(w_lim - object@params@w))] <- 1
            NoW_mat <- which(S == 1)
            if (NoW_mat > 1) {
                biomass[iTime, iPhen, 1:(NoW_mat - 1)] <- 0
            }
        }
    }

    abundanceM <- apply(biomass, c(1, 2), sum)
    colnames(abundanceM) <- TT$Species
    abundanceM[is.na(abundanceM)] <- 0
    abundanceNormal <- matrix(0, nrow = dim(abundanceM)[1], ncol = dim(abundanceM)[2])

    for (iSpecies in SpIdx) {
        abundanceSp <- abundanceM
        abundanceSp[, which(colnames(abundanceM) != iSpecies)] <- 0
        abundanceSp <- sweep(abundanceSp, 1, apply(abundanceSp, 1, sum), "/")
        abundanceSp[is.nan(abundanceSp)] <- 0
        abundanceNormal <- abundanceNormal + abundanceSp
    }

    plotList <- vector(mode = "list", length = length(traitID))
    dataList <- vector(mode = "list", length = length(traitID))

    for (iTrait in traitID) {
        plotStore <- list()
        dataStore <- NULL
        abundanceT <- sweep(abundanceNormal, 2, as.matrix(SumPar[iTrait]), "*")

        max_sp_id <- max(as.numeric(unique(TT$Species)))
        TotMean <- matrix(0, nrow = dim(abundanceT)[1], ncol = max_sp_id,
                          dimnames = list(rownames(abundanceT), as.character(seq_len(max_sp_id))))
        names(dimnames(TotMean)) <- list("time", "species")

        for (i in SpIdx) {
            AMean <- abundanceT[, which(colnames(abundanceT) == i), drop = FALSE]
            if (!is.null(dim(AMean))) AMean <- apply(AMean, 1, sum)
            TotMean[, as.character(i)] <- AMean
        }

        statData <- list()

        for (i in SpIdx) {
            meanSp <- TotMean[, as.character(i)]
            traitSp <- dplyr::filter(TT, Species == i)[[iTrait]]
            weightSp <- abundanceNormal[, which(colnames(abundanceNormal) == i), drop = FALSE]
            statMat <- matrix(cbind(as.numeric(rownames(abundanceT)), meanSp, 0, 0, 0),
                              ncol = 5, nrow = length(meanSp),
                              dimnames = list(NULL, c("time", "mean", "sd", "percentMean", "percentSd")))
            for (itime in seq_along(meanSp)) {
                if (is.null(dim(weightSp)) || ncol(weightSp) == 1) {
                    variance <- sum(((traitSp - meanSp[itime])^2) * weightSp[itime, ])
                } else {
                    variance <- sum(((traitSp - meanSp[itime])^2) * weightSp[itime, ])
                }
                statMat[itime, 3] <- variance
                statMat[itime, 4] <- if (length(traitSp) > 0 && traitSp[1] != 0) (meanSp[itime] - traitSp[1]) / traitSp[1] else 0
                statMat[itime, 5] <- if (length(traitSp) > 0 && traitSp[1] != 0) statMat[itime, 3] / traitSp[1] else 0
            }
            statData[[as.character(i)]] <- as.data.frame(statMat)
        }

        for (i in SpIdx) {
            stat <- statData[[as.character(i)]]
            stat$species <- i
            dataStore <- rbind(dataStore, stat)
            plotStore[[as.character(i)]] <- NULL
        }

        plotList[[which(iTrait == traitID)]] <- plotStore
        dataList[[which(iTrait == traitID)]] <- dataStore
    }

    plot_dat <- dataList[[1]]
    p <- ggplot(plot_dat) +
        geom_smooth(aes(x = time, y = mean, group = as.factor(species), color = as.factor(species)), formula = y ~ x, method = "loess") +
        scale_x_continuous(name = "Time in years", limits = c(NA, NA)) +
        scale_y_continuous(trans = "log10") +
        scale_colour_manual(values = cbPalette, name = "Species") +
        theme(legend.title = element_text(),
              panel.background = element_rect(fill = "white", color = "black"),
              panel.grid.minor = element_line(colour = "grey92"),
              legend.position = "bottom",
              strip.background = element_blank(),
              legend.key = element_rect(fill = "white")) +
        guides(color = guide_legend(nrow = 1)) +
        ggtitle(NULL)

    if (returnData) return(dataList) else return(p)
}
