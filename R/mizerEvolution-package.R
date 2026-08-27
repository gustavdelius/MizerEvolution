#' mizerEvolution: Extends the mizer package with evolution and species invasion dynamics
#'
#' @description
#' This is an extension package for the
#' [mizer package](https://sizespectrum.org/mizer/)
#' that contains additional features
#' contributed by the mizer community. These features can then be further
#' improved while being used by the community. Once matured, frequently-used
#' features can be moved to the core mizer package.
#'
#' @import mizer dplyr ggplot2
#' @importFrom reshape2 melt
#' @importFrom stats pnorm rlnorm runif time
#' @importFrom methods is
#' @importFrom grDevices colorRampPalette
#' @importFrom plyr aaply
"_PACKAGE"

#' @export
reshape2::melt

globalVariables(c("inter","value","sp","phen","Species","species","xint","grp","lineage","size","critical"))

.onLoad <- function(libname, pkgname) {
  mizer::registerExtension(
    pkgname,
    requirement = "gustavdelius/MizerEvolution"
  )
}

record_mizer_evolution <- function(params) {
  version <- if ("mizerEvolution" %in% names(params@extensions)) {
    NULL
  } else {
    as.character(utils::packageVersion("mizerEvolution"))
  }
  mizer::recordExtension(params, "mizerEvolution", version = version)
}
