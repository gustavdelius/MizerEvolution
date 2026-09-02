#' mizerEvolution: Simulating Evolutionary Processes with Mizer
#'
#' With mizerEvolution you can simulate trait evolution and species invasion in
#' the mizer environment.
#'
#' @import mizer
#' @import ggplot2
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom stats pnorm rnorm runif time rlnorm
#' @importFrom grDevices colorRampPalette
#' @importFrom plyr aaply
#' @importFrom lifecycle deprecated is_present deprecate_warn
#' @keywords internal
"_PACKAGE"

globalVariables(c("inter", "value", "sp", "phen", "Species", "species",
                  "xint", "grp", "lineage", "size", "critical", "w",
                  "mean", "time", "trait", "w_mat", "mAmplitude", "eta",
                  "z0pre", "h", "f0", "resident_params"))
