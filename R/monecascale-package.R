#' monecascale: Scalable Backends for MONECA Mobility Clustering
#'
#' Sibling package to \pkg{moneca} providing scalable clustering
#' backends whose output is consumed unchanged by moneca's analysis
#' and plotting stack. The load-bearing theoretical observation is
#' that moneca's relative-risk matrix \code{RR = O / E} is the
#' degree-corrected stochastic block model (DC-SBM) residual under
#' identity block interaction; SBM inference is therefore a
#' principled drop-in for moneca's clique-enumeration step, and it
#' scales to settings (10^6+ nodes) where clique enumeration is not
#' viable.
#'
#' @keywords internal
"_PACKAGE"
NULL
