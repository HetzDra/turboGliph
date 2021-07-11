#' turboGliph: Grouping of Lymphocyte Interactions by Paratope Hotspots
#'
#'
#' @description GLIPH and GLIPH2 clusters TCRs that are predicted to bind
#'     the same HLA-restricted peptide antigen.
#'
#' @section Algorithm:
#' The scripts used in this package are written based on the original scripts of GLIPH which
#' is in perl. The output files produced by this package are exactly similar to the
#' results of perl script. GLIPH in R is parallelized and optimized to run faster. For
#' large input, \eqn{10^5} sequences, it is about 500 times faster than the perl code. For moderate
#' input, \eqn{10^4} sequences, it is about 50 times faster.
#'
#' Recently a new version of the package, GLIPH2, has been introduced where the main change in the
#' algorithm is the speed up of the analysis by substituting repeated sampling by Fisher's exact test
#' to assess the statistical significance of a given motif.
#' However in turboGliph we follow the original GLIPH but accelerate it by using appropriate R tools, efficient codes and parallelization.
#'
#' @source Glanville, Jacob, et al. "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @source Huang, Huang, et al. "Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening." Nature Biotechnology 38.10 (2020): 1194-1202.
#' @source Perl scripts from : https://github.com/immunoengineer/gliph (Aug. 2018)
#'
#' @docType package
#' @name turboGliph
NULL

