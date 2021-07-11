#' Naive reference repertoires
#'
#' A list of naive reference repertoires. Each element of the list consists of
#' another list with a data frame containing the sequences and V-gene information, a data frame containing the V-gene
#' distribution and a data frame containing the CDR3 length distribution. 
#' The following database is available:
#' \itemize{
#' \item "gliph_reference": Reference database of 162,165 CDR3b sequences of naive human CD4+ or CD8+ T cells of two individuals used for the GLIPH paper.
# \item "human_v1.0_CD4": Reference database of 889,808 CDR3b sequences of naive human CD4+ T cells provided for GLIPH2.
# \item "human_v1.0_CD8": Reference database of 336,510 CDR3b sequences of naive human CD8+ T cells provided for GLIPH2.
# \item "human_v1.0_CD48": Reference database of 1,226,318 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD4": Reference database of 772,312 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD8": Reference database of 573,211 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD48": Reference database of 1,329,314 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD4": Reference database of 22,268 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD8": Reference database of 55,847 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD48": Reference database of 74,809 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
#' }
#' Additional reference databases were provided for download by the developers of GLIPH2 in the web tool (http://50.255.35.37:8080/tools).
#'
#' @docType data
#' @keywords datasets
#' @name reference_list
#' @usage data(reference_list)
#' @source Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
"reference_list"

#' GLIPH input data
#'
#' gliph_input_data, consists of around 2000 unique TCRs of known specificity
#' @docType data
#' @keywords datasets
#' @name gliph_input_data
#' @usage data(gliph_input_data)
#' @source Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
"gliph_input_data"

#' Cluster size probabilities in naive reference repertoire
#'
#' \code{"ref_cluster_sizes"} is a list containing two matrices with probabilities for cluster sizes in a naive reference
#' repertoire. The element named \code{"original"} contains the probabilities used in the original GLIPH algorithm. Here,
#' the same probabilities are used for all sample sizes.
#' The element names \code{"simulated"} contains probabilities determined by us. Since the distribution of cluster sizes
#' depend on the sample size, we estimated with GLIPH the probabilities for different sample sizes in a 500-step simulation
#' using different numbers of random sequences from the reference database (125, 250, 500, 1000, 2000, 4000, 6000, 8000,
#' 10000). For scoring, the probabilities from the dataset corresponding the closest to the sample size are used.
#'
#' @docType data
#' @keywords datasets
#' @name ref_cluster_sizes
#' @usage data(ref_cluster_sizes)
#' @source Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
"ref_cluster_sizes"

#' Germline TCR-beta CDR3 fragments
#'
#' \code{"gTRB"} is a list containing three data frames for fragments of germline encoded V genes (\code{"gTRV"}),
#' D genes (\code{"gTRD"}) and J genes (\code{"gTRJ"}) that might be located in the CDR3 region. These fragments are
#' used by GLIPH2 to determine germline-encoded sequence segments.
#'
#' @docType data
#' @keywords datasets
#' @name gTRB
#' @usage data(gTRB)
#' @source Lefranc, M.-P.
#' IMGT, the international ImMunoGeneTics database.
#' Nucl. Acids Res., 29(1):207-209 (2001). DOI:10.1093/nar/29.1.207. PMID:11125093.
"gTRB"

#' Amino acid pairs with non-negative BLOSUM62 score
#'
#' \code{"BlosumVec"} is a vector whose elements represent amino acid pairs in the one-letter code that have a
#' BLOSUM62 score >= 0. This information is used by GLIPH2 to cluster only global similarities whose different
#' amino acids are interpreted as interchangeable according to the BLOSUM62 matrix.
#'
#' @docType data
#' @keywords datasets
#' @name BlosumVec
#' @usage data(BlosumVec)
"BlosumVec"
