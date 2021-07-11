#'  Grouping of Lymphocyte Interactions by Paratope Hotspots version 2
#'
#' @description Identifying specificity groups in the T cell receptor
#'  repertoire. Implementation of GLIPH2 following the instructions of the publication of Huang et al.
#' @param cdr3_sequences vector or dataframe. This dataframe must contain the cdr3 sequences and optional additional information.
#' The columns must be named as specified in the following list in arbitrary order.
#' \itemize{
#' \item "CDR3b": cdr3 sequences of beta chains
#' \item "TRBV": optional. V-genes of beta chains
#' \item "TRBV": optional. V-genes of beta chains
#' \item "patient": optional. Index of donor the appropriate sequence is obtained from. The value is composed of the index of the donor and
#' an optional experimental condition separated by a colon (example: 09/0410:MtbLys). For the calculation of the HLA-scores only the index
#' before the colon is used.
#' \item "HLA": optional. HLA alleles of the appropriate donor. The HLA alleles of a patient are separated by commas.
#' The standard notation of the HLA alleles is expected (example: DPA1*01:03).
#' For the calculation of HLA scores, information after the colon is neglected.
#' \item "counts": optional. Frequency of occurrence of the appropriate clone.
#' }
#' @param result_folder character. By default \code{""}. Path to the folder in which the output files should be stored.
#' If the value is \code{""} the results will not be saved in files.
#' @param refdb_beta character or data frame. By default \code{"gliph_reference"}. Specifies the reference database to be used.
#' For an individual reference database,  a data frame is expected as input. In its first column, the CDR3b sequences
#' must be specified and, if required, the V genes must be specified in the second column. Additional reference databases were provided for download
#' by the developers of GLIPH2 in the web tool (http://50.255.35.37:8080/tools).
#' To use the predefined database, the following keyword must be specified:
#' \itemize{
#' \item "gliph_reference": Reference database of 162,165 CDR3b sequences of naive human CD4+ or CD8+ T cells of two individuals used for the GLIPH paper.
# \item "human_v1.0_CD4": Reference database of 889,808 CDR3b sequences of naive human CD4+ T cells provided for GLIPH2.
# \item "human_v1.0_CD8": Reference database of 336,510 CDR3b sequences of naive human CD8+ T cells provided for GLIPH2.
# \item "human_v1.0_CD48": Reference database of 1,226,318 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD4": Reference database of 772,312 CDR3b sequences of naive human CD4+ T cells provided for GLIPH2.
# \item "human_v2.0_CD8": Reference database of 573,211 CDR3b sequences of naive human CD8+ T cells provided for GLIPH2.
# \item "human_v2.0_CD48": Reference database of 1,329,314 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD4": Reference database of 22,268 CDR3b sequences of naive murine CD4+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD8": Reference database of 55,847 CDR3b sequences of naive murine CD8+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD48": Reference database of 74,809 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
#' }
#' @param v_usage_freq data frame. By default \code{NULL}. This data frame contains the frequency of V-genes in a naive T cell repertoire required for
#' scoring. The first column provides the V-gene alleles and the second column the frequencies.  If the value is \code{NULL}, default frequencies are used.
#' @param cdr3_length_freq data frame. By default \code{NULL}. This data frame contains the frequency of CDR3 lengths in a naive T cell repertoire required for
#' scoring. The first column provides the CDR3 lengths and the second column the frequencies.  If the value is \code{NULL}, default frequencies are used.
#' @param ref_cluster_size character. Either \code{"original"} or \code{"simulated"}, by default \code{"original"}. Defines
#' the probabilities used for calculating the cluster size score. In the case of \code{"original"}, the standard probabilities of
#' the original algorithm which are constant for all sample sizes are used. However, since the distribution of cluster sizes depends on the sample size, we estimated
#' the probabilities for different sample sizes in a 500-step simulation using random sequences from the reference database.
#' To use these probabilities, the keyword \code{"simulated"} must be specified.
#' @param sim_depth numeric. By default 1000. Simulated resampling depth for assessing V gene and CDR3 length enrichment scores of clusters.
#' @param lcminp numeric. By default 0.01. Local convergence maximum probability score cutoff.
#' The score reports the probability that a random sample of the
#' same size as the sample set but of the reference set
#' (i.e. naive repertoire) would generate an enrichment of the
#' given motif at least as high as has been observed in the
#' sample set.
#' @param lcminove numeric. Local convergence minimum observed vs expected fold change.
#' This is a cutoff for the minimum fold enrichment over a
#' reference distribution that a given motif should have in
#' the sample set in order to be considered for further evaluation. By default, the minimum fold enrichment (1000,100,10) is
#' dependent on the motif length (2,3,4 amino acids). \code{lcminove} has to be either a single numeric value or a numeric
#' vector with equal length as \code{motif_length} representing the minimum fold enrichment depending on the respective motif_length.
#' @param motif_distance_cutoff numeric. By default 3. Defines the number of positions between which motifs for a local connection are allowed to vary.
#' @param kmer_mindepth numeric. By default 3. Minimum observations of kmer for it to be evaluated. This is
#' the minimum number of times a kmer should be observed in
#' the sample set in order for it to be considered for further
#' evaluation. The number can be set higher to provide less
#' motif-based clusters with higher confidence. This could be
#' recommended if the sample set is greater than 5000 reads. Lowering
#' the value to 2 will identify more groups but likely at a cost
#' of an increased False Discovery Rate.
#' @param accept_sequences_with_C_F_start_end logical. This logical flag
#' if \code{TRUE}, by default, only accepts sequences with
#' amino-acid C at the start position and amino-acid F at the
#' end position. This flag should be set to \code{FALSE} if
#' you wish to analyze sequences of different origin for example B-cells.
#' @param min_seq_length numeric. By default 8. All the sequences with a length less than this
#' value will be filtered out in input and reference database. If structboundaries
#' is \code{TRUE}, it is recommended not to go below the default. In this case,
#' the min_seq_length will be set to the maximum of 2*boundary_size+2 and
#' min_seq_length.
#' @param structboundaries logical. By default \code{TRUE}. By setting this flag to \code{TRUE}
#' the first boundary_size and the last boundary_size amino acids of each sequence will
#' not be considered in the analysis for computing the Hamming
#' distance and motif enrichment in input and reference database.
#' @param boundary_size numeric. By default 3. Specifies the boundary size if structboundaries is
#' active.
#' @param motif_length accepts a numeric vector of motif lengths
#' you want GLIPH2 to find and study. By default it searches for
#' motifs of size 2, 3 and 4 amino acids.
#' @param discontinuous_motifs logical. By default \code{FALSE}. Determines whether discontinuous_motifs motifs are to be considered.
#' @param local_similarities logical. By default \code{TRUE}. Determines whether the sequences should be analyzed for local similarity.
#' @param global_similarities logical. By default \code{TRUE}. Determines whether the sequences should be analyzed for global similarity.
#' @param global_vgene logical. By default \code{FALSE}. If \code{TRUE} global similarities are restricted to TCRs with shared V-gene.
#' Requires V-gene information in \code{cdr3_sequences}.
#' @param all_aa_interchangeable logical. By default \code{TRUE}.In the case of TRUE, all sequences with a Hamming distance of 1 are evaluated
#' as global similarities. In the case of FALSE, all sequences with a Hamming distance of 1 whose different amino acid has a BLOSUM62 score >= 0
#' are evaluated as global similarities.
#' @param boost_local_significance logical. By default \code{TRUE}. If set to \code{TRUE}, fisher scores of local clusters are repeatedly divided 
#' by 2 for every unique CDR3 sequence in the cluster in which the motif overlaps with non-germline encoded N- or P-nucleotides.
#' @param cluster_min_size numeric. By default 2. Minimal size of a cluster required to be considered for scoring.
#' @param hla_cutoff numeric. By default 0.1. Defines the threshold of HLA probability scores below which HLA alleles are considered significant.
#' @param n_cores numeric. Number of cores to use, by default 1. In case of \code{NULL} it will be set to number of cores in your machine minus 2.
#'
#' @return This function returns a list of six elements whose contents are explained below. If a file path is specified under \code{result_folder},
#' the results are additionally stored there. The individual file names are also specified below (italic name parts indicate the given value of the
#' corresponding parameter).
#'
#' @return $motif_enrichment:
#' A list of two data frames. \code{selected_motifs} contains only the motifs that pass the filtering criterion (ove and p-value),
#' whereas \code{all_motifs} contains p-value and ove of all motifs.\cr
#' File name of \code{selected_motifs}: local_similarities_minp_\code{lcminp}_ove \code{lcminove}_kmer_mindepth \code{kmer_mindepth}.txt
#' File name of \code{all_motifs}: all_motifs.txt
#'
#' @return $global_enrichment:
#' A data frame containing all identified global structures and their corresponding information.\cr
#' File name: global_similarities.txt
#'
#' @return $connections:
#' Contains the edge list. Each row consists of two nodes (cdr3 sequences) and a
#' third column which shows whether they are similar based on global or local similarity.
#' An additional fourth columns contains the cluster tag (motif or sequence structure), by which the sequences are clustered.\cr
#' File name: clone_network.txt
#'
#' @return $cluster_properties:
#' A data frame summarising the following information for each cluster:
#' \itemize{
#' \item "type": Indicates the type of similarity in the cluster (either global or local).
#' \item "tag": In the case of local similarities, the motif is indicated as well as the range of positions where the motif is positioned in the sequences. In the case of global similarities, the basic structure of the sequence is given as well as all amino acids separated by spaces that occur in the sample at the position marked by the % sign.
#' \item "cluster_size": Number of all sample sequences in the cluster.
#' \item "unique_cdr3_sample": Number of all unique CDR3b sequences of the sample in the cluster.
#' \item "unique_cdr3_ref": Number of all unique CDR3b sequences of the reference database matching the cluster properties.
#' \item "OvE": Factor of enrichment of the local or global motif in the sample compared to the reference database.
#' \item "fisher.score": The p-value obtained by performing the Fisher's exact test with a contingency table containing unique_cdr3_sample, unique_cdr3_ref, the number of remaining sample sequences and the number of remaining reference sequences.The score reports the probability that a random sample of the same size as the sample set but into the reference set (i.e. naive repertoire) would generate an enrichment of the given motif at least as high as has been observed in the sample set.
#' \item "members": All unique CDR3b sequences of the cluster separated by spaces.
#' \item "total.score": The product of all following scores.
#' \item "network.size.score": Probability of obtaining a cluster with this size in a naive repertoire.
#' \item "cdr3.length.score": enrichment of CDR3b lengths within the cluster.
#' \item "vgene.score": enrichment of V-genes within the cluster.
#' \item "clonal.expansion.score": enrichment of clonal expansion within the cluster.
#' \item "hla.score": enrichment of common HLA among donor TCR contributors in cluster.
#' \item "lowest.hlas": Enriched HLA alleles within the cluster.
#' } \cr
#' File name: convergence_groups.txt
#'
#'
#' @return $cluster_list:
#' A list containing the members and their additional information of each cluster. The elements of the list are named according to the appropriate cluster tag.\cr
#' File name: cluster_member_details.txt
#' 
#' @return $parameters:
#' A data frame containing all given input parameter values.\cr
#' File name: parameter.txt
#'
#' @examples
#' utils::data("gliph_input_data")
#' res <- gliph2(cdr3_sequences = gliph_input_data[base::seq_len(200),],
#' sim_depth = 50,
#' n_cores = 1)
#'
#' @references Huang, Huang, et al.
#' "Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening." Nature Biotechnology 38.10 (2020): 1194-1202.
#' @references \url{http://50.255.35.37:8080/}
#' @import foreach
#' @export
gliph2 <- function(cdr3_sequences,
                   result_folder = "",
                   # refdb_beta = "human_v1.0_CD4",
                   refdb_beta = "gliph_reference",
                   v_usage_freq = NULL,
                   cdr3_length_freq = NULL,
                   ref_cluster_size = "original",
                   sim_depth = 1000,
                   lcminp = 0.01,
                   lcminove = c(1000,100,10),
                   motif_distance_cutoff = 3,
                   kmer_mindepth = 3,
                   accept_sequences_with_C_F_start_end = TRUE,
                   min_seq_length = 0,
                   structboundaries = TRUE,
                   boundary_size = 3,
                   motif_length = base::c(2,3,4),
                   discontinuous_motifs = FALSE,
                   local_similarities = TRUE, #Draius13042020
                   global_similarities = TRUE, #Draius13042020,
                   global_vgene = FALSE, #Draius14042020
                   all_aa_interchangeable = FALSE,
                   boost_local_significance = TRUE,
                   cluster_min_size = 2,
                   hla_cutoff = 0.1,
                   n_cores = 1){
  t1 <- base::Sys.time()
  
  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################
  
  ### result_folder
  if(!base::is.character(result_folder))base::stop("result_folder has to be a character object")
  if(base::length(result_folder) > 1)base::stop("result_folder has to be a single path")
  save_results <- FALSE
  if(result_folder != ""){
    if(base::substr(result_folder,base::nchar(result_folder),base::nchar(result_folder)) != "/") result_folder <- base::paste0(result_folder,"/")
    if(!base::dir.exists(result_folder)) base::dir.create(result_folder)
    save_results <- TRUE
    
    if(base::file.exists(base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", base::paste(lcminove, collapse = "_"), "_kmer_mindepth_", kmer_mindepth, ".txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", base::paste(lcminove, collapse = "_"), "_kmer_mindepth_", kmer_mindepth, ".txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"all_motifs.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"all_motifs.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder, "global_similarities_minp_",lcminp, ".txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "global_similarities_minp_",lcminp, ".txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"all_global_similarities.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"all_global_similarities.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder, "clone_network.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "clone_network.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"convergence_groups.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"convergence_groups.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"parameter.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"parameter.txt"),"\n"))
    }
  }
  
  ### refdb_beta
  # if(!(refdb_beta %in% base::c("gliph_reference", "human_v1.0_CD4", "human_v1.0_CD8", "human_v1.0_CD48", "human_v2.0_CD4",
  #                        "human_v2.0_CD8", "human_v2.0_CD48", "mouse_v1.0_CD4", "mouse_v1.0_CD8", "mouse_v1.0_CD48")) &&
  #    !base::is.data.frame(refdb_beta)){
  if(!base::is.data.frame(refdb_beta)){
    if(base::length(refdb_beta) != 1 || !is.character(refdb_beta)){
      base::stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    } else if(!(refdb_beta %in% base::c("gliph_reference"))){
      base::stop("refdb_beta has to be a data frame (containing CDR3b sequences in the first column and optional V-gene information in the second column) or the value 'gliph_reference'")
    }
  }
  
  ### v_usage_freq
  if(!base::is.null(v_usage_freq)){
    if(base::is.data.frame(v_usage_freq)){
      if(base::ncol(v_usage_freq) < 2) base::stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(base::nrow(v_usage_freq) < 1) base::stop("v_usage_freq has to contain at least one row.")
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(v_usage_freq[,2])))) == TRUE){
        base::stop("The second column of v_usage_freq must contain the frequency of the corresponding V-gene in the first column in a naive T-cell repertoire.")
      } else v_usage_freq[,2] <- as.numeric(v_usage_freq[,2])
      
    } else {base::stop("v_usage_freq has to be a data frame containing V-gene information in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }
  
  ### cdr3_length_freq
  if(!base::is.null(cdr3_length_freq)){
    if(base::is.data.frame(cdr3_length_freq)){
      if(base::ncol(cdr3_length_freq) < 2) base::stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive  T-cell repertoire in the second column.")
      if(base::nrow(cdr3_length_freq) < 1) base::stop("cdr3_length_freq has to contain at least one row.")
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(cdr3_length_freq[,2])))) == TRUE){
        base::stop("The second column of cdr3_length_freq must contain the frequency of the corresponding CDR3 length in the first column in a naive T-cell repertoire.")
      } else cdr3_length_freq[,2] <- as.numeric(cdr3_length_freq[,2])
      
    } else {base::stop("cdr3_length_freq has to be a data frame containing CDR3 lengths in the first column and the corresponding frequency in a naive T-cell repertoire in the second column.")}
  }
  
  ### ref_cluster_size
  if(!(ref_cluster_size %in% base::c("original", "simulated") ||
       !base::is.character(ref_cluster_size) ||
       base::length(ref_cluster_size) > 1)){
    base::stop("ref_cluster_size has to be either 'original' or 'simulated'.")
  }
  
  ### sim_depth
  if(!base::is.numeric(sim_depth))base::stop("sim_depth has to be numeric")
  if(base::length(sim_depth) > 1)base::stop("sim_depth has to be a single number")
  if(sim_depth < 1)base::stop("sim_depth must be at least 1")
  sim_depth <- base::round(sim_depth)
  
  ### lcminp
  if(!base::is.numeric(lcminp))base::stop("lcminp has to be numeric")
  if(base::length(lcminp) > 1)base::stop("lcminp has to be a single number")
  if(lcminp <= 0)base::stop("lcminp must be greater than 0")
  
  ### lcminove
  if(!base::is.numeric(lcminove)) base::stop("lcminove has to be numeric")
  if(base::length(lcminove) > 1 && base::length(lcminove) != base::length(motif_length)) base::stop("lcminove has to be a single number or of same length as motif_length")
  if(base::any(lcminove < 1)) base::stop("lcminove must be at least 1")
  
  ### motif_distance_cutoff
  if(!base::is.numeric(motif_distance_cutoff)) base::stop("motif_distance_cutoff has to be numeric")
  if(base::length(motif_distance_cutoff) > 1) base::stop("motif_distance_cutoff has to be a single number")
  motif_diffs <- max(0, motif_distance_cutoff-1)
  
  ### kmer_mindepth
  if(!base::is.numeric(kmer_mindepth))base::stop("kmer_mindepth has to be numeric")
  if(base::length(kmer_mindepth) > 1)base::stop("kmer_mindepth has to be a single number")
  if(kmer_mindepth < 1)base::stop("kmer_mindepth must be at least 1")
  kmer_mindepth <- base::round(kmer_mindepth)
  
  ### accept_sequences_with_C_F_start_end
  if(!base::is.logical(accept_sequences_with_C_F_start_end))base::stop("accept_sequences_with_C_F_start_end has to be logical")
  
  ### min_seq_length
  if(!base::is.numeric(min_seq_length))base::stop("min_seq_length has to be numeric")
  if(base::length(min_seq_length) > 1)base::stop("min_seq_length has to be a single number")
  if(min_seq_length < 0)base::stop("min_seq_length must be at least 0")
  min_seq_length <- base::round(min_seq_length)
  
  ### structboundaries
  if(!base::is.logical(structboundaries))base::stop("structboundaries has to be logical")
  
  ### boundary_size
  if(!base::is.numeric(boundary_size))base::stop("boundary_size has to be numeric")
  if(base::length(boundary_size) > 1)base::stop("boundary_size has to be a single number")
  if(boundary_size < 0)base::stop("boundary_size must be at least 0")
  boundary_size <- base::round(boundary_size)
  if(structboundaries == TRUE) min_seq_length <- base::max(min_seq_length, ((boundary_size * 2)+1))
  
  ### motif_length
  if(!base::is.numeric(motif_length))base::stop("motif_length has to be numeric")
  if(base::any(motif_length < 1))base::stop("values of motif_length must be at least 1")
  motif_length <- base::round(motif_length)
  
  ### discontinuous_motifs
  if(!base::is.logical(discontinuous_motifs))base::stop("discontinuous_motifs has to be logical")
  
  ### local_similarities
  if(!base::is.logical(local_similarities))base::stop("local_similarities has to be logical")
  
  ### global_similarities
  if(!base::is.logical(global_similarities))base::stop("global_similarities has to be logical")
  if(local_similarities == FALSE && global_similarities == FALSE)base::stop("Either local_similarities or global_similarities have to be TRUE")
  
  ### global_vgene
  if(!base::is.logical(global_vgene))base::stop("global_vgene has to be logical")
  
  ### cluster_min_size
  if(!base::is.numeric(cluster_min_size))base::stop("cluster_min_size has to be numeric")
  if(base::length(cluster_min_size) > 1)base::stop("cluster_min_size has to be a single number")
  if(cluster_min_size < 1)base::stop("cluster_min_size must be at least 1")
  cluster_min_size <- base::round(cluster_min_size)
  
  ### all_aa_interchangeable
  if(!base::is.logical(all_aa_interchangeable))base::stop("all_aa_interchangeable has to be logical")
  
  ### hla_cutoff
  if(!base::is.numeric(hla_cutoff)) base::stop("hla_cutoff has to be numeric")
  if(base::length(hla_cutoff) > 1) base::stop("hla_cutoff has to be a single number")
  if(hla_cutoff > 1 || hla_cutoff < 0) base::stop("hla_cutoff must be between 0 and 1")
  
  ### n_cores
  if(!base::is.null(n_cores)) {
    if(!base::is.numeric(n_cores))base::stop("n_cores has to be numeric")
    if(base::length(n_cores) > 1)base::stop("n_cores has to be a single number")
    if(n_cores < 1)base::stop("n_cores must be at least 1")
    n_cores <- base::round(n_cores)
  }
  
  ### boost_local_significance
  if(!base::is.logical(boost_local_significance))base::stop("boost_local_significance has to be logical")
  
  ### min_seq_length second round
  if(structboundaries == TRUE) min_seq_length <- base::max(min_seq_length, 2*boundary_size+1)
  
  #################################################################
  ##                      Input preparation                      ##
  #################################################################
  
  vgene.info <- FALSE
  patient.info <- FALSE
  hla.info <- FALSE
  count.info <- FALSE
  sequences <- NULL
  
  ### Check cdr3_sequences for desired structure and information and convert it into a data frame
  if(base::is.atomic(cdr3_sequences)) {
    sequences <- base::data.frame(CDR3b = cdr3_sequences, stringsAsFactors = FALSE)
    if(global_vgene == TRUE){
      base::stop("To restrict global similarities to v-gene sharing, v-gene information is required as column of cdr3_sequences named \"TRBV\"")
    }
  } else if (base::is.data.frame(cdr3_sequences)) {
    if(base::ncol(cdr3_sequences) == 1){
      sequences <- base::data.frame(CDR3b = cdr3_sequences[,1],
                                    stringsAsFactors = FALSE)
      base::cat("Notification: cdr3_sequences is a data frame and the first column is considered as cdr3 sequences.\n")
    } else {
      if(!("CDR3b" %in% base::colnames(cdr3_sequences)))base::stop("cdr3_sequences has to contain a column named \"CDR3b\" comprising cdr3 beta sequences")
      sequences <- base::data.frame(CDR3b = cdr3_sequences$CDR3b,
                                    stringsAsFactors = FALSE)
      base::cat("Notification: cdr3_sequences is a data frame and the column named \"CDR3b\" is considered as cdr3 beta sequences.\n")
    }
    if("TRBV" %in% base::colnames(cdr3_sequences)){
      vgene.info <- TRUE
      sequences$TRBV <- cdr3_sequences$TRBV
    } else if(global_vgene == TRUE){
      base::stop("To restrict global similarities to v-gene sharing, v-gene information is required as column of cdr3_sequences named \"TRBV\"")
    }
    
    if("patient" %in% base::colnames(cdr3_sequences)){
      patient.info <- TRUE
      sequences$patient <- cdr3_sequences$patient
    }
    if("HLA" %in% base::colnames(cdr3_sequences)){
      hla.info <- TRUE
      sequences$HLA <- cdr3_sequences$HLA
    }
    if("counts" %in% base::colnames(cdr3_sequences)){
      count.info <- TRUE
      cdr3_sequences$counts <- base::as.numeric(cdr3_sequences$counts)
      cdr3_sequences$counts[base::is.na(cdr3_sequences$counts)] <- 1
      sequences$counts <- cdr3_sequences$counts
    }
    
    if(base::any(!(base::colnames(cdr3_sequences) %in% base::colnames(sequences)))) sequences <- base::cbind(sequences,
                                                                                                             cdr3_sequences[,!(base::colnames(cdr3_sequences) %in% base::colnames(sequences))])
    
    sequences[] <- base::lapply(sequences, base::as.character)
  } else {
    base::stop("cdr3_sequences must be a character vector or a data frame with sequences in a column named 'CDR3b' .\n")
  }
  if(base::ncol(sequences) == 1){
    sequences <- base::data.frame(CDR3b = sequences[base::grep("^[ACDEFGHIKLMNOPQRSTUVWY]*$", sequences$CDR3b),])
    if(base::nrow(sequences) == 0) base::stop("No CDR3b amino acid sequences found. Check if only letters of the one-letter code are used.")
  } else {
    sequences <- sequences[base::grep("^[ACDEFGHIKLMNPQRSTVWY]*$", sequences$CDR3b),]
    if(base::nrow(sequences) == 0) base::stop("No CDR3b amino acid sequences found. Check if only letters of the one-letter code are used.")
  }
  sequences <- base::cbind(base::data.frame(seq_ID = base::seq_len(base::nrow(sequences))), sequences)
  
  ### Reduce non-redundant sample sequences to the motif region
  seqs <- sequences$CDR3b
  if(accept_sequences_with_C_F_start_end) seqs <- base::grep(pattern = "^C.*F$",x = sequences$CDR3b,perl = TRUE,value = TRUE) else seqs <- base::unlist(sequences$CDR3b)
  if(base::length(seqs) == 0) base::stop("No CDR3b sequences starting with C and ending with F found. Adjust the parameter 'accept_sequences_with_C_F_start_end' for further procedure.")
  seqs <- base::unique(seqs)
  seqs <- seqs[base::nchar(seqs) >= min_seq_length]
  if(base::length(seqs) == 0) base::stop("No CDR3b sequences with a minimum length of 'min_seq_length' found. Adjust this parameter for further procedure.")
  if(structboundaries) motif_region <- base::substr(x = seqs,start = boundary_size + 1 ,stop = base::nchar(seqs) - boundary_size) else motif_region <- seqs
  
  #################################################################
  ##              Preparation of reference database              ##
  #################################################################
  
  ### load non-redundant reference repertoire
  if(base::is.data.frame(refdb_beta)) {
    refseqs <- refdb_beta
    refseqs[] <- base::lapply(refseqs, base::as.character)
    
    if(base::nrow(refseqs) <= 0) base::stop("User-created reference repertoire is of length 0. Without reference sequences the algorithm is not able to work.")
    
    if("CDR3b" %in% base::colnames(refseqs)){
      base::cat("Notification: Column of reference database named 'CDR3b' is considered as cdr3 sequences.\n")
    } else if(base::ncol(refseqs) > 1){
      base::cat("Notification: First column of reference database is considered as cdr3 sequences.\n")
      base::colnames(refseqs)[1] <- "CDR3b"
    }
    if(base::ncol(refseqs) > 1 && global_vgene == TRUE){
      if("TRBV" %in% base::colnames(refseqs)){
        base::cat("Notification: Column of reference database named 'TRBV' is considered as V gene information.\n")
      } else {
        base::cat("Notification: Second column of reference database is considered as V gene information.\n")
        base::colnames(refseqs)[2] <- "TRBV"
      }
    } else if(global_vgene == TRUE){
      base::stop("V-gene restriction for global similarities ('global_vgene' == TRUE) requires V-gene information in second column of reference database.\n")
    } else {
      refseqs$TRBV <- base::rep("", base::nrow(refseqs))
    }
    if(base::ncol(refseqs) == 1) refseqs <- base::cbind(refseqs, base::rep("", base::nrow(refseqs)))
    
    refseqs <- refseqs[, base::c("CDR3b", "TRBV")]
    
    if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[base::grep(pattern = "^C.*F$",x = refseqs[,1],perl = TRUE),]
    if(base::nrow(refseqs) == 0) base::stop("No reference CDR3b sequences starting with C and ending with F found. Adjust the parameter 'accept_sequences_with_C_F_start_end' for further procedure.")
    refseqs <- base::unique(refseqs)
    refseqs <- refseqs[base::which(base::nchar(refseqs[,1]) >= min_seq_length),]
    if(base::nrow(refseqs) == 0) base::stop("No reference CDR3b sequences with a minimum length of 'min_seq_length' found. Adjust this parameter for further procedure.")
    
    refseqs <- refseqs[base::grep("^[ACDEFGHIKLMNPQRSTVWY]*$", refseqs$CDR3b),]
    if(base::nrow(refseqs) == 0) base::stop("No reference CDR3b amino acid sequences found. Check if only letters of the one-letter code are used.")
  } else {
    reference_list <- NULL
    utils::data("reference_list",envir = base::environment(), package = "turboGliph")
    refseqs <- base::as.data.frame(reference_list[[refdb_beta]]$refseqs)
    refseqs[] <- base::lapply(refseqs, base::as.character)
    
    if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[base::grep(pattern = "^C.*F$",x = refseqs[,1],perl = TRUE),]
    refseqs <- base::unique(refseqs)
    refseqs <- refseqs[base::which(base::nchar(refseqs[,1]) >= min_seq_length),]
  }
  
  ### Reduce non-redundant reference sequences to the motif region
  if(structboundaries) refseqs_motif_region <- base::substr(x = base::unique(refseqs[,1]),start = boundary_size + 1,stop = base::nchar(base::unique(refseqs[,1]))-boundary_size) else refseqs_motif_region <- base::unique(refseqs[,1])
  
  if(base::length(refseqs_motif_region) < base::length(motif_region)) {
    base::stop("Analysis failed. Reference database must have more sequences than input data for the significance analysis to make sense and work.\n")
  }
  
  ### Initiate parallel execution
  if(!base::is.null(n_cores)){
    no_cores <- n_cores
  }else{
    if(parallel::detectCores() > 2) {no_cores <- parallel::detectCores() - 2}else{no_cores <- 1}
  }
  
  doParallel::registerDoParallel(no_cores)
  
  ##################################################################
  ##                  Part 1: Local similarities                  ##
  ##################################################################
  
  if(local_similarities == TRUE){
    base::cat(base::paste("Part 1: Searching for local similarities.\n"))
    
    ### Identify all motifs in the reference database
    
    # Divide the sequences equally among all cores
    overhang <- base::length(refseqs_motif_region) %% no_cores
    id_list <- base::list()
    last_id <- 0
    next_id <- 0
    steps <- (base::length(refseqs_motif_region)-overhang)/no_cores
    for(i in base::seq_len(no_cores)){
      next_id <- last_id+steps
      if(overhang > 0){
        next_id <- next_id+1
        overhang <- overhang-1
      }
      id_list[[i]] <- (last_id+1):next_id
      last_id <- next_id
    }
    
    # Receive all motifs in the reference sequences
    ref_motifs_list <- foreach::foreach(i = base::seq_len(no_cores)) %dopar% {
      return(find_motifs(seqs = refseqs_motif_region[id_list[[i]]],
                                     q = motif_length, discontinuous = discontinuous_motifs))
    }
    
    # Convert the list in a more manageable data frame
    ref_motifs_df <- c()
    for(i in base::seq_len(no_cores)){
      if(i == 1) ref_motifs_df <- ref_motifs_list[[i]] else ref_motifs_df <- dplyr::full_join(x = ref_motifs_df, y = ref_motifs_list[[i]], by = "motif")
    }
    ref_motifs_df[base::is.na(ref_motifs_df)] <- 0
    ref_motifs_df <- data.frame(motif = ref_motifs_df$motif, count = base::rowSums(base::data.frame(ref_motifs_df[,-1])))
    
    ### Identify all motifs in the sample set
    
    # Divide the sequences equally among all cores
    overhang <- base::length(motif_region) %% no_cores
    id_list <- base::list()
    last_id <- 0
    next_id <- 0
    steps <- (base::length(motif_region)-overhang)/no_cores
    for(i in base::seq_len(no_cores)){
      next_id <- last_id+steps
      if(overhang > 0){
        next_id <- next_id+1
        overhang <- overhang-1
      }
      id_list[[i]] <- (last_id+1):next_id
      last_id <- next_id
    }
    
    # Receive all motifs in the sample sequences
    motifs_list <- foreach::foreach(i = base::seq_len(no_cores)) %dopar% {
      return(find_motifs(seqs = motif_region[id_list[[i]]],
                                     q = motif_length, discontinuous = discontinuous_motifs))
    }
    
    # Convert the list in a more manageable data frame
    motifs_df <- c()
    for(i in base::seq_len(no_cores)){
      if(i == 1) motifs_df <- motifs_list[[i]] else motifs_df <- dplyr::full_join(x = motifs_df, y = motifs_list[[i]], by = "motif")
    }
    motifs_df[base::is.na(motifs_df)] <- 0
    motifs_df <- base::data.frame(motif = motifs_df$motif, count = base::rowSums(base::data.frame(motifs_df[,-1])))
    
    ### Comparison between motifs found in sample and reference sequences
    
    # summarize motifs found in the sample sequences and their counts in the sample and reference set
    motifs_df <- dplyr::left_join(x = motifs_df, y = ref_motifs_df, by = "motif")
    motifs_df[base::is.na(motifs_df)] <- 0
    base::colnames(motifs_df) <- base::c("motif", "num_in_sample", "num_in_ref")
    
    # assess the significance of motifs for being enriched in the sample set by Fisher's Exact Test (alternative = "greater")
    # Contigency table:
    # ---------------------------------------------------------------------------------------------------------
    #   Number sample sequences containing the motif   |   Number of reference sequences containing the motif
    # ---------------------------------------------------------------------------------------------------------
    # Number sample sequences not containing the motif | Number of reference sequences not containing the motif
    # ---------------------------------------------------------------------------------------------------------
    motifs_df$fisher.score <- stats::phyper(q=motifs_df$num_in_sample-1,
                                            m=motifs_df$num_in_ref+motifs_df$num_in_sample,
                                            n=base::length(base::unique(refseqs$CDR3b))+base::length(base::unique(seqs))-motifs_df$num_in_ref-motifs_df$num_in_sample,
                                            k=base::length(base::unique(seqs)), lower.tail = FALSE)
    motifs_df$fisher.score <- base::as.numeric(base::formatC(motifs_df$fisher.score, digits = 1, format = "e"))
    
    # assess the fold change enrichement ot a motif in the sample set compared to the reference set
    motifs_df$num_fold <- base::round(motifs_df$num_in_sample/(motifs_df$num_in_ref+0.01)/base::length(base::unique(seqs))*base::length(base::unique(refseqs$CDR3b)), digits = 1)
    
    ### Get the minimum fold change for every motif depending on the motif length
    temp_minove <- base::c()
    if(base::length(lcminove) == 1){
      temp_minove <- base::rep(lcminove, base::nrow(motifs_df))
    } else {
      temp_minove <- base::rep(0, base::nrow(motifs_df))
      
      # the number of informative positions in discontinuous motifs is reduced by one compared to continuous motifs with same length
      nam_len <- base::nchar(motifs_df$motif)
      nam_len[base::grep(".", motifs_df$motif, fixed = TRUE)] <- nam_len[base::grep(".", motifs_df$motif, fixed = TRUE)]-1
      
      for(i in base::seq_along(motif_length)){
        temp_minove[nam_len == motif_length[i]] <- lcminove[i]
      }
    }
    
    ### Filter significantly enriched motifs by the following criteria:
    # p-value < lcminp (p-value of the Fisher's Exact Test carried out for the contigency table shown above)
    # fold change >= lcminove (lcminove can be dependent on the motif length; fold change of sample set occurence compared to repeated random subsample occurence)
    # sample set occurence >= kmer_mindepth
    sel_mot_df <- motifs_df[motifs_df$num_fold >= temp_minove &
                            motifs_df$fisher.score <= lcminp &
                            motifs_df$num_in_sample >= kmer_mindepth,]
    
    if(base::nrow(sel_mot_df) > 0){
      ### Preclustering
      # Motifs may only be within a certain distance to be assigned to a cluster (determined by the parameter 'motif_distance_cutoff').
      # Find the range(s) for every motif.
      selected_motifs_pos <- foreach::foreach(i = base::seq_len(base::nrow(sel_mot_df)), .combine = base::c) %dopar% {
        
        # get all N-terminal positions of the motifs in the sample set (starting and stop position)
        all_pos <- stringr::str_locate_all(string = motif_region, pattern = sel_mot_df$motif[i])
        all_pos_ids <- base::which(base::lengths(all_pos) > 0)
        all_pos_df <- base::matrix(base::unlist(all_pos), ncol = 2, byrow = TRUE)
        all_pos_df <- base::data.frame(id = rep(all_pos_ids, times = base::lengths(all_pos[all_pos_ids])/2),
                                       start = all_pos_df[,1],
                                       stop = all_pos_df[,2])
        
        # Find chains of motif positions falling below the maximum distance by sorting all positions.
        
        unq_starts <- base::sort(base::unique(all_pos_df$start))
        unq_starts_rotated <- base::c(unq_starts[-1], base::max(unq_starts))
        
        # Then compare every position with the next heighest position.
        breakpoints <- base::which((unq_starts_rotated-unq_starts) > motif_diffs)
        
        # If the maximum distance from one position to the next is exceeded, the range is split at this point.
        if(motif_distance_cutoff < 1 || base::length(breakpoints) == 0){
          temp_seqs <- base::sort(seqs[base::unique(all_pos_df$id)])
          return(base::c(i,
                         0,
                         0,
                         base::length(temp_seqs),
                         paste(temp_seqs, collapse = " ")))
        } else {
          ret <- base::c()
          min_start <- base::min(all_pos_df$start)
          breakpoints <- base::unique(base::c(breakpoints, base::length(unq_starts)))
          for(j in breakpoints){
            # receive all sequences within a certain position range
            temp_seqs <- base::sort(seqs[base::unique(all_pos_df$id[all_pos_df$start >= min_start & all_pos_df$start <= unq_starts[j]])])
            ret <- base::c(ret, i,
                           min_start,
                           unq_starts[j],
                           base::length(temp_seqs),
                           paste(temp_seqs, collapse = " "))
            min_start <- unq_starts_rotated[j]
          }
          # return for every range the motif ID, start position, stop position, number of sequences with the actual motif in this range and the sequences themself
          return(ret)
        }
      }
      
      ### merge significantly enriched motifs with the division into position ranges
      selected_motifs_pos <- base::data.frame(base::matrix(selected_motifs_pos, ncol = 5, byrow = TRUE))
      base::colnames(selected_motifs_pos) <- base::c("motifID", "start", "stop", "counts", "members")
      selected_motifs_pos$motifID <- base::as.numeric(selected_motifs_pos$motifID)
      selected_motifs_pos$start <- base::as.numeric(selected_motifs_pos$start)
      selected_motifs_pos$stop <- base::as.numeric(selected_motifs_pos$stop)
      selected_motifs_pos$counts <- base::as.numeric(selected_motifs_pos$counts)
      
      selected_clusters <- sel_mot_df[selected_motifs_pos$motifID,]
      selected_clusters$start <- selected_motifs_pos$start
      selected_clusters$stop <- selected_motifs_pos$stop
      selected_clusters$num_in_sample <- selected_motifs_pos$counts
      selected_clusters$members <- selected_motifs_pos$members
      
      selected_clusters$start[selected_clusters$start == 0] <- 1
      selected_clusters$stop[selected_clusters$stop == 0] <- base::max(base::nchar(motif_region))
      
      positions_to_add <- structboundaries*boundary_size
      selected_clusters$start <- selected_clusters$start + positions_to_add
      selected_clusters$stop <- selected_clusters$stop + positions_to_add
      
      # carry out Fisher's Exact Test with the new cluster size number
      selected_clusters$fisher.score <- stats::phyper(q=selected_clusters$num_in_sample-1,
                                                      m=selected_clusters$num_in_ref+selected_clusters$num_in_sample,
                                                      n=base::length(base::unique(refseqs$CDR3b))+base::length(base::unique(seqs))-selected_clusters$num_in_ref-selected_clusters$num_in_sample,
                                                      k=base::length(base::unique(seqs)), lower.tail = FALSE)
      selected_clusters$fisher.score <- base::as.numeric(base::formatC(selected_clusters$fisher.score, digits = 1, format = "e"))
      
      base::cat(base::nrow(sel_mot_df),"significantly enriched motifs found in sample set.\n")
    } else {
      selected_clusters <- base::data.frame(motif = 1,
                                            num_in_sample = 1,
                                            num_in_ref = 1,
                                            fisher.score = 1,
                                            num_fold = 1,
                                            start = 1,
                                            stop = 1,
                                            members = 1)[-1,]
      local_similarities <- FALSE
      
      base::cat("No significantly enriched motifs found in sample set.\n")
    }
    
    
    
    ### save results
    part_1_res <- selected_clusters
    part_1_res_all_clusters <- motifs_df
    
    part_1_res_file_name <- base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", paste(lcminove, collapse = "_"), "_kmer_mindepth_", kmer_mindepth, ".txt")
    part_1_res_all_clusters_file_name <- base::paste0(result_folder, "all_motifs.txt")
    if(save_results == TRUE){
      utils::write.table(x = part_1_res,file = part_1_res_file_name,quote = FALSE,sep = "\t",row.names = FALSE)
      utils::write.table(x = part_1_res_all_clusters,file = part_1_res_all_clusters_file_name,quote = FALSE,sep = "\t",row.names = FALSE)
    }
    
    part_1_res <- base::data.frame(part_1_res, stringsAsFactors = FALSE)
    part_1_res_all_clusters <- base::data.frame(part_1_res_all_clusters, stringsAsFactors = FALSE)
    
    t2_e <- base::Sys.time()
    base::cat("Part 1 cpu time:",t2_e-t1,base::units(t2_e-t1),"\n\n")
  } else {
    part_1_res <- base::c()
    part_1_res_all_clusters <- base::c()
    t2_e <- base::Sys.time()
    base::cat(base::paste("\nLocal similarity search is skipped. \n\n"))
  }
  
  #################################################################
  ##                 Part 2: Global similarities                 ##
  #################################################################
  if(global_similarities == TRUE){
    base::cat(base::paste("Part 2: Searching for global similarities.\n"))
    
    ### load all amino acid pairs with a BLOSUM62-value >= 0
    BlosumVec <- NULL
    utils::data("BlosumVec",envir = base::environment(), package = "turboGliph")
    
    part_2_res_file_name <- base::paste0(result_folder, "global_similarities.txt")
    
    ### Filter sample sequences, extract the motif region (struc) and add some columns for the next steps
    sample_seqs <- sequences$CDR3b
    if(accept_sequences_with_C_F_start_end) sample_seqs <- base::grep(pattern = "^C.*F$",x = sample_seqs,perl = TRUE, value = TRUE)
    sample_seqs <- base::unique(sample_seqs)
    sample_seqs <- sample_seqs[base::nchar(sample_seqs) >= min_seq_length]
    if(structboundaries){
      sample_seqs <- base::data.frame(seq = sample_seqs,
                                      struct = base::substr(x = sample_seqs,start = boundary_size + 1 ,stop = base::nchar(sample_seqs) - boundary_size))
    } else sample_seqs <- sample_seqs <- base::data.frame(seq = sample_seqs,
                                                           struct = sample_seqs)
    sample_seqs$nchar <- base::nchar(sample_seqs$struct)
    sample_seqs$pos <- rep(0, base::nrow(sample_seqs))
    sample_seqs$aa_at_pos <- rep("", base::nrow(sample_seqs))
    sample_seqs$tag <- sample_seqs$struct
    
    ### Receive all structs (motif region with one variable position) that occurs at least twice for every sample sequence
    exp_sample_seqs <- foreach::foreach(i = base::seq_len(base::max(sample_seqs$nchar)), .combine = base::rbind) %dopar% {
      temp_part <- sample_seqs[sample_seqs$nchar >= i,]
      temp_part$pos <- i
      temp_part$aa_at_pos <- base::substr(temp_part$struct,i,i)
      base::substr(temp_part$tag,i,i) <- "%"
      
      # this filters all structs that occur at least twice
      dup_check <- base::duplicated(temp_part$tag) | base::duplicated(temp_part$tag, fromLast = TRUE)
      
      return(temp_part[dup_check,])
    }
    
    # get all unique sample structs
    unq_smaple_struct <- base::unique(exp_sample_seqs$tag)
    
    if(base::length(unq_smaple_struct) > 0){
      ### Extract the motif region (struc) of reference sequences and add some columns for the next steps
      reference_seqs <- base::unique(refseqs$CDR3b)
      if(structboundaries){
        reference_seqs <- base::data.frame(seq = reference_seqs,
                                           struct = base::substr(x = reference_seqs,start = boundary_size + 1 ,stop = base::nchar(reference_seqs) - boundary_size))
      } else sample_seqs <- reference_seqs <- base::data.frame(seq = reference_seqs,
                                                               struct = reference_seqs)
      reference_seqs$nchar <- base::nchar(reference_seqs$struct)
      reference_seqs$pos <- rep(0, base::nrow(reference_seqs))
      reference_seqs$aa_at_pos <- rep("", base::nrow(reference_seqs))
      reference_seqs$tag <- reference_seqs$struct
      
      ### Receive all structs (motif region with one variable position) of reference sequences that occur in the sample set
      exp_reference_seqs <- foreach::foreach(i = base::seq_len(base::min(base::max(reference_seqs$nchar), base::max(sample_seqs$nchar))), .combine = base::rbind) %dopar% {
        temp_part <- reference_seqs[reference_seqs$nchar >= i,]
        temp_part$pos <- i
        temp_part$aa_at_pos <- base::substr(temp_part$struct,i,i)
        base::substr(temp_part$tag,i,i) <- "%"
        
        # this filters all structs that occur in the sample set
        temp_part <- temp_part[temp_part$tag %in% unq_smaple_struct,]
        
        return(temp_part)
      }
      
      ### Summarize the frequencies of structs in sample and reference set
      sample_stats <- base::data.frame(base::table(exp_sample_seqs$tag), stringsAsFactors = FALSE)
      base::colnames(sample_stats) <- c("tag", "Freq")
      if(base::nrow(exp_reference_seqs) == 0){
        ref_stats <- base::data.frame(1, 1)[-1,]
      } else {
        ref_stats <- base::data.frame(base::table(exp_reference_seqs$tag), stringsAsFactors = FALSE)
      }
      base::colnames(ref_stats) <- c("tag", "Freq")
      sample_stats <- dplyr::left_join(sample_stats, ref_stats, by = "tag")
      sample_stats$tag <- base::as.character(sample_stats$tag)
      sample_stats[base::is.na(sample_stats)] <- 0
      base::colnames(sample_stats) <- c("tag", "num_in_sample", "num_in_ref")
      # sample_stats$num_in_ref[sample_stats$num_in_ref == 0] <- 1
      
      ### Add V gene information to all sequences and their corresponding structs
      seqs_w_vgenes <- sequences
      if(global_vgene == FALSE) seqs_w_vgenes$TRBV <- base::rep(" ", base::nrow(seqs_w_vgenes))
      exp_sample_seqs <- dplyr::left_join(exp_sample_seqs, seqs_w_vgenes[, base::c("CDR3b", "TRBV")], by = base::c("seq" = "CDR3b"))
      
      ### Iterate over all structs and filter (if requested) for identical V genes or/and amino acid pairs with a BLOSUM62 value >= 0
      edges <- foreach::foreach(i = base::seq_len(base::nrow(sample_stats)), .combine = base::rbind) %dopar% {
        # get all sequences of current struct
        act_members <- base::unique(exp_sample_seqs[exp_sample_seqs$tag == sample_stats$tag[i],])
        act_members$summary <- base::paste(act_members$seq, act_members$tag, act_members$aa_at_pos, act_members$pos, act_members$TRBV, sep = "#")
        
        # first assume a connection between all this sequences
        if(base::nrow(act_members) >= 2){
          combn_ids <- base::t(utils::combn(x = base::seq_len(base::nrow(act_members)), m = 2))
        } else {
          combn_ids <- base::t(utils::combn(base::rep(1, 2), m = 2))
        }
        ret <- base::data.frame(V1 = act_members$summary[combn_ids[,1]], V2 = act_members$summary[combn_ids[,2]], stringsAsFactors = FALSE)
        
        # filter for identical V genes or replaceable amino acids allowed by BLOSUM62
        keep <- rep(TRUE, base::nrow(ret))
        if(global_vgene == TRUE) keep <- keep & (act_members$TRBV[combn_ids[,1]] == act_members$TRBV[combn_ids[,2]])
        if(all_aa_interchangeable == FALSE) keep <- keep & (base::paste0(act_members$aa_at_pos[combn_ids[,1]], act_members$aa_at_pos[combn_ids[,2]]) %in% BlosumVec)
        return(ret[keep,])
      }
      
      ### Build a cluster network based on all connections
      gr <- igraph::graph_from_edgelist(el = base::as.matrix(edges),directed = FALSE)
      cm <- igraph::components(gr)
      
      if(cm$no > 0){
        ### Extract the clusters, cluster members and cluster sizes (number of unique CDR3b sequences) from the cluster network
        cm_splitted <- base::data.frame(base::matrix(base::unlist(base::strsplit(base::names(cm$membership), split = "#")), ncol = 5, byrow = TRUE), stringsAsFactors = FALSE)
        base::colnames(cm_splitted) <- base::c("CDR3b", "tag", "aa_at_position", "position", "TRBV")
        cm_splitted$position <- base::as.numeric(cm_splitted$position)
        
        cluster_list <- foreach::foreach(i = base::seq_len(cm$no)) %dopar% {
          csize <- cm$csize[i]
          member_df <- cm_splitted[base::which(cm$membership == i),]
          member_df <- base::unique(member_df)
          num_cdr3s <- base::length(base::unique(member_df$CDR3b))
          
          base::return(base::c(member_df$tag[1], csize, num_cdr3s,
                               base::paste(member_df$CDR3b, collapse = " "),
                               base::paste(base::sort(base::unique(member_df$aa_at_position)), collapse = ""), member_df$TRBV[1]))
        }
        cluster_list <- base::data.frame(base::matrix(base::unlist(cluster_list), ncol = 6, byrow = TRUE), stringsAsFactors = FALSE)
        base::colnames(cluster_list) <- base::c("cluster_tag", "cluster_size", "unique_CDR3b", "CDR3b", "aa_at_position", "TRBV")
        cluster_list$cluster_size <- base::as.numeric(cluster_list$cluster_size)
        cluster_list$unique_CDR3b <- base::as.numeric(cluster_list$unique_CDR3b)
        
        ### Add reference frequency to the clusters
        cluster_list <- dplyr::left_join(cluster_list, sample_stats[, c("tag", "num_in_ref")], by = base::c("cluster_tag" = "tag"))
        
        ### Asses significantly enriched structs by applying Fisher's Exact Test with the following contigency table:
        # --------------------------------------------------------------------------------------------------------
        #                     Cluster size               |   Number of reference sequences containing the struct
        # --------------------------------------------------------------------------------------------------------
        # Number sample sequences not within the cluster | Number of reference sequences not containing the struct
        # --------------------------------------------------------------------------------------------------------
        cluster_list$fisher.score <- stats::phyper(q=cluster_list$unique_CDR3b-1,
                                                   m=cluster_list$num_in_ref + cluster_list$unique_CDR3b,
                                                   n=base::length(base::unique(refseqs$CDR3b))+base::length(base::unique(sample_seqs$seq))-cluster_list$num_in_ref-cluster_list$unique_CDR3b,
                                                   k=base::length(base::unique(sample_seqs$seq)), lower.tail = FALSE)
        
        cluster_list <- cluster_list[,c("cluster_tag", "cluster_size", "unique_CDR3b", "num_in_ref", "fisher.score", "aa_at_position", "TRBV","CDR3b")]
        
        base::cat(base::nrow(cluster_list), "clusters based on global similarity found in sample set.\n")
      } else {
        cluster_list <- base::data.frame(cluster_tag = 1,
                                         cluster_size = 1,
                                         unique_CDR3b = 1,
                                         num_in_ref = 1,
                                         fisher.score = 1,
                                         aa_at_position = 1,
                                         TRBV = 1,
                                         CDR3b = 1)[-1,]
        global_similarities <- FALSE
        base::cat("No global similarities found in sample set.\n")
      }
    } else {
      cluster_list <- base::data.frame(cluster_tag = 1,
                                       cluster_size = 1,
                                       unique_CDR3b = 1,
                                       num_in_ref = 1,
                                       fisher.score = 1,
                                       aa_at_position = 1,
                                       TRBV = 1,
                                       CDR3b = 1)[-1,]
      global_similarities <- FALSE
      base::cat("No global similarities found in sample set.\n")
    }
    
    ### save
    part_2_res <- cluster_list
    
    if(save_results == TRUE){
      utils::write.table(x = part_2_res,file = part_2_res_file_name,quote = FALSE,sep = "\t",row.names = FALSE)
    }
    
    t3_e <- base::Sys.time()
    base::cat("Part 2 cpu time:",t3_e-t2_e,base::units(t3_e-t2_e),"\n\n")
  } else {
    part_2_res <- base::c()
    t3_e <- base::Sys.time()
    base::cat(base::paste("\nGlobal similarity search is skipped. \n\n"))
  }
  
  ##################################################################
  ##                      Part 3: Clustering                      ##
  ##################################################################
  
  base::cat(base::paste("Part 3: Clustering sequences.\n"))
  
  ### Identification of N-nucleotides in CDR3b-sequences
  if(boost_local_significance == TRUE){
    # load all CDR3b sequence fragments encoded in the germline
    gTRB <- NULL
    utils::data("gTRB", package = "turboGliph", envir = base::environment())
    
    # framework for the following procedure. Determines up to two regions per sequence that are (probably) not encoded by germline
    range_df <- base::data.frame(start_1 = base::rep(0, base::nrow(sequences)),
                                 stop_1 = base::rep(0, base::nrow(sequences)),
                                 start_2 = base::rep(0, base::nrow(sequences)),
                                 stop_2 = base::rep(0, base::nrow(sequences)))
    
    # For each sequence, find the maximum number of N-terminal positions encoded by a germline V gene 
    remaining_ids <- base::seq_len(base::nrow(sequences))
    for(i in base::max(gTRB$gTRBV[,2]):base::min(gTRB$gTRBV[,2])){
      temp_bool <- (base::substr(sequences$CDR3b[remaining_ids], 1, i) %in% gTRB$gTRBV[,1][gTRB$gTRBV[,2] == i])
      range_df$start_1[remaining_ids][temp_bool] <- i+1
      
      remaining_ids <- remaining_ids[!temp_bool]
      if(base::length(remaining_ids) == 0) break()
    }
    
    # For each sequence without V gene part, find the maximum number of N-terminal positions encoded by a germline J gene
    remaining_ids <- base::seq_len(base::nrow(sequences))
    for(i in base::max(gTRB$gTRBJ[,2]):base::min(gTRB$gTRBJ[,2])){
      temp_bool <- (base::substr(sequences$CDR3b[remaining_ids], base::nchar(sequences$CDR3b[remaining_ids])-i+1, base::nchar(sequences$CDR3b[remaining_ids])) %in% gTRB$gTRBJ[,1][gTRB$gTRBJ[,2] == i])
      range_df$stop_2[remaining_ids][temp_bool] <- i
      
      remaining_ids <- remaining_ids[!temp_bool]
      if(base::length(remaining_ids) == 0) break()
    }
    range_df$stop_2 <- base::nchar(sequences$CDR3b)-range_df$stop_2
    
    # For each sequence without V gene and J gene part, find the maximum number of N-terminal positions encoded by a germline D gene
    remaining_ids <- base::seq_len(base::nrow(sequences))
    if(base::length(remaining_ids) > 0){
      
      # get sequence parts without V & J gene 
      act_subseqs <- base::substr(sequences$CDR3b[remaining_ids], range_df$start_1[remaining_ids], range_df$stop_2[remaining_ids])
      
      # iterate over all D gene fragements lengths (from high to low)
      for(i in base::max(gTRB$gTRBD[,2]):base::min(gTRB$gTRBD[,2])){
        act_subseqs_matrix <- base::c()
        
        if((base::max(base::nchar(act_subseqs))-i+1) >= 1){
          # build a matrix with sequence fragments with length i starting at N-terminal position j
          for(j in base::seq_len(base::max(base::nchar(act_subseqs))-i+1)){
            if(j == 1) act_subseqs_matrix <- base::matrix(data = base::substr(act_subseqs, 1, i), ncol = 1) else act_subseqs_matrix <- base::cbind(act_subseqs_matrix, base::substr(act_subseqs, j, j+i-1))
          }
          
          # Identify all D gene fragments in the sequence fragments, prioritize more N-terminal fragments
          for(j in gTRB$gTRBD[gTRB$gTRBD[,2] == i,1]){
            
            # in which most N-terminal position occurs a D gene fragment
            act_values <- base::do.call(base::pmax, base::data.frame(base::t(base::t(act_subseqs_matrix == j)*(base::ncol(act_subseqs_matrix):1))))
            act_values[act_values > 0] <- base::ncol(act_subseqs_matrix)-act_values[act_values > 0]+1
            
            # save the starting and stop positions of the N-nucleotide positions considering the V, J and D gene fragments
            range_df$start_1[remaining_ids][act_values == 1] <- range_df$start_1[remaining_ids][act_values == 1]+i
            range_df$stop_1[remaining_ids][act_values > 1] <- range_df$start_1[remaining_ids][act_values > 1]+act_values[act_values > 1]-2
            range_df$start_2[remaining_ids][act_values > 1] <- range_df$stop_1[remaining_ids][act_values > 1]+i+1
            temp_bool <- (act_values > 0)
            
            remaining_ids <- remaining_ids[!temp_bool]
            act_subseqs <- act_subseqs[!temp_bool]
            act_subseqs_matrix <- act_subseqs_matrix[!temp_bool,]
            if(base::length(remaining_ids) == 0) break()
          }
          if(base::length(remaining_ids) == 0) break()
        }
      }
    }
  }
  
  sequences$ultCDR3b <- rep("", nrow(sequences))
  
  ### Represent the non-germline coded sequence parts by small letters
  if(boost_local_significance == TRUE){
    sequences$ultCDR3b <- sequences$CDR3b
    sequences$ultCDR3b <- base::toupper(sequences$ultCDR3b)
    base::substr(sequences$ultCDR3b, range_df$start_1, range_df$stop_2) <- base::tolower(base::substr(sequences$ultCDR3b, range_df$start_1, range_df$stop_2))
    base::substr(sequences$ultCDR3b, range_df$stop_1+1, range_df$start_2-1) <- base::toupper(base::substr(sequences$ultCDR3b, range_df$stop_1+1, range_df$start_2-1))
  }
  
  
  ### Combine all information of global and local clusters in one data frame
  merged_clusters <- base::c()
  if(local_similarities == TRUE && base::nrow(part_1_res) > 0){
    local_clusters <- base::data.frame(type = base::rep("local", base::nrow(part_1_res)), stringsAsFactors = FALSE)
    local_clusters$tag <- base::paste(part_1_res$motif, part_1_res$start, part_1_res$stop, sep = "_")
    local_clusters$cluster_size <- base::rep(0, base::nrow(part_1_res))
    local_clusters$unique_cdr3_sample <- part_1_res$num_in_sample
    local_clusters$unique_cdr3_ref <- part_1_res$num_in_ref
    local_clusters$OvE <- part_1_res$num_fold
    local_clusters$fisher.score <- part_1_res$fisher.score
    local_clusters$members <- part_1_res$members
    
    merged_clusters <- local_clusters
  }
  if(global_similarities == TRUE && base::nrow(part_2_res) > 0){
    global_clusters <- base::data.frame(type = base::rep("global", base::nrow(part_2_res)), stringsAsFactors = FALSE)
    if(global_vgene == TRUE){
      global_clusters$tag <- base::paste(part_2_res$cluster_tag, part_2_res$TRBV, part_2_res$aa_at_position, sep = "_")
    } else {
      global_clusters$tag <- base::paste(part_2_res$cluster_tag, part_2_res$aa_at_position, sep = "_")
    }
    
    global_clusters$cluster_size <- part_2_res$cluster_size
    global_clusters$unique_cdr3_sample <- part_2_res$unique_CDR3b
    global_clusters$unique_cdr3_ref <- part_2_res$num_in_ref
    global_clusters$OvE <- base::rep(0, base::nrow(part_2_res))
    global_clusters$fisher.score <- part_2_res$fisher.score
    global_clusters$members <- part_2_res$CDR3b
    
    if(local_similarities == TRUE && base::nrow(part_1_res) > 0) merged_clusters <- base::rbind(merged_clusters, global_clusters) else merged_clusters <- global_clusters
  }
  
  
  cluster_list <- base::list()
  if(!base::is.null(merged_clusters)){
    ### Build a list with elements representing clusters and containing a data frame with all information about the cluster members
    cluster_list <- foreach::foreach(i = base::seq_len(base::nrow(merged_clusters))) %dopar% {
      if(merged_clusters$type[i] == "local"){
        act_seqs <- base::unlist(base::strsplit(merged_clusters$members[i], split = " "))
        base::return(sequences[sequences$CDR3b %in% act_seqs,])
      }
      
      if(merged_clusters$type[i] == "global"){
        act_details <- base::unlist(base::strsplit(merged_clusters$tag[i], split = "_"))
        
        act_seqs <- base::unlist(base::strsplit(merged_clusters$members[i], split = " "))
        
        if(global_vgene == FALSE){
          base::return(base::data.frame(sequences[sequences$CDR3b %in% act_seqs,], stringsAsFactors = FALSE))
        } else {
          base::return(base::data.frame(sequences[sequences$CDR3b %in% act_seqs & sequences$TRBV == act_details[2],], stringsAsFactors = FALSE))
        }
      }
      # })
    }
    base::names(cluster_list) <- merged_clusters$tag
    
    ### Update some cluster information for local clusters
    if(local_similarities == TRUE && base::nrow(part_1_res) > 0 && !base::is.null(merged_clusters)){
      for(i in base::which(merged_clusters$type == "local")){
        
        # cluster size
        merged_clusters$cluster_size[i] <- base::length(base::unique(cluster_list[[i]]$CDR3b))
        
        if(boost_local_significance == TRUE){
          # check if motif positions contain N-nucleotide encoded amino acids
          start_motif <- part_1_res$start[i]
          stop_motif <- part_1_res$stop[i]+base::nchar(part_1_res$motif[i])-1
          act_ids <- cluster_list[[i]]$seq_ID
          motif_in_range <- (range_df$start_1[act_ids] <= stop_motif & range_df$stop_1[act_ids] >= start_motif) | (range_df$start_2[act_ids] <= stop_motif & range_df$stop_2[act_ids] >= start_motif)
        
          # divide Fisher's p-value by 2^n (n = number of unique CDR3b sequences with N-nucleotide encoded amino acids in the motif)
          merged_clusters$fisher.score[i] <- merged_clusters$fisher.score[i]/(2^base::sum(motif_in_range))
        }
      }
    }
    
    ### exclude clusters with less size than the minimum size
    eliminate_ids <- base::which(merged_clusters$cluster_size < cluster_min_size)
    if(base::length(eliminate_ids) > 0){
      merged_clusters <- merged_clusters[-eliminate_ids,]
      cluster_list <- cluster_list[-eliminate_ids]
    }
    if(base::nrow(merged_clusters) == 0){
      merged_clusters <- base::c()
      cluster_list <- base::list()
    }
  }
  
  clone_network <- base::c()
  if(!base::is.null(merged_clusters)){
    
    ### Identify/save all connections between cluster members within a cluster
    
    # first for local similarities
    if(local_similarities == TRUE){
      local_clone_network <- foreach::foreach(i = base::which(merged_clusters$type == "local")) %dopar% {
        
        # get all members of the cluster
        temp_members <- cluster_list[[i]]$CDR3b
        if(structboundaries) temp_members_frags <- base::substr(x = temp_members,start = boundary_size + 1 ,stop = base::nchar(temp_members) - boundary_size) else temp_members_frags <- temp_members
        
        # get the position of the motif in the sequences
        temp_pos <- stringr::str_locate(string = temp_members_frags, pattern = base::strsplit(base::names(cluster_list)[i], split = "_")[[1]][[1]])
        
        # first assume connection between all members
        if(base::length(temp_members) >= 2){
          combn_ids <- base::t(utils::combn(base::seq_along(temp_members), m = 2))
        } else {
          combn_ids <- base::t(utils::combn(base::rep(1, 2), m = 2))
        }
        temp_df <- base::data.frame(V1 = temp_members[combn_ids[,1]], V2 = temp_members[combn_ids[,2]],
                                    type = base::rep("local", base::nrow(combn_ids)),
                                    cluster_tag = base::rep(base::names(cluster_list)[i], base::nrow(combn_ids)))
        
        # Only allow connections between neighbored motif positions (determined by the motif_distance_cutoff) 
        temp_df <- base::unique(temp_df[base::abs(temp_pos[combn_ids[,1],1]-temp_pos[combn_ids[,2],1]) < motif_distance_cutoff ,])
        
        temp_df <- base::t(temp_df)
        base::return(base::unlist(temp_df))
      }
      
      # convert the results into a data frame
      local_clone_network <- base::data.frame(base::matrix(base::unlist(local_clone_network), ncol = 4, byrow = TRUE),
                                              stringsAsFactors = FALSE)
      base::colnames(local_clone_network) <- base::c("V1", "V2", "type", "cluster_tag")
      
      clone_network <- local_clone_network
    }
    
    # second for global similarities
    if(global_similarities == TRUE){
      global_clone_network <- foreach::foreach(i = base::which(merged_clusters$type == "global")) %dopar% {
        
        # get all members of the cluster
        temp_members <- cluster_list[[i]]$CDR3b
        
        # first assume connection between all members
        if(base::length(temp_members) >= 2){
          combn_ids <- base::t(utils::combn(base::seq_along(temp_members), m = 2))
        } else {
          combn_ids <- base::t(utils::combn(base::rep(1, 2), m = 2))
        }
        temp_df <- base::data.frame(V1 = temp_members[combn_ids[,1]], V2 = temp_members[combn_ids[,2]],
                                    type = base::rep("global", base::nrow(combn_ids)),
                                    cluster_tag = base::rep(base::names(cluster_list)[i], base::nrow(combn_ids)))
        
        # Filter (if requested) for connections with replaceable amino acids allowed by BLOSUM62 (value >= 0), filtering for V gene identity already done
        if(all_aa_interchangeable == FALSE){
          temp_pos <- stringr::str_locate(string = base::names(cluster_list)[i], pattern = "%")
          if(structboundaries == TRUE) temp_pos <- temp_pos+boundary_size
          temp_df <- base::unique(temp_df[base::paste0(base::substr(temp_df$V1, temp_pos, temp_pos),
                                                       base::substr(temp_df$V2, temp_pos, temp_pos)) %in% BlosumVec,])
        }
        
        temp_df <- base::t(temp_df)
        base::return(base::unlist(temp_df))
      }
      
      # convert the results into a data frame
      global_clone_network <- base::data.frame(base::matrix(base::unlist(global_clone_network), ncol = 4, byrow = TRUE),
                                               stringsAsFactors = FALSE)
      base::colnames(global_clone_network) <- base::c("V1", "V2", "type", "cluster_tag")
      
      if(local_similarities == FALSE){
        clone_network <- global_clone_network
      } else {
        clone_network <- base::rbind(clone_network, global_clone_network)
      }
    }
    clone_network[] <- base::lapply(clone_network, base::as.character)
    
    # all singletons (sequences without any connetion) to the connection data frame
    not_in_network <- sequences$CDR3b[!(sequences$CDR3b %in% base::c(clone_network$V1, clone_network$V2))]
    if(base::length(not_in_network) > 0){
      clone_network <- base::rbind(clone_network,
                                   base::data.frame(V1 = not_in_network, V2 = not_in_network,
                                                    type = base::rep("singleton", base::length(not_in_network)),
                                                    cluster_tag = base::paste0("singleton_", base::seq_along(not_in_network)),
                                                    stringsAsFactors = FALSE))
    }
  }
  
  ### save
  if(save_results == TRUE) utils::write.table(x = clone_network,file = base::paste0(result_folder, "clone_network.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
  save_cluster_list_df <- base::c()
  if(!base::is.null(merged_clusters)){
  
    ### too many numbers are confusing
    #merged_clusters$fisher.score <- base::round(merged_clusters$fisher.score,
    #                                            digits = base::ceiling(base::abs(base::log10(merged_clusters$fisher.score)))+1)
    merged_clusters$fisher.score <- base::formatC(merged_clusters$fisher.score, digits = 1, format = "e")
    
    merged_clusters$OvE[base::is.infinite(merged_clusters$OvE)] <- 0
    
    ### save cluster_list (must be done as a data frame...)
    save_cluster_list_df <- foreach::foreach(i = base::seq_along(cluster_list), .combine = "rbind") %dopar% {
      temp <- cluster_list[[i]]
      temp <- base::cbind(base::data.frame(tag = base::rep(base::names(cluster_list)[i], base::nrow(temp))), temp)
      return(temp)
    }
  }
  
  if(save_results == TRUE) utils::write.table(x = save_cluster_list_df,file = base::paste0(result_folder, "cluster_member_details.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  ### stop parallelization
  doParallel::stopImplicitCluster()
  
  t4_e <- base::Sys.time()
  base::cat("Part 3 cpu time:",t4_e-t3_e,base::units(t4_e-t3_e),"\n\n")
  
  #################################################################
  ##                       Part 4: Scoring                       ##
  #################################################################
  
  base::cat("Part 4: Scoring convergence groups\n")
  
  if(!base::is.null(merged_clusters)){
    scoring_res <- cluster_scoring(cluster_list = cluster_list,
                                   cdr3_sequences = cdr3_sequences,
                                   refdb_beta = refdb_beta,
                                   v_usage_freq = v_usage_freq,
                                   cdr3_length_freq = cdr3_length_freq,
                                   ref_cluster_size = ref_cluster_size,
                                   sim_depth = sim_depth,
                                   gliph_version = 2,
                                   hla_cutoff = hla_cutoff,
                                   n_cores = n_cores)
    
    merged_clusters <- base::cbind(merged_clusters, scoring_res)
    # some reordering of columns
    merged_clusters <- merged_clusters[, base::c(base::colnames(merged_clusters)[base::colnames(merged_clusters) != "members"], "members")]
  }
  
  if(save_results == TRUE) utils::write.table(x = merged_clusters,file = base::paste0(result_folder, "convergence_groups.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  base::cat("\n")
  t5_e <- base::Sys.time()
  base::cat("Part 4 cpu time:",t5_e-t4_e, base::units(t5_e-t4_e),"\n\n\n")
  
  #################################################################
  ##                        Create output                        ##
  #################################################################
  
  # set all theoretical numeric values (actually characters) to numeric values
  if(base::is.data.frame(part_1_res)){
    for(i in base::seq_len(base::ncol(part_1_res))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(part_1_res[,i])))) == FALSE) part_1_res[,i] <- base::as.numeric(part_1_res[,i])
    }
  }
  if(base::is.data.frame(part_1_res_all_clusters)){
    for(i in base::seq_len(base::ncol(part_1_res_all_clusters))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(part_1_res_all_clusters[,i])))) == FALSE) part_1_res_all_clusters[,i] <- base::as.numeric(part_1_res_all_clusters[,i])
    }
  }
  if(base::is.data.frame(part_2_res)){
    for(i in base::seq_len(base::ncol(part_2_res))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(part_2_res[,i])))) == FALSE) part_2_res[,i] <- base::as.numeric(part_2_res[,i])
    }
  }
  if(base::is.data.frame(merged_clusters)){
    for(i in base::seq_len(base::ncol(merged_clusters))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(merged_clusters[,i])))) == FALSE) merged_clusters[,i] <- base::as.numeric(merged_clusters[,i])
    }
  }
  if(base::is.list(cluster_list)){
    for(i in base::seq_len(base::length(cluster_list))){
      if(base::is.data.frame(cluster_list[[i]])){
        for(j in base::seq_len(base::ncol(cluster_list[[i]]))){
          if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(cluster_list[[i]][,j])))) == FALSE) cluster_list[[i]][,j] <- base::as.numeric(cluster_list[[i]][,j])
        }
      }
      
    }
  }
  
  ### the output list
  output <- base::list(motif_enrichment = NULL,
                       global_enrichment = NULL,
                       connections = NULL,
                       cluster_properties = NULL,
                       cluster_list = NULL,
                       parameters = NULL)
  
  output$motif_enrichment <- base::list(selected_motifs = part_1_res,all_motifs = part_1_res_all_clusters)
  output$global_enrichment <- part_2_res
  output$connections <- clone_network
  output$cluster_properties <- merged_clusters
  output$cluster_list <- cluster_list
  
  output$parameters <- base::list(gliph_version = 2,
                                  result_folder = result_folder,
                                  ref_cluster_size = ref_cluster_size,
                                  sim_depth = sim_depth,
                                  lcminp = lcminp,
                                  lcminove = lcminove,
                                  motif_distance_cutoff = motif_distance_cutoff,
                                  kmer_mindepth = kmer_mindepth,
                                  accept_sequences_with_C_F_start_end = accept_sequences_with_C_F_start_end,
                                  min_seq_length = min_seq_length,
                                  structboundaries = structboundaries,
                                  boundary_size = boundary_size,
                                  motif_length = motif_length,
                                  discontinuous_motifs = discontinuous_motifs,
                                  local_similarities = local_similarities,
                                  global_similarities = global_similarities,
                                  global_vgene = global_vgene,
                                  all_aa_interchangeable = all_aa_interchangeable,
                                  cluster_min_size = cluster_min_size,
                                  hla_cutoff = hla_cutoff,
                                  n_cores = n_cores)
  
  ### save input parameters
  paras <- base::c()
  for(i in base::seq_along(output$parameters)){
    temp_paras <- base::data.frame(base::names(output$parameters)[i], base::paste0(output$parameters[[i]], collapse = ","), stringsAsFactors = FALSE)
    if(i == 1) paras <- temp_paras else paras <- base::rbind(paras, temp_paras)
  }
  if(save_results == TRUE) utils::write.table(x = paras,file = base::paste0(result_folder, "parameter.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
  t2 <- base::Sys.time()
  dt <- (t2-t1)
  base::cat("Total time = ",dt, base::units(dt),"\n")
  
  if(save_results == TRUE) base::cat("Output: results are stored in ", result_folder, "\n")
  
  ### Closing time!
  base::return(output)
}
