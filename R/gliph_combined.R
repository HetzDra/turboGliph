#'  Grouping of Lymphocyte Interactions by Paratope Hotspots Customizable
#'
#' @description Identifying specificity groups in the T cell receptor
#'  repertoire. Detailed customizable implementation of GLIPH/GLIPH2.
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
#' @param min_seq_length numeric. By default 8. All the sequences with a length less than this
#' value will be filtered out in input and reference database. If structboundaries
#' is \code{TRUE}, it is recommended not to go below the default. In this case,
#' the min_seq_length will be set to the maximum of 2*boundary_size+2 and
#' min_seq_length.
#' @param accept_sequences_with_C_F_start_end logical. This logical flag
#' if \code{TRUE}, by default, only accepts sequences with
#' amino-acid C at the start position and amino-acid F at the
#' end position. This flag should be set to \code{FALSE} if
#' you wish to analyze sequences of different origin for example B-cells.
#' @param structboundaries logical. By default \code{TRUE}. By setting this flag to \code{TRUE}
#' the first boundary_size and the last boundary_size amino acids of each sequence will
#' not be considered in the analysis for computing the Hamming
#' distance and motif enrichment in input and reference database.
#' @param boundary_size numeric. By default 3. Specifies the boundary size if structboundaries is
#' active.
#' @param local_similarities logical. By default \code{TRUE}. Determines whether the sequences should be analyzed for local similarity.
#' @param local_method character. Either 'fisher' or 'rrs' (default). Determines the method for searching local similarities.
#' In the case of 'rrs', repeated random sampling is performed as in the GLIPH algorithm.
#' In the case of 'fisher', the long runtime of repeated. random sampling is shortened by approximation using Fisher's Exact Test as established in GLIPH2.
#' @param lcsim_depth numeric. By default 1000. Number of iterations for repeated random sampling, if local_method is set to 'rrs'.
#' @param motif_length accepts a numeric vector of motif lengths
#' you want GLIPH2 to find and study. By default it searches for
#' motifs of size 2, 3 and 4 amino acids.
#' @param motif_distance_cutoff numeric. By default 3. Defines the number of positions between which motifs for a local connection are allowed to vary.
#' @param lcminp numeric. By default 0.01. Local convergence maximum probability score cutoff.
#' The score reports the probability that a random sample of the
#' same size as the sample set but of the reference set
#' (i.e. naive repertoire) would generate an enrichment of the
#' given motif at least as high as has been observed in the
#' sample set.
#' @param lcminove numeric. By default 10. Local convergence minimum observed vs expected fold change.
#' This is a cutoff for the minimum fold enrichment over a
#' reference distribution that a given motif should have in
#' the sample set in order to be considered for further evaluation. By default, the minimum fold enrichment (1000,100,10) is
#' dependent on the motif length (2,3,4 amino acids). \code{lcminove} has to be either a single numeric value or a numeric
#' vector with equal length as \code{motif_length} representing the minimum fold enrichment depending on the respective motif_length.
#' @param lckmer_mindepth numeric. By default 3. Minimum observations of kmer for it to be evaluated. This is
#' the minimum number of times a kmer should be observed in
#' the sample set in order for it to be considered for further
#' evaluation. The number can be set higher to provide less
#' motif-based clusters with higher confidence. This could be
#' recommended if the sample set is greater than 5000 reads. Lowering
#' the value to 2 will identify more groups but likely at a cost
#' of an increased False Discovery Rate.
#' @param discontinuous_motifs logical. By default \code{FALSE}. Determines whether discontinuous motifs are to be considered.
#' @param boost_local_significance logical. By default \code{FALSE}. If set to \code{TRUE}, fisher scores of local clusters are repeatedly divided 
#' by 2 for every unique CDR3 sequence in the cluster in which the motif overlaps with non-germline encoded N- or P-nucleotides.
#' @param cdr3_len_stratify logical. By default \code{FALSE}.
#' Specifies whether the distribution of the cdr3 lengths in the sample should be retained during repeat random sampling.
#' @param vgene_stratify logical. By default \code{FALSE}.
#' Specifies whether the distribution of V-genes in the sample should be retained during repeat random sampling.
#' Requires V-gene information in \code{cdr3_sequences}.
#' @param global_similarities logical. By default \code{TRUE}. Determines whether the sequences should be analyzed for global similarity.
#' @param global_method character. Either 'fisher' or 'cutoff' (default). Determines the method for searching global similarities.
#' In the case of 'cutoff', global similarity is defined as falling below a cutoff (gccutoff) of the Hamming distance between two sequences,
#' as in the GLIPH algorithm. THis option requires clustering_method to be set to 'GLIPH1.0'.
#' In the case of 'fisher', as in GLIPH2, Fisher's Exact Test is used to test for a significant enrichment of global
#' structures in the sample set relative to the reference set.
#' @param gccutoff numeric. Global convergence distance cutoff. Only considered, if global_method is set to 'cutoff'. This is the maximum
#' CDR3 Hamming mutation distance between two clones sharing
#' the same CDR3 length in order for them
#' to be considered to be likely binding the same antigen.
#'
#' This number will change depending on sample depth, as
#' with more reads, the odds of finding a similar sequence
#' increases even in a naive repertoire.
#'
#' This number will also change depending on the species
#' evaluated and even the choice of reference database (memory
#' TCRs will be more likely to have similar TCRs than naive
#' TCR repertoires). Thus, by default this is calculated at
#' runtime if not specified. If the sample depth is less
#' than 125 it will be set to 2, otherwise it will be set to 1.
#' @param gcminp numeric. By default 1. Global convergence maximum probability score cutoff in case of global_method = 'fisher'.
#' The score reports the probability that a random sample of the
#' same size as the sample set but of the reference set
#' (i.e. naive repertoire) would generate an enrichment of the
#' given global structure at least as high as has been observed in the
#' sample set.
#' @param gcminove numeric. By default 0. Global convergence minimum observed vs expected fold change in case of global_method = 'fisher'.
#' This is a cutoff for the minimum fold enrichment over a
#' reference distribution that a given global enrichment should have in
#' the sample set in order to be considered for further evaluation.
#' @param gckmer_mindepth numeric. By default 2. Minimum observations of a global structure for it to be evaluated in case of global_method = 'fisher'. This is
#' the minimum number of times a global structure should be observed in
#' the sample set in order for it to be considered for further
#' evaluation. The number can be set higher to provide less
#' motif-based clusters with higher confidence. Lowering
#' the value will identify more groups but likely at a cost
#' of an increased False Discovery Rate.
#' @param all_aa_interchangeable logical. Only used if global_method = 'fisher'. By default \code{FALSE}.In the case of TRUE, all sequences with a Hamming distance of 1 are evaluated
#' as global similarities. In the case of FALSE, all sequences with a Hamming distance of 1 whose different amino acid has a BLOSUM62 score >= 0
#' are evaluated as global similarities.
#' @param clustering_method character. Either 'GLIPH1.0' (default) or 'GLIPH2.0'. Determines the method for clustering.
#' In the case of 'GLIPH1.0' clusters are generated by all locally and globally connected sequences.
#' In the case of 'GLIPH2.0' clusters only contain sequences with one specific enriched local motif or one specific enriched global structure.
#' @param vgene_match character. Specifies which connections are restricted by shared V genes. Can be set to the following values:
#' \itemize{
#' \item "none": Default. Connections between sequences are not restricted to a shared V gene.
#' \item "local": Local connections between sequences are restricted to a shared V gene.
#' \item "global": Global connections between sequences are restricted to a shared V gene.
#' \item "all": All connections between sequences are restricted to a shared V gene.
#' }
#' @param public_tcrs character. Specifies which connections are restricted by isolation from the same donor. Can be set to the following values:
#' \itemize{
#' \item "none": Connections between sequences are unrestricted from the donor.
#' \item "local": Local connections between sequences are restricted from being obtained from the same donor.
#' \item "global": Global connections between sequences are restricted from being obtained from the same donor.
#' \item "all": Default. All connections between sequences are restricted from being obtained from the same donor.
#' }
#' @param cluster_min_size numeric. By default 2. Minimal size of a cluster required to be considered for scoring.
#' @param scoring_method character. Either 'GLIPH1.0' (default) or 'GLIPH2.0'. Determines from which GLIPH version the scoring algorithm should
#' be used for cluster scoring. The differences are mainly in the final multiplication of all subscores.
#' @param scoring_sim_depth numeric. By default 1000. Simulated resampling depth for assessing V gene and CDR3 length enrichment scores of clusters.
#' @param hla_cutoff numeric. By default 0.1. Defines the threshold of HLA probability scores below which HLA alleles are considered significant.
#' @param n_cores numeric. Number of cores to use, by default 1. In case of \code{NULL} it will be set to number of cores in your machine minus 2.
#'
#' @return This function returns a list of seven elements whose contents are explained below. If a file path is specified under \code{result_folder},
#' the results are additionally stored there. The individual file names are also specified below (italic name parts indicate the given value of the
#' corresponding parameter).
#'
#' @return $sample_log:
#' Only generated if \code{local_method} = 'rrs'. Contains a data frame with \code{1 + lcsim_depth} rows representing observations for all the possible
#' k-mer motifs. The first observation, \code{Discovery} is the actual
#' observation counts in the input sample and the rest shows the observation
#' counts in a subsample from reference database.\cr
#' File name: kmer_resample_\code{lcsim_depth}_log.txt
#'
#' @return $motif_enrichment:
#' A list of two data frames. \code{selected_motifs} contains only the motifs that pass the filtering criterion (ove and p-value),
#' whereas \code{all_motifs} contains p-value and ove of all motifs.\cr
#' File name of \code{selected_motifs}: local_similarities_minp_\code{lcminp}_ove \code{lcminove}_kmer_mindepth \code{lckmer_mindepth}.txt
#' File name of \code{all_motifs}: all_motifs.txt
#'
#' @return $global_enrichment:
#' Only generated if \code{global_method} = 'fisher'. Contains a list of two data frames. \code{selected_structs} contains only the sequence structures that pass the filtering criterion (p-value),
#' whereas \code{all_structs} contains p-value and ove of all sequence structures\cr
#' File name of \code{selected_structs}: global_similarities_minp_\code{gcminp}_ove \code{gcminove}_kmer_mindepth \code{gckmer_mindepth}.txt
#' File name of \code{all_structs}: all_global_similarities.txt
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
#' \item "p.value": The p-value of the motif or global structure, if clustering is performed as in GLIPH2.
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
#' res <- gliph_combined(cdr3_sequences = gliph_input_data[base::seq_len(200),],
#' lcsim_depth = 50,
#' scoring_sim_depth = 50,
#' n_cores = 1)
#'
#' @references Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @references \url{https://github.com/immunoengineer/gliph}
#' @references Huang, Huang, et al.
#' "Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening." Nature Biotechnology 38.10 (2020): 1194-1202.
#' @references \url{http://50.255.35.37:8080/}
#' @import foreach
#' @export
gliph_combined <- function(cdr3_sequences,
                   result_folder = "",
                   # refdb_beta = "human_v1.0_CD4",
                   refdb_beta = "gliph_reference",
                   v_usage_freq = NULL,
                   cdr3_length_freq = NULL,
                   ref_cluster_size = "original",
                   
                   min_seq_length = 8,
                   accept_sequences_with_C_F_start_end = TRUE,
                   structboundaries = TRUE,
                   boundary_size = 3,
                   
                   local_similarities = TRUE, #Draius13042020
                   local_method = "rrs", # or "fisher"
                   lcsim_depth = 1000, # if rrs
                   motif_length = base::c(2,3,4),
                   motif_distance_cutoff = 3,
                   lcminp = 0.01,
                   lcminove = c(1000,100,10),
                   lckmer_mindepth = 3,
                   discontinuous_motifs = FALSE,
                   boost_local_significance = FALSE,
                   cdr3_len_stratify = FALSE, #Draius17042020
                   vgene_stratify = FALSE, #Draius17042020
                   
                   global_similarities = TRUE, #Draius13042020,
                   global_method = "cutoff", # or "fisher"
                   gccutoff = NULL,
                   gcminp = 1,
                   gcminove = 0,
                   gckmer_mindepth = 2,
                   all_aa_interchangeable = FALSE,
                   
                   clustering_method = "GLIPH1.0", # or "GLIPH2.0"
                   vgene_match = "none",
                   public_tcrs = "all",
                   cluster_min_size = 2,
                   
                   scoring_method = "GLIPH1.0", # or "GLIPH2.0"
                   scoring_sim_depth = 1000,
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
    
    if(base::file.exists(base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_log.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_log.txt"),"\n"))
    }
    
    if(base::file.exists(base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", base::paste(lcminove, collapse = "_"), "_kmer_mindepth_", lckmer_mindepth, ".txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", base::paste(lcminove, collapse = "_"), "_kmer_mindepth_", lckmer_mindepth, ".txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"all_motifs.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"all_motifs.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder, "global_similarities_minp_",gcminp, "_minove_", base::paste(gcminove, collapse = "_"), "_kmer_mindepth_", gckmer_mindepth, ".txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "global_similarities_minp_",gcminp, "_minove_", base::paste(gcminove, collapse = "_"), "_kmer_mindepth_", gckmer_mindepth, ".txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder,"all_global_similarities.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder,"all_global_similarities.txt"),"\n"))
    }
    if(base::file.exists(base::paste0(result_folder, "cluster_member_details.txt"))){
      save_results <- FALSE
      base::cat(base::paste("WARNING: Saving the results is cancelled. The following file already exists in the result folder:\n",
                            base::paste0(result_folder, "cluster_member_details.txt"),"\n"))
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
  
  ### local_similarities
  if(!base::is.logical(local_similarities))base::stop("local_similarities has to be logical")
  
  ### local_method
  if(!(local_method %in% base::c("fisher", "rrs")))base::stop("local_method has to be either 'fisher' or 'rrs'.")
  
  ### lcsim_depth
  if(!base::is.numeric(lcsim_depth))base::stop("lcsim_depth has to be numeric")
  if(base::length(lcsim_depth) > 1)base::stop("lcsim_depth has to be a single number")
  if(lcsim_depth < 1)base::stop("lcsim_depth must be at least 1")
  lcsim_depth <- base::round(lcsim_depth)
  
  ### motif_length
  if(!base::is.numeric(motif_length))base::stop("motif_length has to be numeric")
  if(base::any(motif_length < 1))base::stop("values of motif_length must be at least 1")
  motif_length <- base::round(motif_length)
  
  ### motif_distance_cutoff
  if(!base::is.numeric(motif_distance_cutoff)) base::stop("motif_distance_cutoff has to be numeric")
  if(base::length(motif_distance_cutoff) > 1) base::stop("motif_distance_cutoff has to be a single number")
  motif_diffs <- max(0, motif_distance_cutoff-1)
  
  ### lcminp
  if(!base::is.numeric(lcminp))base::stop("lcminp has to be numeric")
  if(base::length(lcminp) > 1)base::stop("lcminp has to be a single number")
  if(lcminp <= 0)base::stop("lcminp must be greater than 0")
  
  ### lcminove
  if(!base::is.numeric(lcminove)) base::stop("lcminove has to be numeric")
  if(base::length(lcminove) > 1 && base::length(lcminove) != base::length(motif_length)) base::stop("lcminove has to be a single number or of same length as motif_length")
  if(base::any(lcminove < 1)) base::stop("lcminove must be at least 1")
  
  ### lckmer_mindepth
  if(!base::is.numeric(lckmer_mindepth))base::stop("lckmer_mindepth has to be numeric")
  if(base::length(lckmer_mindepth) > 1)base::stop("lckmer_mindepth has to be a single number")
  if(lckmer_mindepth < 1)base::stop("lckmer_mindepth must be at least 1")
  lckmer_mindepth <- base::round(lckmer_mindepth)
  
  ### discontinuous_motifs
  if(!base::is.logical(discontinuous_motifs))base::stop("discontinuous_motifs has to be logical")
  
  ### boost_local_significance
  if(!base::is.logical(boost_local_significance))base::stop("boost_local_significance has to be logical")
  
  ### cdr3_len_stratify
  if(!base::is.logical(cdr3_len_stratify)) base::stop("cdr3_len_stratify has to be logical")
  
  ### vgene_stratify
  if(!base::is.logical(vgene_stratify)) base::stop("vgene_stratify has to be logical")
  
  ### global_similarities
  if(!base::is.logical(global_similarities))base::stop("global_similarities has to be logical")
  
  ### global_method
  if(!(global_method %in% base::c("fisher", "cutoff")))base::stop("global_method has to be either 'fisher' or 'cutoff'.")
  
  ### gccutoff
  if(!base::is.null(gccutoff) && !base::is.numeric(gccutoff)) base::stop("gccutoff has to be NULL or numeric")
  if(!base::is.null(gccutoff) && base::length(gccutoff)>1) base::stop("gccutoff has to be NULL or a single number")
  if(!base::is.null(gccutoff) && gccutoff < 0) base::stop("gccutoff must be at least 0")
  
  ### gcminp
  if(!base::is.numeric(gcminp))base::stop("gcminp has to be numeric")
  if(base::length(gcminp) > 1)base::stop("gcminp has to be a single number")
  if(gcminp <= 0)base::stop("gcminp must be greater than 0")
  
  ### gcminove
  if(!base::is.numeric(gcminove)) base::stop("gcminove has to be numeric")
  if(base::length(gcminove) > 1) base::stop("gcminove has to be a single number")
  if(base::any(gcminove < 0)) base::stop("gcminove must be at least 0")
  
  ### gckmer_mindepth
  if(!base::is.numeric(gckmer_mindepth))base::stop("gckmer_mindepth has to be numeric")
  if(base::length(gckmer_mindepth) > 1)base::stop("gckmer_mindepth has to be a single number")
  if(gckmer_mindepth < 1)base::stop("gckmer_mindepth must be at least 1")
  gckmer_mindepth <- base::round(gckmer_mindepth)
  
  ### all_aa_interchangeable
  if(!base::is.logical(all_aa_interchangeable))base::stop("all_aa_interchangeable has to be logical")
  
  ### clustering_method
  if(!(clustering_method %in% base::c("GLIPH1.0", "GLIPH2.0")))base::stop("clustering_method has to be either 'GLIPH1.0' or 'GLIPH2.0'.")
  if(global_method == "cutoff") clustering_method <- "GLIPH1.0"
  
  ### vgene_match
  if(!(vgene_match %in% base::c("none", "local", "global", "all")))base::stop("vgene_match has to be either 'none', 'local', 'global' or 'all'.")
  
  ### public_tcrs
  if(!(public_tcrs %in% base::c("none", "local", "global", "all")))base::stop("public_tcrs has to be either 'none', 'local', 'global' or 'all'.")
  
  ### cluster_min_size
  if(!base::is.numeric(cluster_min_size))base::stop("cluster_min_size has to be numeric")
  if(base::length(cluster_min_size) > 1)base::stop("cluster_min_size has to be a single number")
  if(cluster_min_size < 1)base::stop("cluster_min_size must be at least 1")
  cluster_min_size <- base::round(cluster_min_size)
  
  ### scoring_method
  if(!(scoring_method %in% base::c("GLIPH1.0", "GLIPH2.0")))base::stop("scoring_method has to be either 'GLIPH1.0' or 'GLIPH2.0'.")
  
  ### scoring_sim_depth
  if(!base::is.numeric(scoring_sim_depth))base::stop("scoring_sim_depth has to be numeric")
  if(base::length(scoring_sim_depth) > 1)base::stop("scoring_sim_depth has to be a single number")
  if(scoring_sim_depth < 1)base::stop("scoring_sim_depth must be at least 1")
  scoring_sim_depth <- base::round(scoring_sim_depth)
  
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
    if(vgene_match != "none"){
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
    } else if(vgene_match != "none"){
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
  if(structboundaries) all_motif_region <- base::substr(x = sequences$CDR3b,start = boundary_size + 1 ,stop = base::nchar(sequences$CDR3b) - boundary_size) else all_motif_region <- sequences$CDR3b
  
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
    if(base::ncol(refseqs) > 1 && vgene_match != "none"){
      if("TRBV" %in% base::colnames(refseqs)){
        base::cat("Notification: Column of reference database named 'TRBV' is considered as V gene information.\n")
      } else {
        base::cat("Notification: Second column of reference database is considered as V gene information.\n")
        base::colnames(refseqs)[2] <- "TRBV"
      }
    } else if(vgene_match != "none"){
      base::stop("V-gene restriction for global similarities ('vgene_match' != 'none') requires V-gene information in second column of reference database.\n")
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
  SAMPLE_LOG <- NULL
  CLONE_NETWORK <- NULL
  SELECTED_MOTIFS <- NULL
  ALL_MOTIFS <- NULL
  
  if(local_similarities == TRUE){
    base::cat(base::paste("Part 1: Searching for local similarities.\n"))
    
    ### Method repeated random sampling
    if(local_method == "rrs"){
      
      ### Prepare biased repeated random sampling
      # get v gene distribution for biased repeat random sampling
      motif_region_vgenes_list <- base::list()
      ref_motif_vgenes_id_list <- base::list()
      if(vgene_stratify == TRUE){
        for(act_vgene in base::sort(base::as.character(base::unique(sequences$TRBV)))){
          motif_region_vgenes_ids[[act_vgene]] <- base::sum(sequences$TRBV == act_vgene)
          ref_motif_vgenes_id_list[[act_vgene]] <- base::which(base::unique(refseqs$CDR3b) %in% refseqs$CDR3b[refseqs$TRBV == act_vgene])
        }
      }
      # get CDR3b length distribution for biased repeat random sampling
      motif_lengths_list <- base::list()
      ref_motif_lengths_id_list <- base::list()
      if(cdr3_len_stratify == TRUE){
        motif_region_lengths <- base::nchar(motif_region)
        ref_motif_region_lengths <- base::nchar(refseqs_motif_region)
        
        for(cdr3_length in base::sort(base::unique(motif_region_lengths))){
          motif_lengths_list[[base::as.character(cdr3_length)]] <- base::sum(motif_region_lengths == cdr3_length)
          ref_motif_lengths_id_list[[base::as.character(cdr3_length)]] <- base::which(ref_motif_region_lengths == cdr3_length)
        }
      }
      # get correlated v gene and CDR3b length distribution for biased repeat random sampling
      lengths_vgenes_list <- base::list()
      ref_lengths_vgenes_list <- base::list()
      if((vgene_stratify == TRUE) && (cdr3_len_stratify == TRUE)){
        for(cdr3_length in base::names(ref_motif_lengths_id_list)){
          for(act_vgene in base::names(ref_motif_vgenes_id_list)){
            lengths_vgenes_list[[cdr3_length]][[act_vgene]] <- base::sum(sequences$TRBV == act_vgene & motif_region_lengths == base::as.numeric(cdr3_length))
            ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]] <- ref_motif_lengths_id_list[[cdr3_length]][ref_motif_lengths_id_list[[cdr3_length]] %in% ref_motif_vgenes_id_list[[act_vgene]]]
          }
        }
      }
      
      ### Compute motif frequency in sample set
      discovery <- find_motifs(seqs = motif_region,q = motif_length, discontinuous = discontinuous_motifs)
      
      ### Compute motif frequency in reference subsamples at sample set depth
      res <- foreach::foreach(i = base::seq_len(lcsim_depth), .inorder = FALSE) %dopar% {
        motif_sample <- getRandomSubsample(cdr3_len_stratify = cdr3_len_stratify,
                                                       vgene_stratify = vgene_stratify,
                                                       refseqs_motif_region = refseqs_motif_region,
                                                       motif_region = motif_region,
                                                       motif_lengths_list = motif_lengths_list,
                                                       ref_motif_lengths_id_list = ref_motif_lengths_id_list,
                                                       motif_region_vgenes_list = motif_region_vgenes_list,
                                                       ref_motif_vgenes_id_list = ref_motif_vgenes_id_list,
                                                       ref_lengths_vgenes_list = ref_lengths_vgenes_list,
                                                       lengths_vgenes_list = lengths_vgenes_list)
        
        sim <- find_motifs(seqs = motif_sample,q = motif_length, discontinuous = discontinuous_motifs)
        sim <- dplyr::left_join(discovery,sim,"motif")
        sim$V1.y[base::is.na(sim$V1.y)] <- 0
        simy <- sim$V1.y
        simy
      }
      
      ### Adjust SAMPLE_LOG results
      # convert results into a data frame with iteration represented by the row and motif represented by the column
      res <- base::unlist(res)
      copies <- base::matrix(base::as.numeric(res), ncol = lcsim_depth, byrow = FALSE)
      SAMPLE_LOG <- base::cbind(discovery$motif,discovery$V1,copies)
      SAMPLE_LOG <- base::t(SAMPLE_LOG)
      mot <- SAMPLE_LOG[1,]
      SAMPLE_LOG <- SAMPLE_LOG[-1,]
      SAMPLE_LOG <- base::matrix(base::as.numeric(SAMPLE_LOG),nrow = lcsim_depth+1)
      base::colnames(SAMPLE_LOG) <- mot
      SAMPLE_LOG <- base::as.data.frame(SAMPLE_LOG)
      base::rownames(SAMPLE_LOG) <- base::c("Discovery",base::paste("sim",0:(lcsim_depth-1),sep = "-"))
      
      sample_log_filename <- base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_log.txt")
      if(save_results == TRUE) utils::write.table(x = SAMPLE_LOG,file = sample_log_filename,quote = FALSE,sep = "\t",row.names = TRUE)
      
      ### Evaluate enrichment of motifs
      base::cat(base::paste("Part 2: Selecting motifs significantly more represented in the data compared to the reference DB. \n"))
      
      part_1_res_file_name <- base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_log.txt")
      part_2_res_file_name <- base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_minp",lcminp,"_ove",paste(lcminove, collapse = "_"),".txt")
      all_part_2_res_file_name <- base::paste0(result_folder,"kmer_resample_",lcsim_depth,"_all_motifs.txt")
      
      ### Take along the SAMPLE_LOG
      # actual_reads = motif occurences in the sample set
      # sample_reads = motif occurences in the repeated random sampling
      sample_reads <- SAMPLE_LOG
      actual_reads <- base::as.data.frame(sample_reads[1,])
      base::colnames(actual_reads) <- base::colnames(sample_reads)
      sample_reads <- base::as.data.frame(sample_reads[-1,])
      nam <- base::colnames(actual_reads)
      
      ### Prepare framework
      motifs_df <- base::data.frame(base::matrix(base::rep(0, base::ncol(sample_reads)*7), ncol = 7), stringsAsFactors = FALSE)
      base::colnames(motifs_df) <- base::c("motif",	"counts",	"num_in_ref","avgRef",	"topRef",	"OvE",	"p.value")
      motifs_df$motif <- nam
      motifs_df$counts <- base::unlist(actual_reads)
      motifs_df$avgRef <- base::colMeans(sample_reads)
      motifs_df$topRef <- base::vapply(sample_reads, base::max, FUN.VALUE = base::c(1))
      motifs_df$OvE[motifs_df$avgRef > 0] <- motifs_df$counts[motifs_df$avgRef > 0]/motifs_df$avgRef[motifs_df$avgRef > 0]
      motifs_df$OvE[motifs_df$avgRef == 0] <- 1/(lcsim_depth*base::length(motif_region))
      motifs_df$p.value <- base::vapply(base::seq_len(base::ncol(sample_reads)), function(x) return(base::sum(sample_reads[,x] >= actual_reads[,x])/lcsim_depth),
                                        FUN.VALUE = base::c(1))
      
      motifs_df$avgRef <- base::round(motifs_df$avgRef, digits = 2)
      motifs_df$OvE <- base::round(motifs_df$OvE, digits = 3)
      motifs_df$p.value <- base::round(motifs_df$p.value, digits = 6)
      motifs_df$p.value[motifs_df$p.value == 0] <- 1/lcsim_depth
      
      ### Get the minimum fold change for every motif depending on the motif length
      if(base::length(lcminove) == 1){
        temp_minove <- base::rep(lcminove, base::ncol(actual_reads))
      } else {
        temp_minove <- base::rep(0, base::ncol(actual_reads))
        
        # the number of informative positions in discontinuous motifs is reduced by one 
        nam_len <- base::nchar(nam)
        nam_len[base::grep(".", nam, fixed = TRUE)] <- nam_len[base::grep(".", nam, fixed = TRUE)]-1
        
        for(i in base::seq_along(motif_length)){
          temp_minove[nam_len == motif_length[i]] <- lcminove[i]
        }
      }
      
      ### Filter the significant motifs by the following criteria:
      # p-value < lcminp (p-value = probability of observing the motif in a random reference subsample at sample set depth at least as often as in the sample set)
      # fold change >= lcminove (lcminove can be dependent on the motif length; fold change of sample set occurence compared to repeated random subsample occurence)
      # sample set occurence >= lckmer_mindepth
      sel_mot_df <- motifs_df[motifs_df$OvE >= temp_minove &
                                motifs_df$p.value < lcminp &
                                motifs_df$counts >= lckmer_mindepth,]
    }
    
    ### Method with Fisher's Exact Test
    if(local_method == "fisher"){
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
      base::colnames(motifs_df) <- base::c("motif", "counts", "num_in_ref")
      motifs_df$avgRef <- motifs_df$num_in_ref
      motifs_df$topRef <- base::rep(0, base::nrow(motifs_df))
      motifs_df$OvE <- base::rep(0, base::nrow(motifs_df))
      
      # assess the significance of motifs for being enriched in the sample set by Fisher's Exact Test (alternative = "greater")
      # Contigency table:
      # ---------------------------------------------------------------------------------------------------------
      #   Number sample sequences containing the motif   |   Number of reference sequences containing the motif
      # ---------------------------------------------------------------------------------------------------------
      # Number sample sequences not containing the motif | Number of reference sequences not containing the motif
      # ---------------------------------------------------------------------------------------------------------
      motifs_df$p.value <- stats::phyper(q=motifs_df$counts-1,
                                              m=motifs_df$num_in_ref+motifs_df$counts,
                                              n=base::length(base::unique(refseqs$CDR3b))+base::length(base::unique(seqs))-motifs_df$num_in_ref-motifs_df$counts,
                                              k=base::length(base::unique(seqs)), lower.tail = FALSE)
      motifs_df$p.value <- base::as.numeric(base::formatC(motifs_df$p.value, digits = 1, format = "e"))
      
      # assess the fold change enrichement ot a motif in the sample set compared to the reference set
      motifs_df$OvE <- base::round(motifs_df$counts/(motifs_df$num_in_ref+0.01)/base::length(base::unique(seqs))*base::length(base::unique(refseqs$CDR3b)), digits = 1)
      
      # for comparability with rrs normalize number in reference database to size of sample set
      motifs_df$avgRef <- base::round(motifs_df$num_in_ref/base::length(refseqs_motif_region)*base::length(motif_region), digits = 3)
      
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
      # sample set occurence >= lckmer_mindepth
      sel_mot_df <- motifs_df[motifs_df$OvE >= temp_minove &
                                motifs_df$p.value <= lcminp &
                                motifs_df$counts >= lckmer_mindepth,]
    }

    if(base::nrow(motifs_df) > 0) ALL_MOTIFS <- motifs_df
    if(base::nrow(sel_mot_df) > 0){
      SELECTED_MOTIFS <- sel_mot_df
      
      ### Remember all connections
      
      # Motifs may only be within a certain distance to be connected (determined by the parameter 'motif_distance_cutoff')
      CLONE_NETWORK <- foreach::foreach(i = base::seq_len(base::nrow(sel_mot_df)), .combine = base::rbind) %dopar% {
        # get all sample sequences with the motif
        all_ids <- base::grep(pattern = sel_mot_df$motif[i], x = all_motif_region, value = FALSE)
        
        # first assume connection between all ids
        if(base::length(all_ids) >= 2){
          combn_ids <- base::t(utils::combn(all_ids, m = 2))
        } else {
          combn_ids <- base::t(utils::combn(base::rep(1, 2), m = 2))
        }
        temp_df <- base::data.frame(V1 = combn_ids[,1], V2 = combn_ids[,2],
                                    type = base::rep("local", base::nrow(combn_ids)),
                                    tag = base::rep(sel_mot_df$motif[i], base::nrow(combn_ids)), stringsAsFactors = FALSE)
        
        base::return(temp_df)
      }
      
      ### Build clone network, but restrict to predefined motif distance cutoff
      if(!base::is.null(CLONE_NETWORK)){
        diffs <- base::abs(stringr::str_locate(string = sequences$CDR3b[CLONE_NETWORK$V1], pattern = CLONE_NETWORK$tag)[,1]-stringr::str_locate(string = sequences$CDR3b[CLONE_NETWORK$V2], pattern = CLONE_NETWORK$tag)[,1])
        CLONE_NETWORK <- CLONE_NETWORK[diffs < motif_distance_cutoff,]
      }
      
      base::cat(base::nrow(sel_mot_df),"significantly enriched motifs found in sample set.\n")
    } else {
      sel_mot_df <- base::data.frame(motif = 1,
                                     counts = 1,
                                     num_in_ref = 1,
                                     avgRef = 1,
                                     topRef = 1,
                                     OvE = 1,
                                     p.value = 1)[-1,]
      local_similarities <- FALSE
      
      base::cat("No significantly enriched motifs found in sample set.\n")
    }
    
    ### save results
    sel_motifs_fname <- base::paste0(result_folder, "local_similarities_minp_",lcminp, "_minove_", base::paste(lcminove, collapse = "_"), "_kmer_mindepth_", lckmer_mindepth, ".txt")
    all_motifs_fname <- base::paste0(result_folder, "all_motifs.txt")
    if(save_results == TRUE){
      utils::write.table(x = SELECTED_MOTIFS,file = sel_motifs_fname,quote = FALSE,sep = "\t",row.names = FALSE)
      utils::write.table(x = ALL_MOTIFS,file = all_motifs_fname,quote = FALSE,sep = "\t",row.names = FALSE)
    }
    
    t2_e <- base::Sys.time()
    base::cat("Part 1 cpu time:",t2_e-t1,base::units(t2_e-t1),"\n\n")
  } else {
    t2_e <- base::Sys.time()
    base::cat(base::paste("\nLocal similarity search is skipped. \n\n"))
  }
  
  #################################################################
  ##                 Part 2: Global similarities                 ##
  #################################################################
  SELECTED_GLOBALS <- NULL
  ALL_GLOBALS <- NULL
  
  motif_region_length <- base::nchar(motif_region)
  
  if(global_similarities == TRUE && base::any(motif_region_length > 1)){
    base::cat(base::paste("Part 2: Searching for global similarities.\n"))
    
    ### method as in GLIPH1.0
    if(global_method == "cutoff"){
      
      ### Set global convergency cutoff
      if(base::is.null(gccutoff)) if(base::length(all_motif_region) < 125) gccutoff <- 2 else gccutoff <- 1
      
      if(base::length(all_motif_region) > 1){
        res <- foreach::foreach(i = base::seq_len(base::length(all_motif_region)-1), .combine = rbind) %dopar% {
          global_con <- base::c()
          
          # Compute hamming distance between all sequences and select distances lower or equal to gccutoff
          dis <- stringdist::stringdist(a = all_motif_region[i],b = all_motif_region[(i+1):base::length(all_motif_region)],method = "hamming",nthread = 1) #F: hamming Inf problem
          global_ids <- ((i+1):base::length(all_motif_region))[dis <= gccutoff] #source of error : for short sequences there is a high chance that their hamming distance <=1 when we remove the boundries 3 at start and 3 at the end
          
          # create clone network
          if(base::length(global_ids)>0){
            global_con <- base::cbind(base::rep(i,base::length(global_ids)),global_ids,base::rep("global",base::length(global_ids)), base::rep("",base::length(global_ids)))
          }
          return(global_con)
        }
      } else res <- NULL
      
      ### append to clone network
      if(!base::is.null(res)){
        res <- base::data.frame(res, stringsAsFactors = FALSE) 
        base::colnames(res) <- c("V1", "V2", "type", "tag")
      } else {
        global_similarities <- FALSE
      }
    }
    
    ### (expandable) method of GLIPH2.0
    if(global_method == "fisher"){
      
      motif_region_length <- base::nchar(motif_region)
      
      ### Receive all structs (motif region with one variable position)
      exp_sample_seqs <- foreach::foreach(i = base::unique(motif_region_length), .combine = base::cbind) %dopar% {
        if(i == 1 ){
          return(NULL)
        }
        temp_seqs <- motif_region[motif_region_length == i]
        
        ret <- base::c()
        for(j in base::seq_len(i)){
          temp_struct <- temp_seqs
          substr(temp_struct, j,j) <- "%"
          ret <- cbind(ret, stringdist::qgrams(temp_struct, q = i))
        }
        
        return(ret)
      }
      
      # get all unique sample structs
      unq_sample_struct <- base::colnames(exp_sample_seqs)
      
      ### Receive all reference structs (motif region with one variable position)
      ref_motif_region_length <- base::nchar(refseqs_motif_region)
      exp_reference_seqs <- foreach::foreach(i = base::unique(motif_region_length)[base::unique(motif_region_length) %in% base::unique(ref_motif_region_length)], .combine = base::cbind) %dopar% {
        temp_seqs <- refseqs_motif_region[ref_motif_region_length == i]
        
        ret <- base::c()
        for(j in base::seq_len(i)){
          temp_struct <- temp_seqs
          substr(temp_struct, j,j) <- "%"
          ret <- cbind(ret, stringdist::qgrams(temp_struct, q = i))
        }
        
        return(ret)
      }
      # add columns of structs in sample set but not in reference set
      if(base::is.null(exp_reference_seqs)){
        exp_reference_seqs <- base::matrix(base::rep(0, base::ncol(exp_sample_seqs)), nrow = 1)
        base::colnames(exp_reference_seqs) <- base::colnames(exp_sample_seqs)
      } else {
        not_in_ref <- base::colnames(exp_sample_seqs)[!(base::colnames(exp_sample_seqs) %in% base::colnames(exp_reference_seqs))]
        if(base::length(not_in_ref) > 0){
          temp <- base::matrix(base::rep(0, base::length(not_in_ref)), nrow = 1)
          base::colnames(temp) <- not_in_ref
          exp_reference_seqs <- base::cbind(exp_reference_seqs, temp)
        }
      }
      exp_sample_seqs <- base::rbind(exp_sample_seqs, exp_reference_seqs[,base::colnames(exp_sample_seqs)])
      
      ### Summarize the frequencies of structs in sample and reference set
      structs_df <- base::data.frame(base::matrix(base::rep(0, base::ncol(exp_sample_seqs)*7), ncol = 7), stringsAsFactors = FALSE)
      base::colnames(structs_df) <- base::c("struct",	"counts",	"num_in_ref","avgRef",	"topRef",	"OvE",	"p.value")
      structs_df$struct <- unq_sample_struct
      structs_df$counts <- base::unlist(exp_sample_seqs[1,])
      structs_df$num_in_ref <- base::unlist(exp_sample_seqs[2,])
      structs_df$avgRef <- base::round(structs_df$num_in_ref/base::length(refseqs_motif_region)*base::length(motif_region), digits = 3)
      structs_df$OvE <- base::round(structs_df$counts/(structs_df$num_in_ref+0.01)/base::length(motif_region)*base::length(refseqs_motif_region), digits = 1)
      
      ### Asses significantly enriched structs by applying Fisher's Exact Test with the following contigency table:
      # --------------------------------------------------------------------------------------------------------------
      #    Number sample sequences contaninig the struct    |   Number of reference sequences containing the struct
      # --------------------------------------------------------------------------------------------------------------
      #   Number sample sequences not contaninig the struct | Number of reference sequences not containing the struct
      # --------------------------------------------------------------------------------------------------------------
      structs_df$p.value <- stats::phyper(q=structs_df$counts-1,
                                         m=structs_df$num_in_ref+structs_df$counts,
                                         n=base::length(base::unique(refseqs$CDR3b))+base::length(base::unique(seqs))-structs_df$num_in_ref-structs_df$counts,
                                         k=base::length(base::unique(seqs)), lower.tail = FALSE)
      structs_df$p.value <- base::as.numeric(base::formatC(structs_df$p.value, digits = 1, format = "e"))
      
      ### Filter significantly enriched structs by the following criteria:
      # p-value < gcminp (p-value of the Fisher's Exact Test carried out for the contigency table shown above)
      # fold change >= gcminove (fold change of sample set occurence compared to reference occurence)
      # sample set occurence >= gckmer_mindepth
      sel_struct_df <- structs_df[structs_df$OvE >= gcminove &
                                    structs_df$p.value <= gcminp &
                                    structs_df$counts >= gckmer_mindepth,]
      
      if(base::nrow(sel_struct_df) > 0) SELECTED_GLOBALS <- sel_struct_df
      if(base::nrow(structs_df) > 0) ALL_GLOBALS <- structs_df
      
      ### save results
      sel_structs_fname <- base::paste0(result_folder, "global_similarities_minp_",gcminp, "_minove_", base::paste(gcminove, collapse = "_"), "_kmer_mindepth_", gckmer_mindepth, ".txt")
      all_structs_fname <- base::paste0(result_folder, "all_global_similarities.txt")
      if(save_results == TRUE){
        utils::write.table(x = SELECTED_GLOBALS,file = sel_motifs_fname,quote = FALSE,sep = "\t",row.names = FALSE)
        utils::write.table(x = ALL_GLOBALS,file = all_motifs_fname,quote = FALSE,sep = "\t",row.names = FALSE)
      }
      
      ### Get all connections
      if(base::nrow(sel_struct_df) > 0){
        
        # Motifs may only be within a certain distance to be connected (determined by the parameter 'motif_distance_cutoff')
        sel_structs <- sel_struct_df$struct
        sel_structs <- base::gsub("%", ".", sel_structs, fixed = TRUE)
        
        res <- foreach::foreach(i = base::seq_len(base::nrow(sel_struct_df)), .combine = base::rbind) %dopar% {
          # get all sample sequences with the motif
          all_ids <- base::grep(pattern = paste0("^", sel_structs[i], "$"), x = all_motif_region, value = FALSE)
          
          # first assume connection between all ids
          if(base::length(all_ids) >= 2){
            combn_ids <- base::t(utils::combn(all_ids, m = 2))
          } else {
            combn_ids <- base::t(utils::combn(base::rep(1, 2), m = 2))
          }
          temp_df <- base::data.frame(V1 = combn_ids[,1], V2 = combn_ids[,2],
                                      type = base::rep("global", base::nrow(combn_ids)),
                                      tag = base::rep(sel_structs[i], base::nrow(combn_ids)), stringsAsFactors = FALSE)
          
          base::return(temp_df)
        }
        
        ### Filter connections for interchangeable amino acid pairs (BLOSUM62 >= 0)
        
        # load all amino acid pairs with a BLOSUM62-value >= 0
        BlosumVec <- NULL
        utils::data("BlosumVec",envir = base::environment(), package = "turboGliph")
        
        # filter
        if(!base::is.null(res) && all_aa_interchangeable == FALSE){
          res$tag <- base::gsub(".", "%", res$tag, fixed = TRUE)
          aa_pos <- stringr::str_locate(string = res$tag, pattern = "%")[,1]
          if(structboundaries == TRUE) aa_pos <- aa_pos + boundary_size
          aa_pairs <- paste0(base::substr(sequences$CDR3b[res$V1], aa_pos, aa_pos), base::substr(sequences$CDR3b[res$V2], aa_pos, aa_pos))
          res <- res[aa_pairs %in% BlosumVec,]
          if(base::nrow(res) == 0) res <- NULL
        }
        
        base::cat(base::nrow(sel_struct_df),"significantly enriched global structures found in sample set.\n")
      } else {
        sel_struct_df <- base::data.frame(motif = 1,
                                       counts = 1,
                                       num_in_ref = 1,
                                       avgRef = 1,
                                       topRef = 1,
                                       OvE = 1,
                                       p.value = 1)[-1,]
        glocal_similarities <- FALSE
        
        base::cat("No significantly enriched global structures  found in sample set.\n")
      }
    }
    
    
    ### Build/append clone network
    if(base::is.null(CLONE_NETWORK)){
      CLONE_NETWORK <- res
    } else {
      CLONE_NETWORK <- base::rbind(CLONE_NETWORK, res)
    }
    
    t3_e <- base::Sys.time()
    base::cat("Part 2 cpu time:",t3_e-t2_e,base::units(t3_e-t2_e),"\n\n")
  } else {
    part_2_res <- base::c()
    global_similarities <- FALSE
    t3_e <- base::Sys.time()
    base::cat(base::paste("\nGlobal similarity search is skipped. \n\n"))
  }
  
  ##################################################################
  ##                      Part 3: Clustering                      ##
  ##################################################################
  CLUSTER_LIST <- base::list()
  CLUSTER_INFOS <- NULL
  
  base::cat(base::paste("Part 3: Clustering sequences.\n"))
  
  ### Identification of N-nucleotides in CDR3b-sequences
  if(boost_local_significance == TRUE){
    # load all CDR3b sequence fragments encoded in the germline
    gTRB <- NULL
    utils::data("gTRB", package = "turboGliph", envir = environment())
    
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
  
  ### Clustering
  if(local_similarities == TRUE || global_similarities == TRUE){
    ### Restrict connections to V genes and public tcrs if desired
    if(vgene.info == TRUE){
      exclude <- base::c()
      if(vgene_match == "local" && local_similarities == TRUE){
        exclude <- base::which(sequences$TRBV[CLONE_NETWORK$V1] != sequences$TRBV[CLONE_NETWORK$V2] & CLONE_NETWORK$type == "local")
      }
      if(vgene_match == "global" && global_similarities == TRUE){
        exclude <- base::which(sequences$TRBV[CLONE_NETWORK$V1] != sequences$TRBV[CLONE_NETWORK$V2] & CLONE_NETWORK$type == "global")
      } 
      if(vgene_match == "all" && local_similarities == TRUE){
        exclude <- base::which(sequences$TRBV[CLONE_NETWORK$V1] != sequences$TRBV[CLONE_NETWORK$V2])
      }
      
      if(base::length(exclude) > 0) CLONE_NETWORK <- CLONE_NETWORK[-exclude,]
    }
    if(patient.info == TRUE){
      exclude <- base::c()
      if(public_tcrs == "local" && local_similarities == TRUE){
        exclude <- base::which(sequences$patient[CLONE_NETWORK$V1] != sequences$patient[CLONE_NETWORK$V2] & CLONE_NETWORK$type == "local")
      }
      if(public_tcrs == "global" && global_similarities == TRUE){
        exclude <- base::which(sequences$patient[CLONE_NETWORK$V1] != sequences$patient[CLONE_NETWORK$V2] & CLONE_NETWORK$type == "global")
      } 
      if(public_tcrs == "none" && local_similarities == TRUE){
        exclude <- base::which(sequences$patient[CLONE_NETWORK$V1] != sequences$patient[CLONE_NETWORK$V2])
      }
      
      if(base::length(exclude) > 0) CLONE_NETWORK <- CLONE_NETWORK[-exclude,]
    }
    
    CLONE_NETWORK <- base::unique(CLONE_NETWORK)
    
    # temporary save ID clone network
    ID_CLONE_NETWORK <- CLONE_NETWORK
    
    if(base::nrow(CLONE_NETWORK) > 0){
      
      ### flicked clustering
      if(clustering_method == "GLIPH1.0"){
        clone_network_for_clustering <- base::as.matrix(CLONE_NETWORK[,base::c(1,2)])
        clone_network_for_clustering[,1] <- base::as.character(clone_network_for_clustering[,1])
        clone_network_for_clustering[,2] <- base::as.character(clone_network_for_clustering[,2])
        
        ### Receive the network from the clone network using the igraph-package
        gr <- igraph::graph_from_edgelist(el = clone_network_for_clustering,directed = FALSE)
        cm <- igraph::components(gr)
        
        ### Extract the cluster compositions of the created network
        all_leaders <- base::c()
        for(i in base::seq_len(cm$no)){
          # extract the member tuples of the current cluster
          members_id <- base::as.numeric(base::names(cm$membership)[cm$membership == i])
          members_id <- base::sort(base::unique(members_id))
          members_details <- sequences[members_id,]
          
          # determine the cluster size (= number of redundant (!) members)
          members <- base::sort(base::unique(members_details$CDR3b))
          csize <- base::nrow(members_details)
          
          # get the representative leader tag of the cluster
          leader <- base::paste("CRG",members[1],sep = "-")
          while_counter <- 1
          temp_leader <- leader
          while(temp_leader %in% all_leaders){
            temp_leader <- base::paste0(leader, "-", while_counter)
            while_counter <- while_counter + 1
          }
          leader <- temp_leader
          all_leaders <- base::c(all_leaders, leader)
          
          if(i == 1){
            CLUSTER_INFOS <- base::data.frame(cluster_size = csize, tag = leader,members = base::paste(members,collapse = " "),
                                           stringsAsFactors = FALSE)
          } else {
            CLUSTER_INFOS <- base::rbind(CLUSTER_INFOS, base::data.frame(cluster_size = csize, tag = leader,members = base::paste(members,collapse = " "),
                                                                   stringsAsFactors = FALSE))
          }
          CLUSTER_LIST[[leader]] <- members_details
        }
      }
      
      ### isolated clustering
      if(clustering_method == "GLIPH2.0"){
        
        # prepare network for clustering
        clone_network_for_clustering <- CLONE_NETWORK
        clone_network_for_clustering$V1 <- paste(clone_network_for_clustering$V1, clone_network_for_clustering$tag, clone_network_for_clustering$type, sep = "$%$")
        clone_network_for_clustering$V2 <- paste(clone_network_for_clustering$V2, clone_network_for_clustering$tag, clone_network_for_clustering$type, sep = "$%$")
        
        ### Receive the network from the clone network using the igraph-package
        gr <- igraph::graph_from_edgelist(el = base::as.matrix(clone_network_for_clustering[,base::c(1,2)]),directed = FALSE)
        cm <- igraph::components(gr)
        
        ### Extract the cluster compositions of the created network
        for(i in base::seq_len(cm$no)){
          
          # extract the member tuples of the current cluster
          cluster_id <- base::names(cm$membership)[cm$membership == i]
          cluster_id <- base::data.frame(base::matrix(base::unlist(base::strsplit(cluster_id, split = "$%$", fixed = TRUE)), ncol = 3, byrow = TRUE),
                                         stringsAsFactors = FALSE)
          base::colnames(cluster_id) <- c("ID", "tag", "type")
          cluster_id$ID <- base::as.numeric(cluster_id$ID)
          
          members_id <- cluster_id$ID
          members_id <- base::sort(base::unique(members_id))
          members_details <- sequences[members_id,]
          
          # determine the cluster size (= number of redundant (!) members)
          members <- base::sort(base::unique(members_details$CDR3b))
          csize <- base::nrow(members_details)
          
          # get unique cluster tag (motif/struct_position/amino acids_vgene_patient)
          cluster_name <- ""
          cluster_stats <- NULL
          if(cluster_id$type[1] == "local"){
            # local
            all_pos <- stringr::str_locate(string = base::substr(members,
                                                                 boundary_size*structboundaries+1,
                                                                 base::nchar(members)-boundary_size*structboundaries),
                                           pattern = cluster_id$tag[1])+boundary_size*structboundaries
            
            cluster_name <- cluster_id$tag[1]
            cluster_name <- base::paste(cluster_name, base::paste0(base::min(all_pos[,1]), "-", base::max(all_pos[,1])), sep = "_")
            if(vgene_match == "all" || vgene_match == "local" && vgene.info == TRUE){
              cluster_name <- base::paste(cluster_name, members_details$TRBV[1], sep = "_")
            }
            if(public_tcrs == "none" || public_tcrs == "local" && patient.info == TRUE){
              cluster_name <- base::paste(cluster_name, members_details$patient[1], sep = "_")
            }
            
            cluster_stats <- SELECTED_MOTIFS[SELECTED_MOTIFS$motif == cluster_id$tag[1],]
          } else {
            # gobal
            temp_aa_pos <- stringr::str_locate(string = cluster_id$tag[1], pattern = "%")[,1] + boundary_size*structboundaries
            
            cluster_name <- cluster_id$tag[1]
            cluster_name <- base::paste(cluster_name,
                                        base::paste(base::sort(base::unique(base::substr(members, temp_aa_pos, temp_aa_pos))),
                                                    collapse = ""),
                                        sep = "_")
            if(vgene_match == "all" || vgene_match == "global" && vgene.info == TRUE){
              cluster_name <- base::paste(cluster_name, members_details$TRBV[1], sep = "_")
            }
            if(public_tcrs == "none" || public_tcrs == "global" && patient.info == TRUE){
              cluster_name <- base::paste(cluster_name, members_details$patient[1], sep = "_")
            }
            
            cluster_stats <- SELECTED_GLOBALS[SELECTED_GLOBALS$struct == cluster_id$tag[1],]
          }
          
          # recalculate p.value
          p.value <- stats::phyper(q=base::length(members)-1,
                                   m=cluster_stats$num_in_ref+base::length(members),
                                   n=base::length(refseqs_motif_region)+base::length(motif_region)-cluster_stats$num_in_ref-base::length(members),
                                   k=base::length(motif_region), lower.tail = FALSE)
          
          # boost local motifs for containing conserved N-Nucleotides
          if(boost_local_significance == TRUE && cluster_id$type[1] == "local"){
            # get first position where the motif starts
            start_motif <- base::min(all_pos)
            
            # get last position where the motif ends
            stop_motif <- base::max(all_pos)+base::nchar(cluster_id$tag[1])-1
            
            all_motif_ranges <- base::substr(members, all_pos[,1], all_pos[,2])
            
            motif_in_range <- (range_df$start_1[members_id] <= stop_motif & range_df$stop_1[members_id] >= start_motif) |
                              (range_df$start_2[members_id] <= stop_motif & range_df$stop_2[members_id] >= start_motif)
            
            p.value <- p.value/(2^(base::length(base::grep(".*[A-Z].*", all_motif_ranges))))
          }
          
          # return
          if(i == 1){
            CLUSTER_INFOS <- base::data.frame(type = base::rep(cluster_id$type[1], 1),
                                              tag = base::rep(cluster_id$tag[1], 1),
                                              name = base::rep(cluster_name, 1),
                                              cluster_size = base::rep(csize, 1), 
                                              unique_cdr3_sample = base::rep(base::length(members), 1),
                                              unique_cdr3_ref = base::rep(cluster_stats$num_in_ref, 1),
                                              OvE = base::rep(cluster_stats$OvE, 1),
                                              p.value = base::rep(p.value, 1),
                                              members = base::paste(members,collapse = " "),
                                              stringsAsFactors = FALSE)
          } else {
            CLUSTER_INFOS <- base::rbind(CLUSTER_INFOS,
                                         base::data.frame(type = base::rep(cluster_id$type[1], 1),
                                                          tag = base::rep(cluster_id$tag[1], 1),
                                                          name = base::rep(cluster_name, 1),
                                                          cluster_size = base::rep(csize, 1), 
                                                          unique_cdr3_sample = base::rep(base::length(members), 1),
                                                          unique_cdr3_ref = base::rep(cluster_stats$num_in_ref, 1),
                                                          OvE = base::rep(cluster_stats$OvE, 1),
                                                          p.value = base::rep(p.value, 1),
                                                          members = base::paste(members,collapse = " "),
                                                          stringsAsFactors = FALSE))
          }
          CLUSTER_LIST[[cluster_name]] <- members_details
        }
      }
      
    }
  }
  
  ### Add Singlets to CLONE_NETWORK
  singlets_ids <- (base::seq_along(all_motif_region))[-base::as.numeric(base::unique(base::c(CLONE_NETWORK$V1, CLONE_NETWORK$V2)))]
  if(base::length(singlets_ids) > 0){
    singlet_df <- base::data.frame(V1 = singlets_ids, V2 = singlets_ids, type = rep("singleton", base::length(singlets_ids)),
                                   tag = rep("", base::length(singlets_ids)))
    
    CLONE_NETWORK <- base::rbind(CLONE_NETWORK, singlet_df)
  }
  
  ### regain cdr3b sequences represented by the IDs
  if(base::nrow(CLONE_NETWORK) > 0){
    CLONE_NETWORK$V1 <- sequences$CDR3b[base::as.numeric(CLONE_NETWORK$V1)]
    CLONE_NETWORK$V2 <- sequences$CDR3b[base::as.numeric(CLONE_NETWORK$V2)]
    CLONE_NETWORK <- CLONE_NETWORK[CLONE_NETWORK$V1 != CLONE_NETWORK$V2,]
  }
  
  ### exclude clusters with less size than cluster_min_size
  if(!base::is.null(CLUSTER_INFOS)){
    eliminate_ids <- base::which(CLUSTER_INFOS$cluster_size < cluster_min_size)
    if(base::length(eliminate_ids) > 0){
      CLUSTER_INFOS <- CLUSTER_INFOS[-eliminate_ids,]
      CLUSTER_LIST <- CLUSTER_LIST[-eliminate_ids]
    }
    
    if(base::nrow(CLUSTER_INFOS) == 0){
      CLUSTER_INFOS <- NULL
      CLUSTER_LIST <- base::list()
    }
  }
  
  ### Save clustering results
  if(!base::is.null(CLUSTER_INFOS)){
    ### save part_4_res_list
    save_cluster_list_df <- foreach::foreach(i = base::seq_along(CLUSTER_LIST), .combine = "rbind") %dopar% {
      temp <- CLUSTER_LIST[[i]]
      temp <- base::cbind(base::data.frame(tag = base::rep(base::names(CLUSTER_LIST)[i], base::nrow(temp))), temp)
      return(temp)
    }
  } else {
    save_cluster_list_df <- base::c()
  }
  
  if(save_results == TRUE) utils::write.table(x = save_cluster_list_df,file = base::paste0(result_folder, "cluster_member_details.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  if(save_results == TRUE) utils::write.table(x = CLONE_NETWORK,file = base::paste0(result_folder, "clone_network.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = FALSE)
  
  t4_e <- base::Sys.time()
  base::cat("Part 3 cpu time:",t4_e-t3_e,base::units(t4_e-t3_e), "\n\n\n")
  
  doParallel::stopImplicitCluster()
  
  #################################################################
  ##                       Part 4: Scoring                       ##
  #################################################################
  
  base::cat("Part 4: Scoring convergence groups\n")
  
  if(!base::is.null(save_cluster_list_df)){
    sc_gliph_version <- 0
    if(scoring_method == "GLIPH1.0") sc_gliph_version <- 1
    if(scoring_method == "GLIPH2.0") sc_gliph_version <- 2
    
    scoring_res <- cluster_scoring(cluster_list = CLUSTER_LIST,
                                               cdr3_sequences = cdr3_sequences,
                                               refdb_beta = refdb_beta,
                                               v_usage_freq = v_usage_freq,
                                               cdr3_length_freq = cdr3_length_freq,
                                               ref_cluster_size = ref_cluster_size,
                                               sim_depth = scoring_sim_depth,
                                               gliph_version = sc_gliph_version,
                                               hla_cutoff = hla_cutoff,
                                               n_cores = n_cores)
    
    CLUSTER_INFOS <- base::cbind(CLUSTER_INFOS, scoring_res)
    # some reordering of columns
    CLUSTER_INFOS <- CLUSTER_INFOS[, base::c(base::colnames(CLUSTER_INFOS)[base::colnames(CLUSTER_INFOS) != "members"], "members")]
  }
  
  if(save_results == TRUE) utils::write.table(x = CLUSTER_INFOS,file = base::paste0(result_folder, "convergence_groups.txt"),quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
  
  base::cat("\n")
  t5_e <- base::Sys.time()
  base::cat("Part 4 cpu time:",t5_e-t4_e, base::units(t5_e-t4_e),"\n\n\n")
  
  #################################################################
  ##                        Create output                        ##
  #################################################################
  
  # set all theoretical numeric values (actually characters) to numeric values
  if(base::is.data.frame(SAMPLE_LOG)){
    for(i in base::seq_len(base::ncol(SAMPLE_LOG))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(SAMPLE_LOG[,i])))) == FALSE) SAMPLE_LOG[,i] <- base::as.numeric(SAMPLE_LOG[,i])
    }
  }
  if(base::is.data.frame(SELECTED_MOTIFS)){
    for(i in base::seq_len(base::ncol(SELECTED_MOTIFS))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(SELECTED_MOTIFS[,i])))) == FALSE) SELECTED_MOTIFS[,i] <- base::as.numeric(SELECTED_MOTIFS[,i])
    }
  }
  if(base::is.data.frame(ALL_MOTIFS)){
    for(i in base::seq_len(base::ncol(ALL_MOTIFS))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(ALL_MOTIFS[,i])))) == FALSE) ALL_MOTIFS[,i] <- base::as.numeric(ALL_MOTIFS[,i])
    }
  }
  if(base::is.data.frame(SELECTED_GLOBALS)){
    for(i in base::seq_len(base::ncol(SELECTED_GLOBALS))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(SELECTED_GLOBALS[,i])))) == FALSE) SELECTED_GLOBALS[,i] <- base::as.numeric(SELECTED_GLOBALS[,i])
    }
  }
  if(base::is.data.frame(ALL_GLOBALS)){
    for(i in base::seq_len(base::ncol(ALL_GLOBALS))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(ALL_GLOBALS[,i])))) == FALSE) ALL_GLOBALS[,i] <- base::as.numeric(ALL_GLOBALS[,i])
    }
  }
  if(base::is.data.frame(CLUSTER_INFOS)){
    for(i in base::seq_len(base::ncol(CLUSTER_INFOS))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(CLUSTER_INFOS[,i])))) == FALSE) CLUSTER_INFOS[,i] <- base::as.numeric(CLUSTER_INFOS[,i])
    }
  }
  if(base::is.list(CLUSTER_LIST)){
    for(i in base::seq_along(CLUSTER_LIST)){
      if(base::is.data.frame(CLUSTER_LIST[[i]])){
        for(j in base::seq_len(base::ncol(CLUSTER_LIST[[i]]))){
          if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(CLUSTER_LIST[[i]][,j])))) == FALSE) CLUSTER_LIST[[i]][,j] <- base::as.numeric(CLUSTER_LIST[[i]][,j])
        }
      }
    }
  }
  
  ### the output list
  output <- base::list(sample_log = NULL,
                       motif_enrichment = NULL,
                       global_enrichment = NULL,
                       connections = NULL,
                       cluster_properties = NULL,
                       cluster_list = NULL,
                       parameters = NULL)
  
  output$sample_log <- SAMPLE_LOG
  output$motif_enrichment <- base::list(selected_motifs = SELECTED_MOTIFS,all_motifs = ALL_MOTIFS)
  output$global_enrichment <- base::list(selected_structs = SELECTED_GLOBALS, all_structs = ALL_GLOBALS)
  output$connections <- CLONE_NETWORK
  output$cluster_properties <- CLUSTER_INFOS
  output$cluster_list <- CLUSTER_LIST
  
  output$parameters <- base::list(result_folder = result_folder,
                                  ref_cluster_size = ref_cluster_size,
                                  
                                  min_seq_length = min_seq_length,
                                  accept_sequences_with_C_F_start_end = accept_sequences_with_C_F_start_end,
                                  structboundaries = structboundaries,
                                  boundary_size = boundary_size,
                                  
                                  local_similarities = local_similarities,
                                  local_method = local_method,
                                  lcsim_depth = lcsim_depth,
                                  motif_length = motif_length,
                                  motif_distance_cutoff = motif_distance_cutoff,
                                  lcminp = lcminp,
                                  lcminove = lcminove,
                                  lckmer_mindepth = lckmer_mindepth,
                                  discontinuous_motifs = discontinuous_motifs,
                                  boost_local_significance = boost_local_significance,
                                  cdr3_len_stratify = cdr3_len_stratify, 
                                  vgene_stratify = vgene_stratify, 
                                  
                                  global_similarities = global_similarities, 
                                  global_method = global_method, 
                                  gccutoff = gccutoff,
                                  gcminp = gcminp,
                                  gcminove = gcminove,
                                  gckmer_mindepth = gckmer_mindepth,
                                  all_aa_interchangeable = all_aa_interchangeable,
                                  
                                  clustering_method = clustering_method,
                                  vgene_match = vgene_match,
                                  public_tcrs = public_tcrs,
                                  cluster_min_size = cluster_min_size,
                                  
                                  scoring_method = scoring_method,
                                  scoring_sim_depth = scoring_sim_depth,
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