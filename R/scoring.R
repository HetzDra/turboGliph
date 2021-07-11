#' Scoring cdr3 clusters as performed in GLIPH and GLIPH2 algorithm
#'
#' @description With this method the scores of cdr3 clusters are calculated as in the GLIPH and GLIPH2 algorithm.
#' Depending on the information provided, a final score is calculated based on up to five cluster properties:
#' cluster size, enrichment of cdr3 lengths, enrichment of V genes, enrichment of clonal expansions and enrichment of a common HLA alleles.
#' @param cluster_list list. Each element of this list contains a data frame in which the CDR3b sequences and additional information necessary for
#' scoring are provided. Corresponds to the \code{$cluster_list} element of the output of the functions \code{turbo_gliph} and \code{gliph2}.
#' @param cdr3_sequences vector or dataframe. This dataframe must contain the cdr3 sequences and optional additional information.
# If \code{type = "beta"} (default) column named "CDR3b" is required.
# If \code{type = "alpha"} column named "CDR3a" is required.
# If \code{type = "paired chains"} columns named "CDR3a" (alpha sequences) and "CDR3b" (beta sequences) are required.
#' The columns must be named as specified in the following list in arbitrary order.
#' \itemize{
# \item "CDR3a": cdr3 sequences of alpha chains
#' \item "CDR3b": cdr3 sequences of beta chains
# \item "TRAV": optional. V-genes of alpha chains
#' \item "TRBV": optional. V-genes of beta chains
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
# \item "human_v2.0_CD4": Reference database of 772,312 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD8": Reference database of 573,211 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "human_v2.0_CD48": Reference database of 1,329,314 CDR3b sequences of naive human CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD4": Reference database of 22,268 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
# \item "mouse_v1.0_CD8": Reference database of 55,847 CDR3b sequences of naive murine CD48+ T cells provided for GLIPH2.
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
#' @param gliph_version numeric. Either \code{1} for GLIPH or \code{2} for GLIPH2 algorithm. The scoring of the
#' algorithms differs only in the calculation of the total score. GLIPH2 calculates the product of all individual scores,
#' while GLIPH multiplies this product additionally by 0.064.
#' @param sim_depth numeric. By default 1000. Simulated resampling depth for non-parametric convergence significance tests.
#' A higher number will take longer to run but will produce more reproducible results.
#' @param hla_cutoff numeric. By default 0.1. Defines the threshold of HLA probability scores below which HLA alleles are considered significant.
#' @param n_cores numeric. Number of cores to use, by default 1. In case of \code{NULL} it will be set to number of cores in your machine minus 2.
#'
#' @return This function produces one file in the \code{result_folder} named \code{"GLIPH_scoring_results.txt"}
#' containing the same information as in the returned data frame.
#' The data frame contains the cluster scoring results. The first columns provides the \code{representative_seq} for any evaluated cluster.
#' In the second column the total scores are stored. Additional columns contain up to five scores (cluster size, cdr3 length enrichment, V-gene enrichment,
#' enrichment of clonal expansion and enrichment of common HLA) used to evaluate the total score.
#'
#' @examples
#' utils::data("gliph_input_data")
#' 
#' res <- turbo_gliph(cdr3_sequences = gliph_input_data[base::seq_len(200),],
#' sim_depth = 100,
#' n_cores = 1)
#'
#' scoring_results <- cluster_scoring(cluster_list = res$cluster_list,
#' cdr3_sequences = gliph_input_data[base::seq_len(200),],
#' refdb_beta = "gliph_reference",
#' gliph_version = 1,
#' sim_depth = 100,
#' n_cores = 1)
#'
#' @references Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @references https://github.com/immunoengineer/gliph
#' @import foreach
#' @export
cluster_scoring <- function(cluster_list,
                            cdr3_sequences,
                            refdb_beta = "gliph_reference",
                            v_usage_freq = NULL,
                            cdr3_length_freq = NULL,
                            ref_cluster_size = "original",
                            gliph_version = 1,
                            sim_depth = 1000,
                            hla_cutoff = 0.1,
                            n_cores = 1){
  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################
  
  ### cluster_list
  if(!base::is.list(cluster_list)) base::stop("parameter 'cluster_list' has to be an object of class 'list'.")
  
  ### cdr3_sequences
  if(base::is.atomic(cdr3_sequences)) cdr3_sequences <- base::data.frame("CDR3b" = cdr3_sequences)
  if(!base::is.data.frame(cdr3_sequences)) base::stop("parameter 'cdr3_sequences' has to be an object of class 'data.frame'.")
  cdr3_sequences[] <- base::lapply(cdr3_sequences, as.character)
  
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
  if(!(ref_cluster_size %in% base::c("original", "simulated") || !base::is.character(ref_cluster_size) || base::length(ref_cluster_size) > 1)){
    base::stop("ref_cluster_size has to be either 'original' or 'simulated'.")
  }
  
  ### gliph_version
  if(!(gliph_version %in% base::c(1,2))) base::stop("gliph_version has to be either 1 or 2.")
  
  ### sim_depth
  if(!base::is.numeric(sim_depth)) base::stop("sim_depth has to be numeric")
  if(base::length(sim_depth) > 1) base::stop("sim_depth has to be a single number")
  if(sim_depth < 1) base::stop("sim_depth must be at least 1")
  sim_depth <- base::round(sim_depth)
  
  ### hla_cutoff
  if(!base::is.numeric(hla_cutoff)) base::stop("hla_cutoff has to be numeric")
  if(base::length(hla_cutoff) > 1) base::stop("hla_cutoff has to be a single number")
  if(hla_cutoff > 1 || hla_cutoff < 0) base::stop("hla_cutoff must be between 0 and 1")
  
  ### n_cores
  if(base::is.null(n_cores)) n_cores <- base::max(1, parallel::detectCores()-2) else {
    if(!base::is.numeric(n_cores)) base::stop("n_cores has to be numeric")
    if(base::length(n_cores) > 1) base::stop("n_cores has to be a single number")
    if(n_cores < 1) base::stop("n_cores must be at least 1")

    n_cores <- base::min(n_cores, parallel::detectCores()-2)
  }
  
  #################################################################
  ##                         Preparation                         ##
  #################################################################

  ### Which scores can be calculated from the dataset?
  score_names <- base::c("network.size.score","cdr3.length.score")
  if("TRBV" %in% base::colnames(cdr3_sequences)){
    vgene_info <- TRUE
    score_names <- base::c(score_names, "vgene.score")
  } else vgene_info <- FALSE
  if("counts" %in% base::colnames(cdr3_sequences)) {
    counts_info <- TRUE
    cdr3_sequences$counts <- base::as.numeric(cdr3_sequences$counts)
    cdr3_sequences$counts[base::is.na(cdr3_sequences$counts)] <- 1
    score_names <- base::c(score_names, "clonal.expansion.score")
  } else counts_info <- FALSE
  if("patient" %in% base::colnames(cdr3_sequences)) patient_info <- TRUE else patient_info <- FALSE
  if("HLA" %in% base::colnames(cdr3_sequences)) hla_info <- TRUE else hla_info <- FALSE
  if(hla_info == TRUE && patient_info == TRUE){
    if("patient" %in% base::colnames(cdr3_sequences) && "HLA" %in% base::colnames(cdr3_sequences)){
      cdr3_sequences <- cdr3_sequences[cdr3_sequences$HLA != "",]
      if(nrow(cdr3_sequences > 0)) score_names <- base::c(score_names, "hla.score", "lowest.hlas") else hla_info <- FALSE
    }
  }

  ### load or generate reference tables from reference database
  # ref_cluster_sizes:            data frame containing the cluster size in the first and the probability of observing a cluster with this size
  #                               in a sample from the reference database in the second column
  # vgene_ref_frequencies:        vector containing the (relative) frequencies of v gene usage
  # cdr3_length_ref_frequencies:  vector containing the (relative) frequencies of CDR3b lengths
  utils::data(ref_cluster_sizes, envir = base::environment(), package = "turboGliph")
  ref_cluster_sizes <- ref_cluster_sizes[[ref_cluster_size]]

  # maybe changed in future
  reference_list <- NULL
  refdb_name <- "gliph_reference"
  utils::data("reference_list",envir = base::environment(), package = "turboGliph")
  vgene_ref_frequencies <- reference_list[[refdb_name]][[2]]$freq
  cdr3_length_ref_frequencies <- reference_list[[refdb_name]][[3]]$freq
  
  if(!base::is.null(v_usage_freq)) vgene_ref_frequencies <- base::as.numeric(v_usage_freq[,2])
  if(!base::is.null(cdr3_length_freq)) cdr3_length_ref_frequencies <- base::as.numeric(cdr3_length_freq[,2])

  ### Obtain the distribution of all patients and HLA alleles in the sample
  # all_patients:     vector containing all unique patient indices
  # all_hlas:         vector containing all unique HLA alleles
  # all_patient_hlas: list whose elements are named after the patients and contain the patient's HLA alleles in a vector.
  all_patient_hlas <- base::c()
  if(hla_info == TRUE && patient_info == TRUE){
    cdr3_sequences$patient <- gsub(":.*", "",cdr3_sequences$patient)
    all_patients <- base::sort(base::unique(cdr3_sequences$patient))
    all_hlas <- base::unlist(base::strsplit(base::unique(cdr3_sequences$HLA), ","))
    all_hlas <- base::sort(base::unique(base::gsub(":.*", "", all_hlas, perl = TRUE)))
    num_patients <- base::length(all_patients)
    num_HLAs <- base::length(all_hlas)

    all_patient_hlas <- base::lapply(base::unique(cdr3_sequences$patient), function(x){
      base::sort(base::unique(base::gsub(":.*", "", base::unlist(base::strsplit(cdr3_sequences$HLA[cdr3_sequences$patient == x][1], ",")), perl = TRUE)))
    })
    base::names(all_patient_hlas) <- base::unique(cdr3_sequences$patient)

    all_hlas <- base::data.frame(HLA = all_hlas)
    all_hlas$counts <- base::apply(all_hlas, 1, function(x){
      val <- 0
      for(pat in all_patients){
        if(x %in% all_patient_hlas[[pat]]) val <- val+1
      }
      val
    })
  }

  #################################################################
  ##                           Scoring                           ##
  #################################################################
  
  doParallel::registerDoParallel(n_cores)
  
  actCluster <- NULL
  res <- foreach::foreach(actCluster = base::seq_along(cluster_list)) %dopar% {
    
    ### Get sequences and information of current cluster
    act_seq_infos <- cluster_list[[actCluster]]
    num_members <- base::nrow(act_seq_infos) # number of ALL members
    ori_num_members <- base::length(base::unique(act_seq_infos$CDR3b)) # number of all unique CDR3b sequences
    all_scores <- base::c()
    
    ### Get network size score from lookup table
    score_network_size <- 1
    nearest_sample_size <- base::order(base::abs(1-base::as.numeric(base::colnames(ref_cluster_sizes)[-1])/base::nrow(cdr3_sequences)))[1]
    if(ori_num_members > 100){
      score_network_size <-  ref_cluster_sizes[100,nearest_sample_size+1]
    } else {
      score_network_size <-  ref_cluster_sizes[ori_num_members,nearest_sample_size+1]
    }
    all_scores <- base::c(all_scores, score_network_size)
    
    ### Enrichment of CDR3 length (spectratype) within cluster
    score_cdr3_length <- base::c()
    # calculate score of sample (product of all frequencies)
    pick_freqs <- base::data.frame(base::table(base::nchar(base::unique(act_seq_infos$CDR3b))))
    base::colnames(pick_freqs) <- base::c("object", "probs")
    pick_freqs$probs <- pick_freqs$probs/ori_num_members
    sample_score <- base::round(base::prod(pick_freqs$probs), digits = 3)
    
    # generate random subsamples
    random_subsample <- base::list()
    for(i in base::seq_len(sim_depth)){
      random_subsample[[i]] <- base::sample.int(n = base::length(cdr3_length_ref_frequencies), size = ori_num_members,
                                                prob = cdr3_length_ref_frequencies, replace = TRUE)
    }
    
    # calculate score of subsamples (product of all frequencies)
    pick_freqs <- stringdist::seq_qgrams(.list = random_subsample)[,-1]/ori_num_members
    pick_freqs[pick_freqs == 0] <- 1
    pick_freqs <- base::round(base::exp(base::colSums(base::log(pick_freqs))), digits = 3) # vectorized way to calculate the product of each column
    if(gliph_version == 1){
      score_cdr3_length <- base::sum(pick_freqs >= sample_score)/sim_depth
    } else {
      score_cdr3_length <- base::sum(pick_freqs > sample_score)/sim_depth
    }
    if(score_cdr3_length == 0) score_cdr3_length <- 1/sim_depth # minimum score of 1/sim_depth
    
    all_scores <- base::c(all_scores, score_cdr3_length)
    
    ### Enrichment of v genes within cluster
    score_vgene <- base::c()
    if(vgene_info == TRUE){
      
      # calculate score of sample (product of all frequencies)
      pick_freqs <- base::data.frame(base::table(act_seq_infos$TRBV))
      base::colnames(pick_freqs) <- base::c("object", "probs")
      pick_freqs$probs <- pick_freqs$probs/num_members
      sample_score <- base::round(base::prod(pick_freqs$probs), digits = 3)
      
      # generate random subsamples
      random_subsample <- base::list()
      for(i in base::seq_len(sim_depth)){
        random_subsample[[i]] <- base::sample.int(n = base::length(vgene_ref_frequencies), size = num_members,
                                                  prob = vgene_ref_frequencies, replace = TRUE)
      }
      
      # calculate score of subsamples (product of all frequencies)
      pick_freqs <- stringdist::seq_qgrams(.list = random_subsample)[,-1]/num_members
      pick_freqs[pick_freqs == 0] <- 1
      pick_freqs <- base::round(base::exp(base::colSums(base::log(pick_freqs))), digits = 3) # vectorized way to calculate the product of each column
      if(gliph_version == 1){
        score_vgene <- base::sum(pick_freqs >= sample_score)/sim_depth
      } else {
        score_vgene <- base::sum(pick_freqs > sample_score)/sim_depth
      }
      
      if(score_vgene == 0) score_vgene <- 1/sim_depth # minimum score of 1/sim_depth
      all_scores <- base::c(all_scores, score_vgene)
    }
    
    ### Enrichment of clonal expansion within cluster
    score_clonal_expansion <- base::c()
    if(counts_info == TRUE){
      sample_score <- base::sum(base::as.numeric(act_seq_infos$counts))/num_members
      counter <- 0
      for(i in base::seq_len(sim_depth)){
        random_subsample <- base::sample(x = cdr3_sequences$counts, size = num_members, replace = FALSE)
        test_score <- base::sum(base::as.numeric(random_subsample))/num_members
        if(test_score>=sample_score) counter <- counter+1
      }
      if(counter == 0) score_clonal_expansion <- 1/sim_depth else score_clonal_expansion <- counter/sim_depth
      score_clonal_expansion <- base::round(score_clonal_expansion, digits = 3)
      
      all_scores <- base::c(all_scores, score_clonal_expansion)
    }
    
    ### Enrichment of common HLA among donor TCR contributors in cluster
    score_hla <- base::c()
    lowest_hla <- ""
    if(hla_info == TRUE && patient_info == TRUE){
      act_seq_infos <- act_seq_infos[act_seq_infos$HLA != "",]
      act_seq_infos$patient <- gsub(":.*", "",act_seq_infos$patient)
      
      score_hla <- 1
      for(act_hla in base::seq_len(num_HLAs)){
        crg_patient_count <- base::length(base::unique(act_seq_infos$patient))
        crg_patient_hla_count <- base::sum(base::unlist(base::lapply(all_patient_hlas[base::unique(act_seq_infos$patient)], function(x){
          if(all_hlas$HLA[act_hla] %in% x) 1 else 0
        })))
        if(crg_patient_hla_count > 1){
          act_Prob <- base::sum(base::choose(all_hlas$counts[act_hla], crg_patient_hla_count:crg_patient_count)*base::choose(num_patients-all_hlas$counts[act_hla], crg_patient_count-crg_patient_hla_count:crg_patient_count)/base::choose(num_patients, crg_patient_count))
          if(act_Prob<score_hla) score_hla <- act_Prob
          if(act_Prob < hla_cutoff){
            if(lowest_hla == ""){
              lowest_hla <- base::paste(all_hlas$HLA[act_hla],
                                        " [(", crg_patient_hla_count, "/", crg_patient_count, ") vs (",
                                        all_hlas$counts[act_hla], "/", num_patients,
                                        ") = ",
                                        base::formatC(act_Prob, digits = 1, format = "e"),
                                        "]",
                                        sep = "" )
            } else {
              lowest_hla <- base::paste(lowest_hla, ", ", all_hlas$HLA[act_hla],
                                        " [(", crg_patient_hla_count, "/", crg_patient_count, ") vs (",
                                        all_hlas$counts[act_hla], "/", num_patients,
                                        ") = ",
                                        base::formatC(act_Prob, digits = 1, format = "e"),
                                        "]",
                                        sep = "" )
            }
          }
        }
        
      }
      all_scores <- base::c(all_scores, score_hla)
    }
    
    ### Total score
    if(gliph_version == 1) score_final <- base::prod(all_scores)*0.001*64 else if(gliph_version == 2) score_final <- base::prod(all_scores)
    
    ### Output
    all_scores <- c(score_final, all_scores)
    all_scores <- base::formatC(all_scores, digits = 1, format = "e")
    output <- base::c(base::names(cluster_list)[actCluster], all_scores)
    if(hla_info == TRUE && patient_info == TRUE){
      output <- base::c(output, lowest_hla)
    }
    output
  }
  
  doParallel::stopImplicitCluster()
  
  res <- base::data.frame(base::matrix(base::unlist(res), ncol = 2+base::length(score_names), byrow = TRUE))
  base::colnames(res) <- base::c("leader.tag", "total.score", score_names)
  
  for(i in base::c("total.score", score_names)) if(i != "lowest.hlas") res[,i] <- base::as.numeric(res[,i])

  # set all theoretical numeric values (actually characters) to numeric values
  if(base::is.data.frame(res)){
    for(i in base::seq_len(base::ncol(res))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(res[,i])))) == FALSE) res[,i] <- base::as.numeric(res[,i])
    }
  }
  
  # Closing time!
  base::return(res[,-1])
}

