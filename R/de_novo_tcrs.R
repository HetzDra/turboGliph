#' De novo generation of cdr3 sequences based on GLIPH or GLIPH2
#'
#' @description De novo generation of cdr3 sequences based on GLIPH or GLIPH2. Based on the position-specific abundance of amino acids in the CDR3
#' region of the sequences of a GLIPH or GLIPH2 cluster, artificial sequences are simulated as established in Glanville et al.
#' @param convergence_group_tag character. Tag of the convergence group that shall be used for prediction.
#' @param result_folder character. By default \code{""}. Path to the folder in which the output files of the clustering are stored and the output of this method will be stored.
#' If the value is \code{""} the results are not saved and the output list of the function \code{turboGliph} or \code{gliph2} must be entered under the parameter \code{turboGliph_output}.
#' @param clustering_output list. By default \code{NULL}. If this parameter is specified, the clustering results are loaded directly from the list and not from the files in the \code{result_folder}
#' If the value of \code{result_folder} is \code{""}, the output list of the function \code{turboGliph} or \code{gliph2} must be entered here.
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
#' @param accept_sequences_with_C_F_start_end logical. This logical flag
#' if \code{TRUE}, by default, only accepts sequences with
#' amino-acid C at the start position and amino-acid F at the
#' end position.
#' @param normalization logical. By default \code{FALSE}. If \code{TRUE} the calculated scores are normalized to a reference database
#' and the probability that a reference sequence has a score greater than or equal to the sample sequence is returned.
#' If V gene information is available, only sequences with identical V gene are compared.
#' @param sims numeric. By default 1,000,000. Value of how many de novo cdr3 sequences shall be created.
#' @param num_tops numeric. By default 1000. The \code{num_tops} best scoring de novo created cdr3 sequences are returned
#' @param min_length Numeric value determining the number of N-terminal positions used for scoring. By default it is set to 10.
#' @param make_figure Logical value whether a graph of the \code{num_tops} best scoring de novo created cdr3 sequences in dependence of the rank shall be displayed.
#' @param n_cores numeric. Number of cores to use, by default 1. In case of \code{NULL} it will be set to number of cores in your machine minus 2.
#'
#' @return This function produces one file in the \code{result_folder}, if specified, named \code{convergence_group_tag} followed by \code{_de_novo.txt)}
#' containing the \code{num_tops} best scoring generated sequences and their corresponding scores.
#' A list containing this file and additional information will also be returned as follows:
#'
#' @return $de_novo_sequences
#' A data frame containing the \code{num_tops} best scoring generated sequences and their corresponding scores.
#'
#' @return $sample_sequences_scores
#' A data frame containing the sequences of the used convergence group and their corresponding scores.
#'
#' @return $cdr3_length_probability
#' A data frame with any considered cdr3 length and the probability of occurrence in the convergence group. The distribution of the cdr3 length of all
#' generated sequences resembles this distribution.
#'
#' @return $PWM_Scoring
#' A data frame containing the positional weight matrix used for scoring. The columns represent the different amino acids and the rows represent the
#' position relative to the N-terminus.
#'
#' @return $PWM_Prediction
#' A list of data frames containing the positional weight matrix for any considered cdr3 length used for generation of new sequences.
#' The columns represent the different amino acids and the rows represent the position relative to the N-terminus.
#'
#' @examples
#' utils::data("gliph_input_data")
#' res <- turbo_gliph(cdr3_sequences = gliph_input_data[base::seq_len(200),],
#' sim_depth = 100,
#' n_cores = 1)
#'
#' new_seqs <- de_novo_TCRs(convergence_group_tag = res$cluster_properties$tag[1],
#' clustering_output = res,
#' sims = 10000,
#' make_figure = TRUE,
#' n_cores = 1)
#'
#' @references Glanville, Jacob, et al.
#' "Identifying specificity groups in the T cell receptor repertoire." Nature 547.7661 (2017): 94.
#' @references https://github.com/immunoengineer/gliph
#' @import foreach
#' @export
de_novo_TCRs <- function(convergence_group_tag,
                         result_folder = "",
                         clustering_output = NULL,
                         refdb_beta = "gliph_reference",
                         normalization = FALSE,
                         accept_sequences_with_C_F_start_end = TRUE,
                         sims = 100000,
                         num_tops = 1000,
                         min_length = 10,
                         make_figure = FALSE,
                         n_cores = 1){
  t1 <- base::Sys.time()

  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################
  
  ### convergence_group_tag
  if(!base::is.character(convergence_group_tag)) base::stop("convergence_group_tag has to be a character object")
  if(base::length(convergence_group_tag) > 1) base::stop("convergence_group_tag has to be a single character string")
  
  ### result_folder and clustering_output
  if(!base::is.character(result_folder)) base::stop("result_folder has to be a character object")
  if(base::length(result_folder) > 1) base::stop("result_folder has to be a single path")
  save_results <- FALSE
  if(result_folder != ""){
    if(base::substr(result_folder,base::nchar(result_folder),base::nchar(result_folder)) != "/") result_folder <- base::paste0(result_folder,"/")
    if (!base::dir.exists(result_folder)) base::dir.create(result_folder)
    save_results <- TRUE
  } else {
    if(!base::is.list(clustering_output)) base::stop("If 'result_folder' = \"\" the output list of clustering must be given by 'clustering_output'.")
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
  
  ### accept_sequences_with_C_F_start_end
  if(!base::is.logical(accept_sequences_with_C_F_start_end)) base::stop("accept_sequences_with_C_F_start_end has to be logical")
  
  ### normalization
  if(!base::is.logical(normalization)) base::stop("normalization has to be logical")
  
  ### sims
  if(!base::is.numeric(sims)) base::stop("sims has to be numeric")
  if(base::length(sims) > 1) base::stop("sims has to be a single number")
  if(sims < 1) base::stop("sims must be at least 1")
  sims <- base::round(sims)
  
  ### num_tops
  if(!base::is.numeric(num_tops)) base::stop("num_tops has to be numeric")
  if(base::length(num_tops) > 1) base::stop("num_tops has to be a single number")
  if(num_tops < 1) base::stop("num_tops must be at least 1")
  num_tops <- base::round(num_tops)
  
  ### min_length
  if(!base::is.numeric(min_length)) base::stop("min_length has to be numeric")
  if(base::length(min_length) > 1) base::stop("min_length has to be a single number")
  if(min_length < 1) base::stop("min_length must be at least 1")
  min_length <- base::round(min_length)
  
  ### make_figure
  if(!base::is.logical(make_figure)) base::stop("make_figure has to be logical")
  
  ### n_cores
  if(base::is.null(n_cores)) n_cores <- base::max(1, parallel::detectCores()-2) else {
    if(!base::is.numeric(n_cores)) base::stop("n_cores has to be numeric")
    if(base::length(n_cores) > 1) base::stop("n_cores has to be a single number")
    if(n_cores < 1) base::stop("n_cores must be at least 1")
    
    n_cores <- base::min(n_cores, parallel::detectCores()-2)
  }

  # Amino acids one letter code
  aa_code <- LETTERS[-base::c(2, 10, 15, 21, 24, 26)]

  #################################################################
  ##                      Input preparation                      ##
  #################################################################
  
  ### load convergence groups from file or from input and only use specified convergence group
  base::cat("Loading convergence group with tag", convergence_group_tag, ".\n")
  if(base::is.null(clustering_output)){
    clustering_output <- load_gliph_output(result_folder = result_folder)
    crg <- clustering_output$cluster_list
  } else {
    clustering_output <- clustering_output
    crg <- clustering_output$cluster_list
  }
  if(!(convergence_group_tag %in% base::names(crg))) base::stop("Could not find convergence group with tag ",convergence_group_tag, " in clustering results stored in", result_folder, ".")
  crg <- crg[[convergence_group_tag]]
  all_crg_cdr3_seqs <- crg$CDR3b

  ### Filter sequences for a minimun sequenc length
  excluded <- base::which(base::nchar(all_crg_cdr3_seqs) < min_length)
  if(base::length(excluded) > 0){
    all_crg_cdr3_seqs <- all_crg_cdr3_seqs[-excluded]
    if(base::length(all_crg_cdr3_seqs) > 0){
      base::cat(paste0("Warning: ", base::length(excluded), " sequences of the convergence group were excluded from the further procedure due to falling below a minimum length of ", min_length, "."))
    } else {
      base::stop("No sequences of the convergence group are of minimum length of ", min_length, ". For further procedure, adjust the parameter 'min_length'")
    }
  }
  
  #################################################################
  ##                          Main part                          ##
  #################################################################
  
  ### Initiate parallelization
  doParallel::registerDoParallel(n_cores)
  
  # get all sequences in cluster
  crg_cdr3_seqs <- all_crg_cdr3_seqs
  
  ### load non-redundant reference repertoire and V-genes of sequences in convergence group
  v_genes <- base::c()
  ref_vgenes <- base::c()
  refseqs <- base::c()
  v_gene_norm <- normalization
  if(normalization == TRUE){
    
    if(base::is.data.frame(refdb_beta)) {
      refseqs <- refdb_beta
      refseqs[] <- base::lapply(refseqs, base::as.character)
      
      if(base::ncol(refseqs) > 1){
        base::cat("Notification: First column of reference database is considered as cdr3 sequences.\n")
      }
      if(base::ncol(refseqs) > 1 && v_gene_norm == TRUE){
        base::cat("Notification: Second column of reference database is considered as V-gene information.\n")
      } else if(v_gene_norm == TRUE){
        base::cat("Warning: Beta sequence reference database is missing column containing V-genes. Without V-gene information normalization may be inaccurate.\n")
        v_gene_norm <- FALSE
      }
      if(base::ncol(refseqs) == 1) refseqs <- base::cbind(refseqs, base::rep("", base::nrow(refseqs)))
      
      refseqs <- refseqs[, base::c(1,2)]
      base::colnames(refseqs) <- base::c("CDR3b", "TRBV")
      refseqs <- base::unique(refseqs)
      if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[base::grep(pattern = "^C.*F$",x = refseqs$CDR3b,perl = TRUE),]
      refseqs <- refseqs[base::which(base::nchar(refseqs$CDR3b) >= min_length),]
      refseqs <- refseqs[base::grep("^[ACDEFGHIKLMNOPQRSTUVWY]*$", refseqs$CDR3b),]
      
      if(base::nrow(refseqs) == 0){
        normalization <- FALSE
        v_gene_norm <- FALSE
        base::cat(base::paste0("Warning: No reference sequences with a minimum length of ", min_length, " given. Normalization therefore not possible. Adjust min_length to enable normalization.\n"))
      } else {
        ref_vgenes <- base::as.character(refseqs$TRBV)
        refseqs <- base::as.character(refseqs$CDR3b)
      }
    } else {
      reference_list <- NULL
      utils::data("reference_list",envir = base::environment(), package = "turboGliph")
      refseqs <- base::as.data.frame(reference_list[[refdb_beta]]$refseqs)
      refseqs[] <- base::lapply(refseqs, base::as.character)
      
      if(base::ncol(refseqs) > 1){
        base::cat("Notification: First column of reference database is considered as cdr3 sequences.\n")
      }
      if(base::ncol(refseqs) > 1 && v_gene_norm == TRUE){
        base::cat("Notification: Second column of reference database is considered as V-gene information.\n")
      } else if(v_gene_norm == TRUE){
        base::cat("Warning: Beta sequence reference database is missing column containing V-genes. Without V-gene information normalization may be inaccurate.\n")
        v_gene_norm <- FALSE
      }
      if(base::ncol(refseqs) == 1) refseqs <- base::cbind(refseqs, base::rep("", base::nrow(refseqs)))
      
      refseqs <- refseqs[, base::c(1,2)]
      base::colnames(refseqs) <- base::c("CDR3b", "TRBV")
      refseqs <- base::unique(refseqs)
      if(accept_sequences_with_C_F_start_end) refseqs <- refseqs[base::grep(pattern = "^C.*F$",x = refseqs$CDR3b,perl = TRUE),]
      refseqs <- refseqs[base::which(base::nchar(refseqs$CDR3b) >= min_length),]
      refseqs <- refseqs[base::grep("^[ACDEFGHIKLMNPQRSTVWY]*$", refseqs$CDR3b),]
      
      if(base::nrow(refseqs) == 0){
        normalization <- FALSE
        v_gene_norm <- FALSE
        base::cat(base::paste0("Warning: No reference sequences with a minimum length of ", min_length, " given. Normalization therefore not possible. Adjust min_length to enable normalization.\n"))
      } else {
        ref_vgenes <- base::as.character(refseqs$TRBV)
        refseqs <- base::as.character(refseqs$CDR3b)
      }
    }
    
    # load V genes of cluster members
    if("TRBV" %in% base::colnames(crg) && v_gene_norm == TRUE){
      v_genes <- crg$TRBV
    } else if(v_gene_norm == TRUE){
      v_gene_norm <- FALSE
      base::cat("Warning: Without V-gene information of sample sequences normalization may be inaccurate.\n")
    } else base::cat("Warning: Without V-gene restriction normalization may be inaccurate.\n")
  }
  
  ### Initialization
  crg_num_seqs <- base::length(crg_cdr3_seqs)
  crg_cdr3_scores <- base::rep(1, crg_num_seqs)
  crg_cdr3_norm_scores <- base::rep(0, crg_num_seqs)
  crg_cdr3_lens <- base::nchar(crg_cdr3_seqs)
  max_cdr3_length <- base::max(crg_cdr3_lens)
  
  ### Build Positional Weight Matrix (PWM) of convergence group for min_length N-terminal positions
  base::cat("Calculating positional weight matrix of convergence group.\n")
  
  # initialize the matrix
  crg_pwm_scoring <- base::as.data.frame(base::matrix(base::rep(0, min_length*base::length(aa_code)), ncol = base::length(aa_code)))
  base::colnames(crg_pwm_scoring) <- aa_code
  
  # for every position determine the amino acid frequency
  for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
    aa_freqs <- base::rep(0, base::length(aa_code))
    act_letters <- base::substr(crg_cdr3_seqs, i, i)
    for(j in base::seq_along(aa_freqs)){
      aa_freqs[j] <- base::sum(act_letters == aa_code[j])
    }
    
    # Pseudocounts of 0.5% per aa
    zeroFreqs <- base::which(aa_freqs == 0)
    aa_freqs <- aa_freqs/base::sum(aa_freqs)*(1-base::length(zeroFreqs)*0.005)
    aa_freqs[zeroFreqs] <- 0.005
    crg_pwm_scoring[i,] <- aa_freqs
  }
  
  ### Calculate scores of convergence group members, only use first min_length N-terminal positions
  # score = product of amino acid frequencies from min_length N-terminal positions
  base::cat("Calculating scores of convergence group members.\n")
  
  # is parallelization necessary?
  if(crg_num_seqs > 10000){
    # distribute sequences equally to all cores
    distribute <- base::lapply(base::seq_len(n_cores), function(x){base::return(((x-1)*base::floor(base::length(crg_cdr3_seqs)/n_cores)+1):(x*base::floor(base::length(crg_cdr3_seqs)/n_cores)))})
    
    # calculate the score
    crg_cdr3_scores <- foreach::foreach(i = base::seq_along(distribute), .combine = c) %dopar% {
      temp_scores <- base::rep(1, base::length(distribute[[i]]))
      temp_seqs <- crg_cdr3_seqs[distribute[[i]]]
      for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
        temp_scores <- temp_scores*base::unlist(crg_pwm_scoring[i, base::substr(temp_seqs, i, i)])
      }
      
      base::return(temp_scores)
    }
  } else {
    # calculate the score
    for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
      crg_cdr3_scores <- crg_cdr3_scores*base::unlist(crg_pwm_scoring[i, base::substr(crg_cdr3_seqs, i, i)])
    }
  }
  
  ### Normalize scores of convergence group members, only use first 10 N-terminal positions
  # calculate the probability that a score at least this high occurs in the reference database
  if(normalization == TRUE){
    base::cat("Normalizing scores of convergence group members.\n")
    
    # Calculate scores of reference database (analogue to sample sequences)
    refseq_scores <- base::rep(1, base::length(refseqs))
    if(base::length(refseqs) > 10000){
      distribute <- base::lapply(base::seq_len(n_cores), function(x){base::return(((x-1)*base::floor(base::length(refseqs)/n_cores)+1):(x*base::floor(base::length(refseqs)/n_cores)))})
      
      refseq_scores <- foreach::foreach(i = base::seq_along(distribute), .combine = c) %dopar% {
        temp_scores <- base::rep(1, base::length(distribute[[i]]))
        temp_seqs <- refseqs[distribute[[i]]]
        for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
          temp_scores <- temp_scores*base::unlist(crg_pwm_scoring[i, base::substr(temp_seqs, i, i)])
        }
        
        base::return(temp_scores)
      }
    } else {
      for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
        refseq_scores <- refseq_scores*base::unlist(crg_pwm_scoring[i, substr(refseqs, i, i)])
      }
    }
    
    # Calculate normalized scores, if requested restrict score comparison to identical V genes
    crg_cdr3_norm_scores <- foreach::foreach(i = base::seq_len(crg_num_seqs), .combine = c) %dopar% {
      v_gene_penalty <- base::rep(0, base::length(refseqs))
      if(v_gene_norm == TRUE){
        v_gene_penalty[ref_vgenes != v_genes[i]] <- -2
      }
      
      base::return(base::sum((refseq_scores + v_gene_penalty) >= crg_cdr3_scores[i])/base::length(refseqs))
    }
  }
  
  ### Create global PWM of convergence group
  base::cat("Calculating positional weight matrix for de novo CDR3b-sequence prediction.\n")
  
  # create a PWM for every CDR3b length in the cluster
  crg_pwm_predicting_list <- base::list()
  
  # get the probability of this CDR3b length to occur
  crg_len_prob <- base::data.frame(length = base::unique(crg_cdr3_lens),probability = base::rep(0, base::length(base::unique(crg_cdr3_lens))))
  for(i in base::seq_len(base::nrow(crg_len_prob))){
    crg_len_prob$probability[i] <- base::sum(crg_cdr3_lens == crg_len_prob$length[i])/crg_num_seqs
  }
  
  # receive the positional dependent amino acid frequencies for every CDR3b length
  for(len in base::unique(crg_cdr3_lens)){
    crg_pwm_predicting <- base::as.data.frame(base::matrix(base::rep(0,len*base::length(aa_code)), ncol = base::length(aa_code)))
    base::colnames(crg_pwm_predicting) <- aa_code
    
    seqs <- crg_cdr3_seqs[crg_cdr3_lens == len]
    for(i in base::seq_len(len)){
      aa_freqs <- base::rep(0, base::length(aa_code))
      act_letters <- base::substr(seqs, i, i)
      for(j in base::seq_along(aa_freqs)){
        aa_freqs[j] <- base::sum(act_letters == aa_code[j])
      }
      
      # Pseudocounts of 0.5% per aa
      zeroFreqs <- base::which(aa_freqs == 0)
      aa_freqs <- aa_freqs/base::sum(aa_freqs)*(1-base::length(zeroFreqs)*0.005)
      aa_freqs[zeroFreqs] <- 0.005
      crg_pwm_predicting[i,] <- aa_freqs
    }
    
    # if de novo sequences should only start with C and end with F, correct the PWM for this positions
    if(accept_sequences_with_C_F_start_end == TRUE){
      crg_pwm_predicting[1, ] <- base::rep(0, base::ncol(crg_pwm_predicting))
      crg_pwm_predicting[nrow(crg_pwm_predicting), ] <- base::rep(0, base::ncol(crg_pwm_predicting))
      crg_pwm_predicting$'C'[1] <- 1
      crg_pwm_predicting$'F'[base::nrow(crg_pwm_predicting)] <- 1
    } 
    
    crg_pwm_predicting_list[[base::paste("Length", len)]] <- crg_pwm_predicting
  }
  
  
  ### Create a number of sims de novo TCR sequences
  base::cat("Creating", sims, "de novo sequences.\n")
  
  # randomly select the length of the sequences
  de_novo_lens <- base::sample(x = base::c(crg_len_prob$length, 0), size = sims, prob = base::c(crg_len_prob$probability,0), replace = TRUE)
  de_novo_seqs <- base::rep("", sims)
  
  # based on the PWM randomly create new sequences
  for(len in crg_len_prob$length){
    inds <- base::which(de_novo_lens == len)
    for(i in base::seq_len(len)){
      rands <- aa_code[base::sample.int(n = base::length(aa_code), size = base::length(inds), prob = crg_pwm_predicting_list[[base::paste("Length", len)]][i,], replace = TRUE)]
      de_novo_seqs[inds] <- base::paste(de_novo_seqs[inds], rands, sep = "")
    }
  }
  de_novo_seqs <- base::unique(de_novo_seqs)
  if(accept_sequences_with_C_F_start_end) de_novo_seqs <- base::grep(pattern = "^C.*F$",x = de_novo_seqs ,perl = TRUE,value = TRUE)
  de_novo_lens <- base::nchar(de_novo_seqs)
  
  # Score de_novo seqs as performed above
  base::cat("Calculating scores of de novo sequences.\n")
  de_novo_seqs_scores <- base::rep(1, base::length(de_novo_seqs))
  
  if(base::length(de_novo_seqs) > 10000){
    distribute <- base::lapply(base::seq_len(n_cores), function(x){base::return(((x-1)*base::floor(base::length(de_novo_seqs)/n_cores)+1):(x*base::floor(base::length(de_novo_seqs)/n_cores)))})
    
    de_novo_seqs_scores <- foreach::foreach(i = base::seq_along(distribute), .combine = c) %dopar% {
      temp_scores <- base::rep(1, base::length(distribute[[i]]))
      temp_seqs <- de_novo_seqs[distribute[[i]]]
      for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
        temp_scores <- temp_scores*base::unlist(crg_pwm_scoring[i, base::substr(temp_seqs, i, i)])
      }
      
      base::return(temp_scores)
    }
  } else {
    for(i in base::seq_len(base::nrow(crg_pwm_scoring))){
      de_novo_seqs_scores <- de_novo_seqs_scores*base::unlist(crg_pwm_scoring[i, base::substr(de_novo_seqs, i, i)])
    }
  }
  de_novo_seqs_scores <- base::as.numeric(base::formatC(de_novo_seqs_scores, digits = 1, format = "e"))
  
  # Normalize de_novo seqs scores as performed above
  if(normalization == TRUE){
    base::cat("Normalizing scores of de novo sequences.\n")
    
    de_novo_norm_scores <- foreach::foreach(i = base::seq_along(de_novo_seqs), .combine = c) %dopar% {
      return(base::sum(refseq_scores >= de_novo_seqs_scores[i])/base::length(refseqs))
    }
    de_novo_norm_scores <- base::as.numeric(base::formatC(de_novo_norm_scores, digits = 1, format = "e"))
  }
  
  # sort sequences based on score (normalized scores are prioritized)
  if(normalization == TRUE){
    order_ids <- base::order(de_novo_norm_scores, decreasing=FALSE)
    de_novo_seqs <- de_novo_seqs[order_ids]
    de_novo_lens <- de_novo_lens[order_ids]
    de_novo_seqs_scores <- de_novo_seqs_scores[order_ids]
    de_novo_norm_scores <- de_novo_norm_scores[order_ids]
  } else{
    order_ids <- base::order(de_novo_seqs_scores, decreasing=TRUE)
    de_novo_seqs <- de_novo_seqs[order_ids]
    de_novo_lens <- de_novo_lens[order_ids]
    de_novo_seqs_scores <- de_novo_seqs_scores[order_ids]
  }
  
  # get num_tops heighest scored sequences
  if(base::length(de_novo_seqs) < num_tops) num_tops <- base::length(de_novo_seqs)
  de_novo_seqs <- de_novo_seqs[base::seq_len(num_tops)]
  de_novo_lens <- de_novo_lens[base::seq_len(num_tops)]
  de_novo_seqs_scores <- de_novo_seqs_scores[base::seq_len(num_tops)]
  if(normalization == TRUE){
    de_novo_norm_scores <- de_novo_norm_scores[base::seq_len(num_tops)]
    de_novo <- base::data.frame(length = de_novo_lens, seqs = de_novo_seqs, norm_score = de_novo_norm_scores, score = de_novo_seqs_scores)
  } else {
    de_novo <- base::data.frame(length = de_novo_lens, seqs = de_novo_seqs, score = de_novo_seqs_scores)
  }
  
  ### sort cluster sequences based on their score
  connected_inds <- base::seq_along(crg_cdr3_seqs)
  if(normalization == FALSE){
    crg_cdr3_seqs <- crg_cdr3_seqs[base::order(crg_cdr3_scores, decreasing = TRUE)]
    connected_inds <- connected_inds[base::order(crg_cdr3_scores, decreasing = TRUE)]
    crg_cdr3_scores <- crg_cdr3_scores[base::order(crg_cdr3_scores, decreasing = TRUE)]
  } else{
    crg_cdr3_seqs <- crg_cdr3_seqs[base::order(crg_cdr3_norm_scores, decreasing = FALSE)]
    connected_inds <- connected_inds[base::order(crg_cdr3_norm_scores, decreasing = FALSE)]
    crg_cdr3_scores <- crg_cdr3_scores[base::order(crg_cdr3_norm_scores, decreasing = FALSE)]
    crg_cdr3_norm_scores <- crg_cdr3_norm_scores[base::order(crg_cdr3_norm_scores, decreasing = FALSE)]
  }
  
  ### set all theoretical numeric values (actually characters) to numeric values
  if(base::is.data.frame(de_novo)){
    for(i in base::seq_len(base::ncol(de_novo))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(de_novo[,i])))) == FALSE) de_novo[,i] <- base::as.numeric(de_novo[,i])
    }
  }
  if(base::is.data.frame(crg_cdr3_scores)){
    for(i in base::seq_len(base::ncol(crg_cdr3_scores))){
      if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(crg_cdr3_scores[,i])))) == FALSE) crg_cdr3_scores[,i] <- base::as.numeric(crg_cdr3_scores[,i])
    }
  }
  if(normalization == TRUE){
    if(base::is.data.frame(crg_cdr3_norm_scores)){
      for(i in base::seq_len(base::ncol(crg_cdr3_norm_scores))){
        if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(crg_cdr3_norm_scores[,i])))) == FALSE) crg_cdr3_norm_scores[,i] <- base::as.numeric(crg_cdr3_norm_scores[,i])
      }
    }
  }
  
  ### output
  if(normalization == FALSE){
    output <- base::list(de_novo_sequences = de_novo, sample_sequences_scores = base::data.frame(seqs = all_crg_cdr3_seqs[connected_inds], scores = crg_cdr3_scores),
                         cdr3_length_probability = crg_len_prob, PWM_Scoring = crg_pwm_scoring, PWM_Prediction = crg_pwm_predicting_list)
  } else {
    output <- base::list(de_novo_sequences = de_novo, sample_sequences_scores = base::data.frame(seqs = all_crg_cdr3_seqs[connected_inds], norm_scores = crg_cdr3_norm_scores, scores = crg_cdr3_scores),
                         cdr3_length_probability = crg_len_prob, PWM_Scoring = crg_pwm_scoring, PWM_Prediction = crg_pwm_predicting_list)
  }
  
  ### save
  fname <- base::paste0(result_folder, convergence_group_tag, "_de_novo.txt")
  if(save_results == TRUE) utils::write.table(x = de_novo,file = fname,quote = FALSE,sep = "\t",row.names = FALSE, col.names = TRUE)
  if(save_results == TRUE) base::cat("Output: results are stored in ", fname, "\n")
  
  ### Print a graph with num_tops best scoring de novo sequences
  if(make_figure == TRUE){
    if(normalization == FALSE){
      graphics::plot(x = base::seq_len(num_tops), y = de_novo$score*100, xlab = base::paste("Top", num_tops, "predicted TCRs"), ylab = "TCR probability based on PWM in %", type = "p", log = "xy", col = "grey", pch = 19)
      graphics::points(x = base::seq_len(10), y = de_novo$score[base::seq_len(10)]*100, col = "red", pch = 19)
      
      # which cluster sequences are present in de novo created sequences?
      common_ids <- base::which(de_novo_seqs %in% crg_cdr3_seqs)
      graphics::points(x = common_ids, y = de_novo_seqs_scores[common_ids]*100, col = "yellow", pch = 20)
    } else{
      graphics::plot(x = base::seq_len(num_tops), y = de_novo$norm_score *100, xlab = base::paste("Top", num_tops, "predicted TCRs"), ylab = "Normalized TCR probability based on PWM in %",
                    type = "p", log = "xy", col = "grey", pch = 19)
      graphics::points(x = base::seq_len(10), y = de_novo$norm_score[base::seq_len(10)]*100, col = "red", pch = 19)
      
      # which cluster sequences are present in de novo created sequences?
      common_ids <- base::which(de_novo_seqs %in% crg_cdr3_seqs)
      graphics::points(x = common_ids, y = de_novo_norm_scores[common_ids]*100, col = "yellow", pch = 20)
    }

    graphics::legend("bottomleft", legend=base::c(base::paste0("Top ", num_tops," de novo scores"), "Top 10 de novo scores", "Convergence group members\nin de novo sequences"),
           col=base::c("grey", "red", "yellow"), pch = base::c(19,19,20),title="Legend", cex = 0.75)
  }

  t2 <- base::Sys.time()
  dt <- (t2-t1)
  base::cat("Total time = ",dt, base::units(dt),"\n")

  doParallel::stopImplicitCluster()
  
  
  ### Closing time!
  base::return(output)
}
