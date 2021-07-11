#' Generating a random subsample
#'
#' @description Generating a random subsample. This function is used by the \code{turbo_gliph} function to generate random subsamples of naive reference
#' sequences with a similar size as the sample size for random repeat sampling.
#' If specified the function tries to maintain the distribution of cdr3 lengths and/or V-gene usage of the whole sample in the subsample.
#'
#' @param cdr3_len_stratify logical. By default \code{FALSE}.
#' Specifies whether the distribution of the cdr3 lengths in the sample should be retained during repeat random sampling.
#' @param vgene_stratify logical. By default \code{FALSE}.
#' Specifies whether the distribution of V-genes in the sample should be retained during repeat random sampling.
#' @param refseqs_motif_region character vector. Contains the motif regions of reference sequences.
#' @param motif_region character vector. Contains the motif regions of sample sequences.
#' @param motif_lengths_list list. Required if \code{cdr3_len_stratify = TRUE}.
#' The elements are named after the different cdr3 lengths in the \code{motif_region} vector
#' and contain the frequency of occurrence of the corresponding cdr3 length in the \code{motif_region} vector.
#' @param ref_motif_lengths_id_list list. Required if \code{cdr3_len_stratify = TRUE}.
#' The elements are named after the different cdr3 lengths in the \code{motif_region} vector
#' and contain the indices of sequences in the \code{refseqs_motif_region} vector with the corresponding cdr3 length.
#' @param motif_region_vgenes_list list. Required if \code{vgene_stratify = TRUE}.
#' The elements are named after the different V-genes of the sequences in the \code{motif_region} vector
#' and contain the frequency of occurrence of the corresponding V-genes of the sequences in the \code{motif_region} vector.
#' @param ref_motif_vgenes_id_list list. Required if \code{vgene_stratify = TRUE}.
#' The elements are named after the different V-genes of the sequences in the \code{motif_region} vector
#' and contain the indices of sequences in the \code{refseqs_motif_region} vector with the corresponding V-gene.
#' @param lengths_vgenes_list list. Required if \code{cdr3_len_stratify = TRUE} and \code{vgene_stratify = TRUE}.
#' The elements are lists itself and are named after the different cdr3 lengths in the \code{motif_region} vector.
#' The elements of any list are named after the different V-genes of the sequences in the \code{motif_region} vector
#' and contain the frequency of simultaneous occurrence of the corresponding cdr3 length and V-gene of the sequences in the \code{motif_region} vector.
#' @param ref_lengths_vgenes_list list. Required if \code{cdr3_len_stratify = TRUE} and \code{vgene_stratify = TRUE}.
#' The elements are lists itself and are named after the different cdr3 lengths in the \code{motif_region} vector.
#' The elements of any list are named after the different V-genes of the sequences in the \code{motif_region} vector
#' and contain the frequency of simultaneous occurrence of the corresponding cdr3 length and V-gene of the sequences in the \code{refseqs_motif_region} vector.
#'
#' @return \code{getRandomSubsample} returns a character vector containing a subsample of \code{refseqs_motif_region} with the same size as \code{motif_region}
#'
# @export

getRandomSubsample <- function(cdr3_len_stratify = FALSE,
                               vgene_stratify = FALSE,
                               refseqs_motif_region,
                               motif_region,
                               motif_lengths_list,
                               ref_motif_lengths_id_list,
                               motif_region_vgenes_list,
                               ref_motif_vgenes_id_list,
                               ref_lengths_vgenes_list,
                               lengths_vgenes_list){
  random_subsample <- base::c()
  
  ### Return an unbiased reference subsample
  if(cdr3_len_stratify == FALSE && vgene_stratify == FALSE){
    random_subsample <- refseqs_motif_region[base::sample(x = base::seq_along(refseqs_motif_region),size = base::length(motif_region),replace = FALSE)]
  }
  
  ### Return a reference subsample with biased CDR3b length distribution according to the sample
  if(cdr3_len_stratify == TRUE && vgene_stratify == FALSE){
    
    ## if there are more seqs with specific cdr3 length in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs with different cdr3 length
    random_length <- 0
    for(cdr3_length in base::names(motif_lengths_list)){
      if(base::length(ref_motif_lengths_id_list[[cdr3_length]]) < motif_lengths_list[[cdr3_length]]){
        random_length <- random_length + motif_lengths_list[[cdr3_length]] - base::length(ref_motif_lengths_id_list[[cdr3_length]])
        random_subsample <- base::c(random_subsample, ref_motif_lengths_id_list[[cdr3_length]])
      }else {
        random_subsample <- base::c(random_subsample, base::sample(x = ref_motif_lengths_id_list[[cdr3_length]],size = motif_lengths_list[[cdr3_length]],replace = FALSE))
      }
    }
    if(random_length > 0){
      random_subsample <- base::c(random_subsample, base::sample(x = (base::seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }
  
  ### Return a reference subsample with biased V gene distribution according to the sample
  if(vgene_stratify == TRUE && cdr3_len_stratify == FALSE){
    
    ## if there are more seqs with specific v gene in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs with different vgenes
    random_length <- 0
    for(act_vgene in base::names(motif_region_vgenes_list)){
      if(base::length(ref_motif_vgenes_id_list[[act_vgene]]) < motif_region_vgenes_list[[act_vgene]]){
        random_length <- random_length + motif_region_vgenes_list[[act_vgene]] - base::length(ref_motif_vgenes_id_list[[act_vgene]])
        random_subsample <- base::c(random_subsample, ref_motif_vgenes_id_list[[act_vgene]])
      }else {
        random_subsample <- base::c(random_subsample, base::sample(x = ref_motif_vgenes_id_list[[act_vgene]],size = motif_region_vgenes_list[[act_vgene]],replace = FALSE))
      }
    }
    if(random_length > 0){
      random_subsample <- base::c(random_subsample, base::sample(x = (base::seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }
  
  ### Return a reference subsample with biased V gene and CDR3b length distribution according to the sample
  if(vgene_stratify == TRUE && cdr3_len_stratify == TRUE){
    
    ## if there are more seqs with specific v gene or CDR3b length in sample than in reference database, surplus number of seqs will be taken randomly from remaining seqs
    random_length <- 0
    for(cdr3_length in base::names(motif_lengths_list)){
      for(act_vgene in base::names(motif_region_vgenes_list)){
        if(base::length(ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]]) < lengths_vgenes_list[[cdr3_length]][[act_vgene]]){
          random_length <- random_length + lengths_vgenes_list[[cdr3_length]][[act_vgene]] - base::length(ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]])
          random_subsample <- base::c(random_subsample, ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]])
        }else {
          random_subsample <- base::c(random_subsample, base::sample(x = ref_lengths_vgenes_list[[cdr3_length]][[act_vgene]],size = lengths_vgenes_list[[cdr3_length]][[act_vgene]],replace = FALSE))
        }
      }
    }

    if(random_length > 0){
      random_subsample <- base::c(random_subsample, base::sample(x = (base::seq_along(refseqs_motif_region))[-random_subsample],size = random_length,replace = FALSE))
    }
    random_subsample <- refseqs_motif_region[random_subsample]
  }

  ## Closing time!
  base::return(random_subsample)
}
