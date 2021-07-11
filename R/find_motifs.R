#'  Finding motifs in sequences
#'
#' @description Finding motifs in sequences. This function searches and simultaneously counts continuous and discontinuous motifs in a sequence vector.
#' It is used in the \code{turbo_gliph} function to identify local similarities.
#' @param seqs character vector. This vector must contain the sequences whose motifs are to be identified and quantified.
#' @param q accepts a numeric vector of motif lengths
#' you want to find. By default it searches for
#' motifs of size 2, 3 and 4.
#' @param kmer_mindepth numeric. By default 3. Minimum observations of kmer for it to be evaluated. This is
#' the minimum number of times a kmer should be observed in
#' the sample set in order for it to be considered for being returned.
#' @param discontinuous logical. By default \code{FALSE}. Determines whether discontinuous motifs are to be considered.
#'
#' @return \code{find_motifs} returns a data frame with two columns.
#' The first column contains the motifs and the second column the frequency of the motifs.
#'
#' @examples
#' utils::data("gliph_input_data")
#' sample_seqs <- base::as.character(gliph_input_data$CDR3b)
#' res <- find_motifs(seqs = sample_seqs)
#'
#' @export
find_motifs <- function(seqs, q = 2:4,kmer_mindepth = NULL, discontinuous = FALSE){
  all_kmers <- base::c()
  seqs <- base::as.character(seqs)
  all_q <- q
  if(discontinuous == TRUE) all_q <- base::sort(base::unique(base::c(all_q, q+1)))
  
  ### Iterate for all kmer lengths
  for(i in all_q){
    ## Get all continous motifs with length i
    cont_kmer <- stringdist::qgrams(seqs,q = i)
    cont_kmer <- base::as.data.frame(base::t(cont_kmer))
    cont_kmer$motif <- base::rownames(cont_kmer)

    ## Exclude motifs with less counts than kmer_mindepth
    if(!base::is.null(kmer_mindepth)) cont_kmer <- cont_kmer[base::which(cont_kmer$V1 >= kmer_mindepth),] #dplyr::filter(cont_kmer,(V1 >= kmer_mindepth))
    cont_kmer <- cont_kmer[,2:1]
    
    ## Save motifs and corresponding counts
    if(i %in% q) all_kmers <- base::rbind(all_kmers, cont_kmer)

    ## Get all discontinuous motifs with (i-1) fixed positions with i in q
    disc_kmer <- base::c()
    if(discontinuous == TRUE && i %in% (q+1)) {
      ## Find all discontinuous motifs with (i-1) fixed positions and variable position j
      for(j in 2:(i-1)) {
        ## Replace position j by a dot in all continuous motifs  with length i
        disc_kmer <- cont_kmer
        base::substr(disc_kmer$motif, j, j) <- "."
        disc_kmer <- base::as.data.frame(base::table(base::rep(disc_kmer$motif, disc_kmer$V1)))
        base::colnames(disc_kmer) <- base::c("motif", "V1")

        ## Exclude motifs with less counts than kmer_mindepth
        if(!base::is.null(kmer_mindepth)) disc_kmer <- disc_kmer[base::which(disc_kmer$V1 >= kmer_mindepth),] #dplyr::filter(disc_kmer,(V1 >= kmer_mindepth))
        
        ## Save motifs and corresponding counts
        all_kmers <- base::rbind(all_kmers, disc_kmer)
      }
    }
  }

  ## Closing time!
  base::return(all_kmers)
}
