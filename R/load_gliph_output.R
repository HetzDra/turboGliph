#'  Loading turboGLIPH and GLIPH2 output
#'
#' @description  Loading turboGLIPH and GLIPH2 output. This function loads the output of the GLIPH functions stored in txt files.
#' @param result_folder character. Must be specified. Path to the folder in which the output files of the  gliph function are stored.
#' @return \code{load_gliph_output} returns a list equivalent to either the \code{turboGliph} or \code{gliph2} function (see documentation
#' of the functions for more information).
#'
#' @export
load_gliph_output <- function(result_folder = ""){
  
  ##################################################################
  ##                         Unit-testing                         ##
  ##################################################################
  
  ### result_folder
  if(!base::is.character(result_folder))base::stop("result_folder has to be a character object")
  if(base::length(result_folder) > 1)base::stop("result_folder has to be a single path")
  if(result_folder == "") base::stop("To load output files, the path must be specified unter the parameter 'result_folder'.")
  if(base::substr(result_folder,base::nchar(result_folder),base::nchar(result_folder)) != "/") result_folder <- base::paste0(result_folder,"/")
  base::cat(base::paste("Files are loaded from the following folder:", result_folder, "\n"))
  if(!base::dir.exists(result_folder)) base::stop("Specified path under the parameter 'result folder' does not exist.")
  
  ##################################################################
  ##                          Load files                          ##
  ##################################################################
  
  ### load parameter file: parameters
  fname <- base::paste0(result_folder, "parameter.txt")
  if(!base::file.exists(fname)) base::stop("File named 'parameters.txt' is missing in the specified folder.")
  para_df <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
  parameters <- base::lapply(base::seq_len(base::nrow(para_df)), function(x){
    if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(para_df[x,2])))) == FALSE) return(base::as.numeric(para_df[x,2])) else return(para_df[x,2])
    })
  base::names(parameters) <- para_df[,1]
  parameters$motif_length <- base::as.numeric(base::unlist(base::strsplit(parameters$motif_length, split = ",")))
  parameters$lcminove <- base::as.numeric(base::unlist(base::strsplit(parameters$lcminove, split = ",")))
  base::cat(base::paste("Loaded 'parameters' from parameter.txt", "\n"))
  
  if(!("gliph_version" %in% base::names(parameters))){
    base::cat("Output of gliph_combined function is loaded.", "\n")
    
    ### load sample log
    fname <- base::paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_log.txt")
    if(base::file.exists(fname)){
      sample_log <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'sample_log' from ",fname, "\n"))
    } else sample_log <- NULL
    
    ### load selected_motifs
    fname <- base::paste0(result_folder, "local_similarities_minp_",parameters$lcminp, "_minove_", base::paste(parameters$lcminove, collapse = "_"), "_kmer_mindepth_", parameters$lckmer_mindepth, ".txt")
    if(base::file.exists(fname)){
      selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'selected_motifs' from ",fname, "\n"))
    } else selected_motifs <- NULL
    
    ### load all_motifs
    fname <- base::paste0(result_folder,"all_motifs.txt")
    if(base::file.exists(fname)){
      all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'all_motifs' from ",fname, "\n"))
    } else all_motifs <- NULL
    
    ### build motif_enrichment
    motif_enrichment <- base::list(selected_motifs = selected_motifs, all_motifs = all_motifs)
    
    ### load selected_structs
    fname <- base::paste0(result_folder, "global_similarities_minp_",parameters$gcminp, "_minove_", base::paste(parameters$gcminove, collapse = "_"), "_kmer_mindepth_", parameters$gckmer_mindepth, ".txt")
    if(base::file.exists(fname)){
      selected_structs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'selected_structs' from ",fname, "\n"))
    } else selected_structs <- NULL
    
    ### load all_structs
    fname <- base::paste0(result_folder, "all_global_similarities.txt")
    if(base::file.exists(fname)){
      all_structs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'all_structs' from ",fname, "\n"))
    } else all_structs <- NULL
    
    ### build global_enrichment
    global_enrichment <- base::list(selected_structs = selected_structs, all_structs = all_structs)
    
    ### load connections
    fname <- base::paste0(result_folder,"clone_network.txt")
    if(base::file.exists(fname)){
      connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'connections' from ",fname, "\n"))
    } else connections <- NULL
    
    ### load cluster_properties
    fname <- base::paste0(result_folder,"convergence_groups.txt")
    if(base::file.exists(fname)){
      cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      base::cat(base::paste0("Loaded 'cluster_properties' from ",fname, "\n"))
    } else cluster_properties <- NULL
    
    ### load cluster_list
    fname <- base::paste0(result_folder, "cluster_member_details.txt")
    if(base::file.exists(fname)){
      cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
      for(i in base::seq_len(base::ncol(cluster_list))){
        if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- base::as.numeric(cluster_list[,i])
      }
      tag_names <- base::unique(cluster_list$tag)
      cluster_list <- base::lapply(tag_names, function(x){base::return(cluster_list[cluster_list$tag == x,-1])})
      base::names(cluster_list) <- tag_names
      
      base::cat(base::paste0("Loaded 'cluster_list' from ",fname, "\n"))
    } else cluster_list <- NULL
    
    output <- base::list(motif_enrichment = motif_enrichment,
                         global_enrichment = global_enrichment,
                         connections = connections,
                         cluster_properties = cluster_properties,
                         cluster_list = cluster_list,
                         parameters = parameters)
    
  } else {
    if(parameters$gliph_version == 1){
      base::cat("Output of turboGliph function is loaded.", "\n")
      
      ### load sample log
      fname <- base::paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_log.txt")
      if(base::file.exists(fname)){
        sample_log <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, row.names = 1, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'sample_log' from ",fname, "\n"))
      } else sample_log <- NULL
      
      
      ### load selected_motifs
      fname <- base::paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_minp",parameters$lcminp,"_ove", base::paste(parameters$lcminove, collapse = "_"),".txt")
      if(base::file.exists(fname)){
        selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'selected_motifs' from ",fname, "\n"))
      } else selected_motifs <- NULL
      
      ### load all_motifs
      fname <- base::paste0(result_folder,"kmer_resample_",parameters$sim_depth,"_all_motifs.txt")
      if(base::file.exists(fname)){
        all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'all_motifs' from ",fname, "\n"))
      } else all_motifs <- NULL
      
      ### build motif_enrichment
      motif_enrichment <- base::list(selected_motifs = selected_motifs, all_motifs = all_motifs)
      
      ### load connections
      fname <- base::paste0(result_folder,"clone_network.txt")
      if(base::file.exists(fname)){
        connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'connections' from ",fname, "\n"))
      } else connections <- NULL
      
      ### load cluster_properties
      fname <- base::paste0(result_folder,"convergence_groups.txt")
      if(base::file.exists(fname)){
        cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'cluster_properties' from ",fname, "\n"))
      } else cluster_properties <- NULL
      
      ### load cluster_list
      fname <- base::paste0(result_folder, "cluster_member_details.txt")
      if(base::file.exists(fname)){
        cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        for(i in base::seq_len(base::ncol(cluster_list))){
          if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- base::as.numeric(cluster_list[,i])
        }
        tag_names <- base::unique(cluster_list$tag)
        cluster_list <- base::lapply(tag_names, function(x){base::return(cluster_list[cluster_list$tag == x,-1])})
        base::names(cluster_list) <- tag_names
        
        base::cat(base::paste0("Loaded 'cluster_list' from ",fname, "\n"))
      } else cluster_list <- NULL
      
      output <- base::list(sample_log = sample_log,
                           motif_enrichment = motif_enrichment,
                           connections = connections,
                           cluster_properties = cluster_properties,
                           cluster_list = cluster_list,
                           parameters = parameters)
    }
    if(parameters$gliph_version == 2){
      base::cat("Output of gliph2 function is loaded.", "\n")
      
      ### load selected_motifs
      fname <- base::paste0(result_folder, "local_similarities_minp_",parameters$lcminp, "_minove_", base::paste(parameters$lcminove, collapse = "_"), "_kmer_mindepth_", parameters$kmer_mindepth, ".txt")
      if(base::file.exists(fname)){
        selected_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'selected_motifs' from ",fname, "\n"))
      } else selected_motifs <- NULL
      
      ### load all_motifs
      fname <- base::paste0(result_folder, "all_motifs.txt")
      if(base::file.exists(fname)){
        all_motifs <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'all_motifs' from ",fname, "\n"))
      } else all_motifs <- NULL
      
      ### build motif_enrichment
      motif_enrichment <- base::list(selected_motifs = selected_motifs, all_motifs = all_motifs)
      
      ### load global_enrichment
      fname <- base::paste0(result_folder, "global_similarities.txt")
      if(base::file.exists(fname)){
        global_enrichment <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'global_enrichment' from ",fname, "\n"))
      } else global_enrichment <- NULL
      
      ### load connections
      fname <- base::paste0(result_folder,"clone_network.txt")
      if(base::file.exists(fname)){
        connections <- utils::read.table(file = fname,sep = "\t",quote = "", header = FALSE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'connections' from ",fname, "\n"))
      } else connections <- NULL
      
      ### load cluster_properties
      fname <- base::paste0(result_folder,"convergence_groups.txt")
      if(base::file.exists(fname)){
        cluster_properties <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        base::cat(base::paste0("Loaded 'cluster_properties' from ",fname, "\n"))
      } else cluster_properties <- NULL
      
      ### load cluster_list
      fname <- base::paste0(result_folder, "cluster_member_details.txt")
      if(base::file.exists(fname)){
        cluster_list <- utils::read.table(file = fname,sep = "\t",quote = "", header = TRUE, stringsAsFactors = FALSE)
        for(i in base::seq_len(base::ncol(cluster_list))){
          if(base::suppressWarnings(base::any(base::is.na(base::as.numeric(cluster_list[,i])))) == FALSE) cluster_list[,i] <- base::as.numeric(cluster_list[,i])
        }
        tag_names <- base::unique(cluster_list$tag)
        cluster_list <- base::lapply(tag_names, function(x){base::return(cluster_list[cluster_list$tag == x,-1])})
        base::names(cluster_list) <- tag_names
        
        base::cat(base::paste0("Loaded 'cluster_list' from ",fname, "\n"))
      } else cluster_list <- NULL
      
      output <- base::list(motif_enrichment = motif_enrichment,
                           global_enrichment = global_enrichment,
                           connections = connections,
                           cluster_properties = cluster_properties,
                           cluster_list = cluster_list,
                           parameters = parameters)
    }
  }
  
  base::return(output)
}
