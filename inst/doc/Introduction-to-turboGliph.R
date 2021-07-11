## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.height = 5,
  fig.width = 5,
  fig.align = "center"
)
set.seed(100)

## ---- eval=FALSE--------------------------------------------------------------
#  if (!requireNamespace("BiocManager", quietly = TRUE)){
#      install.packages("BiocManager")
#      BiocManager::install()
#  }
#  
#  BiocManager::install("turboGliph")

## -----------------------------------------------------------------------------
library(turboGliph)
data("gliph_input_data")
gliph_input_data[1:3,]

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
res_turbogliph <- turboGliph::turbo_gliph(cdr3_sequences = gliph_input_data,
                                          n_cores = 1)

## -----------------------------------------------------------------------------
res_turbogliph$sample_log[1:5, 1:6]

## -----------------------------------------------------------------------------
res_turbogliph$motif_enrichment$selected_motifs[1:3,]

## -----------------------------------------------------------------------------
res_turbogliph$connections[1:3,]

## -----------------------------------------------------------------------------
res_turbogliph$cluster_properties[2:4,]

## -----------------------------------------------------------------------------
res_turbogliph$cluster_list[[1]][1:2,]

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
res_gliph2 <- turboGliph::gliph2(cdr3_sequences = gliph_input_data,
                                 n_cores = 1)

## -----------------------------------------------------------------------------
res_gliph2$motif_enrichment$selected_motifs[1:3,]

## -----------------------------------------------------------------------------
res_gliph2$global_enrichment[1:3,]

## -----------------------------------------------------------------------------
res_gliph2$connections[1:3,]

## -----------------------------------------------------------------------------
res_gliph2$cluster_properties[1:3,]

## -----------------------------------------------------------------------------
res_gliph2$cluster_list[[1]][1:2,]

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
res_gliph_combined <- turboGliph::gliph_combined(cdr3_sequences = gliph_input_data,
                                                 n_cores = 1)

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
res_gliph_combined_2 <- turboGliph::gliph_combined(cdr3_sequences = gliph_input_data,
                                                 min_seq_length = 0,
                                                 local_method = "fisher",
                                                 boost_local_significance = TRUE,
                                                 global_method = "fisher",
                                                 clustering_method = "GLIPH2.0",
                                                 scoring_method = "GLIPH2.0",
                                                 n_cores = 1)

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
res_gliph_combined_3 <- turboGliph::gliph_combined(cdr3_sequences = gliph_input_data,
                                                 min_seq_length = 0,
                                                 local_method = "fisher",
                                                 boost_local_significance = TRUE,
                                                 global_method = "fisher",
                                                 clustering_method = "GLIPH1.0",
                                                 scoring_method = "GLIPH1.0",
                                                 n_cores = 1)

## -----------------------------------------------------------------------------
res_gliph_combined$sample_log[1:5, 1:6]

## -----------------------------------------------------------------------------
res_gliph_combined$motif_enrichment$selected_motifs[1:3,]

## -----------------------------------------------------------------------------
res_gliph_combined_2$global_enrichment$selected_structs[1:3,]

## -----------------------------------------------------------------------------
res_gliph_combined$connections[1:3,]

## -----------------------------------------------------------------------------
res_gliph_combined$cluster_properties[1:3,]

## -----------------------------------------------------------------------------
res_gliph_combined$cluster_list[[1]][1:2,]

## ----warning=FALSE,message=FALSE,error=FALSE----------------------------------
plot_network(clustering_output = res_turbogliph,
             n_cores = 1)

## ----warning=FALSE,message=FALSE,error=FALSE----------------------------------
plot_network(clustering_output = res_turbogliph,
             color_info = "antigen.species",
             color_palette = grDevices::rainbow,
             local_edge_color = "darkgrey",
             global_edge_color = "orange",
             show_additional_columns = c("antigen.species", "antigen"),
             cluster_min_size = 6,
             n_cores = 1)

## ----warning=FALSE,message=FALSE,error=FALSE, results='hide',fig.keep='all'----
new_seqs <- turboGliph::de_novo_TCRs(convergence_group_tag = res_turbogliph$cluster_properties$tag[1],
                                     clustering_output = res_turbogliph,
                                     sims = 10000,
                                     make_figure = TRUE,
                                     n_cores = 1)

## ---- message = FALSE---------------------------------------------------------
new_seqs$de_novo_sequences[1:3,]

## -----------------------------------------------------------------------------
packageVersion("turboGliph")

sessionInfo()

