# =================================
# get the sPLS-DA selected features
# =================================
get_modality_features <- function(loadings, dim=-1) {
  if (dim < 1) {
    mask <- which(rowSums(loadings != 0) > 0)
  } else {
    mask <- which(loadings[,dim] != 0)
  }
  return(rownames(loadings)[mask])
}


load("EEN_Post_EEN_selection.RData")
splsda_features <- list(
  Metabolome=get_modality_features(splsda$loadings$Metabolome),
  Microbiome=get_modality_features(splsda$loadings$Microbiome)
)
print(splsda_features)


# ===========================================
# get the features in the correlation network
# ===========================================
load("een_posteen_communities.RData")
igraph::V(net$network)$name


# =============================================
# map metabolite feature IDs to the names shown
# =============================================
load("processed_data.RData")
first_name_map[splsda_features$Metabolome]
