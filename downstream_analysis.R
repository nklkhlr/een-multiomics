source("multi_omics.R")

library(argparse)

set.seed(123)

# =========
# Functions
# =========
largest_connected_component <- function(graph) {
  comps <- clusters(graph, mode="weak")
  v_idxs <- V(graph)[comps$membership == which.max(comps$csize)]
  return(induced_subgraph(graph, v_idxs))
}


community_analysis <- function(
  pargs, subg, net, mgs_data, combined, var_een, timepoint, shapes, node_names,
  first_name_map, mic_annotation, i, file_annotation="community"
) {
  ### plot community network
  pdf(
    paste0(pargs$result_path, "/een_posteen_", file_annotation, "_", i, ".pdf"))
  par(mar=c(2, 1, 1, 1))
  # see `correlation_network` for details
  plot_corr_net(subg, edge_cmap=net$edge_cmap,
                node_cmap=net$node_cmap,
                node_shapes=shapes,
                vertex.label=node_names)
  legend(
    "topright", title="Associated With", legend=c("EEN", "PostEEN"),
    fill=c("#1F77B4", "#FF7F0E"), box.lty=0, cex=1
  )
  dev.off()

  ### store info about community features and plot evaluation
  file_base <- paste0(pargs$result_path, "/", file_annotation, i)

  # boxplot feature expression
  expression_plot <- boxplot_community(
    V(subg)$name, combined$normalised, mgs_data, timepoint,
    # NOTE: metadata is loaded from within multi_omics.R
    #       (metabolomics_meta.RData)
    metadata, healing, name_map=first_name_map, nrow=2
  )
  write_csv(
    expression_plot$data, file=paste0(file_base, "_boxplot_data.csv")
  )
  ggsave(
    filename=paste0(file_base, "_boxplots.pdf"),
    width=16, height=18, plot=expression_plot$plot + scale_colour_d3()
  )

  features_by_modality <- extract_community_features(
    subg, combined$feature_info, mic_annotation, var_een$data, file_base)

  # check if either metabolites or Microbiomes are 1D or empty
  nmets <- length(features_by_modality$metabolites) < 2
  nmics <- length(features_by_modality$mics) < 2
  if (nmets | nmics) {
    # compute ROC-AUC evaluation on communities only. In principle this works
    # the same way as the above evaluation
    roc_eval_models <- roc_test_community_features(
      features_by_modality$metabolites, features_by_modality$mics,
      combined$normalised, mgs_data, timepoint
    )
    # plot evaluation results for a single modality (even if we have a single
    # metabolite or Microbiome, we cannot show these plots
    roc_eval <- plot_single_modality_cv_auc(roc_eval_models$roc)
    mean_auc <- round(mean(roc_eval$aucs[,1]), digits=2)
    sd_auc <- round(sd(roc_eval$aucs[,1]), digits=2)

    # plot annotation
    auc_annotation <- annotate(
      "text", x=.75, y=.25,
      label=paste0("Mean AUC = ", mean_auc, " ± ", sd_auc), size=9
    )
    roc_eval$roc_plot <- roc_eval$roc_plot + auc_annotation

    # plot loading scores
    cv_loadings <- lapply(
      roc_eval_models$model,
      extract_cv_loadings, name_map=first_name_map,
      positive_group=pargs$association_greater
    )
    # plot loadings for each feature over the cross-validation
    loading_plot <- plot_cv_loadings(cv_loadings) +
      theme_pubr() +
      scale_colour_d3("category10")
    loading_data <- loading_plot$data
  } else {
    # compute ROC-AUC evaluation on communities only. In principle this works
    # the same way as the above evaluation
    roc_eval_models <- roc_test_community_features(
      features_by_modality$metabolites, features_by_modality$mics,
      combined$normalised, mgs_data, timepoint
    )

    # plot evaluation results
    roc_eval <- plot_cv_auc(roc_eval_models$roc)
    mean_auc <- sapply(colMeans(roc_eval$aucs), round, digits=2)
    sd_auc <- sapply(apply(roc_eval$aucs, 2, sd), round, digits=2)

    # plot annotation
    auc_labels <- data.frame(
      label=c(
        paste0("Mean AUC = ", mean_auc[1], " ± ", sd_auc[1]),
        paste0("Mean AUC = ", mean_auc[2], " ± ", sd_auc[2])
      ),
      Modality=c("Metabolome", "Microbiome"),
      x=.75,
      y=.25
    )
    auc_annotation <- geom_text(
      data=auc_labels, aes(x=x, y=y, label=label, colour=NULL), size=9,
      show.legend=FALSE
    )
    roc_eval$roc_plot <- roc_eval$roc_plot + auc_annotation

    # plot loading scores
    loadings <- lapply(
      roc_eval_models$model, function(x) {
        # NOTE: the model contains information for both modalities => no need
        #       for inner loop
        lapply(x, extract_cv_loadings, name_map=first_name_map)[[1]]
      }
    )
    # plot loadings for each feature over the cross-validation
    met_loading_plot <- plot_cv_loadings(
      lapply(loadings, function(x) x$Metabolome)) +
      theme_pubr() +
      scale_color_d3("category10")
    mic_loading_plot <- plot_cv_loadings(
      lapply(loadings, function(x) x$Microbiome)) +
      theme_pubr() +
      scale_color_d3("category10")
    loading_plot <- ggarrange(
      met_loading_plot, mic_loading_plot)
    loading_data <- rbind(met_loading_plot$data, mic_loading_plot$data)
  }

  # save all plots computed above
  base_path <- paste0(pargs$result_path, "/", file_annotation, i)
  ggsave(
    filename=paste0(base_path, "_cv_aucs.pdf"),
    width=16, height=9, device="pdf", plot=roc_eval$auc_plot + theme_pubr()
  )
  ggsave(
    filename=paste0(base_path, "_cv_rocs.pdf"),
    width=16, height=9, device="pdf",
    plot=roc_eval$roc_plot + theme_pubr() + scale_colour_d3()
  )
  ggsave(
    filename=paste0(base_path, "_loadings.pdf"),
    width=16, height=9, device="pdf", plot=loading_plot
  )

  loading_data$name <- gsub(
    "\\.\\.\\.*", "", rownames(loading_data))
  write.csv(loading_data, file=paste0(base_path, "_loadings.csv"))

  if (!is.null(roc_eval$mean_curve)) {
    # if iterations don't have the same number of thresholds recored for the ROC
    # calculations, so we can't easily average them. There's ways around it, but I
    # never had the time do so and this case only happened once. Let me know if this
    # is a problem with the metagenomics data - we can probably figure this out.
    if (length(mean_auc) == 1) {
      p <- plot_average_roc(roc_eval$mean_curve) + theme_pubr() + auc_annotation
      ggsave(
        filename=paste0(base_path, "_average_roc.pdf"),
        width=16, height=9, device="pdf", plot=p
      )
    } else {
      p <- plot_average_roc(roc_eval$mean_curve) + theme_pubr() + facet_wrap(.~Modality)
      ggsave(
        filename=paste0(base_path, "_average_roc.pdf"),
        width=16, height=9, device="pdf", plot=p
      )
    }
  }

  return(
    random_community_test(
      features_by_modality$metabolites, features_by_modality$mics,
      combined$normalised, mgs_data, timepoint, 500
    )
  )
}

# ==========
# CLI parser
# ==========
parser <- ArgumentParser(description="multi-omics downstream analysis")

parser$add_argument(
  "result_path", type="character",
  help="Path to the folder where the results should be saved"
)

# feature extraction based on weights
parser$add_argument(
  "--x-threshold", type="double", default=.7,
  help=paste0(
    "Minimum correlation with latent dimension 1 when extracting features ",
    "when extracting most relevant features from loadings."
  )
)
parser$add_argument(
  "--y-threshold", type="double", default=1,
  help=paste0(
    "Minimum correlation with latent dimension 2 when extracting features ",
    "when extracting most relevant features from loadings."
  )
)

parser$add_argument(
  "--association-greater", type="character", default="EEN",
  help="Latent dimension to use for determining association"
)
parser$add_argument(
  "--association-dim", type="character", default="x",
  help="Latent dimension to use for determining association"
)

# to use it from the R console directly
pargs <- tryCatch(
  parser$parse_args(),
  error=function(err) {
    args <- list(
      result_path="Results",
      x_threshold=.7,
      y_threshold=1,
      association_greater="EEN",
      association_dim="x"
    )
  }
)

# ============
# Loading data
# ============
load(paste0(pargs$result_path, "/", "EEN_Post_EEN_selection.RData"))
load(paste0(pargs$result_path, "/", "processed_data.RData"))

# ====
# Main
# ====
### Feature Extraction
# This extracts features by their loading values. By setting both thresholds
# to 0 all selected features are retained. In this case, the y-axis (i.e. the
# second latent dimension) did not compute to separating the samples, so features
# with a high correlation with Dim 2 and low correlation with Dim 1 were removed
var_een <- plot_variable_subset(
  splsda, colour="Block",
  label="names", geom="text",
  x_threshold=pargs$x_threshold, y_threshold=pargs$y_threshold
)
var_een$data <- var_een$data[!grepl("pos1$", rownames(var_een$data)),]
ggsave(
  filename=paste0(pargs$result_path, "/een_posteen_loadings.pdf"), width=12,
  height=9, device="pdf",
  plot=var_een$plot + scale_colour_brewer(palette='Set1')
)
# this is a simple "association" based on loadings: in this case, features
# with higher values on the x-axis were positively associated with EEN.
# => use the loadings + latent space plots to determine your associations
if (pargs$association_greater == "EEN") {
  var_een$data$AssociatedWith <- ifelse(
    var_een$data[[pargs$association_dim]] > 0, "EEN", "PostEEN")
} else {
  var_een$data$AssociatedWith <- ifelse(
    var_een$data[[pargs$association_dim]] > 0, "PostEEN", "EEN")
}

# all metabolite features started with "FT", which was never possible for
# 16S features
metabolome_features <- rownames(var_een$data)[grep("^FT", rownames(var_een$data))]

merged_data <- rbind(corrected_data$z_score[,subsamples],
                     zotu_data[,subsamples])

# split labels => annotates samples by whether they belong to the train or test
# set
set <- c(rep("train", length(split$train)), rep("test", length(split$test)))
names(set) <- c(split$train, split$test)

# feature annotations by names
name_map <- c(combined$feature_info$annot_ms1, mic_annotation)
names(name_map) <- c(rownames(combined$feature_info), names(mic_annotation))

# PCA of the combined data, that only contains the selected features. This is
# not an optimal plot, as it does not consider metabolite and 16S variation
# equally
pca <- prcomp(t(merged_data[rownames(var_een$data),]), scale.=TRUE, center=T)
feature_pca <- plot_pca(pca, timepoint, size=5, shape=set)
ggsave(filename=paste0(pargs$result_path, "/een_posteen_selection_subset_pca.pdf"),
       width=16, heigh=9, plot=feature_pca)


### Correlation network
# comm_file <- paste0(pargs$result_path, "/een_posteen_communities.RData")
comm_file <- paste0(pargs$result_path, "/een_posteen_communities_well_annot.RData")
if (!file.exists(comm_file)) {
  # Spearman's rank-based correlation network. See `correlation_network.R`
  net <- correlation_network(
    data=t(merged_data[rownames(var_een$data),]),
    metabo_micro_only=TRUE,
    threshold=.05, adjust='BH'
  )
  # Removes all unannotated vertices. We used this option as we were
  # interested in annotated metabolites
  net$network <- delete.vertices(
    net$network, grepl("^NA.[0-9]*|^NA$", V(net$network)$name))

  # setup for plotting network
  shape_map <- c("Metabolite"="csquare", "Microbiome"="circle")
  feat_type <- ifelse(
    sapply(V(net$network)$name, function(v) startsWith(v, "FT")),
    "Metabolite", "Microbiome"
  )
  shapes <- shape_map[feat_type]
  names(shapes) <- names(feat_type)

  V(net$network)$colour <- sapply(
    V(net$network),
    function(x) {
      if (var_een$data[x, "AssociatedWith"] == "EEN") return("#1F77B4")
      return("#FF7F0E")
    }
  )

  node_names <- name_map[V(net$network)$name]
  pdf(paste0(pargs$result_path, "/een_posteen_correlation_network.pdf"),
      width=12, height=12)
  plot_corr_net(net$network, edge_cmap=net$edge_cmap,
                node_cmap=net$node_cmap,
                node_shapes=shapes,
                vertex.label=node_names)
  legend(
    "topright", title="Associated With", legend=c("EEN", "PostEEN"),
    fill=c("#1F77B4", "#FF7F0E"), box.lty=0, cex=1
  )
  dev.off()

  # extract communities using louvain modularity
  communities <- cluster_louvain(net$network)
  save(net, communities, file=comm_file)
} else {
  # load previously computed communities
  load(comm_file)

  shape_map <- c("Metabolite"="csquare", "Microbiome"="circle")
  node_names <- name_map[V(net$network)$name]

  feat_type <- ifelse(
    sapply(V(net$network)$name, function(v) startsWith(v, "FT")),
    "Metabolite", "Microbiome"
  )
  shapes <- shape_map[feat_type]
  names(shapes) <- names(feat_type)
}
first_name_map <- sapply(name_map, substr, start=1, stop=30)
save(
  zotu_data, mic_annotation, timepoint, subsamples, healing, first_name_map,
  file=paste0(pargs$result_path, "/", "processed_data.RData")
)

# run the analysis on the entire largest connected component of the network
lcomp <- largest_connected_component(net$network)
community_analysis(
  pargs, lcomp, net, zotu_data, combined, var_een, timepoint,
  shapes, node_names, first_name_map, mic_annotation, "total",
  file_annotation="network"
)

# evaluate and plot communities
empirical_pvalues <- list()
for (i in seq_along(communities)) {

  comm <- communities[[i]]
  if (length(comm) < 5) next

  cat("Community ", i, "\n")

  subg <- induced_subgraph(net$network, comm)
  pval <- community_analysis(
    pargs, subg, net, zotu_data, combined, var_een, timepoint,
    shapes, node_names, first_name_map, mic_annotation, i
  )
  empirical_pvalues[[paste0("Community", i)]] <- pval
}
write.csv(
  as.data.frame(empirical_pvalues),
  file=paste0(pargs$result_path, "/", "empirical_pvalues.csv")
)

