library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

load("Data/comparison_samples.RData")
load("Data/metabolomics_data.RData")
load("Data/metabolomics_meta.RData")
load("Data/zotu_data.RData")


# plot knee points in importance plot
plot_importance_knee <- function(loadings, knee_points) {
  loadings %>%
    arrange(desc(abs(.$`Mean Importance`))) %>%
    ggplot(
      aes(
        x=reorder(Metabolite, abs(`Mean Importance`), FUN=desc),
        y=abs(`Mean Importance`)
      )
    ) +
    geom_point() +
    geom_vline(data=data.frame(), aes(xintercept=knee_points + .5))
}


# load community loading data coming from `downstream_analysis.R`
loading_files <- list.files("./Results/", "community[0-9]+_loadings.csv")
community_index <- gsub("[a-z]|_|\\.", "", loading_files)

# compute feature-wise means
mean_loadings <- list()
for (i in seq_along(loading_files)) {
  idx <- paste0("Community ", community_index[i])
  comm_loadings <- read.csv(paste0("./Results/", loading_files[i]))
  comm_loadings$name <- gsub("\\.\\.\\.[0-9]+$", "", comm_loadings$X)
  mean_loadings[[i]] <- comm_loadings %>%
    group_by(name) %>%
    summarise(mean(importance)) %>%
    add_column(Community=idx)
}

mean_loadings <- bind_rows(mean_loadings) %>%
  rename(Metabolite=name, `Mean Importance`=`mean(importance)`) %>%
  arrange(desc(abs(`Mean Importance`)))

# heatmap of mean loadings
p <- ggplot(
  mean_loadings,
  aes(x=Community, y=reorder(Metabolite, `Mean Importance`))
) +
  geom_tile(aes(fill=`Mean Importance`)) +
  facet_wrap(.~Community, scales="free") +
  scale_fill_gsea() +
  xlab("") +
  ylab("Feature") +
  theme_minimal() +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )

ggsave(filename="Results/mean_loadings.pdf", width=12, height=12, plot=p)


# separate metabolite and and microbiome data
metabolite_mask <- grepl("^FT", mean_loadings$Metabolite)
metabolites <- mean_loadings$Metabolite[metabolite_mask]
zotus <- mean_loadings$Metabolite[!metabolite_mask]


metabolite_data <- corrected_data$z_score[
  metabolites, rownames(comparison_samples)]
zotu_data <- zotu_data$clr[zotus, rownames(comparison_samples)]

feature_meta <- column_to_rownames(mean_loadings, "Metabolite")
feature_meta$Association <- ifelse(
  feature_meta$`Mean Importance` > 0, "EEN", "PostEEN")

# select diet sample annotations
sample_meta <- dplyr::select(
  metadata[colnames(metabolite_data),], "Diet (=PreEEN, EEN, PostEEN)") %>%
  rename("Diet"="Diet (=PreEEN, EEN, PostEEN)")
sample_meta$Diet <- gsub("Post-EEN", "PostEEN", sample_meta$Diet)

# set up colors for heatmap
colors <- list(
  Diet=c(EEN="#1F77B4", PostEEN="#FF7F0E"),
  Association=c(EEN="#1F77B4", PostEEN="#FF7F0E"),
  Community=RColorBrewer::brewer.pal(9, "YlGn")[c(8, 7, 6, 5, 4)],
  `Mean Importance`=c("white", "firebrick")
)
names(colors$Community) <- sapply(1:5, function(i) paste0("Community ", i))

# plot heatmaps over all selected features by modality
pheatmap::pheatmap(
  metabolite_data,
  scale="none",
  annotation_col=sample_meta,
  annotation_row=feature_meta[rownames(metabolite_data),],
  annotation_colors=colors,
  color=colorRampPalette(c("blue", "white", "red"))(20),
  width=16,
  height=9,
  filename="Results/metabolite_loadings.pdf"
)

pheatmap::pheatmap(
  zotu_data,
  scale="none",
  annotation_col=sample_meta,
  annotation_row=feature_meta[rownames(zotu_data),],
  annotation_colors=colors,
  color=colorRampPalette(c("blue", "white", "red"))(20),
  width=16,
  height=9,
  filename="Results/zotu_loadings.pdf",
)


# ============================== #
### plotting top-features only ###
# ============================== #
# fixed sample order by diet and patient
sample_order <- rownames(comparison_samples)[
  order(
    comparison_samples$`Diet (=PreEEN, EEN, PostEEN)`, 
    comparison_samples$Patient
  )
]
 
# selection by knee point of absolute loadings
met_loadings <- mean_loadings[grep("^FT", mean_loadings$Metabolite),]
met_knee_points <- c(11, 22)
plot_importance_knee(met_loadings, met_knee_points) +
  xlab("Metabolite") +
  ylab("Absolute Mean Importance") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

zotu_loadings <- mean_loadings[grep("^Zotu", mean_loadings$Metabolite),]
zotu_knee_points <- c(5, 15, 20)
plot_importance_knee(zotu_loadings, zotu_knee_points) +
  xlab("zOTU") +
  ylab("Absolute Mean Importance") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

# plot top metabolites
for (i in seq_along(met_knee_points)) {
  file <- paste0("Results/metabolite_loadings_top_", met_knee_points[i], ".pdf")
  pheatmap::pheatmap(
    metabolite_data[met_loadings$Metabolite[1:met_knee_points[i]], sample_order],
    scale="row",
    color=colorRampPalette(c("blue", "white", "red"))(20),
    cluster_cols=FALSE,
    annotation_col=sample_meta,
    annotation_row=feature_meta[rownames(metabolite_data),],
    annotation_colors=colors,
    width=16,
    height=9,
    filename=file
  ) 
}

# plot top zOTUs
for (i in seq_along(zotu_knee_points)) {
  file <- paste0("Results/zotu_loadings_top_", zotu_knee_points[i], ".pdf")
  pheatmap::pheatmap(
    zotu_data[zotu_loadings$Metabolite[1:zotu_knee_points[i]], sample_order],
    scale="row",
    color=colorRampPalette(c("blue", "white", "red"))(20),
    cluster_cols=FALSE,
    annotation_col=sample_meta,
    annotation_row=feature_meta[rownames(zotu_data),],
    annotation_colors=colors,
    width=16,
    height=9,
    filename=file
  )
}
