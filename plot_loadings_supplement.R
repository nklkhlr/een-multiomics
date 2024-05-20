library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

load("Data/comparison_samples.RData")
load("~/PhD/EEN_Integration/ReAnalysis/metabolomics_data.RData")
load("~/PhD/EEN_Integration/ReAnalysis/metabolomics_meta.RData")
load("~/PhD/EEN_Integration/ReAnalysis/zotu_data.RData")

loading_files <- list.files("./Results/", "community[0-9]+_loadings.csv")
community_index <- gsub("[a-z]|_|\\.", "", loading_files)

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
  rename(Metabolite=name, `Mean Importance`=`mean(importance)`)
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


metabolite_mask <- grepl("^FT", mean_loadings$Metabolite)
metabolites <- mean_loadings$Metabolite[metabolite_mask]
zotus <- mean_loadings$Metabolite[!metabolite_mask]


metabolite_data <- corrected_data$z_score[
  metabolites, rownames(comparison_samples)]
zotu_data <- zotu_data$clr[zotus, rownames(comparison_samples)]

feature_meta <- column_to_rownames(mean_loadings, "Metabolite")
feature_meta$Association <- ifelse(
  feature_meta$`Mean Importance` > 0, "EEN", "PostEEN")

sample_meta <- dplyr::select(
  metadata[colnames(metabolite_data),], "Diet (=PreEEN, EEN, PostEEN)") %>%
  rename("Diet"="Diet (=PreEEN, EEN, PostEEN)")
sample_meta$Diet <- gsub("Post-EEN", "PostEEN", sample_meta$Diet)

colors <- list(
  Diet=c(EEN="#1F77B4", PostEEN="#FF7F0E"),
  Association=c(EEN="#1F77B4", PostEEN="#FF7F0E"),
  # Community=c(RColorBrewer::brewer.pal(9, "BlGn")[5:9]),
  Community=RColorBrewer::brewer.pal(9, "YlGn")[c(8, 7, 6, 5, 4)],
  `Mean Importance`=c("white", "firebrick")
)
names(colors$Community) <- sapply(1:5, function(i) paste0("Community ", i))

pheatmap::pheatmap(
  metabolite_data,
  scale="none",
  annotation_col=sample_meta,
  annotation_row=feature_meta[rownames(metabolite_data),],
  annotation_colors=colors,
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
  width=16,
  height=9,
  filename="Results/zotu_loadings.pdf",
)
