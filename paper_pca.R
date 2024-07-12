library(tidyverse)
library(see)

# provides: plot_pca
source("een_diversity.R")

# provides: corrected_data and combined_data
load("Data/metabolomics_data.RData")
# provides: metadata
load("Data/metabolomics_meta.RData")


column_to_vector <- function(df, col) {
  vec <- df[,col]
  names(vec) <- rownames(df)
  return(vec)
}


#' Compute effect sizes of Diet and Activity over all features individually
#' while correcting for Patient-specific effects with a linear mixed model of
#' the form `Feature ~ intercept + Diet + Activity + Diet * Activity + (1|Patient)`
#'
#' @param num_data Feature data with samples in rows. All columns will be tested
#'                 individually
#' @param meta_data `data.frame` containing Diet, Activity and Patient annotation
#' @param intercept Whether to fit the LMM with an intercept
#' @param ... Optional parameters to pass to `nlme::lmer`
#'
#' @returns list of length m with per-feature model outputs
activity_diet_lmm <- function(num_data, meta_data, intercept=TRUE, ...) {
  samples <- intersect(rownames(num_data), rownames(meta_data)) 
  lm_data <- cbind(
    num_data[samples,], meta_data[samples, c("Diet", "Activity", "Patient")])
  
  if (intercept) 
    base_formula <- " ~ 1 + Diet + Activity + Diet * Activity + (1|Patient)"
  else 
    base_formula <- " ~ Diet + Activity + Diet * Activity + (1|Patient)" 
  
  models <- lapply(
    colnames(num_data),
    function(feature) {
      lm_out <- lmerTest::lmer(
        as.formula(paste0(feature, base_formula)), data=lm_data, ...)
    }
  ) 
  names(models) <- colnames(num_data)
  
  return(models)
}


# ==== #
# Main #
# ==== #
metadata <- add_column(
  metadata,
  Activity=ifelse(str_to_upper(metadata$PGA) == "REMISSION",
                  "inactive", "active"),
) %>%
  mutate(
    `Diet (=PreEEN, EEN, PostEEN)`=gsub("Post-EEN", "PostEEN",
                                        metadata$`Diet (=PreEEN, EEN, PostEEN)`)
  ) %>%
  filter(rownames(metadata) %in% colnames(combined$z_score))


g_patients <- c(
  "008", "011", "014", "017", "028", "030", "038", "041", "046", "157", "173",
  "180", "181", "186", "194", "203", "217", "218"
)
g_samples <- rownames(metadata)[metadata$Patient %in% g_patients]

well_annot <- rownames(combined$feature_info)[
  combined$feature_info$annot_ms2 %in% c("3", "4")]

# metadata vectors
activity <- column_to_vector(metadata, "Activity")
diet <- column_to_vector(metadata, "Diet (=PreEEN, EEN, PostEEN)")


# Fig. 1G
g_samples <- names(activity)[!is.na(activity)]
g_activity <- activity[g_samples][grepl("CD_EEN", metadata[g_samples, "Group"])]
g_samples <- names(g_activity)

total_pca <- prcomp(t(corrected_data$z_score[, g_samples]))
g_plot <- plot_pca(total_pca, g_activity, shape=diet, size=5) +
  scale_colour_manual(
    values=c("active"=rgb(160, 25, 21, maxColorValue=255),
             "inactive"=rgb(21, 91, 42, maxColorValue=255))
  ) +
  scale_fill_manual(
    values=c("active"=rgb(160, 25, 21, maxColorValue=255),
             "inactive"=rgb(21, 91, 42, maxColorValue=255))
  )
ggsave(
  filename="Results/paper_plots/Fig1G.pdf",
  width=12, height=10, plot=g_plot
)


# Fig. 1H
h_samples <- g_samples[diet[g_samples] != "PreEEN"]

diet_pca <- prcomp(t(corrected_data$z_score[, h_samples]))
h_plot <- plot_pca(diet_pca, diet[h_samples], shape=activity, size=5) +
  scale_colour_manual(
    values=c("EEN"=rgb(33, 62, 131, maxColorValue=255),
             "PostEEN"=rgb(238, 121, 4, maxColorValue=255))
  ) +
  scale_fill_manual(
    values=c("EEN"=rgb(33, 62, 131, maxColorValue=255),
             "PostEEN"=rgb(238, 121, 4, maxColorValue=255))
  )
ggsave(
  filename="Results/paper_plots/Fig1H.pdf",
  width=12, height=10, plot=h_plot
)


# LMM testing
meta <- data.frame(
  Diet=diet[h_samples],
  Activity=activity[h_samples],
  Patient=metadata[h_samples, "Patient"]
)

out <- activity_diet_lmm(diet_pca$x, meta)
pvals <- sapply(out, function(x) summary(x)$coefficients[,"Pr(>|t|)"])
qvals <- t(apply(pvals, 1, p.adjust, method="bonferroni"))

as.data.frame(qvals) %>%
  add_column(Covariate=rownames(qvals)) %>%
  pivot_longer(-Covariate, names_to="Feature", values_to="qval") %>%
  ggplot(aes(x=Covariate, y=-log10(qval), color=Covariate)) +
    geom_violindot() +
    geom_hline(yintercept=1.4)
