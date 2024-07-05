library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)

#' Plot the embedding contained in a `prcomp` object
#'
#' @param pca prcomp object with which the PCA to plot was computed
#' @param group vector containing the sample group annotations for choosing
#'        sample colors
#' @param ellipse whether to plot ellipses around data points indicating
#'        confidence intervals per group
#' @param shape optional vector containing sample group annotations for
#'        choosing sample shapes
#' @param ... optional arguments passed to ggplot::geom_point
plot_pca <- function(pca, group, ellipse=TRUE, shape=NULL, ...) {
    expl_var <- pca$sdev^2/sum(pca$sdev^2)
    data.frame(pca$x) %>%
    add_column(
      Sample=rownames(pca$x),
      Group=group[rownames(pca$x)],
    ) -> plot_data
    if (!is.null(shape)) {
      plot_data$Shape <- shape[rownames(pca$x)]
      p <- ggplot(plot_data,
                  aes(x=PC1, y=PC2, colour=Group, label=Sample,
                      shape=Shape))
    } else {
      p <- ggplot(plot_data,
                  aes(x=PC1, y=PC2, colour=Group, label=Sample,
                      # label2=Baseline, label3=DiseaseGroup)) +
      ))
    }
    p +
    geom_point(...) +
    theme_pubr() +
    scale_colour_d3('category20') +
    scale_fill_d3('category20') +
    xlab(paste0('PC1 (', round(expl_var[1], 3) * 100, '%)')) +
    ylab(paste0('PC2 (', round(expl_var[2], 3) * 100, '%)')) +
    theme(legend.position="right",
          axis.text=element_text(size=26),
          axis.title=element_text(size=28),
          legend.title=element_text(size=28),
          legend.text=element_text(size=26),
          legend.key.size=unit(3, "line")) -> p
    if (ellipse) {
      p <- p + stat_ellipse(
        aes(fill=Group, colour=NULL, shape=NULL),
        geom="polygon", alpha=.15
      )
    }
    return(p)
}

