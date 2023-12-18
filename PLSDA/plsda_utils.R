library(mixOmics)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(tidyverse)


Rcpp::sourceCpp('PLSDA/roc_curve.cpp')


# Plot the latent sapce embedding of a trained PLSDA model
# essentially this is just a wrapper around mixOmics::plotIndiv
plot_samples <- function(plsda.obj, group.name="Group",
                         shape=NULL, plot.ellipse=TRUE,
                         size=2, geom="text", ...) {
  pdf(file=NULL)
  sample.plot <- plotIndiv(plsda.obj, ...)
  dev.off()

  if (!is.null(shape)) {
    sample.plot$df$shape <- shape[[2]][rownames(sample.plot$df)]
    p <- ggplot(data=sample.plot$df,
                aes_string(x="x", y="y",
                           colour="group",
                           shape="shape",
                           label="names"),
                           size=size) +
            labs(shape=shape[[1]],
                 colour=group.name)
  } else {
    p <- ggplot(data=sample.plot$df,
                aes_string(x="x", y="y",
                           colour="group",
                           label="names"),
                size=size) +
            labs(colour=group.name)
  }

  if (plot.ellipse) {
    p <- p + stat_ellipse(data=sample.plot$df,
                          geom="polygon",
                          alpha=.15,
                          aes_string(x="x", y="y",
                                     fill="group"),
                          inherit.aes = FALSE,
                          show.legend=FALSE)
  }

  if (geom == "point") p <- p + geom_point(size=size)
  if (geom == "text") p <- p + geom_text(size=size)

  p +
    scale_colour_d3("category10") +
    scale_fill_d3("category10") +
    theme_pubr() +
    xlab(sample.plot$graph$labels$x) +
    ylab(sample.plot$graph$labels$y)
}

# plot a circle => helper for loading plot
circle.path <- function(r=1, n=100){
  p <- seq(0, 2*pi, length.out=n)
  return(data.frame(x=r*cos(p),
                    y=r*sin(p)))
}


# plot the loadings of PLSDA model
plot_variable_subset <- function(plsda.obj, add.columns=NULL,
                                 colour=NULL, shape=NULL,
                                 label=NULL, geom="point",
                                 x_threshold=.975, y_threshold=.975,
                                 quantiles=TRUE, relation="OR",
                                 return.full=FALSE, size=2) {

  pdf(file=NULL)
  var.plot <- plotVar(plsda.obj)
  dev.off()

  if (quantiles) {
    x.cutoff <- quantile(abs(var.plot$x), x_threshold)
    y.cutoff <- quantile(abs(var.plot$y), y_threshold)
  } else {
    x.cutoff <- x_threshold
    y.cutoff <- y_threshold
  }


  if (!is.null(add.columns)) {
    var.plot <- cbind(var.plot, add.columns[rownames(var.plot),])
  }

  if (relation == "AND") {
    select.var.plot <- dplyr::filter(var.plot,
                                     abs(x) >= x.cutoff & abs(y) > y.cutoff)
  } else {
    select.var.plot <- dplyr::filter(var.plot,
                                     abs(x) >= x.cutoff | abs(y) > y.cutoff)
  }

  if (dim(select.var.plot)[1] == 0) stop("No features left after filtering. Choose different cutoffs or ")

  ggplot() +
    geom_path(data=circle.path(), aes(x=x, y=y)) +
    geom_path(data=circle.path(r=.5), aes(x=x, y=y)) -> p

  if (geom == "point") {
   p <- p + geom_point(data=select.var.plot,
                       aes_string(x="x", y="y", colour=colour,
                                  shape=shape),
                       size=size)
  }
  if (geom == "text") {
    p <- p + geom_text(data=select.var.plot,
                        aes_string(x="x", y="y", colour=colour,
                                   label=label),
                        size=size)
  }

  p <- p +
        theme_pubr() +
        theme(legend.position="right")

  if (return.full) return(list(plot=p,
                               filtered.data=select.var.plot,
                               full.data=var.plot))
  return(list(plot=p,
              data=select.var.plot))
}

features_from_subset <- function(features, meta) {
  feat_sub <- features[features %in% rownames(meta)]
  return(meta[feat_sub,])
}


get_signed_loadings <- function(
  model, block=NULL, ncomp=1, positive_group="EEN", ...
) {
  # NOTE: this only works for __two__ sample groups!
  # preparing sample coordinates
  if (is.null(block))
    block <- "X"

  sample_df <- data.frame(
    coord=model$variates[[block]][,ncomp],
    group=model$Y[rownames(model$variates[[block]])]
  )
  med_coord <- group_by(sample_df, group) %>%
    summarise(median=median(coord)) %>%
    column_to_rownames("group")

  if (nrow(med_coord) != 2)
    stop("Exactly two sample groups are required for this function!")

  neg_group <- rownames(med_coord)[rownames(med_coord) != positive_group]
  if ( all(med_coord < 0) | all(med_coord > 0) ) {
    warning("Median sample group coordinates are on the same side")

    if (med_coord[positive_group,] > med_coord[neg_group,])
      sign_mult <- 1
    else
      sign_mult <- -1
  } else {
    if (med_coord[positive_group,] > 0)
      sign_mult <- 1
    else
      sign_mult <- -1
  }

  # extracting feature importances and correcting sign for consistency over
  # different splits
  loadings <- data.frame(
    name=rownames(model$loadings[[block]]),
    importance=model$loadings[[block]][,ncomp] * sign_mult
  )
  loadings$GroupContrib <- ifelse(
    loadings$importance > 0, positive_group, neg_group)

  return(loadings)
}


loading_plot_per_block_ <- function(
  plsda.obj, name_map, block, n_features=c(20, 10), ncomp=1,
  positive_group="EEN"
) {
  if (!positive_group %in% plsda.obj$Y)
    stop("'positive_group' not found in model target variable")

  loadings <- get_signed_loadings(
    plsda.obj, block=block, ncomp=ncomp, positive_group=positive_group)

  if (block == "Microbiome") {
    # => genera are often duplicated => would be aggregated when plotting
    loadings$name <- sapply(
      rownames(loadings),
      function(x) {
        return(
          paste(
            str_replace(x, "Zotu", "Microbiome"), name_map[x],
            sep=" -- "
          )
        )
      }
    )
  } else {
    loadings$name <- name_map[rownames(loadings)]
    dups <- duplicated(loadings$name)
    if (any(dups)) {
      loadings$name[dups] <- paste0(loadings$name[dups], "__")
    }
  }
  loadings <- loadings[1:n_features[1],]

  max_val <- max(abs(loadings$importance))
  loading_plot <- ggplot(
      loadings,
      aes(x=importance, y=reorder(name, abs(importance)), fill=GroupContrib)
    ) +
    geom_col() +
    xlim(-max_val, max_val) +
    xlab("Feature Importance (sPLS-DA loading)") +
    ylab("") +
    theme_pubr() +
    theme(axis.line.y.left=element_blank(), axis.ticks.y=element_blank())

  return(list(loadings=loadings, plot=loading_plot))
}


plot_loadings <- function(
  plsda.obj, name_map, n_features=c(20, 10), ncomp=1, ...
) {
  met_loadings <- loading_plot_per_block_(
    plsda.obj, name_map, "Metabolome", n_features, ncomp, ...)
  mic_loadings <- loading_plot_per_block_(
    plsda.obj, name_map, "Microbiome", n_features, ncomp, ...)

  return(list(Metabolome=met_loadings, Microbiome=mic_loadings))
}


cv_loadings_per_block <- function(
  model, name_map, block, ncomp=1, positive_group="EEN", ...
) {
  loadings <- get_signed_loadings(
    model, block=block, ncomp=ncomp, positive_group=positive_group)

  if (is.null(block)) {
    block <- "None"
  }
  if (block == "Microbiome") {
    # => genera are often duplicated => would be aggregated when plotting
    loadings$name <- sapply(
      rownames(loadings),
      function(x) {
        return(
          paste(
            str_replace(x, "Zotu", "Microbiome"), name_map[x],
            sep=" -- "
          )
        )
      }
    )
  } else {
    loadings$name <- make.unique(name_map[rownames(loadings)], sep="_")
    # dups <- duplicated(loadings$name)
    # if (any(dups)) {
    #   loadings$name[dups] <- paste0(loadings$name[dups], "__")
    # }
  }

  return(loadings)
}


extract_cv_loadings <- function(plsda.obj, name_map, ncomp=1, ...) {
  if ("block.plsda" %in% class(plsda.obj)) {
    met_loadings <- cv_loadings_per_block(
      plsda.obj, name_map, "Metabolome", ncomp, ...)
    mic_loadings <- cv_loadings_per_block(
      plsda.obj, name_map, "Microbiome", ncomp, ...)

    return(list(Metabolome=met_loadings, Microbiome=mic_loadings))
  } else {
    return(cv_loadings_per_block(plsda.obj, name_map, NULL, ncomp, ...))
  }
}


plot_cv_loadings <- function(cv_loadings) {
  loading_df <- bind_rows(cv_loadings, .id="Iteration")

  max_ <- max(abs(loading_df$importance))
  ggplot(
    loading_df, aes(x=reorder(name, abs(importance)), y=importance,
                    color=GroupContrib)
  ) +
    geom_boxplot() +
    geom_jitter() +
    # TODO: use gghalves
    # geom_half_boxplot() +
    # geom_half_violin(side="l") +
    ylim(-max_, max_) +
    coord_flip() +
    ylab("Feature Importance (sPLS-DA loading)") +
    xlab("") +
    theme_pubr() +
    theme(axis.line.y.left=element_blank(), axis.ticks.y=element_blank())
}


plot_cv_auc <- function(roc_aucs) {
  aucs <- as.data.frame(t(sapply(
    roc_aucs,
    function(x) {
      return(c(x$Metabolome$auc, x$Microbiome$auc))
    }
  )))
  colnames(aucs) <- c("Metabolome", "Microbiome")
  pivot_longer(
    aucs, cols=colnames(aucs), names_to="Modality", values_to="AUC") %>%
      ggplot(aes(x=Modality, y=AUC)) +
      geom_boxplot() +
      ylim(0, 1) +
      theme_pubr() -> auc_plot

  rocs <- lapply(
    roc_aucs,
    function(x) {
      modalities <- c(
        rep("Metabolome", length(x$Metabolome$specificities)),
        rep("Microbiome", length(x$Microbiome$specificities))
      )
      return(
        data.frame(
          Sensitivity=c(x$Metabolome$sensitivities, x$Microbiome$sensitivities),
          Specificity=c(x$Metabolome$specificities, x$Microbiome$specificities),
          Modality=modalities
        )
      )
    }
  )

  sensitivities <- list()
  specificities <- list()
  for (i in seq_along(rocs)) {
    sensitivities[[i]] <- rocs[[i]]$Sensitivity
    specificities[[i]] <- rocs[[i]]$Specificity
  }
  mean_curve <- tryCatch(
    data.frame(
      MeanSensitivity=rowMeans(as.data.frame(sensitivities)),
      MeanSpecificity=rowMeans(as.data.frame(specificities)),
      SDSensitivity=apply(as.data.frame(sensitivities), 1, sd),
      SDSpecificity=apply(as.data.frame(specificities), 1, sd),
      Modality=rocs[[1]]$Modality
    ),
    error=function(err){return(NULL)}
  )

  rocs <-  bind_rows(rocs, .id="Iteration")
  ggplot(rocs, aes(x=1-Sensitivity, y=Specificity, colour=Iteration)) +
    geom_line(size=2) +
    geom_line(
      data=data.frame(
        x=c(0, 1, 0, 1), y=c(0, 1, 0, 1),
        Modality=c("Metabolome", "Metabolome", "Microbiome", "Microbiome")),
      aes(x=x, y=y), color="#D62728", linetype=2, size=1
    ) +
    facet_wrap(.~Modality) -> roc_plot

  return(
    list(auc_plot=auc_plot, roc_plot=roc_plot, mean_curve=mean_curve, aucs=aucs))
}

plot_single_modality_cv_auc <- function(roc_aucs) {
  aucs <- data.frame(
    y=sapply(roc_aucs, function(x) x$auc),
    x=1
  )
  ggplot(aucs, aes(x=x, y=y)) +
    geom_boxplot() +
    xlab("") +
    ylab("AUC") +
    theme_pubr() +
    ylim(0, 1) +
    theme(axis.text.x=element_blank()) -> auc_plot

  rocs <- lapply(
    roc_aucs,
    function(x) {
      return(
        data.frame(
          Sensitivity=x$sensitivities,
          Specificity=x$specificities
        )
      )
    }
  )

  sensitivities <- list()
  specificities <- list()
  for (i in seq_along(rocs)) {
    sensitivities[[i]] <- rocs[[i]]$Sensitivity
    specificities[[i]] <- rocs[[i]]$Specificity
  }
  mean_curve <- data.frame(
    MeanSensitivity=rowMeans(as.data.frame(sensitivities)),
    MeanSpecificity=rowMeans(as.data.frame(specificities)),
    SDSensitivity=apply(as.data.frame(sensitivities), 1, sd),
    SDSpecificity=apply(as.data.frame(specificities), 1, sd)
  )

  rocs <-  bind_rows(rocs, .id="Iteration")
  ggplot(rocs, aes(x=1-Sensitivity, y=Specificity, colour=Iteration)) +
    geom_line(size=2) +
    geom_line(
      data=data.frame(x=c(0, 1), y=c(0, 1)), aes(x=x, y=y), color="#D62728",
      linetype=2, size=1
    ) -> roc_plot

  return(
    list(auc_plot=auc_plot, roc_plot=roc_plot, mean_curve=mean_curve, aucs=aucs))
}
