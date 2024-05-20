library(readxl)
library(dplyr)
library(pROC)
library(pbapply)
library(igraph)

source("PLSDA/diablo_from_list.R")
source("PLSDA/plsda_utils.R")
source("een_diversity.R")
source("correlation_network.R")

load("Data/metabolomics_data.RData")
load("Data/metabolomics_meta.RData")


set.seed(123)


# =============== #
#### Functions ####
# =============== #
subset_by_splsda_ <- function(splsda_obj, met, mic) {
  met_feats <- rownames(splsda_obj$loadings$Metabolome)[
    rowSums(splsda$loadings$Metabolome) > 0]
  mic_feats <- rownames(splsda_obj$loadings$zOTU)[
    rowSums(splsda_obj$loadings$zOTU) > 0]
  sample_names <- intersect(colnames(met), colnames(mic))

  return(
    list(
      Metabolome=met[met_feats, sample_names],
      zOTU=mic[mic_feats, sample_names]
    )
  )
}


#' Plot average roc curves and mark the minimum and maximum
#' interval.
#'
#' @param roc_data List containg the mean, min and max values for specificity and sensitivity
plot_average_roc <- function(roc_data) {
  roc_data$MinSpecificity <- pmax(
    roc_data$MeanSpecificity-roc_data$SDSpecificity, 0)
  roc_data$MaxSpecificity <- pmin(
    roc_data$MeanSpecificity+roc_data$SDSpecificity, 1)
  cv_averaged <- ggplot(
    roc_data, aes(x=1-MeanSensitivity, y=MeanSpecificity)) +
    geom_ribbon(
      aes(ymin=MinSpecificity, ymax=MaxSpecificity),
      alpha=.2, color="lightgrey") +
    geom_line(color="#1F77B4", size=2) +
    geom_line(
      data=data.frame(
        x=c(0, 1, 0, 1), y=c(0, 1, 0, 1),
        Modality=c("Metabolome", "Metabolome", "Microbiome", "Microbiome")),
      aes(x=x, y=y), color="#D62728", linetype=2, size=1
    )

  return(cv_averaged)
}


#' Extract the feature information of a community
#'
#' Given a subgraph of the correlation network, metabolite and Microbiome
#' metadata extract certain metadata fields, associations and
#' neighbors in the network
#'
#' @param sugraph network community (igraph graph)
#' @param metabo_info data.frame containing the metabolite metdata
#' @param mic_info data.frame containing the Microbiome metdata
#' @param associations data.frame containing the assocation annotations in a column named "AssociatedWith"
#' @param file_base basic file path to which extracted data will be saved. This should not contain a file ending - "_metabolite.csv" and "_Microbiome.csv" are appended
#' @param save whether to write the data to .csv files just return them
extract_community_features <- function(
  subgraph, metabo_info, mic_info, associations, file_base, save=TRUE
) {
  comm <- V(subgraph)$name

  metabos <- comm[startsWith(comm, "FT")]
  minfo_sub <- metabo_info[metabos,]
  minfo_sub$AssociatedWith <- associations[metabos, "AssociatedWith"]
  rownames(minfo_sub) <- metabos

  mics <- setdiff(comm, metabos)
  mic_sub <- data.frame(
    Microbe=mic_info[mics],
    AssociatedWith=associations[mics, "AssociatedWith"]
  )
  rownames(mic_sub) <- mics

  minfo_sub$Neighbours <- sapply(
    rownames(minfo_sub),
    function(x) {
      neighbs <- neighbors(subgraph, x)
      paste(mic_sub[neighbs$name, "Microbe"], collapse=" || ")
    }
  )
  mic_sub$Neighbours <- sapply(
    rownames(mic_sub),
    function(x) {
      neighbs <- neighbors(subgraph, x)
      paste(minfo_sub[neighbs$name, "annot_ms1"], collapse=" || ")
    }
  )

  if (save) {
    if (length(metabos) > 0) {
      write.csv(
        apply(
          minfo_sub[,c("annot_ms1", "AssociatedWith", "Neighbours")], 2,
          as.character
        ),
        file=paste0(file_base, "_metabolites.csv"),
        row.names=metabos
      )
    }

    if (length(mics) > 0) {
      write.csv(
        apply(mic_sub, 2, as.character),
        file=paste0(file_base, "_Microbiomes.csv")#, row.names=mics
      )
    }
  }

  return(list(metabolites=metabos, mics=mics))
}


#' Compute the ROC curves on the test set n times for a single modality
#'
#' @param data data.frame containing the data to split and train on
#' @param targets vector of sample group annotations
#' @param n_splits number of iterations to perform
#' @param scale_ whether to scale the data. Should be TRUE (default) if metabolomics data is given
cv_roc_single_modality <- function(data, targets, n_splits=10, scale_=TRUE) {
  splits <- createDataPartition(targets, p=.7, times=n_splits)
  lapply(
    splits, function(k_split) {
      train_samples <- names(targets)[k_split]
      test_samples <- names(targets)[-k_split]

      if (scale_) {
        scale_res <- scale_zscore(t(data[,train_samples]))
        k_train_means <- scale_res$means
        k_train_sds <- scale_res$sds
      } else {
        scale_res <- t(data[,train_samples])
      }

      k_model <- mixOmics::plsda(
        X=scale_res$scaled, Y=targets[train_samples], scale=FALSE)

      if (scale_) {
        k_test_data <- scale_zscore(
          t(data[,test_samples]), k_train_means, k_train_sds)
      } else {
        k_test_data <- t(data[,test_samples])
      }

      preds <- predict(k_model, newdata=k_test_data)$predict
      return(
        list(
          roc=roc(predictor=preds[,,1][,1], response=targets[test_samples]),
          model=k_model
        )
      )
    }
  )
}

#' Compute the ROC curves on the test set n times for multi-omics data
#'
#' @param data list containing the data.frames which caary the data to split and train on
#' @param targets vector of sample group annotations
#' @param n_splits number of iterations to perform
cv_roc_evaluation <- function(met_data, mic_data, targets, n_splits=10) {
  splits <- createDataPartition(targets, p=.7, times=n_splits)
  lapply(
    splits, function(k_split) {
      train_samples <- names(targets)[k_split]
      test_samples <- names(targets)[-k_split]

      scale_res <- scale_zscore(t(met_data[,train_samples]))
      k_train_means <- scale_res$means
      k_train_sds <- scale_res$sds

      k_train_data <- list(
        Metabolome=scale_res$scaled,
        Microbiome=t(mic_data[,train_samples])
      )
      k_model <- block.plsda(
        X=k_train_data, Y=targets[train_samples], design=design, scale=FALSE
      )

      test_metabolome <- scale_zscore(
        t(met_data[,test_samples]), k_train_means, k_train_sds)
      k_test_data <- list(
        Metabolome=test_metabolome,
        Microbiome=t(mic_data[,test_samples])
      )

      yhat <- lapply(
        predict(k_model, newdata=k_test_data)$predict,
        function(preds) {
          # we use the first latent dimension only
          return(
            list(
              roc=roc(predictor=preds[,,1][,1], response=targets[test_samples]),
              model=k_model
            )
          )
        }
      )

      return(yhat)
    }
  )
}


#' Function wrapping the two above defined functions
#'
#' @param metabolites vector of metabolomics features that should be part of the model
#' @param microbes vector of metabolomics features that should be part of the model
#' @param met_data metabolite data (rows metabolites, columns samples)
#' @param mic_data microbial data (rows Microbes, columns samples)
#' @param targets vector of sample group annotations
#' @param n_splits number of iterations to perform
roc_test_community_features <- function(
  metabolites, microbes, met_data, mic_data, targets, n_splits=10
) {
  single_mod <- TRUE
  if (length(metabolites) < 2) {
    roc_data <- cv_roc_single_modality(
      mic_data[microbes,], targets, n_splits, FALSE)
  } else if (length(microbes) < 2) {
    roc_data <- cv_roc_single_modality(
      met_data[metabolites,], targets, n_splits, TRUE)
  } else {
    single_mod <- FALSE
    roc_data <- cv_roc_evaluation(
      met_data[metabolites,], mic_data[microbes,], targets, n_splits)
  }

  if (single_mod) {
    rocs <- lapply(roc_data, function(x) x$roc)
    models <- lapply(roc_data, function(x) x$model)
  } else {
    rocs <- lapply(
      roc_data, function(x) {
        lapply(x, function(y) y$roc)
      }
    )
    models <-  lapply(
      roc_data, function(x) {
        lapply(x, function(y) y$model)
      }
    )
  }
  return(list(roc=rocs, model=models))
}


extract_auc_ <- function(test_results) {
  auc_ <- test_results$roc$Resample1$auc
  if (is.null(auc_)) {
    return(sapply(test_results$roc$Resample1, function(x) x$auc))
  }
  return(auc_)
}


#' Performed random label permutation tests on community features
#' Parameters are the same as for roc_test_community_features except for the
#' ones described below
#'
#' @paaram nrepeats The number of random permutations to perform
#' @param cpus The number of threads to use for parallelization
random_community_test <- function(
  metabolites, microbes, met_data, mic_data, targets, nrepeats=200, cpus=1
) {
  true_auc <- roc_test_community_features(
    metabolites, microbes, met_data, mic_data, targets, 1) %>%
    extract_auc_

  if (cpus == 1) {
    loop_fun <- sapply
  } else {
    plan(multisession, workers=cpus)
    loop_fun <- future_sapply
  }

  rand_aucs <- loop_fun(
    1:nrepeats,
    function(i) {
      rand_labels <- sample(targets)
      roc_test_community_features(
        metabolites, microbes, met_data, mic_data, rand_labels, 1) %>%
        extract_auc_
    }
  )

  if (is.vector(true_auc)) return(rowMeans(true_auc > rand_aucs))

  if (length(metabolites) < 2)
    return(c(Metabolome=NaN, zOTU=mean(true_auc > rand_aucs)))
  return(c(Metabolome=mean(true_auc < rand_aucs), zOTU=NaN))
}


# TOOD: add in statistical tests
#' Plot the 'expression' of all features of a community in boxplots
# TODO: proper documentation
boxplot_community <- function(
  community, met_data, mic_data, groups, metadata, healing_group, tests=NULL,
  name_map=NULL, ...
) {
  mic_mask <- startsWith(community, "Zotu")
  comm_mic <- community[mic_mask]
  comm_met <- community[!mic_mask]

  merged <- bind_rows(
    as.data.frame(met_data[comm_met, names(groups)]),
    as.data.frame(mic_data[comm_mic, names(groups)])
  ) %>%
    t %>%
    as.data.frame %>%
    rownames_to_column("Sample") %>%
    add_column(
      Patient=metadata[.$Sample, "Patient"],
      Diet=groups[.$Sample],
      HealingGroup=healing_group[.$Sample],
      PGA=metadata[.$Sample, "PGA"]
    )

  if (!is.null(name_map)) {
    colnames(merged) <- make.unique(sapply(
      colnames(merged),
      function(col) {
        if (col %in% names(name_map)) {
          return(
            strsplit(
              gsub("REF_", "", name_map[col]),
              ";"
            )[[1]][1]
          )
        }
        return(col)
      }
    ))
  }

  col_na <- is.na(colnames(merged))
  colnames(merged)[col_na] <- make.unique(rep("NA", sum(col_na)), sep="_")

  merged %>%
    pivot_longer(
       -c(Sample, Diet, HealingGroup, PGA, Patient),
       names_to="Feature", values_to="Expression"
    ) %>%
    add_column(
      Modality=ifelse(
        .$Feature %in% name_map[comm_mic], "Microbiome", "Metabolite")) %>%
    ggplot(aes(x=Feature, y=Expression, colour=Diet)) +
      geom_boxplot() +
      geom_jitter(position=position_dodge(.75), size=1) +
      facet_wrap("Modality", scales="free", ...) +
      theme_pubr() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) -> plot

  return(list(plot=plot, data=merged))
}


if (sys.nframe() == 0) {
  source("cli_parser.R")

  pargs <- parser$parse_args()

  # load microbiome data
  # contains:
  # * zotu_data (processed 16S data)
  # * zotu_annotation (taxonomic annotation for zOTUs)
  load(pargs$microbiome_file)
  # NOTE: change this if you want to use annotation on a different level
  mic_annotation <- zotu_annotation$`Genus|Species`
  names(mic_annotation) <- rownames(zotu_annotation)

  # make sure results dir is existing or create it
  if (dir.exists(pargs$result_path)) {
    warning(
      paste0(
        "Result path '", pargs$result_path, "' already exists. If results ",
        "have previously been saved to this directory, they will be ",
        "overwritten."
      )
    )
  } else {
    dir.create(pargs$result_path, recursive=TRUE)
  }

  # create search parameters for feature selection
  X.list <- list(
    Metabolome=create_param_ranges(
      pargs$min_metabolome, pargs$max_metabolome,
      pargs$metabolome_step_size
    ),
    # NOTE: this is called Microbiome so I don't have to make sure all
    #       Microbiome names are replaced downstream
    Microbiome=create_param_ranges(
      pargs$min_microbiome, pargs$max_microbiome,
      pargs$microbiome_step_size
    )
  )
  # set number of cpus
  args$cpus <- pargs$cpus

  # ============ #
  #### PLS-DA ####
  # ============ #
  # reducing metabolites to well-annotated only
  well_ann_mets <- rownames(combined$feature_info)[
    combined$feature_info$annot_ms2 %in% c("3", "4")]
  norm_met_data <- combined$normalised[well_ann_mets,]

  COLOUR_PALETTE <- sapply(pal_d3()(2), substr, start=1, stop=7)
  names(COLOUR_PALETTE) <- c("EEN", "PostEEN")

  # start of the actual integration analysis
  # this file contains the samples Kolja and Debbie selected for analysis
  # we do some formatting to get the metadata annotations into the same
  # format as the previous tables had
  comparison_samples <- read_excel(
    "Data/CD_EEN_Relapse_Remission_Sample_Selection.xlsx",
    sheet="LateEEN_vs_PostEEN"
  ) %>%
    mutate(
      Patient=sapply(
        .$Patient,
        function(x) {
          if (x < 10) {
            return(paste0("00", as.character(x)))
          }
          if (x < 100) {
            return(paste0("0", as.character(x)))
          }
          return(as.character(x))
        }
      ),
      # unify PostEEN naming, this is a little incosistent in the table
      `Diet (=PreEEN, EEN, PostEEN)`=gsub("Post-EEN", "PostEEN",
                                          .$`Diet (=PreEEN, EEN, PostEEN)`)
    ) %>%
    add_column(
      # extract the sample IDs in the metadata table
      MetID=apply(
        ., 1,
        function(row) {
          patient_mask <- metadata$Patient == row['Patient']
          date_mask <- metadata$sampleDate == row['Matches16SMetabolome_sampleDate']
          return(rownames(metadata)[patient_mask & date_mask])
        }
      )
    ) %>%
    column_to_rownames("MetID")

  # extract timepoint annotation (=> PreEEN, EEN, and PostEEN)
  timepoint <- comparison_samples$`Diet (=PreEEN, EEN, PostEEN)`
  names(timepoint) <- rownames(comparison_samples)

  # healing groups are for downstream analyses only
  healing <- comparison_samples$HealingGroup
  names(healing) <- rownames(comparison_samples)

  ### NOTE: these are the results of the code in the if-part
  # set these values if you have already pre-compute the best parameters
  # (line 356) and just want to re-run the downstream analysis
  # keepX <- list(Metabolome=c(170, 170), Microbiome=c(70, 70))
  if (is.null(pargs$feature_dimensions)) keepX <- NULL
  else {
    keepX <- list(
      Metabolome=pargs$feature_dimensions[1:2],
      Microbiome=pargs$feature_dimensions[3:4]
    )
  }
  args$ncomp <- 2

  subsamples <- names(timepoint) %>%
    intersect(colnames(norm_met_data)) %>%
    intersect(colnames(zotu_data))

  save(
    zotu_data, mic_annotation, timepoint, subsamples, healing,
    file=paste0(pargs$result_path, "/", "processed_data.RData")
  )
  stop("stop")
  if (is.null(keepX)) {
    # NOTE: this is now set via CLI args (see lines 317 ff)
    #        unfortunately there is no good way to get an approximation
    #        of what values will work a priori. Eseentially mixOmics
    #        performs a grid search, so n * m models are trained
    # X.list <- list(
    #   Metabolome=seq(170, 250, 10), Microbiome=seq(70, 150, 10))
    # X.list <- list(
    #     Metabolome=c(10), Microbiome=c(10))
    model_file <- paste0(pargs$result_path, "/", "EEN_Post_EEN_selection.RData")
    ### Model computation
    if (!file.exists(model_file)) {
      # This function first selects the optimal parameter combination and then returns
      # the best model. It uses 70% of the data (label-balanced) as training data
      # per default, to adapt this change `split`. The splitting is returned so
      # the test samples can be used later on
      run_splsda(
        subsamples, timepoint, model_file, t(norm_met_data), zotu_data,
        scale=FALSE, split=.7
      )
    }
    load(model_file)
  } else {
    # This means we do not compute the optimal parameter combination, but use a
    # pre-defined one. Note that this does NOT use the model computed in the `if`
    # block and uses a (most likely) different train-test split, so results are
    # likely not the same when re-running
    split_idxs <- createDataPartition(timepoint, p=.7, list=FALSE)
    split <- list(
      train=subsamples[split_idxs],
      test=subsamples[-split_idxs]
    )

    scale_res <- scale_zscore(t(norm_met_data[,split$train]))
    train_means <- scale_res$means
    train_sds <- scale_res$sds

    train_data <- list(
      Metabolome=scale_res$scaled,
      Microbiome=t(zotu_data[,split$train])
    )

    test_metabolome <- scale_zscore(
      t(norm_met_data$normalised[,split$test]), train_means, train_sds)
    test_data <- list(
      Metabolome=test_metabolome,
      Microbiome=t(zotu_data[,split$test])
    )

    splsda <- block.splsda(train_data, timepoint[split$train], scale=FALSE)
    yhat <- predict(splsda, newdata=test_data)
  }

  # Evaluate the reproducibility by repeatedly splitting the data in a 70/30
  # manner. This is NOT a k-fold CV, as we randomly resplit at each iteration.
  # Reported are the ROC-AUC values on the test set for each iteration
  splits <- createDataPartition(timepoint, p=.7, times=10)
  test_auc_cv <- lapply(
    splits, function(k_split) {
      train_samples <- subsamples[k_split]
      test_samples <- subsamples[k_split]

      # scale the train samples
      scale_res <- scale_zscore(t(norm_met_data[,train_samples]))
      k_train_means <- scale_res$means
      k_train_sds <- scale_res$sds

      # generate the training data set
      k_train_data <- list(
        Metabolome=scale_res$scaled,
        Microbiome=t(zotu_data[,train_samples])
      )
      # compute the model
      k_model <- block.splsda(
        X=k_train_data, Y=timepoint[train_samples], keepX=splsda$keepX,
        design=design, scale=FALSE
      )

      # scale test set with parameters from training set
      test_metabolome <- scale_zscore(
        t(norm_met_data[,test_samples]), k_train_means, k_train_sds)
      # generate test data
      k_test_data <- list(
        Metabolome=test_metabolome,
        Microbiome=t(zotu_data[,test_samples])
      )

      # compute ROC curve for test set
      yhat <- lapply(
        predict(k_model, newdata=k_test_data)$predict,
        function(preds) {
          # we use the first latent dimension only
          return(
            roc(predictor=preds[,,1][,1], response=timepoint[test_samples]))
        }
      )

      return(yhat)
    }
  )
  names(test_auc_cv) <- seq_along(test_auc_cv)
  roc_aucs <- plot_cv_auc(test_auc_cv)

  ggsave(
    filename=paste0(pargs$result_path, "/cv_aucs.pdf"), width=16, height=9,
    device="pdf", plot=roc_aucs$auc_plot + theme_pubr()
  )
  ggsave(
    filename=paste0(pargs$result_path, "/cv_rocs.pdf"), width=16, height=9,
    device="pdf", plot=roc_aucs$roc_plot + theme_pubr() + scale_colour_d3()
  )

  # compute the mean test ROC-curve
  roc_aucs$mean_curve$MinSpecificity <- pmax(
    roc_aucs$mean_curve$MeanSpecificity-roc_aucs$mean_curve$SDSpecificity, 0)
  roc_aucs$mean_curve$MaxSpecificity <- pmin(
    roc_aucs$mean_curve$MeanSpecificity+roc_aucs$mean_curve$SDSpecificity, 1)
  cv_averaged <- plot_average_roc(roc_aucs$mean_curve)

  ggsave(
    filename=paste0(pargs$result_path, "/cv_average_roc.pdf"), width=16, height=9,
    device="pdf", plot=cv_averaged + theme_pubr()
  )

  ### Plot the model's latent space
  # split labels => annotates samples by whether they belong to the train or test
  # set
  set <- c(rep("train", length(split$train)), rep("test", length(split$test)))
  names(set) <- c(split$train, split$test)

  # actual plotting
  metabolome_plot <- plot_samples(
    splsda, group.name="Timepoint", size=3, block="Metabolome",
    shape=list("Split", set)
  )
  mic_plot <- plot_samples(
    splsda, group.name="Timepoint", size=3, block="Microbiome",
    shape=list("Split", set)
  )
  combined_plot <- ggarrange(metabolome_plot, mic_plot, common.legend=TRUE)
  ggsave(filename=paste0(pargs$result_path, "/een_posteen_selection.pdf"),
         width=16, height=9, plot=combined_plot)


  # Plot the "feature importances" in the form of loading plots
  name_map <- c(combined$feature_info$annot_ms1, mic_annotation)
  names(name_map) <- c(rownames(combined$feature_info), names(mic_annotation))

  first_name_map <- sapply(name_map, substr, start=1, stop=30)

  loading_plots <- plot_loadings(splsda, first_name_map)
  combined_loading_plot <- ggarrange(
    loading_plots$Metabolome$plot + theme_pubr() + scale_fill_d3("category10"),
    loading_plots$Microbiome$plot + theme_pubr() + scale_fill_d3("category10")
  )

  # computing empirical p-value with label permutations
  feat_subset <- subset_by_splsda(splsda, reduced$z_score, zotu_data$clr)
  rand_eval <- eval.rand.block(
    splsda,
    X=list(Metabolome=feat_subset$Metabolome, zOTU=feat_subset$zOTU),
    y=timepoint,
    Xtest=list(
      Metabolome=t(reduced$z_score[,split$test]),
      zOTU=t(zotu_data$clr[,split$test])
    ),
    ytest=timepoint[split$test],
    n=200, plsda.fun=block.plsda, consensus=TRUE, n.comp=1, cpus=4
  )
  cat("Empirical p-value: ", rand_eval$consensus$p.vals$p_value, "\n")

}
