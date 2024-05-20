library(mixOmics)
library(future.apply)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(caret)
library(pROC)

Rcpp::sourceCpp("PLSDA/roc_curve.cpp")


# compute the AUC given a ROC curve
auc <- function(roc_curve) {
  height <- (roc_curve[,2][-1] + roc_curve[,2][-nrow(roc_curve)])/2
  width <- diff(roc_curve[,1])
  return(sum(height * width))
}


# =============== #
#### Functions ####
# =============== #

#' Tune the parameters for a mixOmics sPLSDA model given some training data
#' and fit a model to the data with the optimal parameter set.
#'
#' Parameter evaluation is done using a weighted vote over all modalities
#' using the centroid distance of the balanced error rate. The final model
#' is trained on the given data and returned.
#'
#' @param ... Data and parameters to pass to mixOmics::block.splsda and
#'            mixOmics::tune.block.splsda
#' @param keepXlist list containing the number of latent dimensions to test
#'        for each modality
#' @param cpus Number of threads to use
#' @param validation Type of cross validation to use. See mixOmics::perf
#'        documentation for available options
#' @param folds Number of folds to use in the cross validation
#' @param nrepeat Number of reapeats of the cross validation procedure
#' @param verbose Whether to be verbose about progress
tune.run.block.splsda <- function(..., keepXlist,
                                  cpus=1,
                                  validation="Mfold",
                                  folds=5, nrepeat=50,
                                  verbose=FALSE) {
  if (verbose) cat("Computing initial splsda model...\n\n")
  ### first checking the best ncomp choice
  # train initial model on the maximum number of components
  first.splsda <- block.splsda(...)

  if (verbose) cat("Estimating performance for ncomp choice...\n\n")
  ncomp <- tryCatch(
    {
      # evaluate the performance of each dimension and return the
      # optimal number
      first.perf <- perf(first.splsda,
                         validation=validation,
                         folds=folds, nrepeat=nrepeat,
                         progressBar=verbose,
                         cpus=cpus)
      # TODO: make this dynamic wrt. to ncomp.choice
      first.perf$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
    },
    error=function(error) {
      # I haven't gotten to the bottom of why mixOmics sometimes fails here
      # you might want to double check if this happens
      warning(
        paste0("ncomp choice failed with the following error\n",
               error, "Defaulting to using ncomp=2"
        )
      )
      # default number of dimensions: 2
      2
    }
  )
  if (is.null(ncomp)) stop("Input parameters for cross-validation led to no result for ncomp")
  # minimum require for further analysis
  if (ncomp < 2) {
    warning("ncomp < 2 found - defaulting to ncomp=2")
    ncomp <- 2
  }
  if (verbose) cat("\nncomp choice: ", ncomp, "\n")
  # setup parameters for
  args <- list(...)
  args$ncomp <- ncomp
  args$test.keepX <- keepXlist
  args$progressBar <- verbose
  args$folds <- folds
  args$nrepeat <- nrepeat
  # the way mixOmics does parallelization is not compatible with all version combinations
  # feel free to remove the comment and test if it works for you - make things much
  # quicker if you have a lot of CPUs available
  # args$BPPARAM <- cpus

  # Parameter tuning
  if (verbose) cat("\nTuning keepX...\n\n")
  tune <- do.call(tune.block.splsda, args=args)

  args$keepX <- tune$choice.keepX
  args$progressBar <- NULL
  args$test.keepX <- NULL
  args$nrepeat <- NULL
  args$folds <- NULL
  # args$BPPARAM <- NULL

  # Compute final model final model
  if (verbose) cat("\nComputing final model...\n\n")
  final.model <- do.call(block.splsda, args=args)

  return(final.model)
}


#' Compute label-randomized evaluations
#'
#' Models are trained on a fraction of the total samples (
#' see `ptrain_`) and then evaluates the performance on the
#' test set. All labels are shuffled randomly in each iteration.
randomised_plsda <- function(data, labels, cpus=1,
                             plsda.fun=block.splsda,
                             ncomp=2, n=100,
                             consensus=TRUE,
                             ptrain_=.7, ...) {
  if (cpus > 1) {
    plan(multisession, workers=cpus)
    future_lapply(1:n,
                  function(i){
                    rand.labels <- sample(labels)
                    # names(rand.labels) <- NULL
                    split_ <- train_test_split(names(rand.labels), ptrain_)
                    train_data_ <- lapply(
                      data, function(x_) t(x_[,split_$train])
                    )
                    rand.plsda <- plsda.fun(
                      X=train_data_, Y=rand.labels[split_$train],
                      ncomp=ncomp, ...
                    )
                    test_data_ <- lapply(
                      data, function(x_) t(x_[,split_$test])
                    )
                    if (consensus) {
                      pred <- predict(rand.plsda, test_data_)
                      return(
                        consensus_auc(y=rand.labels[split_$test],
                                      prediction=pred, dims=1:ncomp)
                      )
                    }
                    return(auroc(rand.plsda, newdata=test_data_,
                                 outcome.test=rand.labels[split_$test],
                                 plot=FALSE, print=FALSE))
                  })
  } else {
    lapply(1:n,
           function(i){
             rand.labels <- sample(labels)
             # names(rand.labels) <- NULL
             split_ <- train_test_split(names(rand.labels), ptrain_)
             train_data_ <- lapply(
               data, function(x_) t(x_[,split_$train])
             )
             rand.plsda <- plsda.fun(
               train_data_, rand.labels[split_$train],
               ncomp=ncomp, ...
             )
             test_data_ <- lapply(
               data, function(x_) t(x_[,split_$test])
             )
             if (consensus) {
               pred <- predict(rand.plsda, test_data_)
               return(
                 consensus_auc(y=rand.labels[split_$test],
                               prediction=pred, dims=1:ncomp)
               )
             }
             return(auroc(rand.plsda, newdata=test_data_,
                          outcome.test=rand.labels[split_$test],
                          plot=FALSE, print=FALSE))
           })
  }
}

# TODO: consensus roc
single.block.aucs <- function(true.auc, rand.auc, block.name,
                              groups, n, consensus=TRUE) {
  if (is.null(n.comps)) n.comps <- length(true.auc[[block.name]])
  groups <- rownames(true.auc[[block.name]][[1]])
  n.groups <- length(groups)
  # arrange randomised AUC values
  if (consensus) {
    browser()
    auc.df <- data.frame(t(data.frame(rand.auc[[block.name]])))
    rownames(auc.df) <- 1:nrow(auc.df)
    colnames(auc.df) <- paste("Comp", 1:n.comps, sep="")
    total.rands <- tidyr::pivot_longer(auc.df, colnames(auc.df),
                                       names_to="Component",
                                       values_to="AUC")
  } else {
    auc.df <- lapply(1:n.comps,
                     function(i) {
                       df <- sapply(1:n.groups,
                                    function(j) {
                                      rand.auc.vals <- sapply(1:n,
                                                              function(k) rand.auc[[k]][[block.name]][[i]][j,1])
                                    })
                       df <- as.data.frame(df)
                       colnames(df) <- groups
                       return(df)
                     })
    names(auc.df) <- paste("Comp", 1:n.comps, sep="")
    total.rands <- bind_rows(auc.df, .id = "Component")
  }

  # arrange actual AUC values
  if (length(unique(groups)) > 2) { # Multi-class comparison
    total.true <- as.data.frame(sapply(1:n.comps, function(i) true.auc[[block.name]][[i]][,1]))
    colnames(total.true) <- paste("Comp", 1:n.comps, sep="")
    total.true <- add_column(total.true, Group=groups) %>%
      pivot_longer(-Group, names_to="Component",
                   values_to="AUC")
    # compute empirical p-values
    emp.p.vals <- sapply(groups,
                         function(g){
                           sapply(names(auc.df),
                                  function(c){
                                    true <- filter(total.true, Group == g,
                                                   Component == c) %>% .$AUC
                                    mean(auc.df[[c]][[g]] >= true)
                                  })
                         }) %>% as.data.frame
    emp.p.vals$Component <- rownames(emp.p.vals)
    emp.p.vals <-  pivot_longer(emp.p.vals, -Component,
                                values_to="p_value", names_to="Group") %>%
      mutate(p_value=paste("p", formatC(.$p_value, format="f", digits=2), sep=" = "))

    # plotting
    pivot_longer(as.data.frame(total.rands), -Component,
                 names_to="Group", values_to="AUC") %>%
      ggplot(aes(x=Component, y=AUC)) +
      geom_boxplot() +
      geom_point(data=total.true,
                 colour="#D62728", shape=18, size=5) +
      geom_text(data=emp.p.vals, aes(y=1.05, label=p_value),
                fontface="bold") +
      facet_wrap(Group~.) +
      theme_pubr() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) -> p
  } else { # binary class-comparison
    true_aucs <- sapply(1:n.comps, function(i) true.auc[[block.name]][[i]][1])
    true <- data.frame(Component=paste("Comp", 1:n.comps, sep=""),
                       AUC=true_aucs)
    rownames(true) <- true$Component
    # compute empirical p-values
    emp.pvs <- sapply(names(auc.df), function(x) mean(auc.df[[x]] >= true[x,"AUC"]))
    emp.p.vals <- data.frame(Component=names(emp.pvs),
                             p_value=paste("p", emp.pvs, sep=" = "))

    # plotting
    colnames(total.rands) <- c("Component", "AUC")
    ggplot(total.rands, aes(x=Component, y=AUC)) +
      geom_boxplot() +
      geom_point(data=true,
                 colour="#D62728", shape=18, size=5) +
      geom_text(data=emp.p.vals, aes(y=1.05, label=p_value),
                fontface="bold") +
      theme_pubr() +
      theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5)) -> p

    total.true <- true
  }

  return(list(rand.data=total.rands,
              true.data=total.true,
              p.vals=emp.p.vals,
              plot=p))
}

consensus_auc <- function(y, prediction, dims=c(1, 2), step_size=.001) {
  return(
    sapply(
      dims,
      function(x) {
        roc_ <- roc(
          predictor=prediction$AveragedPredict[,,x][,2],
          response=as.integer(as.factor(y)) - 1
        )

        return(roc_$auc)
      }
    )
  )
}

eval.rand.block <- function(plsda.obj, X, y,
                            Xtest, ytest,
                            n=100, cpus=1,
                            plsda.fun=block.splsda,
                            consensus=TRUE, ...) {
  if (is.null(n.comp)) ncomp <- plsda.obj$ncomp
  else ncomp <- n.comp

  if (is.list(ncomp) | is.vector(ncomp)) ncomp <- ncomp[1]

  rand.auc <- randomised_plsda(X, y,
                               plsda.fun=plsda.fun,
                               ncomp=ncomp,
                               n, # keepX=plsda.obj$keepX,
                               design=plsda.obj$design,
                               scheme=plsda.obj$scheme,
                               scale=plsda.obj$scale,
                               # mode=plsda.obj$mode,
                               cpus=cpus,
                               consensus=consensus,
                               ...)
  # browser()
  if (consensus) {
    pred <- predict(plsda.obj, Xtest)$AveragedPredict
    true.auc <- lapply(
      1:ncomp,
      function(x) {
        roc_ <- roc(
          predictor=pred[,,x][,2],
          response=as.integer(as.factor(ytest)) - 1,
        )
        return(roc_$auc)
      }
    )
  } else {
    true.auc <- auroc(plsda.obj, newdata=Xtest, outcome.test=ytest,
                      plot=FALSE, print=FALSE)
  }

  if (consensus) {
    return(
      list(
        consensus=single.block.aucs(
          list(consensus=true.auc),
          list(consensus=rand.auc),
          "consensus", groups, n,
          consensus=consensus
        )
      )
    )
  }
  return(lapply(blocks,
                function(block) {
                  single.block.aucs(true.auc, rand.auc,
                                    block, groups, n)
                }
  )
  )
}


balanced_train_test_split <- function(groups, ptrain=.7, times=1) {
  if (times > 1) return(createDataPartition(groups, p=ptrain, times=times, ...))
  train_idxs <- createDataPartition(groups, p=ptrain, times=1, list=FALSE)
  return(
    list(
      train=names(groups)[train_idxs],
      test=names(groups)[-train_idxs]
    )
  )
}

# train-test splitting without balancing
train_test_split <- function(samples, ptrain=.7) {
  n <- length(samples)
  train_idxs <- sample(1:n, ptrain * n)
  return(
    list(train=samples[train_idxs],
         test=samples[-train_idxs])
  )
}


# ======================= #
#### Setting variables ####
# ======================= #
# NOTE: these are the default settings
args <- list(folds=4, cpus=4,
             nrepeat=100,
             ncomp=7)

# overriding defaults by command line parameters
if (sys.nframe() == 0L) {
  o.args <- commandArgs(trailingOnly=TRUE)
  if (length(o.args) > 1) {
    arg.splits <- sapply(o.args, function(x) strsplit(x, "=")[[1]])
    in.args <- as.integer(arg.splits[2,])
    names(in.args) <- arg.splits[1,]

    for (arg.name in names(in.args)) {
      if (arg.name %in% names(args)) {
        args[[arg.name]] <- in.args[[arg.name]]
      }
    }
  }
}

# ========== #
### Design ###
# ========== #
# starting with an expected correlation of .1
# => maybe we should reconsider/test this!
design <- matrix(.4, ncol=2, nrow=2,
                 dimnames=list(c("Metabolome", "Microbiome"),
                               c("Metabolome", "Microbiome")))
# TODO: correlations?
diag(design) <- 0
# design <- rbind(design, Y=c(.1,.1))
# design <- cbind(design, Y=c(.1,.1,0))

### Setting variables for optimisation
X.list <- list(Metabolome=seq(200, 210, 10),
               Microbiome=seq(90, 100, 10))
