library(mixOmics)

# TODO: make relative to current file
source('PLSDA/DIABLO.R')


#' Scale data to z-scores
#'
#' Scale features either based on their own mean and standard
#' deviation or based on previously computed mean and SD
#' (for test sets)
#'
#' @param data data.frame with samples in rows and features
#'        in columns
#' @param means vector, option to pass previously computed means
#' @param sds vector, option to pass previously computed SDs
scale_zscore <- function(data, means=NULL, sds=NULL) {
    return_params <- FALSE
    if (!is.null(means)) {
        if (is.null(sds)) {
            # either both should be given or neither
            stop("'sds' must be given, when 'means' is given")
        }
    } else {
        # compute parameters if not given
        means <- colMeans(data, na.rm=TRUE)
        sds <- apply(data, 2, sd, na.rm=TRUE)
        return_params <- TRUE
    }

    # apply scaling feature wise
    scaled <- sapply(
      seq_len(ncol(data)),
      function(i) (data[,i] - means[i]) / sds[i]
    )
    # restore names
    colnames(scaled) <- colnames(data)

    # if parameters are computed they will also be returned
    if (!return_params) return(scaled)
    return(
      list(
        scaled=scaled,
        means=means,
        sds=sds
      )
    )
}


#' Tune a sPLSDA model on a given data set with random train/test splitting
#'
#' The function first splits samples into train and test samples (balanced).
#' Subsequently, the training data is used to evaluate the optimal number of
#' features for each modality and latent dimension as well as the optimal
#' number of latent dimensions. The best settings are then used to train a
#' model that gets returned
#'
#' @param samples vector of sample names
#' @param target vector of sample lables
#' @param file filepath to save results to. Can be empty if `save_` is FALSE
#' @param metabolome_data microbiome data with *features in columns*
#' @param mic_data microbiome data with *features in rows*
#' @param split fraction of training samples
#' @param return_ Whether to return the results
#' @param save_ Whether to save the results. Note that results will be lost
#'        if `save_` and `return_` are both `FALSE`
#' @param ... optional arguments passed to tune.run.block.splsda (see DIABLO.R)
#'            and from there to mixOmics::block.splsda and mixOmics::tune.block.splsda
run_splsda <- function(
    samples, target, file,
    metabolome_data, mic_data, split=.7,
    return_=FALSE, save_=TRUE, ...
) {
    combined_data <- list(Metabolome=metabolome_data[samples,],
                          Microbiome=t(mic_data[,samples]))

    split <- balanced_train_test_split(target, .7)
    # NOTE: for some reason I cannot comprehend R sometimes throws an out-of-bounds
    #       error here, even though all samples in split$train are rownames of
    #       both matrices
    train_data <- lapply(
      combined_data, function(x) {
          idxs <- which(rownames(x) %in% split$train)
          x[idxs,]
      }
    )

    # NOTE: we are not scaling Microbiome data
    scale_res <- scale_zscore(train_data$Metabolome)
    train_data$Metabolome <- scale_res$scaled
    # we store the parameters to apply the same scaling to the test set later
    split$means <- scale_res$means
    split$sds <- scale_res$sds

    # this is just a sanity check since the above is not necessarily guaranteed
    # to have the same samples and order
    assertthat::are_equal(
      rownames(combined_data$Microbiome), rownames(combined_data$Metabolome))

    # perform the actual hyperparameter optimization and training on the optimal
    splsda <- tune.run.block.splsda(keepXlist=X.list,
                                    X=train_data,
                                    Y=target[split$train],
                                    design=design,
                                    ncomp=args$ncomp,
                                    folds=args$folds,
                                    cpus=args$cpus,
                                    nrepeat=3,
                                    verbose=FALSE,
                                    ...)
    perf <- perf(splsda, folds=args$folds, nrepeat=3)

    if (save_) save(splsda, perf, split, file=file)
    if (return_) {
        return(
          list(model=splsda, perf=perf, split=split)
        )
    }
}

#' Turn a data.frame column into a named vector
column_to_vector <- function(df, column) {
    vec <- df[,column]
    names(vec) <- rownames(df)
    return(vec)
}

