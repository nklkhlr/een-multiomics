library(argparse)


extract_taxonomy <- function(mic_data) {
  rownames(mic_data) <- sapply(
    rownames(mic_data),
    function(x) {
      tail(strsplit(x, "\\|")[[1]], n=1)
    }
  )
  return(mic_data)
}


tax_level_from_file <- function(file_path) {
  # TODO: does this work with windows paths too?
  paths <- strsplit(file_path, "/")[[1]]
  return (
    strsplit(
      strsplit(paths[length(paths)], "_")[[1]][2], ".tab")[[1]][1]
  )
}


format_mgs_names <- function(ids) {
  add_zeros <- function(id_) {
    n <- 6 - nchar(id_)
    id_split <- strsplit(id_, "_")[[1]]
    return(
      paste0(id_split[1], "_", paste(rep("0", n), collapse=""), id_split[2]))
  }
  return(sapply(ids, add_zeros))
}


format_mgsmeta_names <- function(ids) {
  sapply(
    ids,
    function(x) {
      if (x < 10) {
        return(paste0("CF_00", x))
      }
      if (x < 100) {
        return(paste0("CF_0", x))
      }
      return(paste0("CF_", x))
    }
  )
}


create_param_ranges <- function(max_, min_, step_) {
  if (step_ < 1) {
    stop("Step size must be positive!")
  }
  if (min_ > max_) {
    return(seq(max_, min_, step_))
  }
  return(seq(min_, max_, step_))
}


parser <- ArgumentParser(description="multi-omics EEN experiments")

# files and paths
parser$add_argument(
  "microbiome_file", type="character",
  help="Path to the microbiome measurements"
)
parser$add_argument(
  "result_path", type="character",
  help="Path to the folder where the results should be saved"
)

## search parameters
# metabolome
parser$add_argument(
  "--min-metabolome", type="integer", default=170,
  help="Minimum number of metabolite features to retain in feature selection"
)
parser$add_argument(
  "--max-metabolome", type="integer", default=250,
  help="Maximum number of metabolite features to retain in feature selection"
)
parser$add_argument(
  "--metabolome-step-size", type="integer", default=10,
  help="Step size when creating the search range for feature selction"
)

# microbiome
parser$add_argument(
  "--min-microbiome", type="integer", default=170,
  help="Minimum number of microbial features to retain in feature selection"
)
parser$add_argument(
  "--max-microbiome", type="integer", default=250,
  help="Maximum number of microbial features to retain in feature selection"
)
parser$add_argument(
  "--microbiome-step-size", type="integer", default=10,
  help="Step size when creating the search range for feature selction"
)

# skipping feature selection and passing (pre-selected) dimensions
parser$add_argument(
  "--feature-dimensions", type="integer", nargs=4,
  help=paste0(
    "Number dimensions to use when training the sPLS-DA model. ",
    "When given, optimal feature selection will be skipped and the given ",
    "parameter settings are used. The form is Met1, Met2, Mic1, Mic2. ",
    "e.g. `--feature-dimensions 10 20 30 40` means 10 metabolites for the ",
    "first and 20 the second latent dimensions. For microbiome data 30 and 40 ",
    "respectively."
  )
)


# number of cpus to use
parser$add_argument(
  "--cpus", type="integer", default=1, help="Number of CPUs to use")
