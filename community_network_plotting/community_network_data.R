library(argparse)

parser <- ArgumentParser(description="Plotting total community network")

# files and paths
parser$add_argument(
  "input_path", type="character",
  help="Path to the multi-omics downstream results"
)
parser$add_argument(
  "result_path", type="character",
  help="Path to the folder where the results should be saved"
)

# to use it from the R console directly
pargs <- tryCatch(
  parser$parse_args(),
  error=function(err) {
    warning("Using console args")
    args <- list(
      input_path="MGSResults/class",
      result_path="community_network_plotting/class"
    )
  }
)

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

load(paste0(pargs$input_path, "/", "processed_data.RData"))
load(paste0(pargs$input_path, "/", "een_posteen_communities.RData"))

# feature name mapping
write.csv(first_name_map, file=paste0(pargs$result_path, "/", "name_map.csv"))
# saving graph to python-readable format
igraph::write.graph(
  net$network,
  file=paste0(pargs$result_path, "/", "correlation_network.graphml"),
  format="graphml"
)
# correlation network
write.csv(
  net$correlations$Correlation,
  file=paste0(pargs$result_path, "/", "corr_mat.csv")
)
# feature community assignments
write.csv(
  communities$membership,
  file=paste0(pargs$result_path, "/", "community_partitions.csv"),
  row.names=communities$names
)

