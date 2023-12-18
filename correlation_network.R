library(igraph)

# compute all pairwise spearmans rank coefficients, multipletes correct their p-values
# and set all correlations below a given threshold (default .05) to zero
correlations <- function(data, threshold=.05, method='spearman', adjust='BH') {
  cors <- psych::corr.test(data, method=method, adjust=adjust)
  cor_vals <- cors$r
  cor_vals[upper.tri(cor_vals, diag=FALSE)][cors$p.adj > threshold] <- 0
  cor_vals[lower.tri(cor_vals, diag=FALSE)] <- t(cor_vals)[lower.tri(cor_vals, diag=FALSE)]
  return(
    list(
      FDR=cors$p.adj,
      Correlation=cor_vals
    )
  )
}


#' Generate a correlation network either from the data directly or
#' from a precomputed correlation matrix
#'
#' @param data data.frame containing the conacatenated metabolome and microbiome data
#' @param corrs data.frame with all pairwise feature correlations
#' @param metabo_micro_only Whether to remove metabolite-metabolite and microbe-microbe
#'        edges
#' @param threshold corrected p-value threshold to set correaltions to 0
#' @param adjust String indicating which method to use for multiple test correction
correlation_network <- function(data=NULL, corrs=NULL,
                                metabo_micro_only=FALSE,
                                threshold=.05, adjust='BH') {
  # optional correlation calculations
  if (is.null(corrs)) {
    message('Computing correlations')
    if (is.null(data)) stop("'data' has to be provided, when 'corrs' is NULL")
    corrs <- correlations(data, threshold=threshold, adjust=adjust)
    diag(corrs$Correlation) <- 0
  }
  message('Computing network\n')
  # Generate the network. Correlations have to be absolute as igraph
  # treats them as edge weights and negatively weighted edges are removed
  net <- igraph::graph_from_adjacency_matrix(abs(corrs$Correlation),
                                             weighted=TRUE) %>%
    # Remove like self-loops
    simplify() %>%
    # Make undirected
    as.undirected()
  if (metabo_micro_only) {
    # remove all metabolite-metabolite and microbe-microbe edges
    edges_to_remove <- apply(
      ends(net, E(net)),
      1,
      function(edge) {
        src <- startsWith(edge[1], 'ID') | startsWith(edge[1], 'FT')
        tgt <- startsWith(edge[2], 'ID') | startsWith(edge[2], 'FT')
        return(src == tgt)
      }
    )
    net <- delete.edges(net, which(edges_to_remove))
  }

  message('Network attributes\n')
  # colour nodes by whether they are metabolites or zOTUs
  node_cmap <- c(Metabolite="#377EB8", Microbiome="#E41A1C")
  node_colour <- sapply(1:length(V(net)),
                        function(i) {
                          if (startsWith(V(net)[i]$name, 'ID')) return(node_cmap['Metabolite'])
                          else return(node_cmap['Microbiome'])
                        })
  names(node_colour) <- V(net)$name
  V(net)$colour <- node_colour

  # colour edges by correlation coefficient
  edge_weight <- apply(
    ends(net, E(net)),
    1,
    function(edge) {
      return(corrs$Correlation[edge[1], edge[2]])
    }
  )
  # color map for correlations
  edge_cmap <- colorRampPalette(colors = c("red", "lightgrey", "blue"))
  breaks <- cut(edge_weight, 100)
  edge_colour <- edge_cmap(100)[breaks]
  E(net)$colour <- edge_colour

  return(
    list(
      correlations=corrs,
      network=net,
      edge_cmap=edge_cmap(100),
      node_cmap=node_cmap
    )
  )
}

# plot a colorbar => just a helper for plotting legends
colour_bar <- function(
  palette, min, max=-min, nticks=5,
  ticks=seq(min, max, len=nticks),
  new_dev=TRUE, title=''
) {
  scale = (length(palette)-1)/(max-min)
  if (new_dev) dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n',
       xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)

  for (i in 1:(length(palette) - 1)) {
    y = (i-1)/scale + min
    rect(0, y, 10, y+1/scale, col=palette[i], border=NA)
  }
}

#' Plotting a correlation network
#'
#' @param network igraph object representing the correlation network with
#'        all required attributes
#' @param edge_cmap colormap for edges
#' @param node_cmap colormap for nodes
#' @param vertex.label named vector to rename nodes
#' @param node_shapes named vector giving the shape for each node
#' @param layout igrpah layout function to use for plotting
plot_corr_net <- function(network, edge_cmap, node_cmap, vertex.label=NULL,
                          node_shapes=NULL, layout=layout.kamada.kawai, ...) {
  zero_degrees <- which(degree(network) == 0)
  network <- delete.vertices(network, zero_degrees)
  if (!is.null(node_shapes)) {
    node_shapes <- node_shapes[V(network)$name]
  }
  if (!is.null(vertex.label)) {
    vertex.label <- vertex.label[V(network)$name]
  }

  if (typeof(layout) == "closure") {
    layout <- layout(network)
  }

  plot(
    network, vertex.color=V(network)$colour,
    edge.color=E(network)$colour,
    layout=layout,
    vertex.size=10, vertex.label.cex=.7,
    vertex.shape=node_shapes,
    vertex.label=vertex.label,
    vertex.label.color="black",
    ...
  )
  # legend(
  #   "topright", legend=names(node_cmap),
  #   fill=node_cmap, box.lty=0, cex=1
  # )
  # colour_bar(
  #   edge_cmap, -1, new_dev=FALSE, title='Correlation'
  # )
}
