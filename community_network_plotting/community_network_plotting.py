from typing import Dict, Tuple, Union
import pathlib
import numpy as np
import pandas as pd
import networkx as nx
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt


def get_community_graph_centers(
        g: nx.Graph, communities: Dict[str, int], n_communities: int, scale: float
) -> Dict[int, np.ndarray]:
    """Get the center positions for each community"""
    inter_community_weight = {}
    # generate a fake graph where each node is a community and communities
    # are connected via edges if at least one edge between any pair of nodes
    # between them exists. The weight is simply the number of edges between
    # the community pair
    for (src, tgt) in g.edges:
        src_comm = communities[src]
        tgt_comm = communities[tgt]
        if src_comm != tgt_comm:
            inter_community_weight[(src_comm, tgt_comm)] = \
                inter_community_weight.get((src_comm, tgt_comm), 0) + 1
    g = nx.Graph()
    g.add_nodes_from(range(1, n_communities + 1))
    for edge, weight in inter_community_weight.items():
        g.add_edge(*edge, weight=weight)

    # layout of the helper graph gives the centers of the community subgraphs
    return nx.spring_layout(g, weight="weight", scale=scale)


def layout_by_communities(
        g: nx.Graph, communities: Dict[str, int], inter_scale: float,
        intra_scale: float
) -> Dict[str, np.ndarray]:
    """Compute the graph layout separating communities"""
    # layout of community centers
    n_communities = max(list(communities.values()))
    circular_centers = get_community_graph_centers(
        g, communities, n_communities, inter_scale)
    # layout of nodes within communities
    total_pos = {}
    for i, center in circular_centers.items():
        community = [
            node for node, comm in communities.items() if comm == i]
        subgraph = g.subgraph(community)
        total_pos.update(
            nx.kamada_kawai_layout(subgraph, scale=intra_scale, center=center))
    return total_pos


def plot_convex_hull(
        community_positions: Dict[str, np.ndarray], ax: plt.axis, color: tuple,
        scaling: float = 1.2, alpha: float = .5
):
    """Plot the convex hull around a set of vertices"""
    # TODO: ideally we would have a smooth hull
    # TODO: add in scaling
    comm_subset = {}
    for comm, coords in community_positions.items():
        if coords.shape[0] > 1 and coords.size > 2:
            comm_subset[comm] = coords
    if comm_subset:
        hull = ConvexHull(np.array(list(community_positions.values())))
        ax.fill(
            hull.points[hull.vertices, 0], hull.points[hull.vertices, 1],
            facecolor=(*color[:3], alpha), edgecolor=color
        )


def plot_community_graph(
        g: nx.Graph, communities: Dict[str, int], label_dict: Dict[str, str],
        inter_community_scale: float = 2, intra_community_scale: float = .8,
        cmap: str = "tab10", hull_scaling: float = 1.2, hull_alpha: float = .5,
        plot: Tuple[plt.Figure, plt.axis] = None,
        plot_individual_communities: bool = True,
        community_figure: Union[str, pathlib.Path] = "communities.pdf",
        layout: Dict[str, np.ndarray] = None
):
    """Plot the correlation network with community-based node positioning"""
    if plot is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    else:
        fig, ax = plot
    if plot_individual_communities:
        comm_fig, comm_ax = plt.subplots(figsize=(fig.get_size_inches()))
    # get node positions
    if layout is None:
        layout = layout_by_communities(
            g, communities, inter_community_scale, intra_community_scale)
    # plot convex hulls first to be in the background
    n_communities = max(list(communities.values()))
    for i in range(n_communities):
        comm_layout = {
            node: layout[node] for node, comm in communities.items()
            if comm == i + 1 and layout.get(node) is not None
        }
        plot_convex_hull(
            comm_layout, ax, plt.get_cmap(cmap)(i), scaling=hull_scaling,
            alpha=hull_alpha
        )
        if plot_individual_communities:
            plot_convex_hull(
                comm_layout, comm_ax, plt.get_cmap(cmap)(i),
                scaling=hull_scaling, alpha=hull_alpha
            )
    # plot the nodes
    # zOTU nodes
    zotu_nodes = {
        node for node in g.nodes if g.nodes[node]["name"].startswith("Zotu")}
    nx.draw_networkx_nodes(
        g, layout, nodelist=list(zotu_nodes), node_shape="o",
        node_color=[g.nodes[node]["colour"] for node in zotu_nodes],
        ax=ax
    )
    # metabolite nodes
    met_nodes = [node for node in g.nodes if node not in zotu_nodes]
    nx.draw_networkx_nodes(
        g, layout, nodelist=met_nodes, node_shape="s",
        node_color=[g.nodes[node]["colour"] for node in met_nodes],
        ax=ax
    )
    # community figure
    if plot_individual_communities:
        nx.draw_networkx_nodes(
            g, layout, nodelist=list(zotu_nodes), node_shape="o",
            node_color=[g.nodes[node]["colour"] for node in zotu_nodes],
            ax=comm_ax
        )
        nx.draw_networkx_nodes(
            g, layout, nodelist=met_nodes, node_shape="s",
            node_color=[g.nodes[node]["colour"] for node in met_nodes],
            ax=comm_ax
        )

    # drawing edges
    nx.draw_networkx_edges(
        g, layout, edge_color=[g.edges[edge]["corr"] for edge in g.edges],
        width=[g.edges[edge]["weight"] * 2 for edge in g.edges],
        edge_vmin=-1, edge_vmax=1, edge_cmap=plt.get_cmap("seismic"), ax=ax
    )

    # node labels
    sub_labels = {
        node: label for node, label in label_dict.items()
        if layout.get(node) is not None
    }
    nx.draw_networkx_labels(g, layout, labels=sub_labels, font_size=10, ax=ax)

    if plot_individual_communities:
        nx.draw_networkx_labels(
            g, layout, labels=sub_labels, font_size=10,
            ax=comm_ax
        )
        comm_fig.savefig(community_figure)
        plt.close(comm_fig)

    # legends
    # node shape legend
    shape_handles = [
        plt.Line2D(
            [0], [0],  markersize=25, marker="o", linewidth=0,
            markerfacecolor="tab:grey", label="zOTU", markeredgewidth=2,
            markeredgecolor="tab:grey"
        ),
        plt.Line2D(
            [0], [0],  markersize=25, marker="s", linewidth=0,
            markerfacecolor="tab:grey", label="Metabolite", markeredgewidth=2,
            markeredgecolor="tab:grey"
        )
    ]
    ax.add_artist(
        plt.legend(
            handles=shape_handles, loc="upper right", bbox_to_anchor=(1.06, .9),
            labelspacing=2, frameon=False, title="Modality"
        )
    )
    # node shape legend
    col_handles = [
        plt.Line2D(
            [0], [0],  markersize=25, marker="o", linewidth=0,
            markerfacecolor="tab:blue", label="EEN"
        ),
        plt.Line2D(
            [0], [0],  markersize=25, marker="o", linewidth=0,
            markerfacecolor="tab:orange", markeredgecolor="tab:orange",
            label="PostEEN"
        )
    ]
    ax.add_artist(
        plt.legend(
            handles=col_handles, loc="upper right", bbox_to_anchor=(1.05, .8),
            labelspacing=2, frameon=False, title="Associated to"
        )
    )
    # communities
    comm_colours = [
        (*plt.get_cmap("tab10")(i)[:3], hull_alpha)
        for i in range(n_communities)
    ]
    comm_handles = [
        plt.Line2D(
            [0], [0],  markersize=25, marker="s", linewidth=0,
            markerfacecolor=comm_colours[i], markeredgecolor=comm_colours[i],
            label=f"Community {i + 1}"
        )
        for i in range(n_communities)
    ]
    ax.add_artist(
        plt.legend(
            handles=comm_handles, loc="upper right", bbox_to_anchor=(1.07, .7),
            labelspacing=2, frameon=False, title="Communities"
        )
    )

    # edge colour bar
    ecbar_map = plt.cm.ScalarMappable(
        cmap="seismic", norm=plt.Normalize(vmin=-1, vmax=1))
    ecbar_map._A = []
    ecbar = fig.colorbar(ecbar_map, ax=ax, location="left", shrink=.8)
    ecbar.set_label("Correlation")

    # format plotting axis
    ax.axis("off")
    ax.grid(False)

    return fig, ax, layout


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--data-path", type=str, help="Path to the folder with the input data")
    parser.add_argument(
        "--result-path", type=str, required=False, default=None,
        help="Path to the folder to save the plots to. If not passed "
             "'data-path' will be used"
    )

    args = parser.parse_args()

    input_path = pathlib.Path(args.data_path)
    if args.result_path is None:
        result_path = input_path
    else:
        result_path = pathlib.Path(args.result_path)

    # load pre-computed data
    node_annotation = pd.read_csv(
        input_path / "name_map.csv", index_col=0).iloc[:, 0].to_dict()
    node_annotation = {
        k: v.split(";")[0].replace("REF_", "") if isinstance(v, str) else np.nan
        for k, v in node_annotation.items()
    }

    graph = nx.read_graphml(input_path / "correlation_network.graphml")
    # relabel nodes since exported node labels are enumeration only
    graph = nx.relabel_nodes(graph, nx.get_node_attributes(graph, "name"))

    # adding correlation attribute as weight attribute is abs(correlation)
    corr_mat = pd.read_csv(input_path / "corr_mat.csv", index_col=0)
    for e in graph.edges:
        graph.edges[e]["corr"] = corr_mat.loc[e[0], e[1]]

    # loading communities computed in R
    community_partition = pd.read_csv(
        input_path / "community_partitions.csv", index_col=0).iloc[:, 0]
    community_partition = community_partition.to_dict()

    # remove all nodes if they are not in any of the communities
    nodes_to_remove = [
        v for v in graph.nodes if community_partition.get(v) is None]
    graph.remove_nodes_from(nodes_to_remove)

    figure, axis = plt.subplots(figsize=(32, 18))
    figure, _, pos = plot_community_graph(
        graph, community_partition, node_annotation, inter_community_scale=3,
        intra_community_scale=1, plot=(figure, axis),
        plot_individual_communities=True,
        community_figure=result_path / "correlation_network_no_edges.pdf"
    )
    plt.savefig(result_path / "correlation_network.pdf")
    plt.close(figure)

    # Numbered node labels
    numbered_annotation = {}
    number_map = {}
    for k, v in node_annotation.items():
        if k in pos.keys():
            idx = len(numbered_annotation) + 1
            numbered_annotation[k] = str(idx)
            number_map[k] = {"Annotation": v, "Number": idx}

    pd.DataFrame.from_dict(number_map, orient="index").to_csv(
        result_path / "number_legend.csv")

    comm_fig_file = result_path / "correlation_network_numbered_no_edges.pdf"
    figure, axis = plt.subplots(figsize=(32, 18))
    figure, _, __ = plot_community_graph(
        graph, community_partition, numbered_annotation,
        inter_community_scale=3, intra_community_scale=1, plot=(figure, axis),
        plot_individual_communities=True, community_figure=comm_fig_file,
        layout=pos
    )
    plt.savefig(result_path / "correlation_network_numbered.pdf")
    plt.close(figure)

