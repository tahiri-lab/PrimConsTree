from typing import Generator, cast

import ete3
import matplotlib.pyplot as plt
import networkx as nx

from .utils import get_root_id, id_to_clade


def draw_graph(
    graph: nx.Graph,
    taxa: list[str],
) -> None:
    """Draw a networkx graph with matplotlib
    Edge are labeled with "l:<lenght> - f:<frequency>"
    Nodes are labeled with "f:<frequency> - clade content"
    Node corresponding to leaves are red

    Args:
        graph: the graph to display
        taxa: the list of taxa the vertices id are mapped against
    """
    pos = nx.spring_layout(graph)

    # Draw nodes and edges
    leaf_nodes = [1 << i for i in range(len(taxa))]
    node_colors = [
        "red" if node in leaf_nodes else "lightblue" for node in graph.nodes()
    ]

    nx.draw_networkx_nodes(graph, pos, node_color=node_colors, node_size=500)
    nx.draw_networkx_edges(graph, pos, edge_color="gray")

    # Draw edge labels
    edge_labels = {
        (u, v): f"l:{data["avglen"]} - f:{data["edge_freq"]}"
        for u, v, data in graph.edges(data=True)
    }
    nx.draw_networkx_edge_labels(graph, pos, edge_labels=edge_labels)

    # Draw node labels
    node_labels = {
        u: f"f:{data["node_freq"]} - {','.join(id_to_clade(u, taxa))}"
        for u, data in graph.nodes(data=True)
    }
    nx.draw_networkx_labels(graph, pos, labels=node_labels, font_weight="bold")

    # Display the graph
    plt.title("Graph")
    plt.axis("off")
    plt.show()


def draw_tree(tree: ete3.Tree, taxa: list[str]) -> None:
    """Draw a tree with matplotlib for better vizualization
    Tree is mapped onto a nx.Graph and drawn with draw_graph()
    node and edge frequency are all set to 1.

    Args:
        tree: the tree to display
        taxa: the list of taxa node ids are mapped against
    """

    def tree_to_graph(tree: ete3.Tree) -> nx.Graph:
        """Map the tree onto a graph instance"""
        root_id = get_root_id(taxa)
        G = nx.Graph()
        for node in cast(Generator[ete3.Tree, None, None], tree.traverse("preorder")):
            node_id = int(node.name) if node.up else root_id
            G.add_node(node_id, node_freq=1)
            if node.up is not None:
                parent_id = int(node.up.name) if node.up.up else root_id
                edge_len = node.dist if node.dist is not None else 0.0
                G.add_edge(parent_id, node_id, avglen=edge_len, edge_freq=1)

        return G

    tree_graph = tree_to_graph(tree)
    draw_graph(tree_graph, taxa)
