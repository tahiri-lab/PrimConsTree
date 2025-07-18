"""Function for each step of the PCT algorithm"""

import heapq
from statistics import fmean
from typing import Generator, cast

import ete3
import networkx as nx

from .utils import clade_to_id, get_root_id


def add_tree_to_graph(t: ete3.Tree, taxa: list[str], graph: nx.Graph) -> None:
    """Incorporate a new tree into a graph by idendifying nodes
    with clade_to_id(), incorporating new node and edges, and
    updating their frequency
    Also for edges that already exists, lenghts is summed
    """
    for node in cast(Generator[ete3.Tree, None, None], t.traverse(strategy="preorder")):
        # get mapped identifier based on clade content
        nid = clade_to_id(node.get_leaf_names(), taxa)
        # add node to the graph
        if not nid in graph.nodes:
            graph.add_node(nid, node_freq=0, is_leaf=node.is_leaf())

        # if the node is not the root
        if node.up:
            # count the node frequency
            graph.nodes[nid]["node_freq"] += 1
            pid = clade_to_id(node.up.get_leaf_names(), taxa)

            # add incident edge to the graph
            if nid not in graph[pid]:
                graph.add_edge(pid, nid, avglen=0, edge_freq=0)

            # count edge length and edge frequency
            graph[pid][nid]["avglen"] += node.dist
            graph[pid][nid]["edge_freq"] += 1


def build_graph(trees: list[ete3.Tree], taxa: list[str]) -> nx.Graph:
    """Build the supergraph by incorporating every trees
    Then compute average edge length
    Return the graph.
    """
    graph = nx.Graph()

    # incorporate trees incrementally
    for t in trees:
        add_tree_to_graph(t, taxa, graph)

    # compute average edge length
    for u, v in graph.edges():
        graph[u][v]["avglen"] /= graph[u][v]["edge_freq"]

    return graph


def get_weights(graph: nx.Graph, u: int, v: int, crits: list[str]) -> tuple:
    """Assign a weight to an edge of the graph,
    This weights are used to evaluate which edge take next in the MST,
    the less the better.
    First weight has the most priority, in case of equality, next weight
    will be used.

    Args:
        graph: the graph where the edge is
        u: the vertex of the edge which is already in the MST
        v: the fringe vertex (which is not already in the MST)
        crits: the criterions to use, in order of priority, valid ones are
            "max_nfreq_out", "max_nfreq_in", "max_edge_freq", "min_avg_len"

    Return:
        A tuple, the same size of <crits> with the value of each criterion
    """
    criterion = {}
    criterion["max_nfreq_out"] = (
        1 / graph.nodes[v]["node_freq"]
        if graph.nodes[v]["node_freq"] != 0
        else float("inf")
    )
    criterion["max_nfreq_in"] = (
        1 / graph.nodes[u]["node_freq"]
        if graph.nodes[u]["node_freq"] != 0
        else float("inf")
    )
    criterion["max_edge_freq"] = 1 / graph[u][v]["edge_freq"]
    criterion["min_avg_len"] = graph[u][v]["avglen"]

    weights = tuple(criterion[c] for c in crits)

    return weights


def attach_leaves(graph: nx.Graph, taxa: list[str], crits: list[str]) -> None:
    """Attach the node corresponding to leaves to the mst,
    Use the criterion <crits> do decide which parent is the better
    """
    # Attach the leaf nodes by choosing the most profitable edge
    for u in [1 << i for i in range(len(taxa))]:
        for v in graph[u]:
            weights = get_weights(graph, v, u, crits)
            if graph.nodes[u]["key"] > weights:
                graph.nodes[u]["key"] = weights
                graph.nodes[u]["parent"] = v


def mst_to_tree(graph: nx.Graph) -> ete3.Tree:
    """Given the graph with the MST mapped onto it (by assiging
    a "parent" variable to each node), build the tree.
    The tree is returned as a ete3.Tree instance, with internal nodes
    and average edge lengths
    """
    # create all nodes
    tree_nodes = {
        node_id: ete3.TreeNode(name=node_id) for node_id in graph.nodes.keys()
    }
    root = None

    # attach nodes to their parents
    for nid, node in graph.nodes(data=True):
        tree_node = tree_nodes[nid]
        pid = node["parent"]
        if pid == -1:
            root = tree_node
        else:
            parent_node = tree_nodes[pid]
            parent_node.add_child(tree_node, dist=graph[pid][nid]["avglen"])

    if root is None:
        raise Exception("No root found")

    return ete3.Tree(newick=root.write(format=3), format=3)


def build_mst(
    graph: nx.Graph, src: int, taxa: list[str], crits: list[str]
) -> ete3.Tree:
    """Build the MST of the graph, starting at the source node src
    (where src is the node identifier based on clade_to_id()).
    <crits> define the criterion to choose the edges in order of priority,
    see get_weights() doc for the valid criterions.
    The MST is returned as a ete3.Tree instance, with internal nodes
    and average edge lengths
    """
    nodes = graph.nodes(data=True)

    for _, node in nodes:
        # wether the node is already in the mst
        node["in_mst"] = False
        # the parent of the node in the mst (-1) for undetermined
        node["parent"] = -1
        # store the criterion values for the most profitable edge leading to node, found so far
        node["key"] = (float("inf"),) * len(crits)

    # Priority queue storing the next edges to add
    queue = []

    # Append the source node to the queue to start fron it
    graph.nodes[src]["key"] = (0,) * len(crits)
    heapq.heappush(queue, (*graph.nodes[src]["key"], src))

    # Loop until the priority queue becomes empty
    while queue:
        # Pop the next node to add out of the queue
        u = heapq.heappop(queue)[-1]
        node = graph.nodes[u]
        if node["in_mst"]:
            continue

        # Add the node to the mst
        node["in_mst"] = True

        # Check every child of u
        for v in graph[u]:
            next_node = graph.nodes[v]
            # Skip v if it is already included or if it is a leaf
            if next_node["in_mst"] or next_node["is_leaf"]:
                continue

            # Check if this new edge leading to v is better than the one we found so far
            # If so, update the key of v and push the edge on the heap
            weights = get_weights(graph, u, v, crits)
            if next_node["key"] > weights:
                heapq.heappush(queue, (*weights, v))
                next_node["key"] = weights
                next_node["parent"] = u

    # Attach the leaf nodes
    attach_leaves(graph, taxa, crits)
    # Build the tree from the attached parent values
    mst = mst_to_tree(graph)

    return mst


def remove_unecessary_nodes(
    tree: ete3.Tree, taxa: list[str], average_on_merge: bool = False
) -> None:
    """Remove unnecessary / redundant internal nodes from a tree.
        Modify the tree in place.

    Args:
        tree: the tree to clean
        taxa: the list of taxa the node id were mapped against
        average_on_merge: If True, average edge length when
            removing redundant node, else sum. Defaults to False.
    """
    accumulated_lengths = dict()
    taxa_ids = [1 << i for i in range(len(taxa))]

    # postorder traversal to check children before parent
    for node in cast(Generator[ete3.Tree, None, None], tree.traverse("postorder")):
        # if node is not a taxa (not supposed to be a leaf)
        nid = int(node.name) if node.up else get_root_id(taxa)
        if nid not in taxa_ids:
            children = node.get_children()

            # if the non-leaf node has 0 child just remove it
            if len(children) == 0:
                node.delete(prevent_nondicotomic=False, preserve_branch_length=True)

            # if the non-leaf node has only one children contract the edge from node to its parent
            elif len(children) == 1:
                c = children[0]
                if c not in accumulated_lengths:
                    accumulated_lengths[c] = [c.dist]
                accumulated_lengths[c].append(node.dist)
                node.delete(prevent_nondicotomic=False, preserve_branch_length=True)

    if average_on_merge:
        for node, lengths in accumulated_lengths.items():
            node.dist = fmean(lengths)
