"""Module in charge of generating the consensus tree using the PrimConsTree algorithm"""

import ete3

from .algorithm import build_graph, build_mst, remove_unecessary_nodes
from .debug import draw_graph, draw_tree
from .utils import get_root_id, id_to_clade


def primconstree(
    trees: list[ete3.Tree],
    crits: list[str],
    avg_on_merge: bool = False,
    debug: bool = False,
) -> ete3.Tree:
    """Generate the consensus tree from a set of phylogenetic trees
       using the PrimConsTree algorithm

    Args:
        inputs: list of input trees
        crits: list of criterion to use for the MST in priority order. Valid crit are "max_edge_freq", "max_nfreq_out" (fringe vertex), "max_nfreq_in" (mst vertex), "min_avg_len".
        avg_on_merge: By default, branch length are summed in remove_unecessary_nodes, if True average is computed instead (see --help for more info). Defaults to False.
        debug: If True, display informations at different steps, including graph, mst and tree plots.

    Returns:
        ete3.Tree: the consensus tree
    """
    if not trees:
        return ete3.Tree()

    taxa = sorted(trees[0].get_leaf_names())
    if debug:
        print("PCT: " + "Generating PrimConsTree")
        print("PCT: " + f"Building consensus on taxa: {str(taxa)}")

    graph = build_graph(trees, taxa)
    if debug:
        print("PCT: " + f"SuperGraph generated")
        draw_graph(graph, taxa)

    root = get_root_id(taxa)
    if debug:
        print("PCT: " + f"Searching MST from root f{root}")
    mst = build_mst(graph, root, taxa, crits)
    if debug:
        print("PCT: " + f"MST found")
        draw_tree(mst, taxa)

    remove_unecessary_nodes(mst, taxa, avg_on_merge)
    if debug:
        print("PCT: " + "Unecessary internal nodes removed")
        draw_tree(mst, taxa)

    for l in mst.get_leaves():
        l.name = id_to_clade(int(l.name), taxa)[0]

    return mst
