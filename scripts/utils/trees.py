"""Utilities to manipulate trees"""

import ete3
from Bio import Phylo

LEAVES_MAP = {
    "A": "1",
    "B": "2",
    "C": "3",
    "D": "4",
    "E": "5",
    "F": "6",
    "G": "7",
    "H": "8",
    "I": "9",
    "J": "10",
    "K": "11",
    "L": "12",
    "M": "13",
    "N": "14",
    "O": "15",
    "P": "16",
    "Q": "17",
    "R": "18",
    "S": "19",
    "T": "20",
    "U": "21",
    "V": "22",
    "W": "23",
    "X": "24",
    "Y": "25",
    "Z": "26",
    "A1": "27",
    "B1": "28",
    "C1": "29",
    "D1": "30",
    "E1": "31",
    "F1": "32",
    "G1": "33",
    "H1": "34",
    "I1": "35",
    "J1": "36",
    "K1": "37",
    "L1": "38",
    "M1": "39",
    "N1": "40",
    "O1": "41",
    "P1": "42",
    "Q1": "43",
    "R1": "44",
    "S1": "45",
    "T1": "46",
    "U1": "47",
    "V1": "48",
    "W1": "49",
    "X1": "50",
    "Y1": "51",
    "Z1": "52",
}


def read_trees(input_file: str, nwk_format: int = 0) -> list[ete3.Tree]:
    """Read the trees from the input file in Newick format using ete3

    Args:
        input_file (str): Path to the input file.
        nwk_format (int, optional): The format of the Newick tree (see ete3 references tutorial). Defaults to 0.

    Returns:
        list[ete3.Tree]: list of the Trees object
    """
    trees = []
    with open(input_file, "r") as file:
        # Read the tree from the file
        for _, tree in enumerate(file.readlines()):
            tree.strip()
            trees.append(ete3.Tree(tree, format=nwk_format))
    return trees


def phylo_to_ete3(tree: Phylo.BaseTree.Tree) -> ete3.Tree:
    """Convert a Biopython tree to an ete3 tree

    Args:
        tree (Phylo.BaseTree.Tree): the Biopython tree to convert

    Returns:
        ete3.Tree: the ete3 tree
    """
    newick = tree.format("newick")
    ete_tree = ete3.Tree(newick)
    return ete_tree


def map_to_fact(tree: str) -> str:
    """Convert a newick tree yielded by HybridSim in a tree usable in FACT package by renaming leaves.

    Args:
        tree (str): the tree in newick format.

    Returns:
        str: the tree with renamed leaves
    """
    mapped = tree
    for k, v in LEAVES_MAP.items():
        sub1, rep1 = f"({k}:", f"({v}:"
        sub2, rep2 = f",{k}:", f",{v}:"
        mapped = mapped.replace(sub1, rep1).replace(sub2, rep2)
    return mapped


def map_from_fact(tree: str) -> str:
    """Convert a newick tree yielded by FACT algorithm by replacing original leaf names.

    Args:
        tree (str): the tree in newick format.

    Returns:
        str: the tree with renamed leaves
    """
    mapped = tree
    for k, v in LEAVES_MAP.items():
        sub1, rep1 = f"({v},", f"({k},"
        sub2, rep2 = f",{v})", f",{k})"
        sub3, rep3 = f",{v},", f",{k},"
        mapped = mapped.replace(sub1, rep1).replace(sub2, rep2).replace(sub3, rep3)
    return mapped


def set_cst_length(tree: ete3.Tree, cst):
    tree_2 = tree.copy()
    for node in tree_2.traverse():
        if node.up is not None:
            node.dist = cst
    return tree_2
