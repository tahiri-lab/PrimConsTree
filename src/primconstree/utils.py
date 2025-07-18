"""Useful funtion to manipulate trees and encode node ids"""

import ete3


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


def clade_to_id(clade: list[str], taxa: list[str]) -> int:
    """Build the node id as the binary representations of
    the clade for example the clade ABD within the taxa
    A,B,C,D,E get the binary representation 11010 which give
    id 26

    Args:
        clade: the list of taxon in the clade
        taxa: the ordered list of all taxa to map against

    Return:
        int: the binary mapping as an integer
    """
    identifier = 0
    for i, t in enumerate(taxa):
        if t in clade:
            identifier += 1 << i
    return identifier


def get_root_id(taxa: list[str]) -> int:
    """Get the id for the root (clade with all taxa)
    one could just use clade_to_id but this is faster
    """
    return sum(1 << i for i in range(len(taxa)))


def id_to_clade(identifier: int, taxa: list[str]) -> list[str]:
    """Get the clade contend as a list of taxon from
    the id infered with the method of clade_to_id
    """
    clade = []
    for i, t in enumerate(taxa):
        if identifier & (1 << i) == (1 << i):
            clade.append(t)
    return clade
