import subprocess
import timeit

import ete3
from Bio import Phylo
from Bio.Phylo.Consensus import majority_consensus
from primconstree import primconstree

from .trees import map_from_fact, phylo_to_ete3, read_trees, set_cst_length

# Path to fdct (from fact repo) compiled binaries
PATH_TO_FACT1 = "scripts/tools/fact"  # FACT compiled binary
PATH_TO_FACT2 = "scripts/tools/fact2"  # FACT2 compiled binary


def fdct(input_file, cst):
    cmd = [PATH_TO_FACT2, "freq", input_file]
    result = subprocess.run(cmd, capture_output=True, text=True)
    cons = ete3.Tree(map_from_fact(result.stdout.replace("\n", ";")))
    return set_cst_length(cons, cst)


def maj_plus(input_file, cst):
    cmd = ["./scripts/tools/fact1.sh", PATH_TO_FACT1, input_file, "100000000"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    cons = ete3.Tree(map_from_fact(result.stdout.replace("\n", ";")))
    return set_cst_length(cons, cst)


def consensus(
    filename: str, alg: str, coal: float, **args
) -> tuple[ete3.Tree, timeit.Timer]:
    """Compute the consensus tree from a list of input trees using the specified algorithm

    Args:
        filename (str): appropriate input file for the consensus method
        alg (str): algorithm to use (pct, old_pct, maj)
        coal (float): coalescence rate, used to attache average edge length
        args: additional parameters relevant to the algorithm

    Returns:
        tuple: consensus, timit timer for benchmark
    """
    if alg == "pct":
        input_trees = read_trees(filename)
        cons = primconstree(input_trees, debug=False, **args)
        tm = timeit.Timer(lambda: primconstree(input_trees, debug=False, **args))
        return cons, tm

    if alg == "maj":
        input_trees = list(Phylo.parse(filename, "newick"))
        bio_cons = majority_consensus(input_trees, 0)
        cons = phylo_to_ete3(bio_cons)
        tm = timeit.Timer(lambda: majority_consensus(input_trees, 0))
        return cons, tm

    if alg == "fdct":
        tm = timeit.Timer(lambda: fdct(filename, 1))
        if args.get("avg") == True:
            return fdct(filename, 1 / coal), tm
        else:
            return fdct(filename, args.get("cst", 1)), tm

    if alg == "maj_plus":
        tm = timeit.Timer(lambda: maj_plus(filename, 1))
        if args.get("avg") == True:
            return maj_plus(filename, 1 / coal), tm
        else:
            return maj_plus(filename, args.get("cst", 1)), tm

    raise ValueError(f"Unknown algorithm {alg}")
