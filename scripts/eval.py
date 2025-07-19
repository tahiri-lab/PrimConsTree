"""Run several instance of primconstree and extended majority rule
on different datasets, compute metrics, and save results in a file.
"""

import json
import logging
from itertools import product

import ete3

from utils.consensus import consensus
from utils.distances import distance
from utils.misc import create_unique_file, get_alg_id
from utils.trees import read_trees

###########################
### BEGIN OF PARAMETERS ###
###########################

##############
### INPUTS ###
##############

INPUT = "datasets/eval"  # directory to take the inputs from
TXT_DIR = "HS"  # directory to take the inputs for most algorithms (take txts)
NEX_DIR = "FACT"  # directory to take the inputs for FACT2 algorithms

# note parameter must be compatible with the generated trees
K = [10, 30]  # 50, 90, 110, 130, 150]  # values for number of trees
N = [10, 20]  # , 30, 40, 50]  # values for number of leaves
C = [1, 2.5, 5, 7.5, 10]  # values for coalescence rate
NB_BATCH = 5  # number of batch per combination of parameters

BENCHMARK = 0  # number of iteration on benchmark execution time (0 for no benchmark)

###############
### OUTPUTS ###
###############

RESULTS_FILE = "outputs/test.json"  # file to output the results

# algorithms to perfoem (maj, pct, old_pct, freq)
ALGS = [
    ("pct", {}),
    ("pct", {"old_prim": True}),
    ("maj", {}),
    ("fdct", {"cst": 1}),  # Assign cst on all branches
    (
        "fdct",
        {"avg": True},
    ),  # Assign the average edge length on all tree (coalescence rate)
    ("maj_plus", {"cst": 1}),  # Assign cst on all branches
    (
        "maj_plus",
        {"avg": True},
    ),  # Assign the average edge length on all tree (coalescence rate)
]

# distances metrics to evaluates
DISTS = [
    ("rf", {}),
    ("bsd", {}),
    ("tdist", {}),
    ("qdist", {}),
    ("kcdist", {"lamb": 0.0}),
    ("kcdist", {"lamb": 0.5}),
    ("kcdist", {"lamb": 1.0}),
]

#########################
### END OF PARAMETERS ###
#########################


def eval_consensus(
    alg: tuple[str, dict],
    dists: list[tuple[str, dict]],
    input_file: str,
    input_trees: list[ete3.Tree],
    benchmark: int,
    coal: float,
) -> dict:
    """Compute consensus trees and metrics for several batches of input trees

    Args:
        alg (tuple[str, dict]): the name of the consensus algorithm and additional parameters
        dists (list[tuple[str, dict]]): the list of metrics to compute (name + additional parameters)
        input_file (str): input file for the consensus
        input_trees (list): list of input trees as ete3 objects
        benchmark (int): number of iterations for benchmark (0 for no benchmark)
        coal (float): the coalescence rate to get average branch length

    Returns:
        dict: input and consensus as newick strings, metrics
    """
    logging.info("Processing algorithm %s", get_alg_id(alg))

    cons, tm = consensus(input_file, alg[0], coal, **alg[1])
    if benchmark > 0 and tm is not None:
        duration = tm.timeit(benchmark)
    else:
        duration = 0

    results: dict = {
        get_alg_id(m): distance(input_trees, cons, m[0], **m[1]) for m in dists
    }
    results.update({"cons": cons.write(), "duration": duration})

    return results


logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)

if __name__ == "__main__":
    # Execute evaluation on each parameters combinations
    combinations = []
    for k, n, c, b in product(K, N, C, range(NB_BATCH)):
        logging.info(
            "Processing combination k=%i n=%i c=%i b=%i batches, with %i benchmark iterations",
            k,
            n,
            c,
            b,
            BENCHMARK,
        )

        file_txt = f"{INPUT}/{TXT_DIR}/k{k}_n{n}_c{c}_b{b}.txt"
        file_nex = f"{INPUT}/{NEX_DIR}/k{k}_n{n}_c{c}_b{b}.nexus"
        input_trees = read_trees(file_txt)

        # Save parameters
        comb = {
            "file": file_txt,
            "k": k,
            "n": n,
            "c": c,
            "batch": b,
            "benchmark": BENCHMARK,
            "inputs": [t.write() for t in input_trees],
        }

        # Evaluate consensus trees
        for a in ALGS:
            input_file = file_nex if a[0] in ["fdct", "maj_plus"] else file_txt
            comb[get_alg_id(a)] = eval_consensus(
                a, DISTS, input_file, input_trees, BENCHMARK, c
            )

        combinations.append(comb)

    result_file = create_unique_file(RESULTS_FILE)
    with open(result_file, "w") as f:
        json.dump(combinations, f, indent=4)

    logging.info("Saved results to %s", result_file)
