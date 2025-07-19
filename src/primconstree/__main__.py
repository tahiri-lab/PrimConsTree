import argparse

from .primconstree import primconstree
from .utils import read_trees


def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "-f",
        "--file",
        type=str,
        help="input file path, input is a text file with a newick tree on each line.",
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        help=("Choose the seed used to break ties in the MST algorithm"),
        default=0,
    )
    parser.add_argument(
        "-v",
        "--version",
        type=int,
        help=(
            "Choose the version of the algorithm : \n"
            "- 0 : first version of PCT where MST criterion are min edge length, max edge frequency, in this order of priority.\n"
            "- 1 (default) : last version where MST criterion are max edge frequency, max fringe vertex frequency, max mst vertex frequency in this order of priority.\n"
        ),
        default=1,
    )
    parser.add_argument(
        "-l",
        "--len_on_merge",
        type=str,
        help=(
            "Choose how branch length are handled removing unecessary internal node in the last step of the algorithm:\n"
            "- 'sum' (default) : branch length are summed, for instance deleting 'B' in A -(1)-> B -(2.5)-> C result in A -(3.5)-> C \n"
            "- 'avg' : branch length are averaged, for instance deleting 'B' in A -(1)-> B -(2.5)-> C result in A -(1.75)-> C \n"
        ),
        nargs=1,
        default="sum",
    )
    parser.add_argument(
        "-d",
        "--debug",
        action="store_true",
        help="If set, the algorithm will output informations between steps. Especially, this include plotting the supergraph, the mst, and the consensus tree",
    )

    args = parser.parse_args()
    filename = args.file

    if args.version == 0:
        crits = ["min_avg_len", "max_edge_freq"]
    elif args.version == 1:
        crits = ["max_edge_freq", "max_nfreq_out", "max_nfreq_in"]
    else:
        raise Exception(f"PCT version {args.version} invalid")

    if args.len_on_merge == "sum":
        avg_on_merge = False
    elif args.len_on_merge == "avg":
        avg_on_merge = True
    else:
        raise Exception(f"PCT len_on_merge {args.len_on_merge} invalid")

    debug = bool(args.debug)

    input_trees = read_trees(filename)
    consensus = primconstree(
        input_trees, crits=crits, avg_on_merge=avg_on_merge, debug=debug, seed=args.seed
    )
    print(consensus.write())


if __name__ == "__main__":
    main()
