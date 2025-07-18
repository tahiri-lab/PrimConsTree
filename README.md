# PrimConsTree

`PrimConsTree` is a project that:

A consensus tree in phylogenetics is a representation that summarizes the relationships among different phylogenetic trees. It helps to identify the most supported or accepted evolutionary relationships. It thus provides a more robust understanding of evolutionary history. We study the challenge of constructing a consensus tree from a set of phylogenetic trees, each with branch-length data encoding the temporal evolution of genetic mutations. Despite the richness of this temporal information, the prevailing practice in the field tends to prioritize topological features while neglecting branch length. In response, we present a new approach that incorporates both topological information and average branch length to build a more complete consensus tree. Our adaptation of the well-known Prim algorithm efficiently identifies the minimum weight of the branch length in the consensus tree. Through empirical evaluation, we show that the proposed algorithm outperforms existing methods. This research highlights the importance of considering branch-length data in consensus tree construction and makes a valuable contribution to the field of phylogenetic analysis.

## Using PrimConsTree

### Installation

PrimConsTree was developpend under `python 3.12`, it can be installed directly with the following command: 
```bash
pip install -e src
```

PrimConsTree has the following strict dependencies (while it might work with other versions of the listed libraries, it has only been tested with the following):
- networkx 3.3
- ete3 3.1.3
- matplotlib 3.9.0 (for debug purpose only)

### How to run

Once installed the algorithm can be run with
```bash
python3 -m primconstree -f <input_file>
```

The input file should be a text file with a newick string on each line.
Trees are expected to have branch length. See datasets/tests/ for compatible files.

The program has optional arguments that are used for developpement and tests purposes.
If you only wish to compute consensus tree using the PCT algorithm we recommend leaving the default argumens.
The following is accepted :
```
-h, --help            show this help message and exit
-f FILE, --file FILE  input file path, input is a text file with a newick tree on each line.
-v [VERSION], --version [VERSION]
                      Choose the version of the algorithm :
                      - 0 : first version of PCT where MST criterion are min edge length, max edge frequency, in this order of priority.
                      - 1 (default) : last version where MST criterion are max edge frequency, max fringe vertex frequency, max mst vertex frequency in this order of priority.
-l LEN_ON_MERGE, --len_on_merge LEN_ON_MERGE
                      Choose how branch length are handled removing unecessary internal node in the last step of the algorithm:
                      - 'sum' (default) : branch length are summed, for instance deleting 'B' in A -(1)-> B -(2.5)-> C result in A -(3.5)-> C
                      - 'avg' : branch length are averaged, for instance deleting 'B' in A -(1)-> B -(2.5)-> C result in A -(1.75)-> C
-d, --debug           If set, the algorithm will output informations between steps. Especially, this include plotting the supergraph, the mst, and the consensus tree
```

## Developpment and Evaluation 

Along with our code are provided all the script used for evaluation of PrimConsTree, these can be found in `scripts/` directory.
This section provide useful informations to run these scripts.

### Requirements

First you will need `python 3.12` with the packages listed in `requirements.txt`.

**HybridSim :** Tree generation is performed with [hybridsim](https://github.com/MichaelWoodhams/HybridSim) tool, the `.jar` archive is included in `scripts/tools/` for convinience.
Java is required to be in the `PATH` to run the generation script.

**FACT :** If you wish to evaluate PCT against the frequency difference consensus tree you will need the [FACT2](https://github.com/Mesh89/FACT2) program.
In another hand, evaluating the majority rule + consensus tree require the [FACT](https://github.com/Mesh89/FACT) program.
Compiled versions are already included in `scripts/tools/`.
If you find one to not work properly, a bash script `scripts/tools/setup_fact.sh` is provided and can be run from the `scripts/tools/` directory to generate the two executables.

**tqDist :** If you wish to evaluate triplet or quartet distance you will need the [tqDist](https://www.birc.au.dk/~cstorm/software/tqdist/) package.
Installation instruction are provided in the linked website.

### Tree generation

The script `scripts/generate.py` is used to generate trees. A dedicated section of the script allow you to modify simulation parameters.

### Consensus evaluation

The script `scripts/eval.py` is used to compute consensus with several algorithms and evaluate distance to the input trees with several metrics.
Every parameters can be modified in a dedicated section of the script.
It takes as input data generated by the generation script mentioned above, and output the results as a json file.

### Visualizing results

The notebook `scripts/results.ipynb` can be used to visualize results outputed from the previous steps.
Detailed instructions are provided in the notebook.

## Contributing to `PrimConsTree`

To contribute to `PrimConsTree`, follow these steps:

1. Fork this repository.
2. Create a branch: `git checkout -b <branch_name>`.
3. Make your changes and commit them: `git commit -m '<commit_message>'`
4. Push to the original branch: `git push origin <project_name>/<location>`
5. Create the pull request.

Alternatively see the GitHub documentation on [creating a pull request](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request).

## Contact

If you want to contact me you can reach me at `<Nadia.Tahiri@USherbrooke.ca>`.

## License

This project uses the following license: [MIT license](https://github.com/tahiri-lab/PrimConsTree?tab=MIT-1-ov-file#readme).
