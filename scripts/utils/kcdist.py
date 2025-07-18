class _KC_node:
    name: str
    depth_unweighted: int
    depth_weighted: float
    son: list[int]
    label: set
    father: int

def _combinations(iterable, r):
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

def _read_newick(newick_str: str, nodes: list[_KC_node],
                 depth_weighted: float, depth_unweighted: int):
    def _split_newick(newick_str: str) -> list[str]:
        n_lbracket = 0
        rtn = []
        last_i = 0
        for i, char in enumerate(newick_str):
            if char == '(':
                n_lbracket += 1
            elif char == ')':
                n_lbracket -= 1
            elif char == ',' and n_lbracket == 0:
                rtn.append(newick_str[last_i:i])
                last_i = i+1
        rtn.append(newick_str[last_i:])
        return rtn
    
    newick_str = newick_str.rstrip(";\n")
    if newick_str[0] == '(': # not a leaf
        is_leaf = False
        node_info = newick_str[newick_str.rfind(')') + 1:]
        son_str = _split_newick(newick_str[1:newick_str.rfind(')')])
    else:
        is_leaf = True
        node_info = newick_str
    new_node = _KC_node()
    new_node.son = []
    if ':' in node_info:
        new_node.name = node_info.split(':')[0]
        branch_length = float(node_info.split(':')[1])
    else:
        new_node.name = node_info
        branch_length = 0
    new_node.depth_unweighted = depth_unweighted + 1
    new_node.depth_weighted = depth_weighted + branch_length
    if not is_leaf:
        for i in range(len(son_str)):
            son_ind= _read_newick(son_str[i], nodes,
                                  new_node.depth_weighted,
                                  new_node.depth_unweighted)
            new_node.son.append(son_ind)
        new_node.label = set()
        for ind in new_node.son:
            new_node.label = new_node.label | nodes[ind].label
            nodes[ind].father = len(nodes)
    else:
        new_node.label = {new_node.name}
    nodes.append(new_node)
    return len(nodes) - 1

def _get_vector(focal_ind: int, nodes: list[_KC_node], lam: float) -> dict:
    focal_node = nodes[focal_ind]
    if focal_node.son == []:
        dist = focal_node.depth_weighted - nodes[focal_node.father].depth_weighted
        return {focal_node.name: (1 + lam * (dist - 1))}
    rtn_vector = {}
    labels = []
    for ind in focal_node.son:
        rlt = _get_vector(ind, nodes, lam)
        # print(rlt, nodes[ind].name, ind)
        rtn_vector = rtn_vector | rlt
        labels.append(nodes[ind].label)
    dist = focal_node.depth_unweighted * (1 - lam) + focal_node.depth_weighted * lam
    for label_i, label_j in _combinations(labels, 2):
        for item_i in label_i:
            for item_j in label_j:
                rtn_vector[frozenset({item_i, item_j})] = dist
    return rtn_vector

def KC_dist(tree1: str, tree2: str, lam: float = 0) -> float:
    if lam > 1 or lam < 0:
        raise ValueError("Invalid lambda! Lambda value should be between 0 and 1.")
    tree_node1 = []
    _read_newick(tree1, tree_node1, 0, 0)
    tree_node2 = []
    _read_newick(tree2, tree_node2, 0, 0)
    if tree_node1[-1].label != tree_node2[-1].label:
        raise ValueError("Invalid tree nodes! Tips of two trees should be the same.")
    vec1 = _get_vector(len(tree_node1)-1, tree_node1, lam)
    vec2 = _get_vector(len(tree_node2)-1, tree_node2, lam)
    rtn = 0
    for item in vec1:
        rtn += (vec1[item] - vec2[item]) ** 2
    return rtn ** 0.5