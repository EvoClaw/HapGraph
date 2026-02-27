"""
Neighbor-Joining tree construction from F2 distance matrix.
Used as the starting topology for greedy admixture edge search.
"""

import numpy as np
from typing import List, Tuple, Optional
import networkx as nx


def nj_tree_from_f2(
    f2_matrix: np.ndarray,
    populations: List[str],
) -> nx.DiGraph:
    """
    Build an unrooted Neighbor-Joining tree from an F2 distance matrix,
    then root it at an arbitrary internal node to produce a directed tree.

    Parameters
    ----------
    f2_matrix : symmetric N×N array of F2 distances
    populations : list of population names (length N)

    Returns
    -------
    tree : nx.DiGraph with nodes as population names and internal nodes
           labeled 'internal_0', 'internal_1', etc.
           Edges have attribute 'weight' = branch length.
    """
    n = len(populations)
    assert f2_matrix.shape == (n, n), "F2 matrix shape mismatch"

    # Standard NJ algorithm
    D = f2_matrix.copy().astype(float)
    np.fill_diagonal(D, 0.0)

    nodes = list(populations)  # current active nodes
    edges = []  # (u, v, length)
    internal_count = 0

    while len(nodes) > 3:
        n_curr = len(nodes)
        r = D.sum(axis=1)  # row sums

        # Q matrix
        Q = np.zeros((n_curr, n_curr))
        for i in range(n_curr):
            for j in range(n_curr):
                if i != j:
                    Q[i, j] = (n_curr - 2) * D[i, j] - r[i] - r[j]

        # find minimum off-diagonal Q
        np.fill_diagonal(Q, np.inf)
        idx = np.unravel_index(Q.argmin(), Q.shape)
        i, j = int(idx[0]), int(idx[1])

        # branch lengths to new node u
        delta_i = 0.5 * D[i, j] + (r[i] - r[j]) / (2 * (n_curr - 2))
        delta_j = D[i, j] - delta_i

        new_node = f"internal_{internal_count}"
        internal_count += 1

        # distances from new node to remaining nodes
        new_row = np.zeros(n_curr)
        for k in range(n_curr):
            if k != i and k != j:
                new_row[k] = 0.5 * (D[i, k] + D[j, k] - D[i, j])

        # record edges
        edges.append((new_node, nodes[i], max(delta_i, 0.0)))
        edges.append((new_node, nodes[j], max(delta_j, 0.0)))

        # reduce distance matrix
        keep = [k for k in range(n_curr) if k != i and k != j]
        D_new = np.zeros((len(keep) + 1, len(keep) + 1))
        for a_idx, a in enumerate(keep):
            for b_idx, b in enumerate(keep):
                D_new[a_idx, b_idx] = D[a, b]
        for a_idx, a in enumerate(keep):
            D_new[a_idx, len(keep)] = new_row[a]
            D_new[len(keep), a_idx] = new_row[a]

        D = D_new
        nodes = [nodes[k] for k in keep] + [new_node]

    # Connect final three nodes
    assert len(nodes) == 3
    i0, i1, i2 = 0, 1, 2
    d01, d02, d12 = D[i0, i1], D[i0, i2], D[i1, i2]
    b0 = 0.5 * (d01 + d02 - d12)
    b1 = 0.5 * (d01 + d12 - d02)
    b2 = 0.5 * (d02 + d12 - d01)
    new_node = f"internal_{internal_count}"
    edges.append((new_node, nodes[i0], max(b0, 0.0)))
    edges.append((new_node, nodes[i1], max(b1, 0.0)))
    edges.append((new_node, nodes[i2], max(b2, 0.0)))

    # Build directed graph rooted at the last internal node
    G = nx.DiGraph()
    for u, v, w in edges:
        G.add_edge(u, v, weight=w)

    return G


def tree_to_newick(tree: nx.DiGraph, root: Optional[str] = None) -> str:
    """
    Convert a directed tree to Newick format.

    Parameters
    ----------
    tree : nx.DiGraph (edges point from parent to child)
    root : root node name; if None, first node with in-degree 0 is used

    Returns
    -------
    newick : str
    """
    if root is None:
        roots = [n for n in tree.nodes if tree.in_degree(n) == 0]
        if not roots:
            raise ValueError("No root found (no node with in-degree 0)")
        root = roots[0]

    def _recurse(node):
        children = list(tree.successors(node))
        if not children:
            return node
        parts = []
        for child in children:
            w = tree[node][child].get("weight", 0.0)
            child_str = _recurse(child)
            parts.append(f"{child_str}:{w:.6f}")
        return "(" + ",".join(parts) + ")"

    return _recurse(root) + ";"


def robinson_foulds_distance(
    tree_a: nx.DiGraph,
    tree_b: nx.DiGraph,
    leaf_nodes: List[str],
) -> int:
    """
    Compute Robinson-Foulds distance between two unrooted trees.

    Converts edge sets to bipartitions and counts symmetric difference.

    Parameters
    ----------
    tree_a, tree_b : nx.DiGraph representing trees
    leaf_nodes : list of leaf node names

    Returns
    -------
    rf : int (0 = identical topology)
    """

    def _bipartitions(tree: nx.DiGraph) -> set:
        undirected = tree.to_undirected()
        parts = set()
        for edge in undirected.edges:
            # Remove edge, find which leaves are on each side
            temp = undirected.copy()
            temp.remove_edge(*edge)
            components = list(nx.connected_components(temp))
            if len(components) == 2:
                side_a = frozenset(n for n in components[0] if n in leaf_nodes)
                side_b = frozenset(n for n in components[1] if n in leaf_nodes)
                if side_a and side_b:
                    parts.add(frozenset([side_a, side_b]))
        return parts

    bp_a = _bipartitions(tree_a)
    bp_b = _bipartitions(tree_b)
    return len(bp_a.symmetric_difference(bp_b))
