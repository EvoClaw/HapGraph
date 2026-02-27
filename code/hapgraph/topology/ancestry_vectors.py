"""
Ancestry vector (branch contribution) approach for admixture graph F-statistics.

Computes expected F2, F3, F4 from the graph topology using ancestry vectors q_A[b]:
the fraction of population A's lineage passing through branch b.

Formulas:
  E[F2(A,B)]     = sum_b (q_A[b] - q_B[b])^2 * branch_length[b]
  E[F3(C;A,B)]   = sum_b (q_C[b] - q_A[b]) * (q_C[b] - q_B[b]) * branch_length[b]
  E[F4(A,B;C,D)] = sum_b (q_A[b] - q_B[b]) * (q_C[b] - q_D[b]) * branch_length[b]
"""

import numpy as np
import networkx as nx
from typing import List, Tuple, Optional, Dict

# Import AdmixGraph from greedy_search to avoid circular imports at package level
# Callers should pass the graph object


def _get_root(tree: nx.DiGraph) -> str:
    """Return the root node (in-degree 0)."""
    roots = [n for n in tree.nodes if tree.in_degree(n) == 0]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root, got {len(roots)}: {roots}")
    return roots[0]


def _tree_parent(tree: nx.DiGraph, node: str) -> Optional[str]:
    """Return the tree parent of node (unique predecessor in directed tree)."""
    preds = list(tree.predecessors(node))
    if len(preds) == 0:
        return None
    if len(preds) > 1:
        raise ValueError(f"Node {node} has multiple tree parents: {preds}")
    return preds[0]


def _path_from_root(tree: nx.DiGraph, node: str) -> List[Tuple[str, str]]:
    """Return list of (parent, child) edges on the path from root to node."""
    path = []
    current = node
    while True:
        parent = _tree_parent(tree, current)
        if parent is None:
            break
        path.append((parent, current))
        current = parent
    path.reverse()  # root -> ... -> node order
    return path


def _topological_order_admixture(
    admixture_edges: List[Tuple[str, str, float]],
) -> List[Tuple[str, str, float]]:
    """
    Order admixture edges so that when (src, tgt, alpha) is processed,
    all edges into src have already been processed (src's q is final).
    """
    if not admixture_edges:
        return []

    targets = {tgt for _, tgt, _ in admixture_edges}
    sources = {src for src, _, _ in admixture_edges}
    # Sources that are never targets are "done" from the start
    done_sources = sources - targets

    # Build dependency: tgt depends on src
    deps = {tgt: set() for _, tgt, _ in admixture_edges}
    for src, tgt, _ in admixture_edges:
        if src != tgt:
            deps.setdefault(tgt, set()).add(src)

    # Kahn's algorithm
    result = []
    remaining = list(admixture_edges)

    while remaining:
        progress = False
        for i, (src, tgt, alpha) in enumerate(remaining):
            if deps.get(tgt, set()).issubset(done_sources):
                result.append((src, tgt, alpha))
                done_sources.add(tgt)
                remaining.pop(i)
                progress = True
                break
        if not progress:
            result.extend(remaining)
            break

    return result


def compute_ancestry_vectors(
    graph,
    branch_lengths_override: Optional[Dict[Tuple[str, str], float]] = None,
) -> Tuple[np.ndarray, List[str], List[Tuple[str, str]]]:
    """
    Compute the ancestry matrix Q for all leaf populations.

    Parameters
    ----------
    graph : AdmixGraph
        Has .tree (nx.DiGraph), .admixture_edges, .populations
    branch_lengths_override : optional dict (parent, child) -> length
        Override branch lengths from tree edge weights

    Returns
    -------
    Q : np.ndarray, shape (n_pops, n_branches)
        Q[i, b] = ancestry fraction of population i through branch b
    pop_names : List[str]
        Population names in same order as Q rows
    branch_list : List[Tuple[str, str]]
        branch_list[b] = (parent_node, child_node) for each branch
    """
    tree = graph.tree
    admixture_edges = graph.admixture_edges
    populations = graph.populations

    # 1. Build branch list: all (parent, child) tree edges in consistent order
    branch_list = list(tree.edges())
    n_branches = len(branch_list)
    branch_to_idx = {e: i for i, e in enumerate(branch_list)}

    # 2. Compute tree-path ancestry for every node
    root = _get_root(tree)
    q_tree: Dict[str, np.ndarray] = {}

    for node in tree.nodes:
        q = np.zeros(n_branches)
        for parent, child in _path_from_root(tree, node):
            if (parent, child) in branch_to_idx:
                q[branch_to_idx[(parent, child)]] = 1.0
        q_tree[node] = q

    # 3. Initialize q for all nodes (start with tree path)
    q_node = {n: q_tree[n].copy() for n in tree.nodes}

    # 4. Apply admixture in dependency order
    for src, tgt, alpha in _topological_order_admixture(admixture_edges):
        if tgt not in q_node:
            continue
        q_node[tgt] = (1.0 - alpha) * q_node[tgt] + alpha * q_node[src]

    # 5. Extract Q for leaf populations only
    n_pops = len(populations)
    Q = np.zeros((n_pops, n_branches))
    for i, pop in enumerate(populations):
        if pop in q_node:
            Q[i, :] = q_node[pop]
        else:
            Q[i, :] = q_tree.get(pop, np.zeros(n_branches))

    return Q, populations, branch_list


def expected_f2_matrix(Q: np.ndarray, branch_lengths: np.ndarray) -> np.ndarray:
    """
    Compute expected F2 matrix from ancestry vectors.

    E[F2(A,B)] = sum_b (q_A[b] - q_B[b])^2 * L[b]

    Parameters
    ----------
    Q : np.ndarray, shape (n_pops, n_branches)
    branch_lengths : np.ndarray, shape (n_branches,)

    Returns
    -------
    f2_mat : np.ndarray, shape (n_pops, n_pops), symmetric
    """
    n = Q.shape[0]
    n_b = Q.shape[1]
    L = np.asarray(branch_lengths, dtype=float).flatten()
    if len(L) != n_b:
        L = np.resize(L, n_b)

    f2_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            dq = Q[i, :] - Q[j, :]
            f2_mat[i, j] = f2_mat[j, i] = np.sum(dq * dq * L)
    return f2_mat


def expected_f3_array(
    Q: np.ndarray,
    branch_lengths: np.ndarray,
    triples: List[Tuple[int, int, int]],
) -> np.ndarray:
    """
    Compute expected F3 for each triple (c_idx, a_idx, b_idx).

    E[F3(C; A, B)] = sum_b (q_C[b] - q_A[b]) * (q_C[b] - q_B[b]) * L[b]

    Parameters
    ----------
    Q : np.ndarray, shape (n_pops, n_branches)
    branch_lengths : np.ndarray, shape (n_branches,)
    triples : list of (c_idx, a_idx, b_idx)

    Returns
    -------
    f3_arr : np.ndarray, shape (len(triples),)
    """
    L = np.asarray(branch_lengths, dtype=float).flatten()
    n_b = Q.shape[1]
    if len(L) != n_b:
        L = np.resize(L, n_b)

    f3_arr = np.zeros(len(triples))
    for k, (c, a, b) in enumerate(triples):
        dca = Q[c, :] - Q[a, :]
        dcb = Q[c, :] - Q[b, :]
        f3_arr[k] = np.sum(dca * dcb * L)
    return f3_arr


def expected_f4_array(
    Q: np.ndarray,
    branch_lengths: np.ndarray,
    quartets: List[Tuple[int, int, int, int]],
) -> np.ndarray:
    """
    Compute expected F4 for each quartet (a_idx, b_idx, c_idx, d_idx).

    E[F4(A,B;C,D)] = sum_b (q_A[b] - q_B[b]) * (q_C[b] - q_D[b]) * L[b]

    Parameters
    ----------
    Q : np.ndarray, shape (n_pops, n_branches)
    branch_lengths : np.ndarray, shape (n_branches,)
    quartets : list of (a_idx, b_idx, c_idx, d_idx)

    Returns
    -------
    f4_arr : np.ndarray, shape (len(quartets),)
    """
    L = np.asarray(branch_lengths, dtype=float).flatten()
    n_b = Q.shape[1]
    if len(L) != n_b:
        L = np.resize(L, n_b)

    f4_arr = np.zeros(len(quartets))
    for k, (a, b, c, d) in enumerate(quartets):
        dab = Q[a, :] - Q[b, :]
        dcd = Q[c, :] - Q[d, :]
        f4_arr[k] = np.sum(dab * dcd * L)
    return f4_arr
