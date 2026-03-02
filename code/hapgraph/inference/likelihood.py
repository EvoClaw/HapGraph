"""
HapGraph likelihood for Bayesian parameter inference (PyMC 5 + NUTS).

Given a fixed graph topology (with fitted branch lengths from greedy search),
this module:
  1. Decomposes expected F2(alpha) as a quadratic polynomial in alpha — exact
     and differentiable in PyTensor.
  2. Adds IBD timing likelihood (truncation-corrected, L_min = 2 cM):
       E[L | L > 2] = 50/T + 2 cM   (both lineages contribute T meioses)
       mean_IBD_len(tgt, pop) ~ N(50/T + 2, sigma^2)

The quadratic decomposition for K=1 admixture edge (source → target, alpha):
  Q(alpha) = Q0 + alpha * dQ
  where Q0 = ancestry vectors without admixture,
        dQ = (q_source - q_tree_target) for the admixed population row only.

  E[F2(A,B)](alpha) = F2_0[i,j] + alpha * C1[i,j] + alpha^2 * C2[i,j]
  where:
    F2_0[i,j] = sum_b (Q0[i,b] - Q0[j,b])^2 * L[b]   (constant)
    C1[i,j]   = 2 * sum_b (Q0[i,b]-Q0[j,b]) * (dQ[i,b]-dQ[j,b]) * L[b]  (linear)
    C2[i,j]   = sum_b (dQ[i,b]-dQ[j,b])^2 * L[b]      (quadratic)

For K>1 admixture edges: cross-terms arise; handled by precomputing polynomial
coefficients per (alpha_k) independently (first-order interaction approximation).
"""

import numpy as np
import networkx as nx
from scipy.optimize import nnls
from typing import List, Tuple, Optional, Set

from ..topology.ancestry_vectors import compute_ancestry_vectors, expected_f2_matrix


def _get_subtree_leaves(tree: nx.DiGraph, root_node: str, leaf_set: List[str]) -> Set[str]:
    """
    Return the set of leaf populations reachable from root_node in the directed tree.
    If root_node is itself a leaf, returns {root_node}.
    """
    if root_node not in tree:
        return {root_node} if root_node in leaf_set else set()
    reachable = set(nx.descendants(tree, root_node)) | {root_node}
    return reachable & set(leaf_set)


def _get_source_clade_leaves(
    tree: nx.DiGraph, src_node: str, leaf_set: List[str]
) -> Set[str]:
    """
    Return leaves in the clade of src_node's tree parent, including src's sisters.

    For a leaf src_node (e.g. pop1), climbing one level gives the parent clade
    that also includes sister populations (e.g. pop0). This ensures that IBD
    with closely related sister populations is used for timing estimation,
    even when the search happened to select a sister node as the source.
    """
    parents = list(tree.predecessors(src_node))
    if not parents:
        # src is root; use all leaves
        return set(leaf_set)
    parent = parents[0]
    return _get_subtree_leaves(tree, parent, leaf_set)


def _refit_branch_lengths_alpha0(graph, fstats) -> np.ndarray:
    """
    Refit branch lengths via NNLS for the K>0 topology evaluated at alpha=0.

    Why this is needed:
    The NJ K=0 tree branch lengths were fitted BEFORE adding admixture edges.
    They absorb admixture signals by placing admixed populations at
    "compromise" positions (e.g., between the two source clades), so F2_0
    (expected F2 at alpha=0) is already close to the observed data. This causes
    the MCMC to estimate alpha≈0 even when true alpha is large, because adding
    alpha>0 would reduce F2 predictions BELOW the observed values.

    The fix: given the K=1 topology (tree + admixture edge structure), refit
    branch lengths under the assumption alpha=0 (pure tree, no admixture).
    These new branch lengths give F2_0 that reflects the TREE-ONLY prediction,
    so increasing alpha correctly improves the likelihood.

    Parameters
    ----------
    graph : AdmixGraph with K>=1 admixture edges
    fstats : FStatistics object with observed F2 values

    Returns
    -------
    L_alpha0 : np.ndarray, branch lengths fitted at alpha=0 for this topology
    """
    base_graph = graph.copy()
    base_graph.admixture_edges = [(s, t, 0.0) for s, t, _ in graph.admixture_edges]

    Q0, pops, _ = compute_ancestry_vectors(base_graph)
    n = len(pops)

    valid_pairs, obs_vals = [], []
    for i, pa in enumerate(pops):
        for j, pb in enumerate(pops):
            if j <= i:
                continue
            val, _ = fstats.f2(pa, pb)
            if not np.isnan(val):
                valid_pairs.append((i, j))
                obs_vals.append(val)

    if not obs_vals:
        return np.ones(Q0.shape[1]) * 0.01

    n_b = Q0.shape[1]
    n_pairs = len(valid_pairs)
    A = np.zeros((n_pairs, n_b))
    for k, (i, j) in enumerate(valid_pairs):
        dq = Q0[i, :] - Q0[j, :]
        A[k, :] = dq * dq

    L_alpha0, _ = nnls(A, np.array(obs_vals))
    return np.maximum(L_alpha0, 1e-6)


def _compute_f2_polynomial_coeffs(
    graph,
    branch_lengths: np.ndarray,
) -> Tuple[np.ndarray, List[np.ndarray], List[np.ndarray]]:
    """
    Precompute F2 polynomial coefficients for each admixture edge k.

    Returns
    -------
    F2_0 : np.ndarray, shape (n_pairs,)
        F2 under tree backbone only (alpha=0).
    C1_list : list of length K, each shape (n_pairs,)
        Linear coefficient of alpha_k in F2.
    C2_list : list of length K, each shape (n_pairs,)
        Quadratic coefficient of alpha_k in F2.
    pairs : list of (i, j) index pairs
    """
    pops = graph.populations
    n = len(pops)
    K = graph.K

    # Baseline: ancestry vectors with all admixture edges at alpha=0
    base_graph = graph.copy()
    base_graph.admixture_edges = [(s, t, 0.0) for s, t, _ in graph.admixture_edges]
    Q0, _, branch_list = compute_ancestry_vectors(base_graph)

    L = np.asarray(branch_lengths, dtype=float)

    # Build pairs
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]

    # F2_0: baseline with alpha=0
    F2_0_mat = expected_f2_matrix(Q0, L)
    F2_0 = np.array([F2_0_mat[i, j] for i, j in pairs])

    C1_list = []
    C2_list = []

    for k, (src, tgt, _) in enumerate(graph.admixture_edges):
        # dQ[k]: how Q changes for population tgt when alpha_k increases by 1
        # dQ_tgt = q_src - q_base_tgt
        if tgt not in pops or src not in list(graph.tree.nodes):
            C1_list.append(np.zeros(len(pairs)))
            C2_list.append(np.zeros(len(pairs)))
            continue

        tgt_idx = pops.index(tgt)

        # Compute ancestry vectors for pure source path
        src_graph = graph.copy()
        src_graph.admixture_edges = [(s, t, 0.0) for s, t, _ in graph.admixture_edges]
        Q_src, _, _ = compute_ancestry_vectors(src_graph)

        # dQ for this population
        # q_tgt(alpha_k) = (1-alpha_k)*Q0[tgt] + alpha_k*Q_src[src_or_tgt_via_src]
        # We need q_src_node: the ancestry vector of the source node
        # Approximate: if src is a leaf, Q_src[src_idx]; if internal, compute it
        if src in pops:
            src_idx = pops.index(src)
            q_src_node = Q_src[src_idx]
        else:
            # Internal node: compute its ancestry vector from base graph
            # by temporarily making it a leaf
            q_src_node = _get_internal_ancestry_vector(base_graph, src, branch_list)

        dQ_tgt = q_src_node - Q0[tgt_idx]  # shape (n_branches,)

        # Build dQ matrix: only tgt row changes
        dQ = np.zeros_like(Q0)
        dQ[tgt_idx, :] = dQ_tgt

        # C1[i,j] = 2 * sum_b (Q0[i,b]-Q0[j,b]) * (dQ[i,b]-dQ[j,b]) * L[b]
        c1 = np.array([
            2.0 * np.sum((Q0[i] - Q0[j]) * (dQ[i] - dQ[j]) * L)
            for i, j in pairs
        ])
        # C2[i,j] = sum_b (dQ[i,b]-dQ[j,b])^2 * L[b]
        c2 = np.array([
            np.sum((dQ[i] - dQ[j]) ** 2 * L)
            for i, j in pairs
        ])

        C1_list.append(c1)
        C2_list.append(c2)

    return F2_0, C1_list, C2_list, pairs


def _get_internal_ancestry_vector(
    base_graph, internal_node: str, branch_list: List[Tuple[str, str]]
) -> np.ndarray:
    """
    Compute the ancestry vector of an internal (non-leaf) node.
    Uses path from root to the internal node.
    """
    from ..topology.ancestry_vectors import _path_from_root, _get_root
    tree = base_graph.tree
    n_b = len(branch_list)
    branch_to_idx = {e: i for i, e in enumerate(branch_list)}
    q = np.zeros(n_b)
    for parent, child in _path_from_root(tree, internal_node):
        if (parent, child) in branch_to_idx:
            q[branch_to_idx[(parent, child)]] = 1.0
    return q


class HapGraphLikelihood:
    """
    Data container and likelihood specification for PyMC inference.

    Precomputes quadratic polynomial coefficients for expected F2 as a
    function of alpha (exact, differentiable via PyTensor).
    """

    def __init__(
        self,
        graph,
        fstats,
        ibd=None,
        w_ibd: float = 0.5,
        min_se: float = 1e-6,
    ):
        self.graph = graph
        self.fstats = fstats
        self.ibd = ibd
        self.w_ibd = w_ibd if ibd is not None else 0.0
        self.min_se = min_se

        pops = graph.populations
        n = len(pops)
        K = graph.K

        # Branch lengths: refit at alpha=0 under the K>0 topology.
        #
        # We CANNOT use the NJ K=0 branch lengths directly: they absorbed the
        # admixture signal during fitting, so F2_0(alpha=0) already explains
        # the observed data and the MCMC incorrectly estimates alpha≈0.
        #
        # Instead, we refit branch lengths for the K>0 topology at alpha=0.
        # This gives a "tree-only prediction" F2_0 that is correctly larger
        # than the observed data, so increasing alpha reduces residuals and
        # the MCMC correctly recovers the true alpha.
        if K > 0:
            self.branch_lengths = _refit_branch_lengths_alpha0(graph, fstats)
        else:
            self.branch_lengths = np.array([
                d.get("weight", 0.01)
                for _, _, d in graph.tree.edges(data=True)
            ])
        self.n_branches = len(self.branch_lengths)
        self.K = K

        # F2 observations
        f2_obs_list, f2_se_list, f2_pairs_idx = [], [], []
        for i, pa in enumerate(pops):
            for j, pb in enumerate(pops):
                if j <= i:
                    continue
                obs, se = fstats.f2(pa, pb)
                if not np.isnan(obs):
                    f2_obs_list.append(obs)
                    f2_se_list.append(max(se, min_se))
                    f2_pairs_idx.append((i, j))

        self.f2_obs = np.array(f2_obs_list)
        self.f2_se = np.array(f2_se_list)
        self.f2_pairs_idx = f2_pairs_idx
        self.n_pops = n

        # Precompute polynomial coefficients
        if K > 0:
            F2_0_all, C1_list, C2_list, all_pairs = _compute_f2_polynomial_coeffs(
                graph, self.branch_lengths
            )
            # Restrict to observed pairs
            obs_set = {(i, j) for i, j in f2_pairs_idx}
            obs_mask = np.array([
                ap in obs_set for ap in all_pairs
            ])
            self.F2_0 = F2_0_all[obs_mask]
            self.C1 = np.array([c[obs_mask] for c in C1_list])  # (K, n_obs_pairs)
            self.C2 = np.array([c[obs_mask] for c in C2_list])
        else:
            # No admixture: compute F2 from tree directly
            Q, _, _ = compute_ancestry_vectors(graph)
            F2_mat = expected_f2_matrix(Q, self.branch_lengths)
            self.F2_0 = np.array([F2_mat[i, j] for i, j in f2_pairs_idx])
            self.C1 = np.zeros((0, len(f2_pairs_idx)))
            self.C2 = np.zeros((0, len(f2_pairs_idx)))

        # IBD observations (for timing T inference)
        self.ibd_obs_L = None
        self.ibd_se_L = None
        self.ibd_admix_map = None

        if ibd is not None and K > 0:
            # Build "source-side" population sets for each admixture edge.
            # For timing estimation, only use IBD between the admixed population
            # and populations in the SOURCE clade (including src's sister populations).
            # This avoids bias from high background IBD due to direct tree ancestry,
            # and correctly handles the case where a sister of the true source was selected.
            src_side_sets = []
            for src, tgt, _ in graph.admixture_edges:
                src_clade = _get_source_clade_leaves(graph.tree, src, pops)
                src_side_sets.append(src_clade)

            ibd_L_list, ibd_se_list, admix_map = [], [], []
            for m, (src, tgt, _) in enumerate(graph.admixture_edges):
                src_side = src_side_sets[m]
                tgt_idx = pops.index(tgt) if tgt in pops else -1
                if tgt_idx == -1:
                    continue
                for j, pb in enumerate(pops):
                    if pb == tgt or pb not in src_side:
                        continue
                    obs_L = ibd.get(tgt, pb, "mean_len")
                    n_segs = ibd.get(tgt, pb, "n_segs")
                    if np.isnan(obs_L) or obs_L <= 2.0 or n_segs < 5:
                        # Skip pairs with no signal above the 2 cM filter
                        continue
                    ibd_L_list.append(obs_L)
                    # Sigma for the sample mean of IBD lengths, derived from the
                    # truncated-exponential model (L_min = 2 cM, rate = 2T/100):
                    #   E[L | L > 2] = 50/T + 2   ⟹   50/T ≈ obs_L − 2
                    #   Var[L | L > 2] = (50/T)^2  ⟹  std(mean) ≈ (obs_L-2)/√n
                    # Add a model-error floor (background IBD, Ne effects): 2 cM.
                    length_excess = max(obs_L - 2.0, 0.5)
                    sigma_samp = length_excess / max(np.sqrt(n_segs), 1.0)
                    sigma_model = 2.0
                    ibd_se_list.append(float(np.sqrt(sigma_samp**2 + sigma_model**2)))
                    admix_map.append(m)

            if ibd_L_list:
                self.ibd_obs_L = np.array(ibd_L_list)
                self.ibd_se_L = np.array(ibd_se_list)
                self.ibd_admix_map = admix_map

    def add_to_pymc_model(self, alpha_var, T_var, branch_len_var=None):
        """
        Add F-statistics + IBD likelihood to the current PyMC model context.

        F2 expected value is exact quadratic in alpha (PyTensor differentiable):
          exp_F2 = F2_0 + sum_k [alpha_k * C1_k + alpha_k^2 * C2_k]

        IBD timing (truncation-corrected, cross-population admixture IBD):
          Both lineages contribute T meioses each back to the admixture event,
          so the expected segment length for L > L_min (= 2 cM) is:
            E[L | L > L_min] = 100/(2T) + L_min = 50/T + 2  cM
          mean_len(tgt, pop) ~ N(50/T_k + 2, sigma)
        """
        import pymc as pm
        import pytensor.tensor as pt

        K = self.K

        # Compute expected F2 as quadratic polynomial in alpha
        F2_0_t = pt.as_tensor_variable(self.F2_0.astype(np.float64))
        exp_F2 = F2_0_t

        if K > 0:
            for k in range(K):
                c1_k = pt.as_tensor_variable(self.C1[k].astype(np.float64))
                c2_k = pt.as_tensor_variable(self.C2[k].astype(np.float64))
                alpha_k = alpha_var[k]
                exp_F2 = exp_F2 + alpha_k * c1_k + alpha_k ** 2 * c2_k

        f2_obs_t = pt.as_tensor_variable(self.f2_obs.astype(np.float64))
        f2_se_t = pt.as_tensor_variable(self.f2_se.astype(np.float64))

        pm.Normal("ll_f2", mu=exp_F2, sigma=f2_se_t, observed=f2_obs_t)

        if self.ibd_obs_L is not None and len(self.ibd_obs_L) > 0 and K > 0:
            T_for_pair = pt.stack(
                [T_var[self.ibd_admix_map[m]] for m in range(len(self.ibd_obs_L))]
            )
            # E[L | L > 2 cM] = 50/T + 2  (truncation-corrected IBD formula)
            exp_L = 50.0 / T_for_pair + 2.0
            ibd_obs_t = pt.as_tensor_variable(self.ibd_obs_L.astype(np.float64))
            ibd_se_t = pt.as_tensor_variable(self.ibd_se_L.astype(np.float64))
            pm.Normal("ll_ibd", mu=exp_L, sigma=ibd_se_t, observed=ibd_obs_t)
