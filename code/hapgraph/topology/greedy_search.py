"""
Greedy admixture edge search for population graph topology inference.

Algorithm:
1. Start from a Neighbor-Joining tree (from F2 distances).
2. Greedily add admixture edges one at a time:
   - Use F3 scanning to identify admixture-candidate populations (F3(C;A,B) < 0).
   - For each (src, tgt, alpha) candidate: fit branch lengths via NNLS,
     compute log-likelihood improvement.
   - Add the edge giving the largest improvement.
3. Repeat up to K_max times, stopping when BIC no longer improves.
4. Return the top-K candidate topologies (by log-likelihood).

Key fix over v1: uses ancestry-vector expected F2 + scipy NNLS branch-length fitting
instead of the broken path-length approximation.
"""

import copy
import numpy as np
import networkx as nx
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass, field

from ..preprocess.f_stats import FStatistics
from ..preprocess.ibd_stats import IBDStatistics
from .ancestry_vectors import (
    compute_ancestry_vectors,
    expected_f2_matrix,
    expected_f3_array,
)


@dataclass
class AdmixGraph:
    """
    Represents a population admixture graph.

    Attributes
    ----------
    tree : nx.DiGraph — the tree backbone (edges parent→child, weight=branch_length)
    admixture_edges : list of (source, target, alpha) tuples
        source: population that contributed admixture
        target: population that received it
        alpha: admixture proportion (0 < alpha < 1)
    populations : list of leaf node names
    log_lik : float — log-likelihood of this graph
    n_params : int — number of free parameters
    """

    tree: nx.DiGraph
    admixture_edges: List[Tuple[str, str, float]] = field(default_factory=list)
    populations: List[str] = field(default_factory=list)
    log_lik: float = -np.inf
    n_params: int = 0

    @property
    def K(self) -> int:
        return len(self.admixture_edges)

    def bic(self, n_stats: int) -> float:
        """BIC = k*ln(n) - 2*ln(L)."""
        return self.n_params * np.log(max(n_stats, 1)) - 2 * self.log_lik

    def copy(self) -> "AdmixGraph":
        return AdmixGraph(
            tree=copy.deepcopy(self.tree),
            admixture_edges=list(self.admixture_edges),
            populations=list(self.populations),
            log_lik=self.log_lik,
            n_params=self.n_params,
        )


# ---------------------------------------------------------------------------
# Branch-length fitting
# ---------------------------------------------------------------------------

def _build_f2_design_matrix(
    Q: np.ndarray,
    pairs: List[Tuple[int, int]],
) -> np.ndarray:
    """
    Build design matrix A such that A @ branch_lengths ≈ obs_f2_flat.

    A[k, b] = (Q[i, b] - Q[j, b])^2  for pairs[k] = (i, j).
    """
    n_pairs = len(pairs)
    n_b = Q.shape[1]
    A = np.zeros((n_pairs, n_b))
    for k, (i, j) in enumerate(pairs):
        dq = Q[i, :] - Q[j, :]
        A[k, :] = dq * dq
    return A


def fit_branch_lengths(
    Q: np.ndarray,
    obs_f2_flat: np.ndarray,
    pairs: List[Tuple[int, int]],
    min_len: float = 1e-6,
) -> np.ndarray:
    """
    Fit non-negative branch lengths via NNLS.

    Solves  min || A @ L - obs_f2_flat ||^2  s.t. L >= 0
    where A[k, b] = (Q[i,b] - Q[j,b])^2  for pair k = (i,j).

    Returns branch_lengths of shape (n_branches,).
    """
    from scipy.optimize import nnls

    A = _build_f2_design_matrix(Q, pairs)
    L, _ = nnls(A, obs_f2_flat)
    return np.maximum(L, min_len)


def _collect_f2_obs(
    fstats: FStatistics,
    pops: List[str],
    min_se: float = 1e-8,
) -> Tuple[List[Tuple[int, int]], np.ndarray, np.ndarray]:
    """Extract observed F2 pairs, values, and SEs for a population list."""
    pairs, obs, ses = [], [], []
    for i, pa in enumerate(pops):
        for j, pb in enumerate(pops):
            if j <= i:
                continue
            val, se = fstats.f2(pa, pb)
            if np.isnan(val):
                continue
            pairs.append((i, j))
            obs.append(val)
            ses.append(max(se, min_se))
    return pairs, np.array(obs), np.array(ses)


# ---------------------------------------------------------------------------
# Log-likelihood
# ---------------------------------------------------------------------------

def f_stats_log_likelihood(
    graph: AdmixGraph,
    fstats: FStatistics,
    return_fitted_lengths: bool = False,
) -> float | Tuple[float, np.ndarray]:
    """
    Log-likelihood of F2 statistics under the admixture graph.

    1. Computes ancestry vectors Q for the current topology + admixture.
    2. Fits branch lengths via NNLS (non-negative least squares).
    3. Evaluates Gaussian log-likelihood: obs_F2 ~ N(exp_F2, se^2).

    Optionally returns fitted branch lengths alongside the log-likelihood.
    """
    Q, pops, branch_list = compute_ancestry_vectors(graph)
    pairs, obs_flat, se_flat = _collect_f2_obs(fstats, pops)

    if len(pairs) == 0:
        result = (0.0, np.ones(len(branch_list)) * 0.01)
        return result if return_fitted_lengths else result[0]

    L = fit_branch_lengths(Q, obs_flat, pairs)

    f2_exp_mat = expected_f2_matrix(Q, L)
    ll = 0.0
    for k, (i, j) in enumerate(pairs):
        ll += -0.5 * ((obs_flat[k] - f2_exp_mat[i, j]) / se_flat[k]) ** 2

    if return_fitted_lengths:
        return ll, L
    return ll


def ibd_log_likelihood(
    graph: AdmixGraph,
    ibd: IBDStatistics,
) -> float:
    """
    IBD-based log-likelihood during greedy search (topology signal only).

    Uses IBD rate as a proxy for genetic closeness: populations with an admixture
    edge should share more IBD than predicted by tree distance alone.
    During greedy search, T is not yet estimated; we use a fixed-scale exponential.
    """
    ll = 0.0
    pops = graph.populations
    undirected = graph.tree.to_undirected()

    for i, pa in enumerate(pops):
        for j, pb in enumerate(pops):
            if j <= i:
                continue
            obs_rate = ibd.get(pa, pb, "total_rate")
            if np.isnan(obs_rate):
                continue
            try:
                tree_dist = nx.shortest_path_length(undirected, pa, pb, weight="weight")
            except nx.NetworkXNoPath:
                continue
            for src, tgt, alpha in graph.admixture_edges:
                if (tgt == pa and src == pb) or (tgt == pb and src == pa):
                    tree_dist *= max(1.0 - alpha * 0.5, 0.1)
            scale = 0.02
            exp_rate = np.exp(-tree_dist / scale)
            sigma = max(obs_rate * 0.3, 1e-6)
            ll += -0.5 * ((obs_rate - exp_rate) / sigma) ** 2
    return ll


def joint_log_likelihood(
    graph: AdmixGraph,
    fstats: FStatistics,
    ibd: Optional[IBDStatistics] = None,
    w_ibd: float = 0.3,
) -> float:
    """Joint log-likelihood combining F-stats and (optionally) IBD signals."""
    ll_f = f_stats_log_likelihood(graph, fstats)
    if ibd is not None:
        ll_ibd = ibd_log_likelihood(graph, ibd)
        return (1.0 - w_ibd) * ll_f + w_ibd * ll_ibd
    return ll_f


# ---------------------------------------------------------------------------
# F3 scanning for admixture candidate detection
# ---------------------------------------------------------------------------

def scan_f3_admixture_signals(
    fstats: FStatistics,
    z_threshold: float = -2.0,
) -> List[Tuple[str, str, str, float]]:
    """
    Identify admixture-candidate populations via F3 statistics.

    F3(C; A, B) < 0 (significantly negative Z-score) is a strong signal
    that population C is admixed between A-like and B-like sources.

    Returns list of (C, A, B, z_score) sorted by most negative z_score.
    """
    pops = fstats.populations
    signals = []
    for c in pops:
        for i, a in enumerate(pops):
            for b in pops[i+1:]:
                # All three must be distinct; C is the test pop, A and B are sources
                if a == c or b == c:
                    continue
                obs, se = fstats.f3(c, a, b)
                if np.isnan(obs) or se <= 0:
                    continue
                z = obs / se
                if z < z_threshold:
                    signals.append((c, a, b, z))
    signals.sort(key=lambda x: x[3])
    return signals


# ---------------------------------------------------------------------------
# Candidate generation
# ---------------------------------------------------------------------------

def _candidate_admixture_edges(
    graph: AdmixGraph,
    f3_signals: Optional[List[Tuple[str, str, str, float]]] = None,
    max_candidates: int = 200,
) -> List[Tuple[str, str]]:
    """
    Generate candidate (source, target) admixture edges.

    If F3 signals are available, prioritize those implied by negative F3:
    for (C, A, B, z), the target is C and sources are A or B or their ancestors.
    Fall back to all-pairs if fewer than 5 F3-guided candidates exist.
    """
    all_nodes = list(graph.tree.nodes)
    leaves = graph.populations
    existing = {(s, t) for s, t, _ in graph.admixture_edges}

    f3_candidates: List[Tuple[str, str]] = []
    if f3_signals:
        # Use the best (most negative) F3 signal for EACH unique target population.
        # This ensures every admixed population gets candidates, not just the one
        # with the globally strongest signal.
        best_per_target: dict = {}
        for c, a, b, z in f3_signals:
            if c not in best_per_target or z < best_per_target[c][2]:
                best_per_target[c] = (a, b, z)

        # Collect candidates: for each target, propose its top sources + all internal nodes
        for c, (a, b, _) in best_per_target.items():
            for src in [a, b] + [n for n in all_nodes if n not in leaves]:
                if src != c and (src, c) not in existing:
                    f3_candidates.append((src, c))

    if len(f3_candidates) >= 5:
        seen = set()
        ordered = []
        for pair in f3_candidates:
            if pair not in seen:
                seen.add(pair)
                ordered.append(pair)
        return ordered[:max_candidates]

    # Fallback: all pairs
    all_pairs = [
        (src, tgt)
        for src in all_nodes
        for tgt in leaves
        if src != tgt and (src, tgt) not in existing
    ]
    return all_pairs[:max_candidates]


# ---------------------------------------------------------------------------
# Main greedy search
# ---------------------------------------------------------------------------

def greedy_admixture_search(
    init_graph: AdmixGraph,
    fstats: FStatistics,
    ibd: Optional[IBDStatistics] = None,
    k_max: int = 4,
    n_alpha_grid: int = 9,
    w_ibd: float = 0.3,
    top_k: int = 5,
    verbose: bool = True,
) -> List[AdmixGraph]:
    """
    Greedy admixture-edge search using ancestry-vector F2 + NNLS branch fitting.

    1. Scan F3 signals to identify candidate admixture targets.
    2. For each candidate (src → tgt, alpha), compute ancestry vectors,
       fit branch lengths via NNLS, evaluate log-likelihood.
    3. Add the edge with largest positive Δlog-lik.
    4. Repeat up to k_max times, stopping when Δlog-lik < 0.5.

    Parameters
    ----------
    init_graph : starting AdmixGraph (NJ tree, K=0)
    fstats : FStatistics object with .f2(), .f3() methods
    ibd : IBDStatistics object (or None for F-only mode)
    k_max : maximum admixture edges to add
    n_alpha_grid : alpha values tested per candidate (grid 0.05–0.95)
    w_ibd : weight for IBD component in joint likelihood
    top_k : number of candidate topologies to return
    verbose : print progress

    Returns
    -------
    List[AdmixGraph] sorted by log-likelihood (best first)
    """
    alpha_grid = np.linspace(0.05, 0.95, n_alpha_grid)
    n_stats = len(fstats.populations) * (len(fstats.populations) - 1) // 2

    # Score the initial tree
    current = init_graph.copy()
    base_ll, base_L = f_stats_log_likelihood(current, fstats, return_fitted_lengths=True)
    # Store fitted branch lengths back onto tree edges
    edges_list = list(current.tree.edges())
    for b, (u, v) in enumerate(edges_list):
        current.tree[u][v]["weight"] = float(base_L[b])

    if ibd is not None:
        base_ll = base_ll * (1 - w_ibd) + ibd_log_likelihood(current, ibd) * w_ibd

    current.log_lik = base_ll
    current.n_params = len(edges_list) + len(current.admixture_edges)

    all_candidates = [current]

    # Pre-scan F3 signals to guide candidate search
    f3_signals = scan_f3_admixture_signals(fstats)
    if verbose and f3_signals:
        print(f"  F3 scan: {len(f3_signals)} negative-F3 signals detected.")
        for c, a, b, z in f3_signals[:5]:
            print(f"    F3({c}; {a}, {b}): Z = {z:.2f}")

    for k in range(1, k_max + 1):
        if verbose:
            print(f"\n  Greedy step {k}/{k_max}: testing candidate edges...")

        best_delta = -np.inf
        best_proposal = None

        candidates = _candidate_admixture_edges(current, f3_signals)
        n_tested = 0

        # Each population may be the target of at most one admixture edge.
        # This prevents the greedy search from over-fitting a single strongly
        # admixed population (e.g. ASW) with multiple edges, at the cost of
        # missing weaker signals in other populations.
        already_tgt = {t for _, t, _ in current.admixture_edges}

        for src, tgt in candidates:
            if tgt in already_tgt:
                continue
            for alpha in alpha_grid:
                proposal = current.copy()
                proposal.admixture_edges.append((src, tgt, alpha))
                proposal.n_params = current.n_params + 2  # +alpha, +T

                ll_f, _ = f_stats_log_likelihood(proposal, fstats, return_fitted_lengths=True)
                if ibd is not None:
                    ll_ibd = ibd_log_likelihood(proposal, ibd)
                    ll = (1 - w_ibd) * ll_f + w_ibd * ll_ibd
                else:
                    ll = ll_f

                delta = ll - current.log_lik
                if delta > best_delta:
                    best_delta = delta
                    best_proposal = proposal
                    best_proposal.log_lik = ll
                n_tested += 1

        if verbose:
            print(f"  Tested {n_tested} (src, tgt, alpha) combinations.")
            print(f"  Best Δlog-lik = {best_delta:.6f}")

        # BIC-based stopping criterion.
        #
        # Adding one admixture edge introduces 1 new free parameter (alpha; branch
        # lengths are refitted by NNLS at every step so they don't count as new params
        # in the BIC sense). We require "strong evidence" for the new edge:
        #
        #   ΔlogL > 0.5 × k × ln(n_stats)   [BIC improvement]
        #
        # With an additional safety margin of +5 to guard against false positives:
        #   min_delta = 5.0 + 0.5 × ln(n_stats)
        #
        # For 7 populations (21 F2 pairs):  threshold ≈ 5 + 1.5 = 6.5
        # For 26 populations (325 F2 pairs): threshold ≈ 5 + 2.9 = 7.9
        #
        # This replaces the old hardcoded threshold of 0.05, which was calibrated
        # to an SE that was ~49× too large (block jackknife bug).
        min_delta = 5.0 + 0.5 * np.log(max(n_stats, 2))
        if best_proposal is None or best_delta < min_delta:
            if verbose:
                print(f"  Stopping: Δlog-lik {best_delta:.3f} below BIC threshold {min_delta:.3f}.")
            break

        # Store fitted branch lengths onto winning proposal
        winning_ll, winning_L = f_stats_log_likelihood(best_proposal, fstats, return_fitted_lengths=True)
        for b, (u, v) in enumerate(list(best_proposal.tree.edges())):
            best_proposal.tree[u][v]["weight"] = float(winning_L[b])

        if verbose:
            src, tgt, alpha = best_proposal.admixture_edges[-1]
            print(f"  Added: {src} → {tgt} (α={alpha:.2f}), Δlog-lik = {best_delta:.4f}")

        current = best_proposal
        all_candidates.append(current.copy())

    all_candidates.sort(key=lambda g: g.log_lik, reverse=True)
    return all_candidates[:top_k]
