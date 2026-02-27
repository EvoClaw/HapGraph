"""
F-statistics computation using scikit-allel.
Computes F2, F3, F4 (Patterson statistics) and block-jackknife standard errors.
No R/rpy2 dependency.
"""

import numpy as np
from itertools import combinations
from typing import Dict, List, Tuple, Optional
import allel


def compute_allele_counts(
    gt: allel.GenotypeArray,
    pop_indices: Dict[str, List[int]],
) -> Dict[str, np.ndarray]:
    """
    Compute per-population allele counts from a genotype array.

    Parameters
    ----------
    gt : allel.GenotypeArray, shape (n_variants, n_samples, 2)
    pop_indices : dict mapping population name -> list of sample indices

    Returns
    -------
    ac : dict mapping population name -> allele_counts array (n_variants, 2)
    """
    ac = {}
    for pop, idx in pop_indices.items():
        ac[pop] = gt.take(idx, axis=1).count_alleles()
    return ac


def _allele_freqs(ac: np.ndarray) -> np.ndarray:
    """Derived allele frequency from allele count array (n_variants, 2)."""
    total = ac.sum(axis=1)
    with np.errstate(invalid="ignore"):
        freq = ac[:, 1] / total
    # Mask monomorphic or missing sites
    freq = np.where(total == 0, np.nan, freq)
    return freq


def _f2_single(p_a: np.ndarray, p_b: np.ndarray) -> float:
    """
    Unbiased F2 estimator between populations A and B.
    F2(A,B) = E[(p_A - p_B)^2] (up to correction for sample size,
    but for large samples the correction is small; we use the simple version).
    """
    mask = ~(np.isnan(p_a) | np.isnan(p_b))
    if mask.sum() == 0:
        return np.nan
    d = (p_a[mask] - p_b[mask]) ** 2
    return float(np.mean(d))


def _f3_single(
    p_c: np.ndarray,
    p_a: np.ndarray,
    p_b: np.ndarray,
    n_c_haps: Optional[int] = None,
) -> float:
    """
    Bias-corrected F3(C; A, B) estimator (Patterson et al. 2012).

    Raw F3 = E[(p_C - p_A)(p_C - p_B)] is positively biased by sampling
    variance in C by approximately mean[p_C*(1-p_C)] / (n_C_haps - 1).
    This correction is applied when n_c_haps is provided.

    Negative F3 indicates admixture.
    """
    mask = ~(np.isnan(p_c) | np.isnan(p_a) | np.isnan(p_b))
    if mask.sum() == 0:
        return np.nan
    raw = float(np.mean((p_c[mask] - p_a[mask]) * (p_c[mask] - p_b[mask])))
    if n_c_haps is not None and n_c_haps > 1:
        het_c = float(np.mean(p_c[mask] * (1.0 - p_c[mask])))
        correction = het_c / (n_c_haps - 1)
        return raw - correction
    return raw


def _f4_single(
    p_a: np.ndarray, p_b: np.ndarray, p_c: np.ndarray, p_d: np.ndarray
) -> float:
    """
    F4(A, B; C, D) = E[(p_A - p_B)(p_C - p_D)]
    Also known as the D-statistic (when standardised).
    """
    mask = ~(
        np.isnan(p_a) | np.isnan(p_b) | np.isnan(p_c) | np.isnan(p_d)
    )
    if mask.sum() == 0:
        return np.nan
    return float(
        np.mean((p_a[mask] - p_b[mask]) * (p_c[mask] - p_d[mask]))
    )


def block_jackknife_se(values: np.ndarray, n_blocks: int = 50) -> float:
    """
    Block-jackknife standard error of the mean (Patterson et al. 2012 formula).

    For the mean estimator θ̂ = mean(values), the jackknife SE is:
        SE = std(block_means) / sqrt(n_blocks)

    This is derived from the pseudovalue formula: since the pseudovalue of a mean
    estimator equals the block mean itself (pseudo_i = block_mean_i), the SE reduces
    to the standard error of the mean of block means.

    This correctly accounts for LD between adjacent SNPs because each block spans
    ~1-5 Mb (typical LD decay distance), making block means approximately independent.

    Parameters
    ----------
    values : 1-D array of per-SNP contributions (e.g. (p_A - p_B)^2 for F2)
    n_blocks : number of jackknife blocks (50 per chromosome is standard)

    Returns
    -------
    se : float, jackknife standard error of the mean
    """
    n = len(values)
    if n < n_blocks:
        return float(np.std(values) / np.sqrt(max(n, 1)))
    block_size = n // n_blocks
    n_use = block_size * n_blocks
    # block_means[i] = mean of per-site values in block i
    block_means = values[:n_use].reshape(n_blocks, block_size).mean(axis=1)
    # SE = std(block_means) / sqrt(n_blocks) — standard SE of the mean of block means
    se = float(np.std(block_means, ddof=1) / np.sqrt(n_blocks))
    return se


class FStatistics:
    """
    Container for all pairwise F-statistics among a set of populations.
    Computed from allele frequency arrays.
    """

    def __init__(
        self,
        freqs: Dict[str, np.ndarray],
        n_haps: Optional[Dict[str, int]] = None,
        n_blocks: int = 50,
    ):
        """
        Parameters
        ----------
        freqs : dict mapping pop_name -> allele_frequency array (n_variants,)
        n_haps : dict mapping pop_name -> number of haplotypes (for F3 bias correction)
        n_blocks : blocks for jackknife SE
        """
        self.freqs = freqs
        self.n_haps = n_haps or {}
        self.populations = list(freqs.keys())
        self.n_blocks = n_blocks
        self._f2: Dict[Tuple[str, str], Tuple[float, float]] = {}
        self._f3: Dict[Tuple[str, str, str], Tuple[float, float]] = {}
        self._f4: Dict[Tuple[str, str, str, str], Tuple[float, float]] = {}
        self._compute_f2()

    def _compute_f2(self):
        for a, b in combinations(self.populations, 2):
            pa = self.freqs[a]
            pb = self.freqs[b]
            mask = ~(np.isnan(pa) | np.isnan(pb))
            vals = (pa[mask] - pb[mask]) ** 2
            est = float(vals.mean())
            se = block_jackknife_se(vals, self.n_blocks)
            self._f2[(a, b)] = (est, se)
            self._f2[(b, a)] = (est, se)

    def f2(self, a: str, b: str) -> Tuple[float, float]:
        """Return (F2_estimate, jackknife_SE) for populations a, b."""
        key = (a, b) if (a, b) in self._f2 else (b, a)
        return self._f2[key]

    def f2_matrix(self) -> np.ndarray:
        """Return symmetric N x N F2 distance matrix."""
        n = len(self.populations)
        mat = np.zeros((n, n))
        for i, a in enumerate(self.populations):
            for j, b in enumerate(self.populations):
                if i != j:
                    mat[i, j] = self._f2[(a, b)][0]
        return mat

    def f3(self, c: str, a: str, b: str) -> Tuple[float, float]:
        """
        Return (F3_estimate, SE) for admixture test F3(C; A, B).

        Uses bias-corrected estimator (Patterson et al. 2012) when n_haps[c]
        is available. Negative F3 indicates C is admixed between A- and B-like
        lineages.
        """
        key = (c, a, b)
        if key not in self._f3:
            pc = self.freqs[c]
            pa = self.freqs[a]
            pb = self.freqs[b]
            mask = ~(np.isnan(pc) | np.isnan(pa) | np.isnan(pb))
            raw_vals = (pc[mask] - pa[mask]) * (pc[mask] - pb[mask])
            raw_est = float(raw_vals.mean())

            # Bias correction for sampling variance in C
            n_c = self.n_haps.get(c, None)
            if n_c is not None and n_c > 1:
                het_c = float(np.mean(pc[mask] * (1.0 - pc[mask])))
                correction = het_c / (n_c - 1)
                est = raw_est - correction
                # Jackknife on the corrected per-site values
                corr_vals = raw_vals - pc[mask] * (1.0 - pc[mask]) / (n_c - 1)
                se = block_jackknife_se(corr_vals, self.n_blocks)
            else:
                est = raw_est
                se = block_jackknife_se(raw_vals, self.n_blocks)

            self._f3[key] = (est, se)
        return self._f3[key]

    def f4(self, a: str, b: str, c: str, d: str) -> Tuple[float, float]:
        """Return (F4_estimate, SE)."""
        key = (a, b, c, d)
        if key not in self._f4:
            pa = self.freqs[a]
            pb = self.freqs[b]
            pc = self.freqs[c]
            pd = self.freqs[d]
            mask = ~(
                np.isnan(pa) | np.isnan(pb) | np.isnan(pc) | np.isnan(pd)
            )
            vals = (pa[mask] - pb[mask]) * (pc[mask] - pd[mask])
            est = float(vals.mean())
            se = block_jackknife_se(vals, self.n_blocks)
            self._f4[key] = (est, se)
        return self._f4[key]


def fstats_from_genotype_array(
    gt: allel.GenotypeArray,
    pop_indices: Dict[str, List[int]],
    n_blocks: int = 50,
) -> FStatistics:
    """
    Convenience function: genotype array → FStatistics object.

    Parameters
    ----------
    gt : allel.GenotypeArray
    pop_indices : dict pop_name -> list of sample column indices
    n_blocks : jackknife blocks

    Returns
    -------
    FStatistics instance
    """
    freqs = {}
    for pop, idx in pop_indices.items():
        ac = gt.take(idx, axis=1).count_alleles()
        total = ac.sum(axis=1)
        with np.errstate(invalid="ignore"):
            freq = np.where(total == 0, np.nan, ac[:, 1] / total)
        freqs[pop] = freq
    return FStatistics(freqs, n_blocks=n_blocks)


def fstats_from_msprime_ts(ts, pop_map: Dict[int, str]) -> FStatistics:
    """
    Compute FStatistics directly from an msprime TreeSequence.
    Uses allele frequencies computed from the tree sequence genotype matrix.

    Parameters
    ----------
    ts : msprime TreeSequence
    pop_map : dict mapping population_id (int) -> population_name (str)

    Returns
    -------
    FStatistics instance
    """
    import msprime  # noqa

    # Build sample index map: pop_name -> list of haplotype indices
    pop_indices_hap: Dict[str, List[int]] = {name: [] for name in pop_map.values()}
    for ind in ts.individuals():
        pop_id = ts.node(ind.nodes[0]).population
        if pop_id in pop_map:
            pop_indices_hap[pop_map[pop_id]].extend(list(ind.nodes))

    # Genotype matrix: (n_variants, n_haplotypes)
    G = ts.genotype_matrix()  # shape (n_variants, n_haplotypes)

    freqs: Dict[str, np.ndarray] = {}
    n_haps: Dict[str, int] = {}
    for pop_name, hap_idx in pop_indices_hap.items():
        if len(hap_idx) == 0:
            continue
        g_pop = G[:, hap_idx]  # (n_variants, n_haps_in_pop)
        freqs[pop_name] = g_pop.mean(axis=1)  # derived allele frequency
        n_haps[pop_name] = len(hap_idx)

    return FStatistics(freqs, n_haps=n_haps)
