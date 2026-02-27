import numpy as np
from typing import Tuple


def f3_mom_alpha(fstats, tgt: str, src_a: str, src_b: str) -> Tuple[float, float, float]:
    """
    F3 method-of-moments alpha estimator (Patterson et al. 2012).

    For tgt = alpha*src_a + (1-alpha)*src_b:
        F3(tgt; src_a, src_b) = -alpha*(1-alpha) * F2(src_a, src_b)
    Solving: alpha = (1 - sqrt(1 - 4*(-F3/F2))) / 2

    Topology-free: does not use the NJ tree at all.
    """
    f3_obs, f3_se = fstats.f3(tgt, src_a, src_b)
    f2_obs, f2_se = fstats.f2(src_a, src_b)

    if np.isnan(f3_obs) or np.isnan(f2_obs) or f2_obs <= 0:
        return np.nan, np.nan, np.nan

    x = -f3_obs / f2_obs
    if x <= 0:
        return 0.0, 0.0, 0.0

    disc = 1.0 - 4.0 * x
    alpha_est = (1.0 - np.sqrt(max(disc, 0.0))) / 2.0

    dg_dx = 1.0 / max(np.sqrt(abs(disc)), 1e-6)
    dx_dF3 = -1.0 / f2_obs
    dx_dF2 = f3_obs / f2_obs ** 2
    x_se = np.sqrt((dx_dF3 * f3_se) ** 2 + (dx_dF2 * f2_se) ** 2)
    alpha_se_sampling = dg_dx * x_se

    # Model uncertainty floor: the 2-source admixture approximation has inherent
    # inaccuracy due to drift, background relatedness, and multi-source ancestry.
    # A minimum of 3% absolute uncertainty is conservative and ensures calibrated CIs
    # across both simulation (where sampling SE is artificially small) and real data.
    alpha_se = max(alpha_se_sampling, 0.03)

    return alpha_est, max(0.0, alpha_est - 1.96 * alpha_se), min(1.0, alpha_est + 1.96 * alpha_se)


def find_best_sources(fstats, tgt: str, z_threshold: float = -2.0):
    """Scan all F3(tgt; A, B) triples, return best (src_a, src_b) by most negative Z."""
    pops = [p for p in fstats.populations if p != tgt]
    best = None
    best_z = 0.0
    for i, pa in enumerate(pops):
        for j, pb in enumerate(pops):
            if j <= i:
                continue
            obs, se = fstats.f3(tgt, pa, pb)
            if np.isnan(obs) or se <= 0:
                continue
            z = obs / se
            if z < best_z:
                best_z = z
                ae, alo, ahi = f3_mom_alpha(fstats, tgt, pa, pb)
                best = (pa, pb, ae, alo, ahi, z)
    return best
