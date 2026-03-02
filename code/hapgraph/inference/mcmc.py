"""
Bayesian inference of alpha and T for fixed topology, using PyMC 5 + NUTS.

Priors
------
alpha ~ Beta(1, 1)    — uniform on [0, 1], shape (K,)
T_raw ~ HalfNormal(σ=20)   — T = T_raw + 5, so T ≥ 5 generations
    (prior mean ≈ 5 + 20*√(2/π) ≈ 21 generations)

IBD likelihood (via HapGraphLikelihood.add_to_pymc_model):
    E[L | L > L_min] = 50/T + 2  cM   (L_min = 2 cM truncation correction)
    Both admixed and source lineages contribute T meioses back to the admixture
    event, so the effective rate parameter is 2T/100 per cM, giving a mean of
    100/(2T) = 50/T cM for an untruncated segment.  With L_min = 2 cM the
    conditional mean shifts to 50/T + 2 cM.
"""

import numpy as np
from typing import Optional, Dict, Any
import warnings


def run_mcmc(likelihood, chains=4, draws=1000, tune=1000,
             target_accept=0.9, random_seed=42, verbose=True):
    """
    NUTS sampling for alpha and T given fixed topology.
    Returns dict with trace, posterior summaries, and diagnostics.
    """
    import pymc as pm
    import arviz as az
    import pytensor.tensor as pt

    K = likelihood.K

    with pm.Model():
        if K > 0:
            alpha = pm.Beta("alpha", alpha=1.0, beta=1.0, shape=(K,))
            T_raw = pm.HalfNormal("T_raw", sigma=20.0, shape=(K,))
            T = pm.Deterministic("T", T_raw + 5.0)
        else:
            alpha = pt.zeros(0)
            T = pt.zeros(0)

        likelihood.add_to_pymc_model(alpha, T)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            trace = pm.sample(
                draws=draws, tune=tune, chains=chains,
                target_accept=target_accept, random_seed=random_seed,
                progressbar=verbose, return_inferencedata=True,
            )

    summary = az.summary(trace, var_names=["alpha", "T"] if K > 0 else [])
    if verbose:
        print("\n=== MCMC Convergence Summary ===")
        if not summary.empty:
            print(summary[["mean", "sd", "hdi_3%", "hdi_97%", "r_hat", "ess_bulk"]])

    rhat_max = float(summary["r_hat"].max()) if not summary.empty else 1.0
    ess_min  = float(summary["ess_bulk"].min()) if not summary.empty else 0.0
    n_div    = int(trace.sample_stats["diverging"].sum())

    if K > 0:
        alpha_flat = trace.posterior["alpha"].values.reshape(-1, K)
        T_flat     = trace.posterior["T"].values.reshape(-1, K)
        alpha_mean = alpha_flat.mean(axis=0)
        alpha_hdi  = np.percentile(alpha_flat, [2.5, 97.5], axis=0).T
        T_mean     = T_flat.mean(axis=0)
        T_hdi      = np.percentile(T_flat, [2.5, 97.5], axis=0).T
    else:
        alpha_mean = np.array([])
        alpha_hdi  = np.zeros((0, 2))
        T_mean     = np.array([])
        T_hdi      = np.zeros((0, 2))

    return dict(trace=trace, alpha_mean=alpha_mean, alpha_hdi=alpha_hdi,
                T_mean=T_mean, T_hdi=T_hdi, rhat_max=rhat_max,
                ess_min=ess_min, n_divergences=n_div)
