"""
Shared simulation utilities for HapGraph benchmarks.
Generates msprime tree sequences for scenarios S1-S7.

Single-simulation architecture (v3 — final):
  Each scenario runs ONE msprime tree sequence per seed, used for BOTH
  F-statistics (Stage I) and IBD detection (Stage II).

  Key design parameters
  ---------------------
  Ne = 1_000          — small effective size; coalescence fast → fewer ARG trees
  split times         — proportional to Ne (same F2 ratios as Ne=10_000 designs)
  seq_len  = 50 Mbp   — physical length
  recomb   = 1e-8 /bp — 0.5 Morgan equivalent

  IBD genome-length:  genome_len_morgan = 50e6 × 1e-8 = 0.5 Morgan
  bp_per_cM          = 50e6 / (0.5 × 100) = 1,000,000 bp/cM  ✓
  min_span_bp (2 cM) = 2,000,000 bp = 2 Mbp  ✓

  At T=20, α=0.3, this simulation produces ~100-200 IBD segments > 2 cM
  between the admixed population and the source clade (n=30 per pop),
  sufficient for stable mean-length estimation.

Benchmarks show ~0.5 s/seed for simulation (vs ~20 s with Ne=10,000).

S3 design: 3 clades (A:3, B:3, C:2 pops), two independent cross-clade events.
  popM1 = alpha1*popA0 + (1-alpha1)*popB0,  T1 ~ DiscreteUniform{10,...,49}
  popM2 = alpha2*popC0 + (1-alpha2)*popB1,  T2 ~ DiscreteUniform{10,...,49}
"""

import numpy as np
import msprime
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

# Simulation parameters  (Ne=1000, proportionally scaled split times)
_NE            = 1_000
_SEQ_LEN       = 50_000_000    # bp
_RECOMB        = 1e-8          # per bp per generation
_MUT           = 1e-8          # per bp per generation

# IBD genome length for ibd_stats_from_msprime_ts
IBD_GENOME_LEN_MORGAN = _SEQ_LEN * _RECOMB   # = 0.5 Morgan


@dataclass
class SimResult:
    ts: object                                           # tree sequence with mutations
    pop_names: List[str]
    pop_map: Dict[int, str]
    admixture_edges: List[Tuple[str, str, float, float]]
    # Each edge: (source, target, alpha_true, T_true)
    n_pops: int
    seed: int
    scenario: str
    ibd_genome_len_morgan: float = field(default=IBD_GENOME_LEN_MORGAN)

    # Convenience alias: same ts used for IBD detection
    @property
    def ts_ibd(self):
        return self.ts


def sim_s1_tree(
    n_pops: int = 6,
    n_samples: int = 30,
    seed: int = 42,
) -> SimResult:
    """S1: Pure tree, no admixture. Balanced 2-clade structure."""
    pop_names = [f"pop{i}" for i in range(n_pops)]
    d = msprime.Demography()
    for name in pop_names:
        d.add_population(name=name, initial_size=_NE)
    for name in ["anc01", "anc23", "anc45", "root"]:
        d.add_population(name=name, initial_size=_NE)

    d.add_population_split(time=250, derived=["pop0", "pop1"], ancestral="anc01")
    d.add_population_split(time=250, derived=["pop2", "pop3"], ancestral="anc23")
    if n_pops >= 6:
        d.add_population_split(time=250, derived=["pop4", "pop5"], ancestral="anc45")
        d.add_population_split(time=500, derived=["anc01", "anc23", "anc45"], ancestral="root")
    else:
        d.add_population_split(time=500, derived=["anc01", "anc23"], ancestral="root")

    samples = {n: n_samples for n in pop_names}

    ts = msprime.sim_ancestry(
        samples=samples, demography=d,
        sequence_length=_SEQ_LEN, recombination_rate=_RECOMB, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=_MUT, random_seed=seed)

    pop_map = {ts.population(i).id: ts.population(i).metadata['name']
               for i in range(len(pop_names))}
    return SimResult(
        ts=ts, pop_names=pop_names, pop_map=pop_map,
        admixture_edges=[], n_pops=n_pops, seed=seed, scenario="S1",
    )


def sim_s2_cross_clade(
    n_pops: int = 7,
    n_samples: int = 30,
    seed: int = 42,
    t_admix_range: Tuple[int, int] = (10, 50),
) -> SimResult:
    """
    S2: One cross-clade admixture event.
    Tree: ((popA0,popA1,popA2)@250, (popB0,popB1,popB2)@250)@500
    popM = alpha*popA0 + (1-alpha)*popB0,  T ~ DiscreteUniform{10,...,49}
    """
    rng = np.random.default_rng(seed)
    alpha_true = float(rng.uniform(0.1, 0.5))
    T_true = float(rng.integers(*t_admix_range))

    clade_a = ["popA0", "popA1", "popA2"]
    clade_b = ["popB0", "popB1", "popB2"]
    admixed = ["popM"]
    all_pops = clade_a + clade_b + admixed

    d = msprime.Demography()
    for name in all_pops:
        d.add_population(name=name, initial_size=_NE)
    for name in ["ancA", "ancB", "root"]:
        d.add_population(name=name, initial_size=_NE)

    d.add_admixture(
        time=T_true, derived="popM",
        ancestral=["popA0", "popB0"],
        proportions=[alpha_true, 1.0 - alpha_true],
    )
    d.add_population_split(time=250, derived=clade_a, ancestral="ancA")
    d.add_population_split(time=250, derived=clade_b, ancestral="ancB")
    d.add_population_split(time=500, derived=["ancA", "ancB"], ancestral="root")

    samples = {n: n_samples for n in all_pops}

    ts = msprime.sim_ancestry(
        samples=samples, demography=d,
        sequence_length=_SEQ_LEN, recombination_rate=_RECOMB, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=_MUT, random_seed=seed)

    pop_map = {ts.population(i).id: ts.population(i).metadata['name']
               for i in range(len(all_pops))}
    return SimResult(
        ts=ts, pop_names=all_pops, pop_map=pop_map,
        admixture_edges=[("popA0", "popM", alpha_true, T_true)],
        n_pops=len(all_pops), seed=seed, scenario="S2",
    )


def sim_s3_two_admixture(
    n_samples: int = 30,
    seed: int = 42,
    t_admix_range: Tuple[int, int] = (10, 50),
) -> SimResult:
    """
    S3: Two independent cross-clade admixture events, 10 populations.

    Tree backbone: 3 clades (A:3 pops, B:3 pops, C:2 pops)
      ((A,B)@500, C)@800
    Admixture:
      popM1 = alpha1*popA0 + (1-alpha1)*popB0, T1~DiscreteUniform{10,...,49}
      popM2 = alpha2*popC0 + (1-alpha2)*popB1, T2~DiscreteUniform{10,...,49}

    Tests K=2 recovery: both distinct admixture events must be found,
    with correct sources and independent alpha/T estimates.
    """
    rng = np.random.default_rng(seed)
    alpha1 = float(rng.uniform(0.1, 0.5))
    T1 = int(rng.integers(*t_admix_range))
    alpha2 = float(rng.uniform(0.1, 0.5))
    T2 = int(rng.integers(*t_admix_range))

    # Ensure T2 != T1 to avoid msprime ordering collision
    if T2 == T1:
        T2 = T1 + 1

    clade_a = ["popA0", "popA1", "popA2"]
    clade_b = ["popB0", "popB1", "popB2"]
    clade_c = ["popC0", "popC1"]
    admixed = ["popM1", "popM2"]
    all_pops = clade_a + clade_b + clade_c + admixed

    d = msprime.Demography()
    for name in all_pops:
        d.add_population(name=name, initial_size=_NE)
    for name in ["ancA", "ancB", "ancC", "ancAB", "root"]:
        d.add_population(name=name, initial_size=_NE)

    d.add_admixture(
        time=T1, derived="popM1",
        ancestral=["popA0", "popB0"],
        proportions=[alpha1, 1.0 - alpha1],
    )
    d.add_admixture(
        time=T2, derived="popM2",
        ancestral=["popC0", "popB1"],
        proportions=[alpha2, 1.0 - alpha2],
    )

    d.add_population_split(time=250, derived=clade_a, ancestral="ancA")
    d.add_population_split(time=250, derived=clade_b, ancestral="ancB")
    d.add_population_split(time=250, derived=clade_c, ancestral="ancC")
    d.add_population_split(time=500, derived=["ancA", "ancB"], ancestral="ancAB")
    d.add_population_split(time=800, derived=["ancAB", "ancC"], ancestral="root")

    d.sort_events()

    samples = {n: n_samples for n in all_pops}

    ts = msprime.sim_ancestry(
        samples=samples, demography=d,
        sequence_length=_SEQ_LEN, recombination_rate=_RECOMB, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=_MUT, random_seed=seed)

    pop_map = {ts.population(i).id: ts.population(i).metadata['name']
               for i in range(len(all_pops))}
    return SimResult(
        ts=ts, pop_names=all_pops, pop_map=pop_map,
        admixture_edges=[
            ("popA0", "popM1", alpha1, float(T1)),
            ("popC0", "popM2", alpha2, float(T2)),
        ],
        n_pops=len(all_pops), seed=seed, scenario="S3",
    )
