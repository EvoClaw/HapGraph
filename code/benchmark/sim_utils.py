"""
Shared simulation utilities for HapGraph benchmarks.
Generates msprime tree sequences for scenarios S1-S7.
"""
# S3 design: 10 populations, K=2 admixture events, two independent cross-clade events.
# Tree backbone: 3 clades (A, B, C), each with 3 populations.
#   Clade A: popA0, popA1, popA2 — split at 2500 gen
#   Clade B: popB0, popB1, popB2 — split at 2500 gen
#   Clade C: popC0, popC1       — split at 2500 gen
#   ((A, B)@5000, C)@8000
# Admixture events:
#   popM1 = alpha1*popA0 + (1-alpha1)*popB0,  T1 ~ Uniform(10,50)
#   popM2 = alpha2*popC0 + (1-alpha2)*popB1,  T2 ~ Uniform(10,50) (independent)
# This tests whether HapGraph can simultaneously recover two independent
# admixture events with different sources and targets.

import numpy as np
import msprime
from typing import Dict, List, Tuple
from dataclasses import dataclass, field


@dataclass
class SimResult:
    ts: object
    pop_names: List[str]
    pop_map: Dict[int, str]
    admixture_edges: List[Tuple[str, str, float, float]]
    # (source, target, alpha_true, T_true)
    n_pops: int
    seed: int
    scenario: str


def sim_s1_tree(
    n_pops: int = 6,
    n_samples: int = 30,
    ne: int = 10_000,
    seq_len: int = 50_000_000,
    seed: int = 42,
) -> SimResult:
    """S1: Pure tree, no admixture. Balanced 2-clade structure."""
    pop_names = [f"pop{i}" for i in range(n_pops)]
    d = msprime.Demography()
    for name in pop_names:
        d.add_population(name=name, initial_size=ne)
    for name in ["anc01", "anc23", "anc45", "root"]:
        d.add_population(name=name, initial_size=ne)

    d.add_population_split(time=2500, derived=["pop0", "pop1"], ancestral="anc01")
    d.add_population_split(time=2500, derived=["pop2", "pop3"], ancestral="anc23")
    if n_pops >= 6:
        d.add_population_split(time=2500, derived=["pop4", "pop5"], ancestral="anc45")
        d.add_population_split(time=5000, derived=["anc01", "anc23", "anc45"], ancestral="root")
    else:
        d.add_population_split(time=5000, derived=["anc01", "anc23"], ancestral="root")

    ts = msprime.sim_ancestry(
        samples={n: n_samples for n in pop_names},
        demography=d, sequence_length=seq_len,
        recombination_rate=1e-8, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)
    pop_map = {ts.population(i).id: ts.population(i).metadata['name'] for i in range(n_pops)}
    return SimResult(ts=ts, pop_names=pop_names, pop_map=pop_map,
                     admixture_edges=[], n_pops=n_pops, seed=seed, scenario="S1")


def sim_s2_cross_clade(
    n_pops: int = 7,
    n_samples: int = 30,
    ne: int = 10_000,
    seq_len: int = 50_000_000,
    seed: int = 42,
    t_admix_range: Tuple[int, int] = (10, 50),
) -> SimResult:
    """
    S2: One cross-clade admixture event.
    Tree: ((popA0,popA1,popA2)@2500, (popB0,popB1,popB2)@2500)@5000
    popM = alpha*popA0 + (1-alpha)*popB0, at T~Uniform(t_admix_range) gen ago.
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
        d.add_population(name=name, initial_size=ne)
    for name in ["ancA", "ancB", "root"]:
        d.add_population(name=name, initial_size=ne)

    d.add_admixture(
        time=T_true, derived="popM",
        ancestral=["popA0", "popB0"],
        proportions=[alpha_true, 1.0 - alpha_true],
    )
    d.add_population_split(time=2500, derived=clade_a, ancestral="ancA")
    d.add_population_split(time=2500, derived=clade_b, ancestral="ancB")
    d.add_population_split(time=5000, derived=["ancA", "ancB"], ancestral="root")

    ts = msprime.sim_ancestry(
        samples={n: n_samples for n in all_pops},
        demography=d, sequence_length=seq_len,
        recombination_rate=1e-8, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)
    pop_map = {ts.population(i).id: ts.population(i).metadata['name']
               for i in range(len(all_pops))}
    return SimResult(
        ts=ts, pop_names=all_pops, pop_map=pop_map,
        admixture_edges=[("popA0", "popM", alpha_true, T_true)],
        n_pops=len(all_pops), seed=seed, scenario="S2",
    )


def sim_s3_two_admixture(
    n_samples: int = 30,
    ne: int = 10_000,
    seq_len: int = 50_000_000,
    seed: int = 42,
    t_admix_range: Tuple[int, int] = (10, 50),
) -> SimResult:
    """
    S3: Two independent cross-clade admixture events, 10 populations.

    Tree backbone: 3 clades (A:3 pops, B:3 pops, C:2 pops)
      ((A,B)@5000, C)@8000
    Admixture:
      popM1 = alpha1*popA0 + (1-alpha1)*popB0, T1~Uniform(10,50)
      popM2 = alpha2*popC0 + (1-alpha2)*popB1, T2~Uniform(10,50)

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
        d.add_population(name=name, initial_size=ne)
    for name in ["ancA", "ancB", "ancC", "ancAB", "root"]:
        d.add_population(name=name, initial_size=ne)

    # Admixture events (must occur before clade splits)
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

    d.add_population_split(time=2500, derived=clade_a, ancestral="ancA")
    d.add_population_split(time=2500, derived=clade_b, ancestral="ancB")
    d.add_population_split(time=2500, derived=clade_c, ancestral="ancC")
    d.add_population_split(time=5000, derived=["ancA", "ancB"], ancestral="ancAB")
    d.add_population_split(time=8000, derived=["ancAB", "ancC"], ancestral="root")

    d.sort_events()

    ts = msprime.sim_ancestry(
        samples={n: n_samples for n in all_pops},
        demography=d, sequence_length=seq_len,
        recombination_rate=1e-8, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)
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
