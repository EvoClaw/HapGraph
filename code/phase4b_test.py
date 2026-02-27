"""
Phase 4b Prototype Test: HapGraph with cross-clade admixture scenario.

Scenario:
  Tree: ((pop0,pop1)@200, (pop2,pop3)@200)@400, with pop4 as outgroup @600
  pop5 = admixed from pop0 (alpha=0.3) + pop2 (0.7) at T=20 gen ago

This cross-clade admixture should give:
  F3(pop5; pop0, pop2) << 0  (clearly negative, strong Z-score)
  NJ tree places pop5 awkwardly between clades -> detectable
"""

import sys
sys.path.insert(0, "/home/yanlin/1KGothers/code")

import time
import numpy as np
import msprime

print("=== Phase 4b: HapGraph Cross-Clade Admixture Test ===\n")

TRUE_ALPHA = 0.3   # pop0 contribution to pop5
TRUE_T     = 20    # generations ago
N_SAMPLES  = 30
Ne         = 10_000   # Standard Ne; deeper divergence gives larger F2
pop_names  = [f"pop{i}" for i in range(6)]

# ----------------------------------------------------------
# 1. Simulate: cross-clade admixture at inter-continental scale
# Tree: ((pop0,pop1)@2500, (pop2,pop3)@2500)@5000, pop4 @7500
# F2(pop0,pop2) ≈ 0.02-0.03 (Africa-Europe scale) → F3 clearly negative
# pop5 = alpha*pop0 + (1-alpha)*pop2 at T=20
# ----------------------------------------------------------
print("[1] Simulating cross-clade demographic history (inter-continental scale)...")
d = msprime.Demography()

for name in pop_names:
    d.add_population(name=name, initial_size=Ne)
for name in ["anc01", "anc23", "anc0123", "root"]:
    d.add_population(name=name, initial_size=Ne)

# pop5 admixed from pop0 (alpha=0.3) and pop2 (0.7) at T=20
d.add_admixture(
    time=TRUE_T,
    derived="pop5",
    ancestral=["pop0", "pop2"],
    proportions=[TRUE_ALPHA, 1.0 - TRUE_ALPHA],
)

# Tree mergers — 5x deeper divergence than before
d.add_population_split(time=2500, derived=["pop0", "pop1"],  ancestral="anc01")
d.add_population_split(time=2500, derived=["pop2", "pop3"],  ancestral="anc23")
d.add_population_split(time=5000, derived=["anc01", "anc23"], ancestral="anc0123")
d.add_population_split(time=7500, derived=["anc0123", "pop4"], ancestral="root")

samples = {name: N_SAMPLES for name in pop_names}

t0 = time.time()
ts = msprime.sim_ancestry(
    samples=samples, demography=d,
    sequence_length=50_000_000, recombination_rate=1e-8, random_seed=42,
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)
t_sim = time.time() - t0
print(f"   Done: {ts.num_sites} variants, {t_sim:.1f}s")

leaf_pop_map = {ts.population(i).id: ts.population(i).metadata['name'] for i in range(6)}

# ----------------------------------------------------------
# 2. F-statistics + sanity checks
# ----------------------------------------------------------
print("\n[2] Computing F-statistics...")
from hapgraph.preprocess.f_stats import fstats_from_msprime_ts

t0 = time.time()
fstats = fstats_from_msprime_ts(ts, leaf_pop_map)
t_fstats = time.time() - t0
print(f"   Done {t_fstats:.2f}s")

print("\n   F2 matrix (diag=0):")
for a in pop_names:
    for b in pop_names:
        if b > a:
            obs, se = fstats.f2(a, b)
            print(f"   F2({a},{b}) = {obs:.5f}")

print("\n   F3 signals for pop5:")
for a in pop_names[:5]:
    for b in pop_names[:5]:
        if a >= b: continue
        obs, se = fstats.f3("pop5", a, b)
        z = obs/se if se>0 else 0
        marker = " *** NEGATIVE ***" if z < -2 else ""
        print(f"   F3(pop5; {a}, {b}) = {obs:.6f}  Z={z:.2f}{marker}")

# ----------------------------------------------------------
# 3. IBD statistics
# ----------------------------------------------------------
print("\n[3] Computing IBD statistics...")
from hapgraph.preprocess.ibd_stats import ibd_stats_from_msprime_ts

t0 = time.time()
ibd = ibd_stats_from_msprime_ts(ts, leaf_pop_map, min_len_cM=2.0)
t_ibd = time.time() - t0
print(f"   Done {t_ibd:.1f}s")

for pa, pb in [("pop5","pop0"), ("pop5","pop2"), ("pop0","pop2"), ("pop5","pop3"), ("pop5","pop4")]:
    mean_L = ibd.get(pa, pb, "mean_len")
    n_s    = ibd.get(pa, pb, "n_segs")
    print(f"   IBD {pa}-{pb}: mean_len={mean_L:.2f}cM  n_segs={n_s:.0f}")

# ----------------------------------------------------------
# 4. NJ tree
# ----------------------------------------------------------
print("\n[4] NJ tree...")
from hapgraph.topology.nj_tree import nj_tree_from_f2, tree_to_newick

f2_mat = fstats.f2_matrix()
t0 = time.time()
nj_tree = nj_tree_from_f2(f2_mat, pop_names)
t_nj = time.time() - t0
print(f"   Done {t_nj:.4f}s")
print(f"   Newick: {tree_to_newick(nj_tree)[:200]}")

# ----------------------------------------------------------
# 5. Greedy admixture search
# ----------------------------------------------------------
print("\n[5] Greedy admixture search (ancestry-vector + NNLS)...")
from hapgraph.topology.greedy_search import (
    AdmixGraph, greedy_admixture_search, f_stats_log_likelihood,
    scan_f3_admixture_signals
)

# Show F3 scan result
signals = scan_f3_admixture_signals(fstats, z_threshold=-2.0)
print(f"   F3 scan found {len(signals)} signals with Z < -2:")
for c, a, b, z in signals[:8]:
    print(f"     F3({c}; {a}, {b}): Z={z:.2f}")

init_graph = AdmixGraph(tree=nj_tree, admixture_edges=[], populations=pop_names)

# Show LL for true edge
ll_base, _ = f_stats_log_likelihood(init_graph, fstats, return_fitted_lengths=True)
true_graph = init_graph.copy()
true_graph.admixture_edges = [("pop0", "pop5", 0.3)]
ll_true, _ = f_stats_log_likelihood(true_graph, fstats, return_fitted_lengths=True)
print(f"\n   LL base: {ll_base:.4f}  LL(true edge pop0->pop5, a=0.3): {ll_true:.4f}  Delta={ll_true-ll_base:.4f}")

t0 = time.time()
top_graphs = greedy_admixture_search(
    init_graph, fstats, ibd, k_max=2, n_alpha_grid=19, w_ibd=0.3, top_k=3, verbose=True
)
t_greedy = time.time() - t0
best = top_graphs[0]
print(f"\n   Done {t_greedy:.1f}s | Best edges: {best.admixture_edges}")

if best.admixture_edges:
    src, tgt, alpha_est = best.admixture_edges[0]
    correct_tgt = tgt == "pop5"
    correct_src = src in ["pop0", "pop2", "anc01", "anc23"]
    print(f"   Detected: {src}->{tgt} alpha={alpha_est:.2f}")
    print(f"   True: pop0->pop5 alpha={TRUE_ALPHA}")
    print(f"   Topology: target={'PASS' if correct_tgt else 'FAIL'}  source={'PASS' if correct_src else 'PARTIAL'}")
else:
    print("   No admixture detected!")

# ----------------------------------------------------------
# 6. MCMC
# ----------------------------------------------------------
print("\n[6] PyMC NUTS (2 chains x 500 draws)...")
from hapgraph.inference.likelihood import HapGraphLikelihood
from hapgraph.inference.mcmc import run_mcmc

# Use best graph if admixture detected, else inject true edge for MCMC test
mcmc_graph = best if best.admixture_edges else true_graph
if not best.admixture_edges:
    print("   [Using true edge for MCMC test since greedy failed]")

lik = HapGraphLikelihood(mcmc_graph, fstats, ibd, w_ibd=0.5)
t0 = time.time()
result = run_mcmc(lik, chains=2, draws=500, tune=500, random_seed=42, verbose=True)
t_mcmc = time.time() - t0

print(f"\n   MCMC done {t_mcmc:.1f}s  R-hat={result['rhat_max']:.3f}  "
      f"ESS={result['ess_min']:.0f}  diverg={result['n_divergences']}")

if len(result["alpha_mean"]) > 0:
    a      = result["alpha_mean"][0]
    a_lo   = result["alpha_hdi"][0, 0]
    a_hi   = result["alpha_hdi"][0, 1]
    T      = result["T_mean"][0]
    T_lo   = result["T_hdi"][0, 0]
    T_hi   = result["T_hdi"][0, 1]
    print(f"\n   === Posterior Estimates ===")
    print(f"   alpha: {a:.3f} [{a_lo:.3f},{a_hi:.3f}]  (true={TRUE_ALPHA})")
    print(f"   T:     {T:.1f}  [{T_lo:.1f},{T_hi:.1f}] gen  (true={TRUE_T})")
    print(f"   alpha MAE = {abs(a - TRUE_ALPHA):.3f}  (target <=0.05)")
    print(f"   T MAE     = {abs(T - TRUE_T):.1f} gen   (target <=10)")
    print(f"   True alpha in CI: {a_lo <= TRUE_ALPHA <= a_hi}")
    print(f"   True T in CI:     {T_lo <= TRUE_T <= T_hi}")

print(f"\n=== Timing: sim={t_sim:.1f}s fstats={t_fstats:.1f}s ibd={t_ibd:.1f}s "
      f"nj={t_nj:.4f}s greedy={t_greedy:.1f}s mcmc={t_mcmc:.1f}s ===")
print("\nPhase 4b basic test complete.")
