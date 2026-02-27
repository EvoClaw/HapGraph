"""Phase 4a Prototype Test: HapGraph on msprime S2 scenario."""

import sys
sys.path.insert(0, "/home/yanlin/1KGothers/code")

import time
import numpy as np
import msprime

print("=== Phase 4a: HapGraph Prototype Test ===\n")

TRUE_ALPHA = 0.3
TRUE_T     = 20
N_SAMPLES  = 20
N_POPS     = 6
pop_names  = [f"pop{i}" for i in range(N_POPS)]
Ne         = 10_000

# ----------------------------------------------------------
# 1. Simulate
# Tree: ((pop0,pop1)@50, pop2)@100 and (pop3,pop4)@80
# pop5 is admixed from pop4 (70%) and pop3 (30%) at T=20 gen ago
# After admixture, pop5 lineages -> pop4 and pop3 (both merge at t=80)
# ----------------------------------------------------------
print("[1] Simulating demographic history...")
d = msprime.Demography()

for name in pop_names:
    d.add_population(name=name, initial_size=Ne)
for name in ["anc01", "anc34", "anc012", "root"]:
    d.add_population(name=name, initial_size=Ne)

# Admixture: at T=20 gen ago, pop5 lineages go 70% to pop4, 30% to pop3
d.add_admixture(
    time=TRUE_T,
    derived="pop5",
    ancestral=["pop4", "pop3"],
    proportions=[1 - TRUE_ALPHA, TRUE_ALPHA],
)

# Tree mergers (pop5 is now empty after admixture, only pop0-pop4 merge)
d.add_population_split(time=50,  derived=["pop0", "pop1"],  ancestral="anc01")
d.add_population_split(time=80,  derived=["pop3", "pop4"],  ancestral="anc34")
d.add_population_split(time=100, derived=["anc01", "pop2"], ancestral="anc012")
d.add_population_split(time=200, derived=["anc012", "anc34"], ancestral="root")

samples = {name: N_SAMPLES for name in pop_names}

t0 = time.time()
ts = msprime.sim_ancestry(
    samples=samples, demography=d,
    sequence_length=30_000_000, recombination_rate=1e-8, random_seed=42,
)
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=42)
t_sim = time.time() - t0
print(f"   Done: {ts.num_sites} variants, {t_sim:.1f}s")

leaf_pop_map = {ts.population(i).id: ts.population(i).metadata['name'] for i in range(N_POPS)}

# ----------------------------------------------------------
# 2. F-statistics
# ----------------------------------------------------------
print("\n[2] Computing F-statistics...")
from hapgraph.preprocess.f_stats import fstats_from_msprime_ts

t0 = time.time()
fstats = fstats_from_msprime_ts(ts, leaf_pop_map)
t_fstats = time.time() - t0
print(f"   Done {t_fstats:.2f}s")

f2_mat = fstats.f2_matrix()
f2_53 = fstats.f2("pop5", "pop3")[0]
f2_54 = fstats.f2("pop5", "pop4")[0]
f2_50 = fstats.f2("pop5", "pop0")[0]
print(f"   F2(pop5,pop3)={f2_53:.5f}  F2(pop5,pop4)={f2_54:.5f}  F2(pop5,pop0)={f2_50:.5f}")
print(f"   Sanity check: pop5 should be closer to pop3/pop4 than pop0  -> {'OK' if f2_53 < f2_50 else 'CHECK'}")

# ----------------------------------------------------------
# 3. IBD statistics
# ----------------------------------------------------------
print("\n[3] Computing IBD statistics...")
from hapgraph.preprocess.ibd_stats import ibd_stats_from_msprime_ts

t0 = time.time()
ibd = ibd_stats_from_msprime_ts(ts, leaf_pop_map, min_len_cM=2.0)
t_ibd = time.time() - t0
print(f"   Done {t_ibd:.1f}s")

for pa, pb in [("pop5","pop3"), ("pop5","pop4"), ("pop3","pop4"), ("pop5","pop0")]:
    mean_L = ibd.get(pa, pb, "mean_len")
    n_s    = ibd.get(pa, pb, "n_segs")
    print(f"   IBD {pa}-{pb}: mean_len={mean_L:.2f}cM  n_segs={n_s:.0f}")

# Expect: pop5-pop3 and pop5-pop4 IBD mean_len >> pop5-pop0
ibd_54 = ibd.get("pop5","pop4","mean_len")
ibd_50 = ibd.get("pop5","pop0","mean_len")
print(f"   IBD sanity: pop5-pop4 ({ibd_54:.2f}cM) >> pop5-pop0 ({ibd_50:.2f}cM)  -> {'OK' if (not np.isnan(ibd_54) and not np.isnan(ibd_50) and ibd_54 > ibd_50) else 'NaN or CHECK'}")

# ----------------------------------------------------------
# 4. NJ tree
# ----------------------------------------------------------
print("\n[4] NJ tree...")
from hapgraph.topology.nj_tree import nj_tree_from_f2, tree_to_newick

t0 = time.time()
nj_tree = nj_tree_from_f2(f2_mat, pop_names)
t_nj = time.time() - t0
try:
    newick = tree_to_newick(nj_tree)
    print(f"   Done {t_nj:.4f}s | {newick[:100]}")
except Exception as e:
    print(f"   Done {t_nj:.4f}s (newick: {e})")

# ----------------------------------------------------------
# 5. Greedy search
# ----------------------------------------------------------
print("\n[5] Greedy admixture search...")
from hapgraph.topology.greedy_search import AdmixGraph, greedy_admixture_search

init_graph = AdmixGraph(tree=nj_tree, admixture_edges=[], populations=pop_names)

t0 = time.time()
top_graphs = greedy_admixture_search(
    init_graph, fstats, ibd, k_max=2, n_alpha_grid=5, w_ibd=0.3, top_k=3, verbose=True
)
t_greedy = time.time() - t0
best = top_graphs[0]
print(f"   Done {t_greedy:.1f}s | Best edges: {best.admixture_edges}")
if best.admixture_edges:
    src, tgt, alpha_est = best.admixture_edges[0]
    correct_tgt = tgt == "pop5"
    correct_src = src in ["pop3", "pop4"]
    print(f"   Detected: {src}->{tgt} alpha={alpha_est:.2f}  (true: pop3->pop5 alpha={TRUE_ALPHA})")
    print(f"   Topology: target={'PASS' if correct_tgt else 'FAIL'} source={'PASS' if correct_src else 'PARTIAL'}")

# ----------------------------------------------------------
# 6. MCMC
# ----------------------------------------------------------
print("\n[6] PyMC NUTS (2 chains x 500 draws)...")
from hapgraph.inference.likelihood import HapGraphLikelihood
from hapgraph.inference.mcmc import run_mcmc

lik = HapGraphLikelihood(best, fstats, ibd, w_ibd=0.5)
t0 = time.time()
result = run_mcmc(lik, chains=2, draws=500, tune=500, random_seed=42, verbose=True)
t_mcmc = time.time() - t0

print(f"\n   MCMC done {t_mcmc:.1f}s  R-hat={result['rhat_max']:.3f}  ESS={result['ess_min']:.0f}  diverg={result['n_divergences']}")

if len(result["alpha_mean"]) > 0:
    a, a_lo, a_hi = result["alpha_mean"][0], result["alpha_hdi"][0,0], result["alpha_hdi"][0,1]
    T, T_lo, T_hi = result["T_mean"][0], result["T_hdi"][0,0], result["T_hdi"][0,1]
    print(f"\n   === Posterior Estimates ===")
    print(f"   alpha: {a:.3f} [{a_lo:.3f},{a_hi:.3f}]  (true={TRUE_ALPHA})")
    print(f"   T:     {T:.1f}  [{T_lo:.1f},{T_hi:.1f}] gen  (true={TRUE_T})")
    print(f"   alpha MAE={abs(a-TRUE_ALPHA):.3f}  T MAE={abs(T-TRUE_T):.1f}gen")
    print(f"   True alpha in CI: {a_lo<=TRUE_ALPHA<=a_hi}")
    print(f"   True T in CI:     {T_lo<=TRUE_T<=T_hi}")

print(f"\n=== Timing: sim={t_sim:.1f}s fstats={t_fstats:.1f}s ibd={t_ibd:.1f}s nj={t_nj:.4f}s greedy={t_greedy:.1f}s mcmc={t_mcmc:.1f}s ===")
print("Phase 4a complete.")
