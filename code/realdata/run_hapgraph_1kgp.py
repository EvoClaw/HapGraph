"""
Phase 4b-C: Run HapGraph on 1kGP chr1+chr22 data.

Pipeline:
  1. Load precomputed F-statistics (26 pops, chr1+chr22)
  2. Build NJ tree from F2 matrix
  3. Greedy admixture search (BIC threshold)
  4. Estimate alpha for each admixture event via F3 MOM
  5. Save results and print summary
"""
import sys, json, pickle, time
import numpy as np

sys.path.insert(0, '/home/yanlin/1KGothers/code')

from hapgraph.preprocess.f_stats import FStatistics
from hapgraph.topology.nj_tree import nj_tree_from_f2, tree_to_newick
from hapgraph.topology.greedy_search import greedy_admixture_search, AdmixGraph, f_stats_log_likelihood
from hapgraph.inference.f3_estimator import f3_mom_alpha, find_best_sources

FSTATS_PKL = '/home/yanlin/1KGothers/results/realdata/fstats_chr1_22_n50.pkl'
OUT_JSON   = '/home/yanlin/1KGothers/results/realdata/hapgraph_1kgp_result.json'

print('=' * 60)
print('HapGraph on 1kGP (chr1 + chr22, 26 populations)')
print('=' * 60)

# ── 1. Load F-statistics ─────────────────────────────────────
print('\n[1] Loading F-statistics...')
with open(FSTATS_PKL, 'rb') as f:
    fstats = pickle.load(f)
pops = fstats.populations
print(f'  {len(pops)} populations: {pops}')

# ── 2. NJ tree ───────────────────────────────────────────────
print('\n[2] Building NJ tree from F2 matrix...')
f2_mat = fstats.f2_matrix()
nj = nj_tree_from_f2(f2_mat, pops)
newick = tree_to_newick(nj)
print(f'  NJ tree:\n  {newick}')

# ── 3. Greedy admixture search ───────────────────────────────
print('\n[3] Greedy admixture search (BIC threshold)...')
t0 = time.time()
init_graph = AdmixGraph(tree=nj, admixture_edges=[], populations=pops)
candidates = greedy_admixture_search(init_graph, fstats, ibd=None, k_max=8,
                                     w_ibd=0.0, verbose=True)
# Select BIC-optimal graph (lowest BIC = best tradeoff of fit vs complexity)
n_stats = len(pops) * (len(pops) - 1) // 2
bic_scores = [(g.bic(n_stats), g) for g in candidates]
best_graph = min(bic_scores, key=lambda x: x[0])[1]
print(f'\n  BIC scores: ' + ', '.join(f'K={g.K}:{bic:.1f}' for bic, g in bic_scores[:6]))
print(f'  BIC-optimal: K={best_graph.K}')
history = [{'K': g.K, 'bic': float(g.bic(n_stats))} for g in candidates]
elapsed = time.time() - t0
print(f'\n  Best graph: K={best_graph.K} admixture edges')
print(f'  Runtime: {elapsed:.1f}s')

# ── 4. F3 MOM alpha for each admixture edge ──────────────────
print('\n[4] F3 MOM alpha estimation for each admixture event...')
admixture_events = []
for k, (src, tgt, alpha_init, *rest) in enumerate(best_graph.admixture_edges):
    print(f'\n  Edge {k+1}: {src} -> {tgt}')

    # Find best source pair by most-negative F3
    best_pair = find_best_sources(fstats, tgt, z_threshold=-2.0)
    if best_pair is not None:
        src_a, src_b, ae, alo, ahi, z = best_pair
        print(f'    Best F3 pair: F3({tgt}; {src_a}, {src_b})  Z={z:.1f}')
        print(f'    Alpha MOM: {ae:.3f}  95%CI [{alo:.3f}, {ahi:.3f}]')
    else:
        src_a, src_b, ae, alo, ahi, z = src, tgt, 0.25, 0.0, 0.5, 0.0
        print(f'    No significant F3 signal found for {tgt}')

    # Note: src may be an internal node; skip direct F3 check with it

    admixture_events.append({
        'edge_idx': k,
        'src_greedy': src,
        'tgt': tgt,
        'alpha_greedy': float(alpha_init),
        'src_a_f3': src_a,
        'src_b_f3': src_b,
        'alpha_mom': float(ae),
        'alpha_ci_lo': float(alo),
        'alpha_ci_hi': float(ahi),
        'f3_z': float(z),
    })

# ── 5. Summary table ─────────────────────────────────────────
print('\n' + '=' * 60)
print('ADMIXTURE SUMMARY (26 populations, chr1+chr22)')
print('=' * 60)
print(f'  NJ tree: {len(pops)} populations')
print(f'  Admixture edges: K = {best_graph.K}')
print()
print(f'  {"Tgt":<6} {"Src_A":<6} {"Src_B":<6} {"Alpha_MOM":>10} {"95% CI":>18} {"F3_Z":>8}')
print(f'  {"-"*6} {"-"*6} {"-"*6} {"-"*10} {"-"*18} {"-"*8}')
for ev in admixture_events:
    ci = f'[{ev["alpha_ci_lo"]:.3f}, {ev["alpha_ci_hi"]:.3f}]'
    print(f'  {ev["tgt"]:<6} {ev["src_a_f3"]:<6} {ev["src_b_f3"]:<6} '
          f'{ev["alpha_mom"]:>10.3f} {ci:>18} {ev["f3_z"]:>8.1f}')

# ── 6. All significant F3 signals ────────────────────────────
print('\n[6] All significant admixture signals (F3 Z < -3.0):')
sig_f3 = []
for i, ta in enumerate(pops):
    for j, sa in enumerate(pops):
        if sa == ta:
            continue
        for k2, sb in enumerate(pops):
            if k2 <= j or sb == ta:
                continue
            obs, se = fstats.f3(ta, sa, sb)
            if np.isnan(obs) or se <= 0:
                continue
            z = obs / se
            if z < -3.0:
                sig_f3.append((z, ta, sa, sb, obs, se))

sig_f3.sort()
print(f'  {"Pop":<6} {"Src_A":<6} {"Src_B":<6} {"F3":>10} {"SE":>8} {"Z":>8}')
for z, ta, sa, sb, obs, se in sig_f3[:30]:
    print(f'  {ta:<6} {sa:<6} {sb:<6} {obs:>10.5f} {se:>8.5f} {z:>8.1f}')

# ── 7. Save ──────────────────────────────────────────────────
    result = {
    'n_pops': len(pops),
    'populations': pops,
    'nj_newick': newick,
    'K': best_graph.K,
    'admixture_edges': [
        {'src': s, 'tgt': t, 'alpha_greedy': float(a)}
        for s, t, a in best_graph.admixture_edges
    ],
    'admixture_events': admixture_events,
    'runtime_s': elapsed,
}
with open(OUT_JSON, 'w') as f:
    json.dump(result, f, indent=2)
print(f'\nSaved: {OUT_JSON}')
