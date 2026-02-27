"""Quick F3 significant signal scan for 1kGP data."""
import sys, pickle
import numpy as np

sys.path.insert(0, '/home/yanlin/1KGothers/code')

FSTATS_PKL = '/home/yanlin/1KGothers/results/realdata/fstats_chr1_22_n50.pkl'

with open(FSTATS_PKL, 'rb') as f:
    fstats = pickle.load(f)
pops = fstats.populations

print('Significant admixture signals F3(C; A, B) Z < -3.0:')
print(f'  {"Pop":<6} {"Src_A":<6} {"Src_B":<6} {"F3":>10} {"SE":>8} {"Z":>8}')

sig_f3 = []
for ta in pops:
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
for z, ta, sa, sb, obs, se in sig_f3[:40]:
    print(f'  {ta:<6} {sa:<6} {sb:<6} {obs:>10.5f} {se:>8.5f} {z:>8.1f}')

# Summary: unique target populations with significant F3
targets = set()
for z, ta, sa, sb, obs, se in sig_f3:
    targets.add(ta)
print(f'\nAdmixed populations detected (F3 Z < -3.0): {sorted(targets)}')
print(f'Total significant F3 triplets: {len(sig_f3)}')
