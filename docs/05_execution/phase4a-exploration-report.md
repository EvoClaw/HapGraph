# Phase 4a Exploration Report (Type C)
# HapGraph Prototype Test
# Date: 2026-02-26

## Quick Benchmark: msprime S2 scenario
Scenario: N=6 populations, K=1 admixture edge
True: pop3 -> pop5, alpha=0.3, T=20 generations

### Component Results

| Component | Status | Time | Notes |
|-----------|--------|------|-------|
| F-statistics (scikit-allel) | PASS | 0.21s | Correct F2 values; sanity check OK |
| IBD statistics (msprime batch) | PASS | 1.9s | mean_len ~ 4.1 cM (theoretical 5 cM for T=20) |
| NJ tree construction | PASS | 0.0002s | Correct topology: pop5 near pop3/pop4 |
| Greedy admixture search | FAIL | 0.5s | F2-based likelihood doesn't detect admixture |
| PyMC NUTS (fixed topology) | PASS | 12.9s | R-hat=1.01, ESS=744, 0 divergences |
| T (timing) recovery | PASS | — | MAE=5.8 gen (target <= 10) ✓ |
| alpha recovery | FAIL | — | MAE=0.20 (target <= 0.05); alpha ~ prior flat |

### Installation Test
- scikit-allel 1.3.13: installed ✓
- PyMC 5.28.0: installed ✓ (pip)
- msprime 1.3.4: installed ✓
- networkx 3.5: installed ✓
- TreeMix: NOT installed (needed for Phase 4b comparison)
- AdmixtureBayes: NOT installed (needed for Phase 4b comparison)

## Claimed Advantage
Claimed: joint F-stats + IBD improves inference over F-only
Actual (Phase 4a): IBD component DOES provide T signal; F-stats component NOT YET providing alpha signal
Gap: PARTIAL — the IBD timing claim holds; alpha estimation needs F3 statistics

## Issues Found

### ISSUE 1 (CRITICAL): F2 likelihood doesn't constrain alpha
- Root cause: Simple tree-path approximation for expected F2 is insensitive to admixture
- F2(A,B) = path_length doesn't change significantly with admixture correction
- Solution: Implement proper admixture graph COVARIANCE MATRIX (e.g., Patterson 2012 / TreeMix approach)
- Formula: E[f2(A,B)] = sum of squared branch-length contributions, accounting for admixture fractions
- This requires computing the "admixture graph matrix" where each population's allele freq
  is a linear combination of ancestral populations

### ISSUE 2 (IMPORTANT): Greedy search uses F2 signal → also needs F3
- Root cause: Greedy search only checks F2 log-likelihood improvement
- F3(admixed; pop1, pop2) < 0 is a much stronger signal for admixture detection
- Solution: Add F3-based admixture detection to greedy search
  - Screen population triples for negative F3
  - Use as a prior for which edges to try in greedy search

### ISSUE 3 (MINOR): NJ tree puts pop5 as sister to pop4 (not fully correct topology)
- In truth, pop5 is admixed from pop3+pop4; it shouldn't be in a clade with either alone
- NJ tree can't represent admixture, so this is expected
- Greedy search should add the admixture edge after the NJ backbone

## Assessment
The core idea is sound: IBD timing works. The F-statistics component needs a proper
admixture graph covariance matrix implementation. This is a known mathematical challenge
(same as used in TreeMix/AdmixtureBayes) and is achievable in Phase 4b.

## Surprises
- msprime ts.ibd_segments(between=...) is extremely fast (1.9s for 6 pops, 30Mb)
- PyMC NUTS converges well (ESS=744, R-hat=1.01) in 12.9s with 2 chains
- F2 path-length approximation is too simplistic for admixture detection
  (known limitation, not unexpected for experts)
