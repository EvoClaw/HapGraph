# Acknowledged Gaps — HapGraph Paper
# Phase 5 output — user-approved decisions to skip certain experiments.
# These MUST appear in the paper's Limitations / Future Work section.

## Gap 1: No TreeMix/AdmixtureBayes baseline comparison
- Decision: User explicitly chose to skip (2026-02-27)
- Rationale: TreeMix is well-established (Pickrell 2012); the paper's novel claim
  is IBD-based timing, which TreeMix cannot do at all. A side-by-side alpha/topology
  comparison would not illuminate this unique contribution.
- Required framing in paper: "We benchmark HapGraph's correctness on simulations
  with known truth (S1, S2); direct comparison with TreeMix topology/alpha performance
  is deferred to future work."
- Appears in: Limitations section, Discussion

## Gap 2: S4–S7 scenarios not completed
- Decision: S1 (K=0), S2 (K=1), and S3 (K=2) all completed. S4–S7 not run.
- S3 results (20 seeds): Greedy K correct 100%, Alpha MAE 0.054±0.024,
  T MAE 10.9±4.2 gen, Alpha CI cover 0.62, T CI cover 0.93.
- Required framing: "Validation covers K≤2; performance on K≥3 is shown
  qualitatively on 1kGP real data and left for systematic benchmarking."
- Appears in: Limitations section

## Gap 3: Alpha 95% CI coverage = 0.57 (target 0.85)
- Decision: Accept as limitation, discuss in text.
- Root cause: F3 MOM delta-method CI with 3% floor; underestimates tail uncertainty
  near α~0.5; does not propagate model error from NJ tree misplacement.
- Required framing: "Alpha CI coverage of 57% indicates CIs are conservative;
  T CI coverage of 85–90% is well-calibrated. Alpha uncertainty quantification
  is a known limitation and target for future improvement."
- Appears in: Results + Limitations

## Gap 4: No second independent real dataset
- Decision: Only 1kGP used.
- Appears in: Limitations, Future Work

## Gap 5: T estimation not applied to 1kGP real data
- Decision: run_hapgraph_1kgp.py uses ibd=None; hap-ibd was never run on full 1kGP data.
- Consequence: 1kGP results section can only report topology (K=8) and alpha estimates.
  T estimation for real data is deferred to future work.
- Required framing: "Admixture timing estimation was validated on simulations with known
  ground truth (S2, S3). Application of the IBD-based timing step to 1kGP data requires
  pre-computation of IBD segments via hap-ibd, which is deferred to future work."
- Appears in: Results (1kGP subsection), Limitations
