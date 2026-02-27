# HapGraph — Tool Design Document
# Phase 3 (Type C Path), C-Steps 1-3

## C-Step 1: Utility Claim

### Primary Claim
"HapGraph enables simultaneous population admixture graph inference
(topology + admixture proportions alpha) AND admixture timing estimation (T, in
generations) from a phased VCF file, with full Bayesian uncertainty
quantification — a capability that was previously impossible with any
single existing tool."

### Current Pain (User Workflow)
Current 2-step workflow:
  Step 1: TreeMix / AdmixtureBayes  ->  graph topology + alpha  (no timing, or no uncertainty)
  Step 2: ALDER / ROLLOFF / IBDNe   ->  admixture timing T  (separate analysis, separate tool)
  Problem: No joint uncertainty; topology from step 1 is fixed assumption in step 2

HapGraph 1-step workflow:
  Phased VCF + pop labels -> preprocessing (ADMIXTOOLS2 + hap-ibd)
                          -> graph topology + alpha + T + credible intervals
  Result: Coherent Bayesian model; alpha and T uncertainty jointly propagated

### Quantitative Targets

| Target | Metric | Reference Tool |
|--------|--------|---------------|
| Topology accuracy | >= 5% better Robinson-Foulds distance | vs TreeMix (complex graphs K=2, N=10 pops) |
| Admixture proportion | MAE <= 0.03 for alpha in [0.05, 0.5] | vs AdmixtureBayes (F-only) |
| Timing accuracy | MAE <= 10 generations for T in [5, 100] | No existing graph-based tool has this |
| CI calibration | 95% CI covers true value >= 85% | On simulations with known truth |
| Runtime | < 2 hours for N=26 populations | Single CPU server |

## C-Step 2: Evaluation Framework

### Simulation Scenarios (msprime)

| Scenario | Populations | Admixture edges | T range | Replicates | Purpose |
|----------|-------------|----------------|---------|------------|---------|
| S1 | 6 | 0 | none | 20 | Baseline tree recovery |
| S2 | 8 | 1 | 10-50 gen | 30 | Simple admixture, known timing |
| S3 | 10 | 2 | 5-100 gen | 30 | Complex admixture, 2 timings |
| S4 | 15 | 3 | 3-80 gen | 20 | Larger graph, multiple events |
| S5 | 26 | 4 | 5-150 gen | 10 | Full scale (1kGP-size) |
| S6 | 10 | 2 | <5 gen | 20 | Very recent admixture |

### Ablation Study (REQUIRED)

| Configuration | F-stats | IBD-stats |
|--------------|---------|----------|
| F-only (Bayesian) | YES | NO |
| IBD-only | NO | YES |
| F+IBD (HapGraph full) | YES | YES |
| TreeMix (ML baseline) | YES (ML) | NO |
| AdmixtureBayes (Bayes baseline) | YES (RJMCMC) | NO |

### Real Data Applications

1. 1kGP 2022: 26 populations, primary application
2. Second dataset: TBD (populations with archaeologically-dated admixture for timing validation)

## C-Step 3: Competing Tools

| Tool | Year | Method | Key Weakness | HapGraph Advantage |
|------|------|--------|-------------|-------------------|
| TreeMix | 2012 | ML + F-stats | No timing, no uncertainty | + Timing + Bayesian UQ |
| qpGraph/ADMIXTOOLS 2 | 2009/2022 | F-stats moment-matching | Manual topology, no timing | + Automated + IBD + Timing |
| AdmixtureBayes | 2023 | RJMCMC + F-stats | No timing; slow for N>15 | + IBD + Timing + Greedy (faster for N=26) |
| gLike | 2025 | Genealogical likelihood | Requires tree sequences | + Standard phased VCF input |
| ALDER/ROLLOFF | 2012 | IBD/LD timing | Only timing, no graph | HapGraph integrates into unified framework |

## Statistical Model

### Joint Likelihood
L(G, alpha, T | data) = L_F(F-stats | G, alpha) x L_IBD(IBD-stats | G, alpha, T, N_e) x priors

F-stats: Observed ~ MVN(expected from graph, Sigma_F), Sigma_F by block jackknife
IBD: L-bar(i,j) | T_k, alpha_k ~ E[100/T_k] (mean segment length ~ 100/T cM)
     r(i,j) | alpha_k, T_k, N_e ~ admixture proportion and timing theory
Covariance: Block-diagonal (F block + IBD block); cross-block correlation = sensitivity analysis

### Topology Search Strategy
1. Neighbor-Joining tree from F2 distance matrix
2. Greedy admixture edge addition (maximize joint log-likelihood)
3. Repeat for K = 1, 2, 3, 4 (collect top-5 topologies)
4. Full Bayesian inference (NUTS) on each top topology
5. Model comparison by WAIC

### Bayesian Inference (PyMC + NUTS)
Parameters: alpha (Beta prior), T (HalfNormal prior), N_e (fixed or HalfNormal)
Sampler: NUTS (gradient-based), 4 chains, 2000 warmup + 2000 draws
Convergence: R-hat < 1.05, ESS > 400, trace plot inspection

## Implementation Timeline (single developer, 8-10 months)

Month 1-2: F-stats + IBD pipeline, likelihood functions
Month 3-4: Greedy topology search, basic PyMC model
Month 5-6: Full pipeline, msprime simulation benchmark (S1-S3)
Month 7:   Ablation + tool comparison, timing validation (S4-S6)
Month 8:   1kGP application + second dataset
Month 9:   Packaging (pip/conda), documentation, tutorial
Month 10:  Paper writing

## Design Updates from Phase 3 Deliberation (Round 1)

### CRIT-1: Replace admixr/rpy2 with scikit-allel
- Remove R dependency entirely
- scikit-allel provides Patterson F2, F3, F4 + block-jackknife SEs
- Install chain: Python + NumPy/SciPy/Pandas + PyMC + scikit-allel + matplotlib/networkx

### CRIT-2: IBD Timing Scope Restriction
- Restrict admixture timing claims to T in [5, 50] generations
- Rationale: hap-ibd detects IBD >= 2 cM; mean length = 100/T cM
  - T=50: mean 2cM (marginal), T=100: mean 1cM (below threshold)
- Add simulation scenario showing degradation at T > 50
- Clearly state limitation in paper and documentation

### CRIT-3: Add ALDER/ROLLOFF as Timing Baseline
- Compare HapGraph T estimates vs ALDER (LD-decay dating) on same populations
- This answers: "Does joint inference improve timing over a dedicated timing tool?"

### CRIT-4: Add Model Misspecification Scenario
- Simulate 2 admixture pulses, fit single-pulse model -> test bias/robustness

### CRIT-5: CLI Design
hapgraph run --vcf data.vcf.gz --popfile pops.txt --out results/ [--max-edges 3] [--seed 42] [--config config.yaml]
hapgraph prep --vcf data.vcf.gz --popfile pops.txt --out preprocessed/
hapgraph run --fstats fstats_dir/ --ibd ibd_summary.tsv --popfile pops.txt --out results/

### CRIT-6: Output Formats
- CSV/TSV: posteriors (mean, SD, 2.5%, 50%, 97.5%) for all parameters
- Newick-like: topology representation
- PDF/SVG: publication-ready graph figure
- WAIC table: model comparison across top-5 topologies
- ArviZ NetCDF: full posterior for downstream analysis
- Convergence report: R-hat, ESS, divergences per parameter
