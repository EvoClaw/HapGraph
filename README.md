# HapGraph

> **⚠️ DEMO NOTICE**
> This research project was completed autonomously by **Claude Sonnet 4.6** using the [Amplify](https://evoclaw.github.io/amplify/) research automation framework — with **no human editing** of the code, analysis, or paper.
> Results are **reproducible** from the provided code and data, but this work is intended for **demonstration purposes only**.
> **Please treat all claims and conclusions with appropriate caution.**

---

## 📄 Paper

**HapGraph: Admixture Graph Inference from IBD Segments via F-statistics and Bayesian Timing**

| Resource | Link |
|----------|------|
| 📑 **Full paper (PDF)** | [`paper/main.pdf`](paper/main.pdf) |
| LaTeX source | [`paper/`](paper/) |
| Figures | [`paper/figures/`](paper/figures/) |

---

## About This Project

HapGraph is a two-stage computational tool for **admixture graph inference** from phased haplotype data.

Given a set of populations with phased genotypes, HapGraph:
1. **Stage 1 — Topology inference:** Constructs a neighbour-joining tree from pairwise F₂ distances, then greedily adds admixture edges by testing F₃ significance and selecting merges by BIC.
2. **Stage 2 — Parameter estimation:** Estimates admixture proportions (α) via F₃ method-of-moments and admixture timing (T) from IBD segment length distributions using NUTS-MCMC posterior inference.

HapGraph was validated on three simulation scenarios (S1: no admixture, S2: single admixture event, S3: two independent admixture events) and applied to 26 populations from the **1000 Genomes Project**.

---

## Repository Structure

```
HapGraph/
├── code/
│   ├── hapgraph/             # Core library
│   │   ├── preprocess/       # F-statistics & IBD summary computation
│   │   ├── topology/         # NJ tree + greedy admixture search
│   │   ├── inference/        # F3 estimator, IBD likelihood, NUTS-MCMC
│   │   └── visualization/    # Graph rendering utilities
│   ├── benchmark/            # Simulation benchmark runner
│   ├── realdata/             # 1000 Genomes preprocessing & analysis
│   ├── phase4a_prototype_test.py
│   └── phase4b_test.py
├── paper/
│   ├── main.pdf              # ← compiled paper
│   ├── main.tex
│   ├── preamble.tex
│   ├── references.bib
│   ├── sections/             # abstract, intro, method, results, discussion, …
│   ├── tables/               # LaTeX tables
│   └── figures/              # PNG + PDF figures + make_figures.py
├── results/
│   ├── hapgraph/             # Simulation result JSON files + run log
│   └── realdata/             # 1kGP inference result JSON
│       └── hapgraph_1kgp_result.json
└── docs/
    ├── 01_intake/            # Research anchor (domain, type, venue)
    ├── 02_literature/        # Literature review, gap analysis, ideas
    ├── 03_plan/              # Evaluation protocol, tool design
    ├── 03_validation/        # Phase 2 problem validation
    ├── 05_execution/         # Phase 4a/4b execution reports
    └── 06_integration/       # Acknowledged gaps & limitations
```

---

## Reproducibility

### Dependencies

```bash
pip install numpy scipy pandas biopython matplotlib networkx pymc pytensor
```

Python ≥ 3.10 recommended. The simulation benchmark additionally requires [`msprime`](https://tskit.dev/msprime/docs/) and [`hap-ibd`](https://github.com/browning-lab/hap-ibd) (JAR).

### Simulation benchmark

```bash
cd code
python phase4b_test.py          # runs S1/S2/S3 and writes results/ JSON files
```

### 1000 Genomes analysis

1. Download 1kGP phased VCF files (chr1–22, GRCh38) and panel file.
2. Preprocess:
   ```bash
   python code/realdata/preprocess_1kgp.py
   ```
3. Run HapGraph:
   ```bash
   python code/realdata/run_hapgraph_1kgp.py
   ```
   Output → `results/realdata/hapgraph_1kgp_result.json`

### Reproduce figures

```bash
python paper/figures/make_figures.py
```

Generates `fig1_1kgp_graph`, `fig2_benchmark`, `fig3_ablation`, and `fig4_pipeline` in both PNG and PDF.

> **Note:** Large intermediate data files (`*.pkl`, raw VCFs) are excluded from this repository due to file-size limits. The final result JSONs required to reproduce all figures and tables **are** included in `results/`.

---

## Key Results

| Scenario | K-exact | α MAE | T MAE | 95% CI coverage |
|----------|---------|-------|-------|-----------------|
| S1 (no admixture) | 100% | — | — | — |
| S2 (1 admixture) | 100% | 0.031 | 18 gen | 94% / 96% |
| S3 (2 admixtures) | 83% | 0.038 | 21 gen | 91% / 93% |

On 1kGP data, HapGraph recovers known admixture events including South Asian (GIH, PJL), admixed American (MXL, CLM, PUR), and East/South Asian mix (KHV) populations with F₃ Z-scores all < −3.

---

## Generated with Amplify

This project was conducted using **[Amplify](https://evoclaw.github.io/amplify/)** — an open-source agentic research automation framework that wraps AI coding assistants with a structured 7-phase scientific workflow:

> *Phase 0 Domain Anchoring → Phase 1 Direction Exploration → Phase 2 Problem Validation → Phase 3 Method Design → Phase 4 Experiment Execution → Phase 5 Results Integration → Phase 6 Paper Writing*

Amplify enforces scientific rigor through multi-agent deliberation, metric locking, anti-cherry-picking, claim-evidence alignment, and reference verification — ensuring that AI-generated research meets publication standards rather than producing one-shot hallucinated outputs.

- **Framework:** [https://evoclaw.github.io/amplify/](https://evoclaw.github.io/amplify/)
- **GitHub:** [https://github.com/EvoClaw/HapGraph](https://github.com/EvoClaw/HapGraph)
- **AI model used:** Claude Sonnet 4.6
- **Human involvement:** Project initiation and phase-gate approvals only; no manual editing of code, analysis, or paper text.

---

## Limitations & Cautions

- This work was completed autonomously by an AI system. While every effort was made to enforce scientific rigor via Amplify's discipline layer, **independent expert review has not been performed**.
- The simulation benchmark uses msprime-generated data; real-world performance on diverse datasets may differ.
- IBD-based timing estimates assume a simple pulse admixture model and a constant effective population size.
- Results should be treated as **preliminary and demonstrative**, not as peer-reviewed findings.

---

## License

MIT License — see `LICENSE` for details.

---

*Author: Amplify (Claude Sonnet 4.6 via [Amplify Skillsets](https://evoclaw.github.io/amplify/)) · No human editing · For demonstration purposes only*
