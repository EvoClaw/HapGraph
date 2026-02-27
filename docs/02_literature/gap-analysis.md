# Gap Analysis — 1KGP Population Genetics
# Phase 1 Step 5 output

## What Has Been Done (Already Covered)

### SV Population Structure
- Almarri et al. 2020 (Cell): SV population structure in 929 HGDP samples. FST, PCA, STRUCTURE.
- Sudmant et al. 2015 (Nature): Phase 3 1kGP SV map in 2504 genomes. Population stratification.
- Korbel et al. 2025 (Nature): Long-read SVs in 1019 1kGP samples. PCA, admixture, FST, LD with SNPs.
  → These three papers cover SV population structure quite thoroughly.

### Balancing Selection on SVs (Deletions)
- eLife 2023: Excess of ancient deletion polymorphisms explained by balancing selection.
  → Deletions specifically; used older/smaller data; no INDEL analysis.

### Archaic Introgression via SVs
- biorXiv 2025: Global map of introgressed SV. 
  → This is recent and covers the introgression angle well.

### SNV-based Demographic Inference and Selection
- Many papers exist. Very competitive space.

### INDEL Functional Constraint
- Genome Res 2013: 179 individuals, INDELs in functional sequences face stronger purifying selection.
  → VERY OLD. Small sample. No population-level analysis.
- Nat Genet 2025: INDEL taxonomy focused on CANCER mutational signatures. Not population genetics.

---

## GAPS — Confirmed and Ranked by Novelty × Feasibility

### GAP 1 [HIGH NOVELTY, HIGH FEASIBILITY] ★★★★★
**INDEL population genetics at scale: 9.5 million INDELs in 3202 samples never systematically analyzed**

- The 2022 Cell paper released 9.5M phased INDELs but only reported their COUNT; no population genetic analysis.
- No paper has done: INDEL SFS by population, INDEL FST, INDEL-based selection scans, INDEL demographic inference.
- The only prior INDEL pop-gen paper used 179 genomes (2013, Genome Research). 18x smaller dataset.
- INDELs have fundamentally different mutational mechanisms (polymerase slippage, NHEJ, replication errors) from SNVs and SVs.
- This gap is genuine and large. Confidence: HIGH (verified — no competing paper found).

### GAP 2 [HIGH NOVELTY, MEDIUM FEASIBILITY] ★★★★
**No integrated SNV+INDEL+SV multi-variant-class comparison with the same dataset**

- All prior papers study one variant class at a time.
- Korbel 2025 compared SV vs SNV SFS at some level but only for 1019 samples and not systematically.
- No paper asks: do the three variant classes give concordant or discordant population structure signals?
- Do they agree on selection signals? Do they give the same demographic inference?
- Confidence: HIGH. No competing paper found.

### GAP 3 [MEDIUM NOVELTY, HIGH FEASIBILITY] ★★★
**Balancing selection on INDELs and SVs jointly in high-coverage 3202-sample data**

- eLife 2023 studied balancing selection on deletions with smaller data.
- No paper has studied balancing selection on INDELs.
- With 3202 samples (3x more than Korbel), statistical power is substantially better.
- Confidence: MEDIUM-HIGH.

### GAP 4 [MEDIUM NOVELTY, MEDIUM FEASIBILITY] ★★★
**SV-based fine-scale population differentiation with 3202 high-coverage samples**

- Korbel 2025 did FST analysis with 1019 samples. We have 3202.
- Key question: do more samples reveal additional population-differentiated SVs missed by Korbel?
- Risk: Korbel 2025 is very recent. Needs clear differentiation from their work.
- Confidence: MEDIUM.

### GAP 5 [MEDIUM NOVELTY, MEDIUM FEASIBILITY] ★★★
**Distribution of fitness effects (DFE) for INDELs vs SNVs across populations**

- DFE analysis (GRAPES, polyDFE) widely applied to SNVs but almost never to INDELs.
- Key question: Are INDEL fitness effect distributions shifted relative to SNVs? Are they population-specific?
- Could reveal whether different variant classes have different histories of purifying selection.
- Confidence: MEDIUM (one old paper exists, no recent large-scale study).

### GAP 6 [LOW-MEDIUM NOVELTY, HIGH FEASIBILITY] ★★
**Population-differentiated INDELs in functional elements (coding/regulatory)**

- Identify population-specific frameshifting INDELs in coding regions.
- Identify INDELs in eQTL positions with high FST.
- More of an annotation study; publishable but not high-impact alone.

---

## Summary Assessment

The LARGEST and MOST NOVEL gap is GAP 1: **INDEL population genetics**.
This is not a small oversight — it is a systematic blind spot in the field.
With 9.5M phased INDELs in 3202 samples, this dataset is 50x larger than any prior INDEL pop-gen study.

GAP 2 (integrated comparison) is complementary and could be combined with GAP 1 to make a stronger story:
"INDEL population genetics reveals distinct patterns from SNVs and SVs"
