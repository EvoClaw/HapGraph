# Candidate Research Ideas — Phase 1 Step 5c
# Pre-brainstorming: 5 idea cards

═══════════════════════════════════════════════════════════════
Idea 1: INDEL Population Genetics at Genome Scale (SAFE)
═══════════════════════════════════════════════════════════════
Core question:   What are the population-genetic properties of 9.5 million phased INDELs 
                 across 26 human populations? Do they reveal the same demographic history 
                 and selection signatures as SNVs, or distinct patterns?

Novelty source:  Trend extrapolation (#5) + Assumption Challenging (#1) + Limitation-to-Opportunity (#3)
                 The dataset paper explicitly did not analyze INDELs. No prior paper has done 
                 genome-wide INDEL population genetics with >200 samples.

Why it matters:  INDELs represent 15% of all variants in 1kGP 2022 and are systematically ignored.
                 If INDEL patterns differ from SNVs, it changes how we interpret population history.
                 If they agree, it validates the robustness of our inference methods.
                 Frameshifting INDELs in coding regions could reveal population-specific functional effects.

Feasibility:     HIGH — Standard tools (PLINK, ANGSD, vcftools) work directly on INDEL VCFs.
                 CPU-only, standard population genetics pipeline. 3-6 month timeline.
                 Data already in hand. No new data collection.

Risk:            Results might look "too similar to SNVs" — publication challenge if no interesting 
                 differences are found. Mitigation: focus on INDEL size classes and repeat context 
                 which are guaranteed to differ.

Competition:     LOW — No competing paper found. This is the lowest-risk gap.

Estimated scope: Full journal article, Genome Research / MBE / PLOS Genetics

Key analyses:
  - INDEL allele frequency spectrum by population × INDEL class (frameshift/in-frame; size 1bp/2-5bp/6-50bp)
  - INDEL-based FST across 26 populations; compare to SNV FST
  - INDEL-based PCA and ADMIXTURE
  - Tajima's D, Fu & Li's D for INDELs by population
  - Selection scan: iHS and XP-EHH on phased INDEL haplotypes
  - Functional annotation: INDELs in coding, UTR, promoter, enhancer regions
  - LD between INDELs and nearby SNVs/SVs

═══════════════════════════════════════════════════════════════
Idea 2: Integrated Three-Variant-Class Population Portrait (AMBITIOUS)
═══════════════════════════════════════════════════════════════
Core question:   Do SNVs, INDELs, and SVs give concordant or discordant signals of population 
                 structure, selection, and demographic history when analyzed with the same 
                 dataset? Where do they agree and where do they diverge?

Novelty source:  Contradiction Mining (#2) + Cross-Domain Transfer (#4) + Counterfactual (#6)
                 Plant genomics routinely uses all variant types jointly. Human genomics never has.

Why it matters:  If the three classes agree: robustness of our demographic/evolutionary inferences.
                 If they disagree in specific regions: those regions are under variant-class-specific 
                 evolutionary forces (e.g., repeat-driven SV formation, slippage-driven INDEL expansion).
                 This provides a complete picture of human genetic diversity.

Feasibility:     MEDIUM — Requires careful normalization across variant classes. SVs have different 
                 properties (phased but lower quality in short reads). Analysis pipeline more complex.
                 CPU feasible. 4-8 months.

Risk:            Results might be "each class tells the same story" — hard to publish if no contrast.
                 Risk of scope creep (doing three separate analyses poorly rather than one well).

Competition:     LOW for the integrated approach. MEDIUM for individual components.

Estimated scope: Ambitious journal article, potentially Genome Biology or Genome Research

Key analyses:
  - Parallel FST analysis: SNV-FST vs INDEL-FST vs SV-FST; correlation and outliers
  - Parallel PCA: do 3 classes produce same clusters?
  - Concordance of demographic inference: do SFS-based Ne estimates agree across classes?
  - Identify "class-discordant" regions (high SV-FST but low SNV-FST) as novel discovery

═══════════════════════════════════════════════════════════════
Idea 3: Distribution of Fitness Effects for INDELs vs SNVs (METHOD-ADJACENT)
═══════════════════════════════════════════════════════════════
Core question:   What is the distribution of fitness effects (DFE) for INDELs in each of the 
                 26 populations, and how does it differ from the DFE of SNVs?

Novelty source:  Assumption Challenging (#1) + Trend Extrapolation (#5)
                 DFE analysis (GRAPES, polyDFE) is standard for SNVs. Never applied to INDELs at scale.

Why it matters:  Understanding the DFE reveals the proportion of new mutations that are strongly/mildly 
                 deleterious vs neutral. If INDEL DFE differs from SNV DFE, it suggests INDELs 
                 face fundamentally different selective constraints — with implications for GWAS 
                 variant interpretation and disease genetics.

Feasibility:     MEDIUM-HIGH — GRAPES and polyDFE are CPU tools, work on folded SFS.
                 Need to separate synonymous (neutral) from nonsynonymous INDELs for calibration.
                 Challenge: defining "synonymous" INDELs is non-trivial (3bp in-frame vs frameshift).

Risk:            DFE methods have strong demographic assumptions; need careful validation.
                 May require simulation-based benchmarking.

Competition:     LOW — No INDEL DFE paper at this scale found.

Estimated scope: Journal article, MBE / Genome Biology

═══════════════════════════════════════════════════════════════
Idea 4: INDEL Functional Divergence and Population-Stratified Disease Risk (TRANSLATIONAL)
═══════════════════════════════════════════════════════════════
Core question:   Which INDELs show extreme population differentiation (high FST), and do they 
                 disproportionately affect functional elements? Are population-stratified 
                 frameshifting INDELs enriched in disease-associated genes?

Novelty source:  Cross-Domain Transfer (#7 — clinical genetics) + Limitation-to-Opportunity (#3)

Why it matters:  Population-stratified functional INDELs could explain population differences in 
                 disease prevalence that SNV-GWAS misses. This bridges pop gen with medical genetics.
                 Directly actionable for precision medicine / pharmacogenomics.

Feasibility:     HIGH — FST calculation, functional annotation (VEP), disease database overlap 
                 (OMIM, ClinVar). All CPU tools. 2-4 months for core analysis.

Risk:            Many groups work on population-stratified disease risk. Risk of "too incremental."
                 Needs strong functional validation story.

Competition:     MEDIUM — Many papers on population-stratified SNVs and SVs. INDEL angle is novel.

Estimated scope: Journal article, Human Molecular Genetics / BMC Genomics / PLOS Genetics

═══════════════════════════════════════════════════════════════
Idea 5: INDEL Haplotype Blocks and Population History (DEMOGRAPHIC)
═══════════════════════════════════════════════════════════════
Core question:   Do phased INDEL haplotypes define population boundaries more or less precisely 
                 than SNV haplotypes? Can INDEL-tagged haplotype blocks serve as markers of 
                 archaic introgression or population admixture events?

Novelty source:  Trend Extrapolation (#5) + Counterfactual (#6)

Why it matters:  Phased data allows reconstruction of INDEL haplotype structure. If INDELs 
                 mark ancient haplotype blocks, they could reveal admixture/introgression events 
                 that SNVs (which recombine freely) miss. Extended INDEL haplotypes in specific 
                 populations could mark selective sweeps or introgression.

Feasibility:     MEDIUM — LD analysis with phased data, iHS/nSL for INDELs. 
                 Haplotype phasing of INDELs is already done in 2022 dataset.
                 Some computational challenges with INDEL LD calculation.

Risk:            INDELs in repeat regions might have inflated LD due to technical artifacts.
                 Need careful filtering of repeat-context INDELs for LD analysis.

Competition:     LOW — No prior paper on INDEL haplotype structure in human populations.

Estimated scope: Journal article, Genome Research / MBE


═══════════════════════════════════════════════════════════════
MULTI-AGENT BRAINSTORMING — Pending
═══════════════════════════════════════════════════════════════
Status: Ready for Round 1 deliberation

═══════════════════════════════════════════════════════════════
ROUND 1 AGENT VERDICTS SYNTHESIS
═══════════════════════════════════════════════════════════════

| Idea | Visionary | Pragmatic | Scout | Consensus |
|------|-----------|-----------|-------|-----------|
| 1. INDEL Pop Gen | STRONG | VIABLE | STRONG | ✅ PROCEED |
| 2. 3-class Integrated | STRONG | WEAK* | STRONG | ⚠️ MODIFY |
| 3. DFE INDELs | VIABLE | WEAK | VIABLE | → Sub-analysis only |
| 4. High-FST INDEL+Disease | VIABLE | STRONG | VIABLE | → Include as component |
| 5. INDEL Haplotypes+Archaic | WEAK | KILL | WEAK | ❌ DROP |

*Pragmatic's concern: SVs are technically complex. Safer version: drop SVs, do SNV vs INDEL only.

CONVERGENCE: 
Ideas 1 + Modified 2 (SNV vs INDEL, no SV) = Visionary's New Idea A
All three agents rate this STRONG/VIABLE. Idea 4 becomes a component subsection.

═══════════════════════════════════════════════════════════════
FINAL RANKED TOP IDEAS (Post Round 1 Deliberation)
═══════════════════════════════════════════════════════════════

RANK 1 [TOP RECOMMENDATION] ★★★★★
════════════════════════════════
Title: INDEL Population Genetics in 26 Human Populations: Patterns of 
       Diversity, Selection, and Comparison with SNVs

Core question: What are the genome-wide population-genetic properties of 9.5 million 
               phased INDELs across 26 human populations? Do INDELs and SNVs show 
               concordant or discordant signals of population structure, selection, 
               and demographic history?

Components:
  A. Allele frequency spectra: INDEL SFS by population × INDEL class (size 1bp/2-5bp/
     6-49bp; frameshift/in-frame/non-coding)
  B. Population structure: INDEL-based PCA, ADMIXTURE, FST heatmap — compare to SNV
  C. Selection signatures: Tajima's D, Fu&Li's D, iHS, XP-EHH on INDEL haplotypes
  D. Concordance/discordance: Where do INDEL and SNV FST/selection signals disagree?
     (These "class-discordant" loci are the main discovery angle)
  E. Functional annotation: High-FST frameshift INDELs in coding/regulatory elements 
     + disease database overlap (translational angle from Idea 4)

Why STRONG:
  - 0 competing papers for INDEL pop gen at this scale
  - Dataset paper explicitly deferred this analysis (published invitation)
  - Methods are all CPU-feasible with standard tools
  - Clear story: "A comprehensive second look at 1kGP through the INDEL lens"
  - Discordance = guaranteed novel finding (either agreement confirms robustness, 
    or disagreement reveals new biology — both publishable)

Target venue: Genome Research (primary), MBE or Genome Biology (backup)
Timeline: 4-6 months for core analysis; 2 months writing = ~8 months total

RANK 2 [OPTIONAL EXTENSION]
════════════════════════════
Title: INDEL Distribution of Fitness Effects across Populations
(only if GRAPES/polyDFE support INDELs; otherwise drop)
Timeline: add 2 months if included

ELIMINATED:
  Idea 3 (DFE) — only as secondary analysis if tools support
  Idea 5 (archaic introgression) — crowded, feasibility concerns
