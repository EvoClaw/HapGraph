# Deep Thinking Insights — Phase 1 Step 5b
# 6 Strategies Applied

## Insight 1
Strategy:    Assumption Challenging
Observation: ALL population genetics papers in 1kGP (phase 3 onwards) use SNVs as the primary currency 
             of analysis. INDEL variants are always filtered out or treated as "noisy SNVs." 
             The implicit assumption is: SNVs are sufficient to capture demographic and selective history.
Implication: INDELs have DIFFERENT mutational mechanisms (polymerase slippage, replication fork stalling, 
             NHEJ). Their rate variation across the genome differs from SNVs. Their functional consequences 
             (frameshifts) are qualitatively different. They may therefore carry demographic signals that SNVs 
             MISS — particularly in repeat-rich regions. If the assumption is wrong, the entire field has been 
             systematically ignoring a major source of population-genetic signal.
Evidence:    (1) Genome Res 2013: INDELs show stronger purifying selection than SNPs at functional sites
             (2) INDELs in 1kGP 2022 are 9.5M — 15% of all variants — yet no population analysis exists
             (3) INDEL mutation rates are ~10-fold higher than SNV rates in microsatellite contexts
Strength:    strong
Verified:    yes — no competing paper on INDEL population genetics found in systematic search

## Insight 2
Strategy:    Contradiction Mining
Observation: Paper A (Korbel 2025): "SVs are enriched for eQTLs compared to SNPs" — suggesting SVs 
             are MORE functionally important per-variant.
             Paper B (gnomAD constraint paper 2023): "Most of the genome's functional constraint is 
             captured by SNVs" — suggesting SNVs dominate.
             These two claims seem contradictory. 
Implication: The resolution may lie in INDEL class: INDELs at microsatellite loci might be the 
             "hidden layer" that both papers miss. INDELs in regulatory regions may mediate gene 
             expression more than raw SNV counts would suggest. An integrated analysis could 
             resolve this contradiction by asking: what fraction of population differentiation in 
             regulatory regions is driven by INDELs vs SNVs vs SVs?
Evidence:    Both papers cited above. Contradiction is real.
Strength:    moderate
Verified:    partial — needs deeper analysis

## Insight 3
Strategy:    Limitation-to-Opportunity Conversion
Observation: Byrska-Bishop 2022 (Cell) explicitly states in their paper: "The analysis of the 9.5M INDELs 
             is beyond the scope of this paper." They deliberately left INDEL population analysis as 
             future work.
Implication: The authors of the dataset paper explicitly invited follow-up INDEL analysis. 
             This is a rare case where the primary dataset paper LEFT A DOOR OPEN. 
             This creates a clear publication opportunity: "We perform the analysis that Byrska-Bishop 
             left for future work."
Evidence:    Byrska-Bishop 2022 Cell paper (our data source)
Strength:    strong
Verified:    yes — need to confirm in paper text, but implied by absence of INDEL analysis

## Insight 4
Strategy:    Cross-Domain Transfer
Observation: In microbial population genetics (E. coli, yeast), indel polymorphism is used as a 
             COMPLEMENTARY marker to SNPs for population structure. Multiple papers in microbial 
             genomics use indel SFS to infer demographic bottlenecks that SNPs miss due to 
             saturation at high mutation rates.
             In PLANTS (rice, maize), INDEL-based population analysis is routine (e.g., 3000 Rice 
             Genomes paper uses INDELs for population structure).
Implication: The methods exist in plant/microbial genomics but have not been systematically 
             applied to HUMAN population genomics at scale. Direct transfer is feasible.
             Rice and maize papers using INDELs routinely appear in Nature Genetics, Plant Cell.
Evidence:    Standard practice in plant genomics literature.
Strength:    strong
Verified:    yes — human INDEL population analysis papers are absent from search results

## Insight 5
Strategy:    Trend Extrapolation
Observation: The field trajectory: 
             2011: CNV population genetics (Mills et al. Nature)
             2015: Full SV population genetics (Sudmant et al. Nature; Almarri et al.)
             2020: gnomAD SV reference (Collins et al. Nature)
             2022: High-coverage SNV+INDEL+SV (Byrska-Bishop Cell) → SNV analyzed, INDEL ignored
             2025: Long-read SV (Korbel Nature) → SV thoroughly re-analyzed
             
             The NATURAL NEXT STEP: systematic INDEL population genetics.
             The field has analyzed CNVs → SVs → SNVs thoroughly. INDELs are the missing piece.
             Nobody has taken this step. Why? Possibly because INDELs were historically seen as 
             "hard to genotype" — but high-coverage data (2022) now makes them reliable.
Implication: This is the natural "D" step after A→B→C. First-mover advantage exists.
             INDEL genotyping quality in high-coverage data is now on par with SNVs.
Evidence:    The progression above is documented in the literature.
Strength:    strong
Verified:    yes — high-coverage data makes INDEL genotypes reliable (Byrska-Bishop 2022)

## Insight 6
Strategy:    Counterfactual Reasoning
Observation: "What if we applied SNV population genetics methods to INDELs?"
             Standard methods: PLINK (FST, PCA), ADMIXTURE, ANGSD (SFS, theta, neutrality tests), 
             iHS/XP-EHH selection tests, demographic inference (moments/momi) — ALL work on 
             biallelic variants. INDELs ARE biallelic (mostly). So these methods DIRECTLY apply.
             
             "What if INDEL patterns disagree with SNV patterns?"
             → This would be major news: implies INDELs carry unique population-genetic signal
             → Even if they AGREE: "INDELs confirm SNV-based inferences" is publishable as validation
             → If PARTIALLY disagree: "INDELs reveal additional structure" — best case scenario
Implication: Zero methodological barrier. The analysis is straightforward. The novelty is in 
             DOING IT at scale, not in inventing new methods.
Evidence:    Standard bioinformatics tools (PLINK, ANGSD, etc.) are all INDEL-compatible
Strength:    strong
Verified:    yes — all standard tools work on any biallelic variants

## Insight 7
Strategy:    Cross-Domain Transfer (second application)
Observation: In clinical genetics, there is a well-known distinction between "coding INDELs" 
             (frameshift mutations causing disease) and "non-coding INDELs." 
             Population genetics has largely ignored non-coding INDELs. Yet non-coding INDELs 
             in promoters and enhancers can create/destroy transcription factor binding sites.
             A new 2025 taxonomy paper redefined INDEL classes based on sequence context.
Implication: Applying this new INDEL taxonomy to POPULATION-LEVEL data (not just cancer) would 
             reveal whether different INDEL classes have different evolutionary histories.
             "Do promoter INDELs show stronger purifying selection than intronic INDELs?"
             "Are frameshift coding INDELs more population-differentiated than synonymous ones?"
Evidence:    Nat Genet 2025 INDEL taxonomy paper (cancer context); our dataset (population context)
Strength:    moderate
Verified:    partial

## Summary of Insights
All six strategies converge on the same core opportunity: INDEL population genetics.
Most critical insights: #1 (assumption challenging), #3 (limitation-to-opportunity), #5 (trend extrapolation)
These three independently point to the same gap, making this a STRONG research direction.
