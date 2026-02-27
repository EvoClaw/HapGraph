"""
1kGP real data preprocessing for HapGraph.

Design decisions:
  - INDEL removal: bcftools pre-filter to biallelic SNPs only.
    F-statistics are only defined on SNPs. INDELs would bias F2.
  - No LD pruning: F-statistics use block jackknife (1-2 Mb blocks) which
    accounts for LD. LD pruning would reduce information without benefit.
  - No SNP downsampling: use all SNPs after quality filtering.
    More SNPs = lower SE. Downsampling to every 10th SNP reduces power.
  - MAF filter >= 0.01: removes rare variants with high sampling noise.
    Standard in ADMIXTOOLS2 (Patterson 2022).

Pipeline:
  1. bcftools: per-chromosome pre-filter → biallelic SNPs, MAF >= 0.01
  2. scikit-allel: load pre-filtered VCF, compute allele frequencies
  3. FStatistics: block jackknife over all chromosomes concatenated
  4. hap-ibd: IBD segment detection for T estimation (optional)

Usage:
  cd /home/yanlin/1KGothers/code
  python3 realdata/preprocess_1kgp.py --chroms 1 22 --n_per_pop 50
  python3 realdata/preprocess_1kgp.py --all_chroms  # full analysis
"""

import sys
import os
import argparse
import subprocess
import pickle
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional

sys.path.insert(0, str(Path(__file__).parent.parent))


# ── Population metadata ─────────────────────────────────────────────────────

SUPERPOP_MAP = {
    "YRI": "AFR", "LWK": "AFR", "GWD": "AFR", "MSL": "AFR",
    "ESN": "AFR", "ACB": "AFR", "ASW": "AFR",
    "CEU": "EUR", "TSI": "EUR", "GBR": "EUR", "FIN": "EUR", "IBS": "EUR",
    "CHB": "EAS", "CHS": "EAS", "JPT": "EAS", "CDX": "EAS", "KHV": "EAS",
    "BEB": "SAS", "GIH": "SAS", "ITU": "SAS", "PJL": "SAS", "STU": "SAS",
    "MXL": "AMR", "PUR": "AMR", "CLM": "AMR", "PEL": "AMR",
}

POPS_26 = sorted(SUPERPOP_MAP.keys())

ADMIXED_POPS = {"ACB", "ASW", "MXL", "PUR", "CLM", "PEL"}


def load_sample_map(
    samples_info_file: str,
    pops: Optional[List[str]] = None,
) -> Dict[str, List[str]]:
    """
    Load population-to-sample mapping from samples.info (pre-filtered unrelated).
    Format: sample_id<TAB>population<TAB>gender
    """
    pop_map = {}
    with open(samples_info_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            sample_id, pop = parts[0], parts[1]
            if pops and pop not in pops:
                continue
            pop_map.setdefault(pop, []).append(sample_id)
    return pop_map


def subsample_populations(
    pop_map: Dict[str, List[str]],
    n_per_pop: int = 50,
    seed: int = 42,
) -> Dict[str, List[str]]:
    """Subsample to at most n_per_pop samples per population."""
    rng = np.random.default_rng(seed)
    return {
        pop: list(rng.choice(samples, min(n_per_pop, len(samples)), replace=False))
        for pop, samples in pop_map.items()
    }


def bcftools_filter_chrom(
    vcf_gz: str,
    out_vcf_gz: str,
    maf: float = 0.01,
    n_threads: int = 4,
) -> bool:
    """
    Use bcftools to extract biallelic SNPs with MAF >= maf.

    Removes INDELs and multiallelic sites. No LD pruning (not needed for F-stats:
    block jackknife accounts for LD). Output as bgzipped VCF for allel compatibility.

    Returns True if successful.
    """
    if Path(out_vcf_gz).exists():
        print(f"(cached) ", end="")
        return True

    cmd = [
        "bcftools", "view",
        "--type", "snps",     # SNPs only, removes INDELs/MNPs/SVs
        "-m2", "-M2",         # biallelic only
        "-q", f"{maf}:minor", # MAF filter (both alleles >= maf)
        "--threads", str(n_threads),
        vcf_gz,
        "-O", "z",            # bgzipped VCF (scikit-allel compatible)
        "-o", out_vcf_gz,
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        subprocess.run(["bcftools", "index", "--tbi",
                        "--threads", str(n_threads), out_vcf_gz], check=True)
        return True
    except subprocess.CalledProcessError as e:
        print(f"  ERROR bcftools: {e.stderr}")
        return False


def extract_allele_freqs_from_bcf(
    bcf_path: str,
    pop_map: Dict[str, List[str]],
    all_samples: Optional[List[str]] = None,
) -> Tuple[np.ndarray, List[str], List[int]]:
    """
    Extract allele frequencies from a pre-filtered BCF file.

    Parameters
    ----------
    bcf_path : str path to pre-filtered biallelic SNP BCF
    pop_map : dict pop -> list of sample IDs
    all_samples : list of all sample IDs in BCF (from VCF header); if None, loaded automatically

    Returns
    -------
    freqs : (n_snps, n_pops) float32 array
    pop_names : list of population names (sorted)
    n_haps : list of haplotype counts per population
    """
    import allel

    callset = allel.read_vcf(
        bcf_path,
        fields=["calldata/GT", "samples"],
        alt_number=1,
    )
    if callset is None:
        return None, None, None

    all_sample_list = list(callset["samples"])
    gt = allel.GenotypeArray(callset["calldata/GT"])
    n_snps = gt.shape[0]

    pop_names = sorted(pop_map.keys())
    n_pops = len(pop_names)
    freqs = np.full((n_snps, n_pops), np.nan, dtype=np.float32)
    n_haps = []

    for k, pop in enumerate(pop_names):
        idxs = [all_sample_list.index(s) for s in pop_map[pop] if s in all_sample_list]
        if not idxs:
            n_haps.append(0)
            continue
        gt_pop = gt.take(idxs, axis=1)
        ac_pop = gt_pop.count_alleles()
        total = ac_pop.sum(axis=1)
        freq = np.where(total > 0, ac_pop[:, 1] / total.clip(min=1), np.nan)
        freqs[:, k] = freq.astype(np.float32)
        n_haps.append(len(idxs) * 2)

    return freqs, pop_names, n_haps


def main():
    parser = argparse.ArgumentParser(
        description="1kGP preprocessing for HapGraph (bcftools + scikit-allel)"
    )
    parser.add_argument(
        "--samples_info",
        default="/home/yanlin/1KGothers/1000GP/samples.info",
    )
    parser.add_argument(
        "--vcf_dir",
        default="/home/yanlin/1KGothers/1000GP/20220422_3202_phased_SNV_INDEL_SV",
    )
    parser.add_argument(
        "--chroms", nargs="+", default=["22"],
        help="Chromosomes (e.g. 1 2 22). Use --all_chroms for 1-22.",
    )
    parser.add_argument("--all_chroms", action="store_true",
                        help="Use all autosomes (chr1-22)")
    parser.add_argument(
        "--pops", nargs="+", default=None,
        help="Populations to include (default: all 26)",
    )
    parser.add_argument("--n_per_pop", type=int, default=50)
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--n_threads", type=int, default=4)
    parser.add_argument(
        "--bcf_cache",
        default="/tmp/hapgraph_bcf_cache",
        help="Directory for pre-filtered BCF cache",
    )
    parser.add_argument(
        "--outdir",
        default="/home/yanlin/1KGothers/results/realdata",
    )
    args = parser.parse_args()

    chroms = [str(c) for c in range(1, 23)] if args.all_chroms else args.chroms
    pops = args.pops or POPS_26
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    bcf_cache = Path(args.bcf_cache)
    bcf_cache.mkdir(parents=True, exist_ok=True)

    # ── Step 1: Load sample map ─────────────────────────────────────────────
    print("Loading sample map...")
    pop_map = load_sample_map(args.samples_info, pops=pops)
    pop_map = subsample_populations(pop_map, n_per_pop=args.n_per_pop)
    pop_names = sorted(pop_map.keys())
    total_samples = sum(len(s) for s in pop_map.values())
    print(f"  {len(pop_names)} populations, {total_samples} total samples")
    print(f"  Populations: {pop_names}")

    # ── Step 2: bcftools pre-filter each chromosome ─────────────────────────
    print(f"\nPre-filtering chromosomes (biallelic SNPs, MAF >= {args.maf})...")
    bcf_files = []
    for chrom in chroms:
        vcf_gz = (
            f"{args.vcf_dir}/1kGP_high_coverage_Illumina.chr{chrom}"
            f".filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        )
        if not Path(vcf_gz).exists():
            print(f"  SKIP: VCF not found for chr{chrom}")
            continue
        bcf_out = bcf_cache / f"chr{chrom}_snps_maf{int(args.maf*100)}.vcf.gz"
        print(f"  chr{chrom}...", end=" ", flush=True)
        ok = bcftools_filter_chrom(vcf_gz, str(bcf_out),
                                   maf=args.maf, n_threads=args.n_threads)
        if ok:
            bcf_files.append((chrom, str(bcf_out)))
            print("ok")

    # ── Step 3: Load allele frequencies ─────────────────────────────────────
    print("\nLoading allele frequencies...")
    all_freqs = []
    n_haps_final = None

    for chrom, bcf_path in bcf_files:
        print(f"  chr{chrom}...", end=" ", flush=True)
        freqs, _pop_names, n_haps = extract_allele_freqs_from_bcf(bcf_path, pop_map)
        if freqs is None:
            print("ERROR (empty)")
            continue
        all_freqs.append(freqs)
        if n_haps_final is None:
            n_haps_final = n_haps
        n_snps = freqs.shape[0]
        print(f"{n_snps:,} SNPs")

    if not all_freqs:
        print("ERROR: No data loaded!")
        return

    freqs_all = np.concatenate(all_freqs, axis=0)
    n_snps_total, n_pops_total = freqs_all.shape
    print(f"\nTotal: {n_snps_total:,} SNPs across {n_pops_total} populations")

    # ── Step 4: Compute F-statistics ────────────────────────────────────────
    print("\nComputing F-statistics (block jackknife, n_blocks=50*n_chroms)...")
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from hapgraph.preprocess.f_stats import FStatistics

    freqs_dict = {pop: freqs_all[:, k].astype(np.float64)
                  for k, pop in enumerate(pop_names)}
    n_haps_dict = {pop: n_haps_final[k] for k, pop in enumerate(pop_names)}

    n_blocks = max(50, 50 * len(bcf_files))  # 50 blocks per chromosome
    fstats = FStatistics(freqs=freqs_dict, n_haps=n_haps_dict, n_blocks=n_blocks)

    # ── Step 5: Save ─────────────────────────────────────────────────────────
    label = f"chr{'_'.join(chroms)}_n{args.n_per_pop}"
    fstats_file = outdir / f"fstats_{label}.pkl"
    with open(fstats_file, "wb") as f:
        pickle.dump(fstats, f)
    print(f"Saved: {fstats_file}")

    # Save metadata
    meta = {
        "chroms": chroms,
        "n_per_pop": args.n_per_pop,
        "populations": pop_names,
        "n_snps": int(n_snps_total),
        "maf_filter": args.maf,
        "n_haps": n_haps_dict,
        "note": "Biallelic SNPs only (INDELs removed by bcftools). No LD pruning."
                " Block jackknife accounts for LD in SE estimates.",
    }
    with open(outdir / f"metadata_{label}.json", "w") as f:
        json.dump(meta, f, indent=2)

    # ── Step 6: Print F2 summary ─────────────────────────────────────────────
    print("\nF2 matrix (representative pairs):")
    test_pairs = [
        ("YRI", "CEU"), ("YRI", "CHB"), ("CEU", "FIN"),
        ("CHB", "JPT"), ("ACB", "CEU"), ("ACB", "YRI"),
        ("MXL", "CEU"), ("MXL", "CHB"),
    ]
    for pa, pb in test_pairs:
        if pa in pop_names and pb in pop_names:
            obs, se = fstats.f2(pa, pb)
            z = obs / se if se > 0 else 0
            print(f"  F2({pa:3s}, {pb:3s}) = {obs:.5f} ± {se:.5f}  (Z={z:.1f})")

    # ── Step 7: F3 admixture signals ─────────────────────────────────────────
    print("\nF3 admixture signals (testing admixed populations):")
    test_trios = [
        ("ACB", "YRI", "CEU"),   # African Caribbean: AFR + EUR
        ("ASW", "YRI", "CEU"),   # African American: AFR + EUR
        ("MXL", "CHB", "CEU"),   # Mexican: EAS + EUR
        ("PUR", "YRI", "CEU"),   # Puerto Rican: AFR + EUR
        ("CLM", "YRI", "CEU"),   # Colombian: AFR + EUR
    ]
    for pop_c, pa, pb in test_trios:
        if all(p in pop_names for p in [pop_c, pa, pb]):
            obs, se = fstats.f3(pop_c, pa, pb)
            z = obs / se if se > 0 else 0
            sig = " *** ADMIXTURE SIGNAL" if z < -2 else (" (weak)" if z < 0 else "")
            print(f"  F3({pop_c}; {pa}, {pb}) = {obs:.5f} ± {se:.5f}  Z={z:.1f}{sig}")

    return fstats


if __name__ == "__main__":
    main()
