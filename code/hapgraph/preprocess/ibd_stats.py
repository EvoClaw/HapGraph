"""
IBD statistics computation.
Two backends:
1. msprime TreeSequence (simulations) - fast via ts.ibd_segments(between=...)
2. hap-ibd output file parser (real data)
"""

import numpy as np
from typing import Dict, List, Tuple, Optional

MORGAN_TO_CM = 100.0


class IBDStatistics:
    """Container for inter-population IBD sharing statistics."""

    def __init__(self):
        self.populations = []
        self.mean_len = {}
        self.total_rate = {}
        self.n_segs = {}

    def get(self, a, b, stat="mean_len"):
        for key in [(a,b),(b,a)]:
            if stat == "mean_len"   and key in self.mean_len:   return self.mean_len[key]
            if stat == "total_rate" and key in self.total_rate: return self.total_rate[key]
            if stat == "n_segs"     and key in self.n_segs:     return float(self.n_segs[key])
        return float("nan")

    def as_matrix(self, stat="mean_len"):
        n = len(self.populations)
        mat = [[float("nan")]*n for _ in range(n)]
        for i,a in enumerate(self.populations):
            for j,b in enumerate(self.populations):
                mat[i][j] = 0.0 if i==j else self.get(a,b,stat)
        import numpy as np
        return np.array(mat), self.populations


def ibd_stats_from_msprime_ts(ts, pop_map, min_len_cM=2.0, genome_len_morgan=35.0):
    """
    Fast inter-population IBD statistics from msprime TreeSequence.
    Uses ts.ibd_segments(between=[A_samples, B_samples]) for each pair.
    """
    import numpy as np
    bp_per_cM = ts.sequence_length / (genome_len_morgan * MORGAN_TO_CM)
    min_span_bp = min_len_cM * bp_per_cM

    pop_samples = {name: [] for name in pop_map.values()}
    for pop_id, pop_name in pop_map.items():
        pop_samples[pop_name] = ts.samples(population=pop_id).tolist()

    pop_names = list(pop_map.values())
    ibd = IBDStatistics()
    ibd.populations = pop_names

    for i, pa in enumerate(pop_names):
        for j, pb in enumerate(pop_names):
            if j <= i:
                continue
            sA = pop_samples[pa]
            sB = pop_samples[pb]
            if not sA or not sB:
                ibd.mean_len[(pa,pb)] = float("nan")
                ibd.total_rate[(pa,pb)] = float("nan")
                ibd.n_segs[(pa,pb)] = 0
                continue
            try:
                ibd_result = ts.ibd_segments(
                    between=[sA, sB],
                    min_span=min_span_bp,
                    store_segments=True,
                )
                lens_bp = [
                    seg.right - seg.left
                    for _, segs in ibd_result.items()
                    for seg in segs
                ]
            except Exception:
                lens_bp = []

            if lens_bp:
                lens_cM = np.array(lens_bp) / bp_per_cM
                n_pairs = len(sA) * len(sB)
                ibd.mean_len[(pa,pb)]   = float(lens_cM.mean())
                ibd.total_rate[(pa,pb)] = float(lens_cM.sum() / n_pairs)
                ibd.n_segs[(pa,pb)]     = len(lens_cM)
            else:
                ibd.mean_len[(pa,pb)]   = float("nan")
                ibd.total_rate[(pa,pb)] = float("nan")
                ibd.n_segs[(pa,pb)]     = 0
    return ibd


def ibd_stats_from_hapibd_file(hapibd_file, pop_of_sample, min_len_cM=2.0):
    """Parse hap-ibd output and aggregate inter-population IBD statistics."""
    import gzip, numpy as np
    opener = gzip.open if hapibd_file.endswith(".gz") else open
    pops = sorted(set(pop_of_sample.values()))
    ibd = IBDStatistics()
    ibd.populations = pops
    seg_lengths = {}
    with opener(hapibd_file, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            s1, _, s2, _, _, _, _, length_str = parts[:8]
            try:
                length = float(length_str)
            except ValueError:
                continue
            if length < min_len_cM:
                continue
            p1 = pop_of_sample.get(s1)
            p2 = pop_of_sample.get(s2)
            if p1 is None or p2 is None or p1 == p2:
                continue
            key = (min(p1,p2), max(p1,p2))
            seg_lengths.setdefault(key, []).append(length)
    hap_count = {}
    for sample, pop in pop_of_sample.items():
        hap_count[pop] = hap_count.get(pop, 0) + 2
    for (pa,pb), lengths in seg_lengths.items():
        arr = np.array(lengths)
        n_pairs = hap_count.get(pa,1) * hap_count.get(pb,1)
        ibd.mean_len[(pa,pb)]   = float(arr.mean())
        ibd.total_rate[(pa,pb)] = float(arr.sum() / n_pairs)
        ibd.n_segs[(pa,pb)]     = len(arr)
    for pa in pops:
        for pb in pops:
            if pa >= pb:
                continue
            key = (pa,pb)
            if key not in ibd.mean_len:
                ibd.mean_len[key]   = float("nan")
                ibd.total_rate[key] = float("nan")
                ibd.n_segs[key]     = 0
    return ibd

