"""
Generate all three paper figures for HapGraph.
Run from: /home/yanlin/1KGothers/paper/figures/
"""
import json, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from pathlib import Path

RESULTS = Path("/home/yanlin/1KGothers/results")
OUTDIR  = Path("/home/yanlin/1KGothers/paper/figures")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.08,
})

# Colorblind-safe Wong palette
BLUE   = "#0072B2"
ORANGE = "#E69F00"
GREEN  = "#009E73"
RED    = "#D55E00"
PURPLE = "#CC79A7"
GRAY   = "#999999"
LBLUE  = "#56B4E9"


# ─────────────────────────────────────────────────────────────────────────────
# Figure 1: 1kGP admixture graph — proper cladogram + admixture arrows
# ─────────────────────────────────────────────────────────────────────────────
def fig1_1kgp_graph():
    from Bio import Phylo
    from io import StringIO

    with open(RESULTS / "realdata/hapgraph_1kgp_result.json") as f:
        d = json.load(f)

    events   = d["admixture_events"]
    nwk      = d["nj_newick"]
    tree     = Phylo.read(StringIO(nwk), "newick")

    # Super-population colors
    sp_color = {}
    for pops, col in [
        (["YRI","ESN","GWD","MSL","LWK","ACB","ASW"], BLUE),
        (["CEU","GBR","FIN","IBS","TSI"],              ORANGE),
        (["CHB","CHS","CDX","KHV","JPT"],              GREEN),
        (["GIH","PJL","BEB","STU","ITU"],              RED),
        (["MXL","PUR","CLM","PEL"],                    PURPLE),
    ]:
        for p in pops:
            sp_color[p] = col

    admixed = {ev["tgt"] for ev in events}

    # ── Phylogram layout: x = cumulative branch length from root ─────────────
    def assign_x(clade, x0=0.0):
        bl = clade.branch_length if clade.branch_length else 0.0
        clade._x = x0 + bl
        for child in clade.clades:
            assign_x(child, clade._x)

    leaf_counter = [0]
    def assign_y(clade):
        if clade.is_terminal():
            clade._y = leaf_counter[0]
            leaf_counter[0] += 1
        else:
            for child in clade.clades:
                assign_y(child)
            clade._y = (clade.clades[0]._y + clade.clades[-1]._y) / 2.0

    assign_x(tree.root, 0.0)
    assign_y(tree.root)

    def all_clades(clade):
        yield clade
        for c in clade.clades:
            yield from all_clades(c)

    leaf_pos = {}
    node_pos = {}
    for cl in all_clades(tree.root):
        node_pos[id(cl)] = (cl._x, cl._y)
        if cl.is_terminal():
            leaf_pos[cl.name] = (cl._x, cl._y)

    n_leaves  = leaf_counter[0]
    x_max_tip = max(x for x, y in leaf_pos.values())

    # ── Draw ─────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 8.5))
    fig.subplots_adjust(left=0.01, right=0.80, top=0.93, bottom=0.10)
    ax.axis("off")

    # Draw tree branches (phylogram style)
    def draw_branches(clade):
        px, py = node_pos[id(clade)]
        if clade.clades:
            ys = [node_pos[id(c)][1] for c in clade.clades]
            ax.plot([px, px], [min(ys), max(ys)], color="#cccccc", lw=1.0, zorder=1)
        for child in clade.clades:
            cx, cy = node_pos[id(child)]
            ax.plot([px, cx], [cy, cy], color="#cccccc", lw=1.0, zorder=1)
            draw_branches(child)

    draw_branches(tree.root)

    # Draw leaf nodes and labels
    # For short-branched leaves (EUR), extend a dotted line to a common x
    x_label = x_max_tip * 1.02   # all labels at this x column
    for name, (x, y) in leaf_pos.items():
        col = sp_color.get(name, GRAY)
        mk  = "D" if name in admixed else "o"
        sz  = 100 if name in admixed else 65
        ax.scatter(x, y, s=sz, c=col, marker=mk, zorder=5,
                   edgecolors="white", linewidths=0.9)
        # Dotted leader line from tip to label column (for short branches)
        if x < x_label - 0.0005:
            ax.plot([x, x_label], [y, y], color=col,
                    lw=0.5, ls=":", alpha=0.45, zorder=1)
        ax.text(x_label + 0.00035, y, name, fontsize=8,
                va="center", ha="left", color=col,
                fontweight="bold" if name in admixed else "normal",
                zorder=6)

    # Scale bar
    sb_len = 0.005
    sb_y   = -0.8
    ax.plot([0, sb_len], [sb_y, sb_y], color="black", lw=1.5)
    ax.text(sb_len / 2, sb_y - 0.35, f"{sb_len:.3f}", fontsize=7,
            ha="center", va="top", color="black")

    # ── Admixture arrows ──────────────────────────────────────────────────────
    # Events sorted: long-range first so they don't crowd near the leaves
    events_sorted = sorted(events, key=lambda e: -abs(
        leaf_pos.get(e["src_a_f3"], (0,0))[1] - leaf_pos.get(e["tgt"], (0,0))[1]
    ))

    arrow_colors = [RED, PURPLE, "#7b2d8b", "#b05000", ORANGE, BLUE, "#006060", GREEN]
    # Pre-compute label y-positions to stagger them
    label_y_used = []

    for ev_idx, ev in enumerate(events_sorted):
        tgt   = ev["tgt"]
        src   = ev["src_a_f3"]
        alpha = ev["alpha_mom"]
        z     = ev["f3_z"]
        if tgt not in leaf_pos or src not in leaf_pos:
            continue
        x0, y0 = leaf_pos[src]
        x1, y1 = leaf_pos[tgt]
        col = arrow_colors[ev_idx % len(arrow_colors)]

        dy = y1 - y0
        # Vary curvature so multiple arrows from same src fan out
        base_rad = 0.38 + 0.06 * (ev_idx % 3)
        rad = base_rad if dy < 0 else -base_rad

        ax.annotate(
            "", xy=(x1, y1), xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>",
                color=col, lw=1.4,
                connectionstyle=f"arc3,rad={rad}",
                linestyle=(0, (5, 3)),
                mutation_scale=11,
            ),
            zorder=4,
        )
        # Place label left of midpoint; stagger vertically to avoid overlap
        lx = (x0 + x1) / 2 - 0.006   # shift left by ~20% of total x range
        ly = (y0 + y1) / 2
        # Nudge if too close to an already-placed label
        for used_y in label_y_used:
            if abs(ly - used_y) < 0.9:
                ly = used_y - 1.1 if ly <= used_y else used_y + 1.1
        label_y_used.append(ly)

        ax.text(lx, ly, f"α={alpha:.2f}  Z={z:.0f}",
                fontsize=6.5, color=col, va="center", ha="right",
                bbox=dict(facecolor="white", edgecolor=col, alpha=0.90,
                          pad=1.4, boxstyle="round,pad=0.28", linewidth=0.7),
                zorder=7)

    # ── Super-population bracket on the right ─────────────────────────────────
    # Only draw contiguous brackets (skip if max-min span > group size + 1)
    sp_groups = [
        ("AFR", ["YRI","ESN","GWD","MSL","LWK","ACB","ASW"], BLUE),
        ("EUR", ["GBR","CEU","TSI","IBS"],                    ORANGE),  # FIN excluded (outlier)
        ("EAS", ["CHB","CHS","CDX","KHV","JPT"],              GREEN),
        ("SAS", ["GIH","STU","ITU","BEB","PJL"],              RED),
        ("AMR", ["PEL","MXL","CLM"],                          PURPLE),  # PUR excluded (outlier)
    ]
    x_bar = x_label + 0.0025
    for sp_label, grp, col in sp_groups:
        ys = [leaf_pos[p][1] for p in grp if p in leaf_pos]
        if not ys:
            continue
        span = max(ys) - min(ys)
        if span > len(grp) + 1:
            continue
        y_lo, y_hi = min(ys) - 0.3, max(ys) + 0.3
        ax.plot([x_bar]*2, [y_lo, y_hi], color=col, lw=3.5,
                solid_capstyle="round", zorder=3)
        ax.text(x_bar + 0.0004, (y_lo + y_hi) / 2, sp_label,
                va="center", ha="left", fontsize=9.5, color=col, fontweight="bold")
    for name, note, col in [("FIN","EUR*",ORANGE), ("PUR","AMR*",PURPLE)]:
        if name in leaf_pos:
            _, y = leaf_pos[name]
            ax.text(x_bar + 0.0004, y, note, va="center", ha="left",
                    fontsize=7.5, color=col, fontstyle="italic")

    ax.set_xlim(-0.003, x_label + 0.006)
    ax.set_ylim(-1.2, n_leaves + 0.5)

    # Legend
    legend_els = [
        mpatches.Patch(color=BLUE,   label="AFR"),
        mpatches.Patch(color=ORANGE, label="EUR"),
        mpatches.Patch(color=GREEN,  label="EAS"),
        mpatches.Patch(color=RED,    label="SAS"),
        mpatches.Patch(color=PURPLE, label="AMR"),
        plt.scatter([], [], s=70, c="#888888", marker="D",
                    edgecolors="white", label="Admixed (target)"),
        plt.scatter([], [], s=50, c="#888888", marker="o",
                    edgecolors="white", label="Non-admixed"),
        plt.Line2D([0],[0], color=RED, lw=1.3, ls="--",
                   label="Admixture edge (α, F3 Z-score)"),
    ]
    ax.legend(handles=legend_els, loc="lower left", ncol=2,
              fontsize=7.5, frameon=True, framealpha=0.9,
              edgecolor="#dddddd", bbox_to_anchor=(0.0, -0.01))

    ax.set_title("HapGraph admixture graph on 1000 Genomes Project  (K = 8 events)\n"
                 "Tree backbone: NJ on genome-wide F2 statistics",
                 fontsize=9.5, pad=8)

    fig.savefig(OUTDIR / "fig1_1kgp_graph.pdf")
    fig.savefig(OUTDIR / "fig1_1kgp_graph.png")
    plt.close(fig)
    print("Figure 1 saved.")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 2: Simulation benchmark scatter plots
# ─────────────────────────────────────────────────────────────────────────────
def fig2_benchmark_scatter():
    def load_oracle(path):
        with open(path) as f:
            return list(json.load(f).values())[0]

    s2 = load_oracle(RESULTS / "hapgraph/S2_oracle_results.json")
    s3 = load_oracle(RESULTS / "hapgraph/S3_oracle_results.json")

    def extract_alpha(rows, max_k=2):
        trues, ests, ci_lo, ci_hi = [], [], [], []
        for r in rows:
            for k in range(max_k):
                at = r.get(f"alpha_true_{k}")
                ae = r.get(f"alpha_hat_{k}")
                ci = r.get(f"alpha_ci_{k}")
                if at is not None and ae is not None:
                    trues.append(at)
                    ests.append(ae)
                    ci_lo.append(ci[0] if isinstance(ci, list) else ae)
                    ci_hi.append(ci[1] if isinstance(ci, list) else ae)
        return map(np.array, [trues, ests, ci_lo, ci_hi])

    def extract_T(rows, max_k=2):
        trues, ests, ci_lo, ci_hi = [], [], [], []
        for r in rows:
            for k in range(max_k):
                tt = r.get(f"T_true_{k}")
                te = r.get(f"T_hat_{k}")
                ci = r.get(f"T_ci_{k}")
                if tt is not None and te is not None:
                    trues.append(tt)
                    ests.append(te)
                    ci_lo.append(ci[0] if isinstance(ci, list) and ci[0] is not None else te)
                    ci_hi.append(ci[1] if isinstance(ci, list) and ci[1] is not None else te)
        return map(np.array, [trues, ests, ci_lo, ci_hi])

    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6.5))
    fig.subplots_adjust(hspace=0.44, wspace=0.40)

    scenarios = [("S2 (K=1)", s2, BLUE, 1), ("S3 (K=2)", s3, ORANGE, 2)]

    # Row 0: Alpha
    for col_idx, (label, rows, col, mk) in enumerate(scenarios):
        ax = axes[0, col_idx]
        ta, ea, clo_a, chi_a = extract_alpha(rows, max_k=mk)
        if len(ta):
            elo = np.clip(ea - clo_a, 0, 1)
            ehi = np.clip(chi_a - ea, 0, 1)
            ax.errorbar(ta, ea, yerr=[elo, ehi],
                        fmt="o", color=col, ms=5, lw=0.8,
                        capsize=2.5, elinewidth=0.8, alpha=0.80)

        lim = [-0.02, 0.72]
        ax.plot(lim, lim, "k--", lw=0.9, alpha=0.45, zorder=0)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("True α", fontsize=8.5)
        ax.set_ylabel("Estimated α", fontsize=8.5)
        ax.set_title(label, fontsize=9.5, fontweight="bold", color=col, pad=6)

        a_maes = [r["alpha_mae"] for r in rows if r.get("alpha_mae") is not None]
        a_covs = [r["alpha_ci_cover"] for r in rows if r.get("alpha_ci_cover") is not None]
        txt = f"MAE = {np.mean(a_maes):.3f}"
        if a_covs:
            txt += f"\n95% CI cover = {np.mean(a_covs):.0%}"
        ax.text(0.04, 0.93, txt, transform=ax.transAxes,
                fontsize=7.5, va="top", color=col, fontweight="bold")

    # Row 1: T (generations)
    for col_idx, (label, rows, col, mk) in enumerate(scenarios):
        ax = axes[1, col_idx]
        tt, te, clo_t, chi_t = extract_T(rows, max_k=mk)
        if len(tt):
            elo = np.clip(te - clo_t, 0, None)
            ehi = np.clip(chi_t - te, 0, None)
            ax.errorbar(tt, te, yerr=[elo, ehi],
                        fmt="s", color=col, ms=5, lw=0.8,
                        capsize=2.5, elinewidth=0.8, alpha=0.80)

        lim = [0, 65]
        ax.plot(lim, lim, "k--", lw=0.9, alpha=0.45, zorder=0)
        ax.set_xlim(lim); ax.set_ylim(lim)
        ax.set_xlabel("True T (generations)", fontsize=8.5)
        ax.set_ylabel("Estimated T (generations)", fontsize=8.5)
        ax.set_title(label, fontsize=9.5, fontweight="bold", color=col, pad=6)

        t_maes = [r["T_mae"] for r in rows if r.get("T_mae") is not None]
        t_covs = [r["T_ci_cover"] for r in rows if r.get("T_ci_cover") is not None]
        txt = f"MAE = {np.mean(t_maes):.1f} gen"
        if t_covs:
            txt += f"\n95% CI cover = {np.mean(t_covs):.0%}"
        ax.text(0.04, 0.93, txt, transform=ax.transAxes,
                fontsize=7.5, va="top", color=col, fontweight="bold")

    # Row labels
    fig.text(0.00, 0.74, "Admixture proportion (α)",
             va="center", fontsize=8.5, rotation=90, transform=fig.transFigure)
    fig.text(0.00, 0.28, "Admixture timing (T)",
             va="center", fontsize=8.5, rotation=90, transform=fig.transFigure)

    fig.suptitle("Simulation benchmark: parameter recovery under oracle topology\n"
                 "(20 seeds per scenario; error bars = 95% credible interval)",
                 fontsize=9.5, y=1.02)

    fig.savefig(OUTDIR / "fig2_benchmark.pdf")
    fig.savefig(OUTDIR / "fig2_benchmark.png")
    plt.close(fig)
    print("Figure 2 saved.")


# ─────────────────────────────────────────────────────────────────────────────
# Figure 3: Oracle vs E2E comparison + CI calibration
# ─────────────────────────────────────────────────────────────────────────────
def fig3_ablation():
    """
    Panel A: Topology detection accuracy across S1 (K=0), S2 (K=1), S3 (K=2).
    Panel B: 95% CI coverage for alpha and T (S2 and S3 oracle).
    """
    def load_json(path):
        with open(path) as f:
            return list(json.load(f).values())[0]

    s1 = load_json(RESULTS / "hapgraph/S1_oracle_results.json")
    s2 = load_json(RESULTS / "hapgraph/S2_oracle_results.json")
    s3 = load_json(RESULTS / "hapgraph/S3_oracle_results.json")

    def stats(rows, key):
        vals = [r[key] for r in rows if r.get(key) is not None]
        return (np.mean(vals), np.std(vals)) if vals else (float("nan"), 0)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0))
    fig.subplots_adjust(wspace=0.48)

    # ── Panel A: Topology detection ──────────────────────────────────────────
    ax = axes[0]

    # S1: false positive rate (detected K>0 when K_true=0)
    s1_fp = np.mean([r["detected_K"] > 0 for r in s1 if r["K_true"] == 0])
    s1_k_correct = np.mean([r["greedy_correct_K"] for r in s1])

    s2_recall_m, s2_recall_s = stats(s2, "greedy_recall")
    s2_prec_m,   s2_prec_s   = stats(s2, "greedy_precision")
    s2_k_m,      _           = stats(s2, "greedy_correct_K")
    s3_recall_m, s3_recall_s = stats(s3, "greedy_recall")
    s3_prec_m,   s3_prec_s   = stats(s3, "greedy_precision")
    s3_k_m,      _           = stats(s3, "greedy_correct_K")

    x = np.arange(3)
    w = 0.28
    # Recall (= 1 - FP rate for S1)
    recalls = [1.0 - s1_fp, s2_recall_m, s3_recall_m]
    precs   = [1.0 - s1_fp, s2_prec_m,   s3_prec_m]
    recall_errs = [0, s2_recall_s, s3_recall_s]
    prec_errs   = [0, s2_prec_s,   s3_prec_s]

    b1 = ax.bar(x - w/2, recalls, w, color=BLUE,   label="Recall",
                yerr=recall_errs, capsize=3, error_kw={"lw": 1.0},
                edgecolor="white", zorder=3)
    b2 = ax.bar(x + w/2, precs,   w, color=ORANGE, label="Precision",
                yerr=prec_errs,   capsize=3, error_kw={"lw": 1.0},
                edgecolor="white", zorder=3)

    ax.axhline(1.0, color=GRAY, lw=0.8, ls="--", alpha=0.5, zorder=2)
    ax.set_ylim(0, 1.22)
    ax.set_xticks(x)
    ax.set_xticklabels(["S1\n(K=0 true)", "S2\n(K=1 true)", "S3\n(K=2 true)"], fontsize=8.5)
    ax.set_ylabel("Admixture detection performance", fontsize=8.5)
    ax.set_title("Greedy topology search accuracy\nacross simulation scenarios", fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right", frameon=False)
    for bar in list(b1) + list(b2):
        h = bar.get_height()
        if not np.isnan(h):
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.01,
                    f"{h:.0%}", ha="center", va="bottom", fontsize=7)
    ax.text(-0.32, -0.08, "(A)", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top")

    # ── Panel B: CI coverage ─────────────────────────────────────────────────
    ax = axes[1]
    s2_acov_m, _ = stats(s2, "alpha_ci_cover")
    s2_tcov_m, _ = stats(s2, "T_ci_cover")
    s3_acov_m, _ = stats(s3, "alpha_ci_cover")
    s3_tcov_m, _ = stats(s3, "T_ci_cover")

    scenario_labels = ["S2 (K=1)", "S3 (K=2)"]
    alpha_covs = [s2_acov_m, s3_acov_m]
    t_covs     = [s2_tcov_m, s3_tcov_m]
    xb = np.arange(len(scenario_labels))
    wb = 0.32
    b3 = ax.bar(xb - wb/2, alpha_covs, wb, color=ORANGE, label="α CI coverage",
                edgecolor="white", zorder=3)
    b4 = ax.bar(xb + wb/2, t_covs,     wb, color=BLUE,   label="T CI coverage",
                edgecolor="white", zorder=3)
    ax.axhline(0.95, color="black", lw=1.1, ls="--", alpha=0.6,
               label="Nominal 95%", zorder=2)
    ax.set_ylim(0, 1.18)
    ax.set_xticks(xb); ax.set_xticklabels(scenario_labels, fontsize=9)
    ax.set_ylabel("Empirical 95% CI coverage", fontsize=9)
    ax.set_title("Credible interval calibration\n(oracle topology)", fontsize=9)
    ax.legend(fontsize=7, loc="lower right", frameon=False)
    for bars in [b3, b4]:
        for bar in bars:
            h = bar.get_height()
            if not np.isnan(h):
                ax.text(bar.get_x() + bar.get_width()/2, h + 0.01,
                        f"{h:.0%}", ha="center", va="bottom", fontsize=7.5)
    ax.text(-0.32, -0.08, "(B)", transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top")

    fig.suptitle("Topology detection accuracy and credible interval calibration",
                 fontsize=10, y=1.03)

    fig.savefig(OUTDIR / "fig3_ablation.pdf")
    fig.savefig(OUTDIR / "fig3_ablation.png")
    plt.close(fig)
    print("Figure 3 saved.")


# ── Run all ───────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig1_1kgp_graph()
    fig2_benchmark_scatter()
    fig3_ablation()
    print("\nAll figures saved to", OUTDIR)


# ─────────────────────────────────────────────────────────────────────────────
# Figure 4: Pipeline overview schematic
# ─────────────────────────────────────────────────────────────────────────────
def fig4_pipeline():
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4)
    ax.axis("off")

    # Color scheme
    C_INPUT  = "#E8F4FD"
    C_STAGE1 = "#D4EDDA"
    C_STAGE2 = "#FFF3CD"
    C_OUTPUT = "#F8D7DA"
    BORDER   = "#555555"

    def box(ax, x, y, w, h, label, sublabel, color, fontsize=9):
        rect = mpatches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.08",
            facecolor=color, edgecolor=BORDER, linewidth=1.2, zorder=3
        )
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2 + 0.08, label,
                ha="center", va="center", fontsize=fontsize,
                fontweight="bold", zorder=4)
        if sublabel:
            ax.text(x + w/2, y + h/2 - 0.22, sublabel,
                    ha="center", va="center", fontsize=7,
                    color="#444444", zorder=4, style="italic")

    def arrow(ax, x0, x1, y, label=""):
        ax.annotate("", xy=(x1, y), xytext=(x0, y),
                    arrowprops=dict(arrowstyle="-|>", color=BORDER,
                                   lw=1.3, mutation_scale=12),
                    zorder=2)
        if label:
            ax.text((x0+x1)/2, y+0.14, label, ha="center", va="bottom",
                    fontsize=6.5, color="#666666")

    # ── Boxes ─────────────────────────────────────────────────────────────────
    # Input
    box(ax, 0.15, 1.2, 1.5, 1.6, "Input",
        "Phased VCF\n(N pops, SNPs)", C_INPUT)

    # Stage I: two sub-boxes stacked
    ax.add_patch(mpatches.FancyBboxPatch(
        (2.0, 0.5), 3.5, 3.0,
        boxstyle="round,pad=0.12",
        facecolor="#F0FBF0", edgecolor=GREEN, linewidth=1.8, zorder=2
    ))
    ax.text(3.75, 3.35, "Stage I — F-statistics topology & proportions",
            ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=GREEN, zorder=4)

    box(ax, 2.15, 1.85, 1.5, 1.1, "F₂/F₃ stats",
        "Block-jackknife\n± bias-correction", C_STAGE1, fontsize=8)
    box(ax, 3.85, 1.85, 1.5, 1.1, "NJ tree +\nGreedy search",
        "BIC-adaptive\nstopping", C_STAGE1, fontsize=8)
    box(ax, 2.15, 0.65, 3.2, 0.95, "F₃ MOM estimator",
        "α̂ = (1−√(1−4x))/2   ·   delta-method CI", C_STAGE1, fontsize=8)

    arrow(ax, 3.65, 3.85, 2.40, "F-stats")

    # Stage II
    ax.add_patch(mpatches.FancyBboxPatch(
        (5.8, 0.5), 2.85, 3.0,
        boxstyle="round,pad=0.12",
        facecolor="#FFFEF0", edgecolor=ORANGE, linewidth=1.8, zorder=2
    ))
    ax.text(7.22, 3.35, "Stage II — IBD timing",
            ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=ORANGE, zorder=4)

    box(ax, 5.95, 1.85, 2.55, 1.1, "IBD segment lengths",
        "hap-ibd  ≥ 2 cM\nE[L̄] = 100/T cM", C_STAGE2, fontsize=8)
    box(ax, 5.95, 0.65, 2.55, 0.95, "NUTS-MCMC",
        "4 chains · 1000 draws\nR̂ < 1.05, ESS > 200", C_STAGE2, fontsize=8)

    # Output
    box(ax, 8.9, 1.2, 1.0, 1.6, "Output",
        "Graph + α CI\n+ T CI", C_OUTPUT)

    # ── Arrows ────────────────────────────────────────────────────────────────
    arrow(ax, 1.65, 2.0, 2.0)
    arrow(ax, 3.65, 3.85, 2.40)
    arrow(ax, 5.35, 5.80, 2.0, "topology, α̂")
    arrow(ax, 8.65, 8.90, 2.0)

    # internal Stage I arrows
    ax.annotate("", xy=(3.85, 2.22), xytext=(3.65, 2.22),
                arrowprops=dict(arrowstyle="-|>", color=BORDER, lw=1, mutation_scale=9), zorder=3)
    ax.annotate("", xy=(3.75, 1.85), xytext=(3.75, 1.60),
                arrowprops=dict(arrowstyle="-|>", color=BORDER, lw=1, mutation_scale=9), zorder=3)
    # internal Stage II
    ax.annotate("", xy=(7.22, 1.85), xytext=(7.22, 1.60),
                arrowprops=dict(arrowstyle="-|>", color=BORDER, lw=1, mutation_scale=9), zorder=3)

    ax.set_title("HapGraph two-stage pipeline", fontsize=11, fontweight="bold", pad=6)

    fig.savefig(OUTDIR / "fig4_pipeline.pdf")
    fig.savefig(OUTDIR / "fig4_pipeline.png")
    plt.close(fig)
    print("Figure 4 saved.")


if __name__ == "__main__":
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig1_1kgp_graph()
    fig2_benchmark_scatter()
    fig3_ablation()
    fig4_pipeline()
    print("\nAll figures saved to", OUTDIR)
