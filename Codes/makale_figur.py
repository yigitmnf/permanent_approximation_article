import csv
import glob
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({"font.family": "serif", "font.size": 8,
                     "axes.titlesize": 8.5, "axes.labelsize": 8,
                     "legend.fontsize": 7, "xtick.labelsize": 6.5,
                     "ytick.labelsize": 7, "pdf.fonttype": 42})

ORDER = ["ras", "once", "light", "md", "md_once", "md_light"]
NAME = {"ras": "Rasmussen", "once": "Once-Scale", "light": "Light-Scale",
        "md": "MD", "md_once": "MD+Once", "md_light": "MD+Light"}
MARKER = {"ras": "x", "once": "s", "light": "^", "md": "o", "md_once": "D", "md_light": "*"}
COLOR = {"ras": "#666666", "once": "#1f77b4", "light": "#17becf",
         "md": "#ff7f0e", "md_once": "#9467bd", "md_light": "#d62728"}
FACE = {"ras": "#666666", "once": "none", "light": "none",
        "md": "#ff7f0e", "md_once": "none", "md_light": "#d62728"}
OUT = "makale/figs"

def rows():
    with open("results/toplu.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def val(r, k):
    v = r.get(k, "")
    return float(v) if v not in ("", None) else None

def short(m):
    m = m.replace("block_diagonal_", "blk_").replace("generated_", "gen_")
    return m

def scatter_grid(ax, mats, data, getter, ylab, logy=False, msz=6):
    for v in ORDER:
        xs, ys = [], []
        for i, m in enumerate(mats):
            r = data.get((m, v))
            if r is None:
                continue
            y = getter(r)
            if y is None:
                continue
            xs.append(i)
            ys.append(y)
        ax.plot(xs, ys, ls="none", marker=MARKER[v], ms=msz if v != "md_light" else msz + 2.5,
                mfc=FACE[v], mec=COLOR[v], mew=1.1, label=NAME[v], alpha=0.9)
    for i in range(len(mats)):
        ax.axvline(i, color="#e3e3e3", lw=0.5, zorder=0)
    ax.set_xticks(range(len(mats)))
    ax.set_xticklabels(mats, rotation=62, ha="right")
    ax.set_ylabel(ylab)
    if logy:
        ax.set_yscale("log")
    ax.tick_params(length=2)

def fig_completion(rs):
    known, real, er = {}, {}, {}
    for r in rs:
        key = (r["matris"], r["varyant"])
        if r["kume"] == "K":
            known[key] = r
        elif r["kume"] == "G":
            real[key] = r
        elif r["kume"] == "E":
            er[key] = r
    order_known = ["band_n100_k2", "band_n100_k3", "band_n500_k3", "band_n1000_k2", "band_n1000_k3",
                   "block_diagonal_n500_p0p1", "block_diagonal_n500_p0p3", "block_diagonal_n500_p0p5",
                   "block_diagonal_n500_p0p9", "block_diagonal_n1000_p0p1", "block_diagonal_n1000_p0p3",
                   "block_diagonal_n1000_p0p5", "block_diagonal_n1000_p0p9",
                   "generated_500_494_1", "generated_500_499_4", "generated_1000_989_4", "generated_1000_995_5"]
    order_real = sorted({k[0] for k in real},
                        key=lambda m: (int(real[(m, "ras")]["n"]) if (m, "ras") in real else 0, m))
    fig, axs = plt.subplots(2, 1, figsize=(7.0, 4.9))
    d1 = {(short(m), v): known[(m, v)] for (m, v) in known}
    scatter_grid(axs[0], [short(m) for m in order_known], d1, lambda r: val(r, "pm_orani"),
                 "perfect-matching rate")
    axs[0].set_title("(a) Synthetic classes with known permanents (band / block-diagonal / generated)")
    axs[0].set_ylim(-0.04, 1.06)
    axs[0].legend(ncol=6, loc="center right", frameon=False, columnspacing=0.9, handletextpad=0.25)
    d2 = {(short(m), v): real[(m, v)] for (m, v) in real}
    scatter_grid(axs[1], [short(m) for m in order_real], d2, lambda r: val(r, "pm_orani"),
                 "perfect-matching rate")
    axs[1].set_title("(b) 33 real SuiteSparse matrices, all with full structural rank")
    axs[1].set_ylim(-0.04, 1.06)
    fig.tight_layout(h_pad=0.7)
    fig.savefig(f"{OUT}/completion.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 2.5))
    xpos, labels = {}, []
    x = 0
    for n in (500, 1000):
        for k in (3, 4, 5):
            for c in range(1, 6):
                xpos[f"er_n{n}_k{k}_{c}"] = x
                labels.append(f"n{n} K{k} \\#{c}".replace("\\#", "#"))
                x += 1
            x += 1
    for v in ORDER:
        xs, ys = [], []
        for m, p in xpos.items():
            r = er.get((m, v))
            if r is None:
                continue
            xs.append(p)
            ys.append(val(r, "pm_orani"))
        ax.plot(xs, ys, ls="none", marker=MARKER[v], ms=5 if v != "md_light" else 7.5,
                mfc=FACE[v], mec=COLOR[v], mew=1.0, label=NAME[v], alpha=0.9)
    for p in xpos.values():
        ax.axvline(p, color="#efefef", lw=0.4, zorder=0)
    ax.set_xticks(list(xpos.values()))
    ax.set_xticklabels(labels, rotation=75, ha="right", fontsize=5.5)
    ax.set_ylabel("perfect-matching rate")
    ax.set_ylim(-0.04, 1.06)
    ax.set_title("Erd\u0151s\u2013R\u00e9nyi instances: 5 random matrices per (n, K) configuration")
    ax.legend(ncol=6, loc="upper right", frameon=False, columnspacing=0.9, handletextpad=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUT}/completion_er.pdf")
    plt.close(fig)

def fig_error(rs):
    data = {}
    for r in rs:
        if r["kume"] == "K" and r["mutlak_bagil_hata"] != "":
            data[(short(r["matris"]), r["varyant"])] = r
    mats = sorted({k[0] for k in data},
                  key=lambda m: (0 if "band" in m else (1 if "blk" in m else 2), m))
    fig, ax = plt.subplots(figsize=(7.0, 2.9))
    scatter_grid(ax, mats, data, lambda r: max(val(r, "mutlak_bagil_hata"), 1e-6) if val(r, "mutlak_bagil_hata") is not None else None,
                 r"$|\hat{P}-\mathrm{per}(A)|\,/\,\mathrm{per}(A)$", logy=True)
    ax.axhline(1.0, color="#999999", ls="--", lw=0.8)
    ax.text(0.1, 1.25, "error = 1 (estimate collapsed to ~0)", fontsize=6.5, color="#555555")
    ax.legend(ncol=6, loc="lower right", frameon=False, columnspacing=0.9, handletextpad=0.25)
    fig.tight_layout()
    fig.savefig(f"{OUT}/error_known.pdf")
    plt.close(fig)

def read_doom(m, v):
    path = f"results/D__{m}__{v}.dogum"
    if not os.path.exists(path):
        return []
    out = []
    for line in open(path):
        a, b, c, d = line.split()
        out.append((int(a), int(b), int(c), float(d)))
    return out

def fig_collapse():
    mats = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1",
            "er_n1000_k3_1", "block_diagonal_n1000_p0p3", "dwt_503", "olm500", "band_n1000_k3"]
    birth = {v: [] for v in ORDER}
    lag = {v: [] for v in ORDER}
    for m in mats:
        for v in ORDER:
            recs = read_doom(m, v)
            if not recs:
                continue
            n = recs[0][2]
            birth[v].append(float(np.mean([a / n for (a, _, _, _) in recs])))
            lag[v].append(float(np.mean([(b - a) / n for (a, b, _, _) in recs])))
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.5))
    for ax, src, ttl, yl in [(axs[0], birth, "(a) Where doom is born", "birth iteration / n"),
                             (axs[1], lag, "(b) Detection lag = wasted work", "(detected \u2212 born) / n")]:
        for i, v in enumerate(ORDER):
            pts = src[v]
            if not pts:
                continue
            ax.plot([i] * len(pts), pts, ls="none", marker=MARKER[v], ms=5.5,
                    mfc=FACE[v], mec=COLOR[v], mew=1.1, alpha=0.85)
            ax.plot([i - 0.28, i + 0.28], [np.mean(pts)] * 2, color=COLOR[v], lw=2)
        ax.set_xticks(range(len(ORDER)))
        ax.set_xticklabels([NAME[v] for v in ORDER], rotation=28, ha="right")
        ax.set_title(ttl)
        ax.set_ylabel(yl)
        ax.tick_params(length=2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/collapse.pdf")
    plt.close(fig)

def fig_doomhist():
    mats = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 2.4))
    for ax, v, ttl in [(axs[0], "once", "(a) Once-Scale"), (axs[1], "light", "(b) Light-Scale")]:
        for m in mats:
            vals = [d for (_, _, _, d) in read_doom(m, v) if d >= 0]
            if not vals:
                continue
            hi = max(0.2, min(1.0, max(vals) * 1.05))
            ax.hist(vals, bins=np.arange(0.0, hi + 0.01, 0.01), histtype="step",
                    lw=1.1, label=short(m), density=True)
        ax.set_title(ttl)
        ax.set_xlabel(r"scaled value $s=r_i\,a_{ij}\,c_j$ of the killer edge")
        ax.set_ylabel("density")
        ax.legend(frameon=False)
        ax.tick_params(length=2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/doomhist.pdf")
    plt.close(fig)

def fig_adaptive(rs):
    modes = ["never", "adaptive", "always"]
    mname = {"never": "no dyn. pruning", "adaptive": "adaptive", "always": "every step"}
    mcol = {"never": "#8c8c8c", "adaptive": "#d62728", "always": "#1f77b4"}
    data = {}
    for r in rs:
        if r["kume"] == "A":
            data[(r["matris"], r["tassa_modu"])] = r
    mats = sorted({k[0] for k in data})
    fig, axs = plt.subplots(1, 3, figsize=(7.0, 2.6))
    x = np.arange(len(mats))
    w = 0.26
    for j, mo in enumerate(modes):
        pm = [val(data[(m, mo)], "pm_orani") if (m, mo) in data else 0 for m in mats]
        axs[0].bar(x + (j - 1) * w, pm, w, label=mname[mo], color=mcol[mo])
        ms = [val(data[(m, mo)], "ornek_ms") if (m, mo) in data else 0 for m in mats]
        axs[1].bar(x + (j - 1) * w, ms, w, color=mcol[mo])
    for j, mo in enumerate(["adaptive", "always"]):
        sh = []
        for m in mats:
            r = data.get((m, mo))
            sh.append(100 * (val(r, "sure_tassa") or 0) / (val(r, "sure_toplam") or 1) if r else 0)
        axs[2].bar(x + (j - 0.5) * w, sh, w, color=mcol[mo])
    for ax, t in zip(axs, ["(a) perfect-matching rate", "(b) time per sample (ms, log)",
                           "(c) pruning share of runtime (%)"]):
        ax.set_xticks(x)
        ax.set_xticklabels([short(m) for m in mats], rotation=62, ha="right", fontsize=5.8)
        ax.set_title(t)
        ax.tick_params(length=2)
    axs[1].set_yscale("log")
    axs[0].legend(frameon=False, fontsize=6.3)
    fig.tight_layout()
    fig.savefig(f"{OUT}/adaptive.pdf")
    plt.close(fig)

def fig_cost(rs):
    data = {}
    for r in rs:
        if r["kume"] in ("K", "C"):
            data[(r["matris"], r["varyant"])] = r
    mats = ["band_n1000_k3", "block_diagonal_n1000_p0p3", "generated_500_494_1", "dwt_503"]
    vs = ["md", "md_once", "md_light", "md_every"]
    nm = {"md": "MD", "md_once": "MD+Once", "md_light": "MD+Light", "md_every": "MD+FullEvery"}
    cols = ["#ff7f0e", "#9467bd", "#d62728", "#1f77b4"]
    fig, ax = plt.subplots(figsize=(3.35, 2.3))
    x = np.arange(len(mats))
    w = 0.2
    for j, v in enumerate(vs):
        ys = [val(data.get((m, v)), "ornek_ms") if (m, v) in data else 0 for m in mats]
        b = ax.bar(x + (j - 1.5) * w, ys, w, label=nm[v], color=cols[j])
        for bb, y in zip(b, ys):
            if y:
                ax.text(bb.get_x() + bb.get_width() / 2, y * 1.08, f"{y:.2f}", ha="center", fontsize=5)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([short(m) for m in mats], fontsize=6, rotation=12)
    ax.set_ylabel("ms per sample (log)")
    ax.legend(frameon=False, fontsize=6, ncol=2)
    ax.tick_params(length=2)
    fig.tight_layout()
    fig.savefig(f"{OUT}/cost.pdf")
    plt.close(fig)

def fig_livennz():
    mats = ["generated_1000_995_5", "er_n1000_k3_1", "block_diagonal_n1000_p0p3", "dwt_503"]
    fig, axs = plt.subplots(1, 4, figsize=(7.0, 1.9), sharey=False)
    for ax, m in zip(axs, mats):
        for v in ORDER:
            path = f"results/D__{m}__{v}.egri"
            if not os.path.exists(path):
                continue
            ts, ys = [], []
            for line in open(path):
                a, b = line.split()
                ts.append(int(a))
                ys.append(float(b))
            n = ts[-1]
            ax.plot([t / n for t in ts], ys, lw=1.0, color=COLOR[v],
                    marker=MARKER[v], markevery=max(1, n // 8), ms=3.6,
                    mfc=FACE[v], mec=COLOR[v], label=NAME[v])
        ax.set_title(short(m), fontsize=7)
        ax.set_xlabel("iteration / n", fontsize=6.5)
        ax.tick_params(length=2, labelsize=6)
    axs[0].set_ylabel("live nnz (zero-padded mean)", fontsize=6.5)
    axs[0].legend(frameon=False, fontsize=5, loc="upper right")
    fig.tight_layout()
    fig.savefig(f"{OUT}/livennz.pdf")
    plt.close(fig)

def main():
    os.makedirs(OUT, exist_ok=True)
    rs = rows()
    fig_completion(rs)
    fig_error(rs)
    fig_collapse()
    fig_doomhist()
    fig_adaptive(rs)
    fig_cost(rs)
    fig_livennz()
    print("paper figures written to", OUT)

if __name__ == "__main__":
    main()
