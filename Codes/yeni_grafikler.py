import glob
import math
import os
import statistics
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

EKSEN_RENK = "#7a1f1f"
CIKTI = "yeni_figs"

def eksen_boya(ax):
    for k in ax.spines.values():
        k.set_color(EKSEN_RENK)
        k.set_linewidth(1.1)
    ax.tick_params(colors=EKSEN_RENK, labelsize=7.5)
    ax.xaxis.label.set_color(EKSEN_RENK)
    ax.yaxis.label.set_color(EKSEN_RENK)
    ax.title.set_color("#222222")

def ozet(yol):
    d = {}
    for satir in open(yol):
        p = satir.rstrip("\n").split("\t")
        if len(p) >= 2:
            d[p[0]] = p[1]
    return d

def bilinen():
    t = {}
    for satir in open("known_log10.tsv"):
        p = satir.rstrip("\n").split("\t")
        t[p[0]] = float(p[3])
    return t

KT = bilinen()

def hata(d, matris):
    est = d.get("tahmin_log10", "YOK")
    if est == "YOK":
        return 1.0
    fark = float(est) - KT[matris]
    if fark > 30:
        return 1e30
    return abs(10.0 ** fark - 1.0)

def kisa(m):
    return m.replace("block_diagonal_", "blok_").replace("generated_", "gen_").replace("erdos_renyi_", "pER_")

def g1_esik_taramasi(makro):
    hard = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]
    safe = ["block_diagonal_n1000_p0p1", "dwt_503", "olm500"]
    taus = [0.10, 0.20, 0.30, 0.50, 0.80]
    marker = {"generated_500_494_1": "o", "generated_1000_995_5": "s", "er_n500_k3_1": "^",
              "er_n1000_k3_1": "D", "block_diagonal_n1000_p0p1": "v", "dwt_503": "P", "olm500": "X"}
    renk = {"generated_500_494_1": "#d62728", "generated_1000_995_5": "#8c1515", "er_n500_k3_1": "#ff7f0e",
            "er_n1000_k3_1": "#b25900", "block_diagonal_n1000_p0p1": "#1f77b4", "dwt_503": "#2ca02c", "olm500": "#17becf"}
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.9), gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    x_inf = 1.15
    for m in hard + safe:
        pmler, msler = [], []
        for t in taus:
            d = ozet(f"results/T__{m}__t{t:.2f}.ozet")
            pmler.append(100 * float(d["tam_eslesme_orani"]))
            msler.append(float(d["ornek_basina_ms"]))
        cagrilar = []
        for t in taus:
            d = ozet(f"results/T__{m}__t{t:.2f}.ozet")
            cagrilar.append(float(d["tassa_cagri"]) / float(d["ornek"]))
        da = ozet(f"results/A__{m}__adaptive.ozet")
        pm_inf = 100 * float(da["tam_eslesme_orani"])
        cagri_inf = float(da["tassa_cagri"]) / float(da["ornek"])
        taban = 0.04
        cagri_ciz = [max(c, taban) for c in cagrilar] + [max(cagri_inf, taban)]
        stil = "-" if m in hard else "--"
        axs[0].plot(taus + [x_inf], pmler + [pm_inf], stil, marker=marker[m], ms=6,
                    color=renk[m], mfc="none" if m in safe else renk[m], mew=1.4, label=kisa(m))
        axs[1].plot(taus + [x_inf], cagri_ciz, stil, marker=marker[m], ms=6,
                    color=renk[m], mfc="none" if m in safe else renk[m], mew=1.4)
        if m == "generated_500_494_1":
            makro["esikGenAPmUcte"] = f"{pmler[1]:.1f}"
            makro["esikGenAPmNoktaUc"] = f"{pmler[2]:.1f}"
            makro["esikGenAPmYarim"] = f"{pmler[3]:.1f}"
            makro["esikGenAPmInf"] = f"{pm_inf:.1f}"
            makro["esikGenACagriYarim"] = f"{cagrilar[3]:.0f}"
            makro["esikGenACagriNoktaUc"] = f"{cagrilar[2]:.0f}"
            makro["esikGenACagriInf"] = f"{cagri_inf:.0f}"
        if m == "block_diagonal_n1000_p0p1":
            makro["esikBlkCagriKucuk"] = f"{cagrilar[0]:.0f}"
            makro["esikBlkCagriNoktaUc"] = f"{cagrilar[2]:.0f}"
            makro["esikBlkCagriInf"] = f"{cagri_inf:.0f}"
    for eksen in axs[:2]:
        eksen.set_xticks(taus + [x_inf])
        eksen.set_xticklabels([f"{t:.1f}" for t in taus] + [r"$\infty$"])
        eksen.set_xlabel(r"esik $\tau$  (yalniz $s_{sec} < \tau$ ise buda)")
        eksen_boya(eksen)
    axs[0].set_ylabel("tam eslesme orani (%)")
    axs[0].set_title("(a) Tamamlanma vs esik", fontsize=9)
    axs[0].legend(fontsize=6, ncol=1, loc="center right")
    axs[1].set_yscale("log")
    axs[1].set_ylabel("budama cagrisi / ornek (log)")
    axs[1].set_title("(b) Maliyet vs esik (deterministik cagri sayisi)", fontsize=9)
    axs[1].axhline(0.04, color="#bbbbbb", lw=0.7)
    axs[1].text(0.105, 0.05, "0 cagri (tetik hic atesle-\nmiyor: bedava)", fontsize=6, color="#555555")
    for eksen in axs[:2]:
        eksen.axvspan(0.30, 0.50, color="#f4c7c7", alpha=0.35, zorder=0)
    for m in ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]:
        yol = f"results/D__{m}__md_light.dogum"
        if not os.path.exists(yol):
            continue
        vals = [float(l.split()[3]) for l in open(yol) if float(l.split()[3]) >= 0]
        if not vals:
            continue
        ust = min(1.2, max(vals) * 1.05)
        axs[2].hist(vals, bins=np.arange(0, ust + 0.02, 0.02), histtype="step", lw=1.3,
                    density=True, label=kisa(m), color=renk[m])
    axs[2].set_title("(c) MD+Light katil kenar $s$ dagilimi", fontsize=9)
    axs[2].set_xlabel(r"$s = r_i a_{ij} c_j$")
    axs[2].set_ylabel("yogunluk")
    axs[2].legend(fontsize=6)
    eksen_boya(axs[2])
    fig.suptitle(r"Deney A — Esikli adaptif budama: tetik kosuluna $s_{sec}<\tau$ eklenirse ne olur? (3.000 ornek, tohum 42; $\infty$ = esiksiz adaptif)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(f"{CIKTI}/G1_esik_taramasi.png", dpi=150)
    plt.close(fig)

def g2_hata_vs_sure(makro):
    butceler = [5, 15, 30, 60]
    modlar = [("never", "Budama yok", "#8c8c8c", "x"), ("adaptive", "Adaptif", "#d62728", "*"),
              ("always", "Her adimda", "#1f77b4", "o")]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for mod, ad, renk, mk in modlar:
        ort, sd = [], []
        for b in butceler:
            hatalar = []
            for s in (42, 43, 44):
                if b == 30:
                    yol = f"results/P__generated_500_494_1__eq30_{mod}_s{s}.ozet"
                else:
                    yol = f"results/Z__generated_500_494_1__b{b}_{mod}_s{s}.ozet"
                hatalar.append(hata(ozet(yol), "generated_500_494_1"))
            ort.append(statistics.mean(hatalar))
            sd.append(statistics.stdev(hatalar))
        ax.errorbar(butceler, ort, yerr=sd, fmt="-" + mk, ms=9 if mk == "*" else 6.5, capsize=3,
                    color=renk, mfc="none" if mk == "o" else renk, mew=1.5, lw=1.4, label=ad)
        if mod == "adaptive":
            makro["zAdaAltmis"] = f"{ort[-1]:.3f}"
        if mod == "never":
            makro["zNevAltmis"] = f"{ort[-1]:.3f}"
        if mod == "always":
            makro["zAlwAltmis"] = f"{ort[-1]:.3f}"
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(butceler)
    ax.set_xticklabels([str(b) for b in butceler])
    ax.set_xlabel("duvar saati butcesi (saniye)")
    ax.set_ylabel(r"$|\hat{P}-\mathrm{per}|/\mathrm{per}$  (3 tohum ort. $\pm$ ss)")
    ax.set_title("Deney B — Esit surede hata: gen_500_494 (MD+Light tabani)", fontsize=10)
    ax.legend(fontsize=8)
    eksen_boya(ax)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/G2_hata_vs_sure.png", dpi=150)
    plt.close(fig)

def g3_hata_vs_ornek(makro):
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for mod, ad, renk, mk in [("never", "gen_500 — budama yok", "#8c8c8c", "x"),
                              ("adaptive", "gen_500 — adaptif", "#d62728", "*")]:
        xs, ys, sd = [], [], []
        for ns in (3000, 10000, 30000):
            hatalar = []
            for s in (42, 43, 44):
                if ns == 10000:
                    yol = (f"results/P__generated_500_494_1__rep_{mod}_s{s}.ozet"
                           if s != 42 else f"results/A__generated_500_494_1__{mod}.ozet")
                else:
                    yol = f"results/S__generated_500_494_1__n{ns}_{mod}_s{s}.ozet"
                hatalar.append(hata(ozet(yol), "generated_500_494_1"))
            xs.append(ns)
            ys.append(statistics.mean(hatalar))
            sd.append(statistics.stdev(hatalar))
        ax.errorbar(xs, ys, yerr=sd, fmt="-" + mk, ms=9 if mk == "*" else 6.5, capsize=3,
                    color=renk, lw=1.4, label=ad)
        if mod == "adaptive":
            makro["sAdaOtuzK"] = f"{ys[-1]:.3f}"
    xs, ys, sd2 = [], [], []
    for ns in (1000, 3000, 30000):
        hatalar = []
        for sd_ in (42, 43, 44):
            yol = f"results/S__band_n1000_k3__n{ns}_never_s{sd_}.ozet"
            hatalar.append(hata(ozet(yol), "band_n1000_k3"))
        xs.append(ns)
        ys.append(statistics.mean(hatalar))
        sd2.append(statistics.stdev(hatalar))
    ax.errorbar(xs, ys, yerr=sd2, fmt="-^", ms=6.5, capsize=3, color="#2ca02c", mfc="none", mew=1.5, lw=1.4,
                label="band_n1000_k3 — saglikli referans")
    makro["sBandOtuzK"] = f"{ys[-1]:.4f}"
    kilavuz_x = np.array([1000, 30000])
    kilavuz_y = ys[0] * np.sqrt(1000.0 / kilavuz_x)
    ax.plot(kilavuz_x, kilavuz_y, ":", color="#999999", lw=1.2, label=r"$1/\sqrt{N}$ kilavuzu")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("ornek sayisi N")
    ax.set_ylabel(r"$|\hat{P}-\mathrm{per}|/\mathrm{per}$")
    ax.set_title("Deney C — Yakinsama: saglikli sinif $1/\\sqrt{N}$ gibi iner, gen sinifi budamasiz tikanir", fontsize=10)
    ax.legend(fontsize=7.5)
    eksen_boya(ax)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/G3_hata_vs_ornek.png", dpi=150)
    plt.close(fig)

def g4_budama_hassasiyeti(makro):
    matrisler = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n500_k5_1",
                 "er_n1000_k3_1", "block_diagonal_n1000_p0p1", "dwt_503", "olm500"]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 3.9))
    x = np.arange(len(matrisler))
    w = 0.36
    for j, (mod, ad, renk) in enumerate([("adaptive", "Adaptif", "#d62728"), ("always", "Her adimda", "#1f77b4")]):
        oranlar, iptaller = [], []
        for m in matrisler:
            d = ozet(f"results/A__{m}__{mod}.ozet")
            cagri = float(d["tassa_cagri"])
            oranlar.append(float(d["tassa_silinen"]) / cagri if cagri > 0 else 0)
            iptaller.append(100 * float(d["tassa_erken_sifir"]) / float(d["ornek"]))
        axs[0].bar(x + (j - 0.5) * w, oranlar, w, label=ad, color=renk)
        axs[1].bar(x + (j - 0.5) * w, iptaller, w, label=ad, color=renk)
        if mod == "adaptive":
            makro["killPerCallAdaGenA"] = f"{oranlar[0]:.2f}"
        else:
            makro["killPerCallAlwGenA"] = f"{oranlar[0]:.2f}"
    for eksen, baslik, yl in [(axs[0], "(a) Cagri basina silinen kenar — budama 'isabeti'", "silinen kenar / cagri"),
                              (axs[1], "(b) Erken iptal edilen ornek orani (%)", "erken sifir / ornek (%)")]:
        eksen.set_xticks(x)
        eksen.set_xticklabels([kisa(m) for m in matrisler], rotation=55, ha="right", fontsize=7)
        eksen.set_title(baslik, fontsize=9)
        eksen.set_ylabel(yl)
        eksen.legend(fontsize=7.5)
        eksen_boya(eksen)
    fig.suptitle("Deney D — Budama hassasiyeti (mevcut kosu sayaclarindan): adaptif cagrilar daha 'dolu', iptaller aninda", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(f"{CIKTI}/G4_budama_hassasiyeti.png", dpi=150)
    plt.close(fig)

def g5_ds_sapmasi(makro):
    matrisler = ["generated_500_494_1", "er_n500_k3_1", "band_n1000_k3"]
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 3.4), sharey=False)
    for eksen, m in zip(axs, matrisler):
        for v, ad, renk, mk in [("md_once", "MD+Once (donuk)", "#9467bd", "D"),
                                ("md_light", "MD+Light (yerel tazeleme)", "#d62728", "*")]:
            yol = f"results/DS__{m}__{v}.dsat"
            xs, ys = [], []
            n = None
            for satir in open(yol):
                a, b, c = satir.split()
                xs.append(int(a))
                ys.append(float(b))
            n = max(xs)
            eksen.plot([t / n for t in xs], ys, "-", marker=mk, ms=6 if mk == "D" else 8.5,
                       color=renk, mfc="none" if mk == "D" else renk, mew=1.3, lw=1.2, label=ad)
            if m == "generated_500_494_1":
                makro[f"dsSonGenA{'Once' if v=='md_once' else 'Light'}"] = f"{ys[-1]:.2f}"
        eksen.set_title(kisa(m), fontsize=9)
        eksen.set_xlabel("kardinalite / n")
        eksen.set_yscale("log")
        eksen_boya(eksen)
    axs[0].set_ylabel(r"ortalama $|\,\Sigma_j\, r_i a_{ij} c_j - 1\,|$ (canli satirlar, log)")
    axs[0].legend(fontsize=7)
    fig.suptitle("Deney E — Cifte-stokastikten sapma: Once donuk kaliyor ve sapma buyuyor, Light yerel tazelemeyle bastiriyor (300 ornek ort.)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(f"{CIKTI}/G5_ds_sapmasi.png", dpi=150)
    plt.close(fig)

def g6_agirlikli_er(makro):
    ayarlar = [(100, 5), (100, 7), (500, 5), (500, 7), (1000, 5)]
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 3.6))
    x = np.arange(len(ayarlar))
    w = 0.36
    for eksen, v, baslik in [(axs[0], "md", "(a) MD"), (axs[1], "md_light", "(b) MD+Light")]:
        for j, (suf, ad, renk) in enumerate([("", "0/1 desen", "#8c8c8c"), ("_weighted", "agirlikli", "#d62728")]):
            ys = []
            for n, k in ayarlar:
                d = ozet(f"results/W__erdos_renyi_n{n}_k{k}{suf}__{v}.ozet")
                ys.append(100 * float(d["tam_eslesme_orani"]))
            eksen.bar(x + (j - 0.5) * w, ys, w, label=ad, color=renk)
            if v == "md_light":
                makro[f"wErSonPm{'W' if suf else 'P'}"] = f"{ys[-1]:.1f}"
        eksen.set_xticks(x)
        eksen.set_xticklabels([f"n{n} K{k}" for n, k in ayarlar], fontsize=8)
        eksen.set_title(baslik, fontsize=9)
        eksen.set_ylabel("tam eslesme orani (%)")
        eksen.set_ylim(0, 105)
        eksen.legend(fontsize=7.5)
        eksen_boya(eksen)
    fig.suptitle("Deney F — Proje Erdos-Renyi matrisleri: agirlik tamamlanmayi degistiriyor mu? (1.000 ornek, tohum 42)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    fig.savefig(f"{CIKTI}/G6_agirlikli_er.png", dpi=150)
    plt.close(fig)

def main():
    os.makedirs(CIKTI, exist_ok=True)
    makro = {}
    g1_esik_taramasi(makro)
    g2_hata_vs_sure(makro)
    g3_hata_vs_ornek(makro)
    g4_budama_hassasiyeti(makro)
    g5_ds_sapmasi(makro)
    g6_agirlikli_er(makro)
    def tr(x):
        return x.replace(".", "{,}")
    with open("yeni_figs/rakamlar2.tex", "w") as f:
        for k, v in sorted(makro.items()):
            f.write(f"\\newcommand{{\\{k}}}{{{tr(v)}}}\n")
    for k, v in sorted(makro.items()):
        print(k, "=", v)

if __name__ == "__main__":
    main()
