import csv
import glob
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

VARYANT_SIRA = ["ras", "once", "light", "md", "md_once", "md_light"]
VARYANT_AD = {"ras": "Rasmussen", "once": "Once-Scale", "light": "Light-Scale",
              "md": "MD", "md_once": "MD+Once", "md_light": "MD+Light"}
MARKER = {"ras": "x", "once": "s", "light": "^", "md": "o", "md_once": "D", "md_light": "*"}
RENK = {"ras": "#666666", "once": "#1f77b4", "light": "#17becf",
        "md": "#ff7f0e", "md_once": "#9467bd", "md_light": "#d62728"}
DOLGU = {"ras": "#666666", "once": "none", "light": "none",
         "md": "#ff7f0e", "md_once": "none", "md_light": "#d62728"}

EKSEN_RENK = "#7a1f1f"
CIKTI = "figs"

def eksenleri_boya(ax):
    for kenar in ax.spines.values():
        kenar.set_color(EKSEN_RENK)
        kenar.set_linewidth(1.2)
    ax.tick_params(colors=EKSEN_RENK, labelsize=8)
    ax.xaxis.label.set_color(EKSEN_RENK)
    ax.yaxis.label.set_color(EKSEN_RENK)

def tabloyu_oku():
    kayitlar = []
    with open("results/toplu.tsv") as f:
        okuyucu = csv.DictReader(f, delimiter="\t")
        for satir in okuyucu:
            kayitlar.append(satir)
    return kayitlar

def deger(satir, alan):
    v = satir.get(alan, "")
    return float(v) if v not in ("", None) else None

def dikey_hizali_ciz(ax, matris_listesi, veri, y_alici, y_etiket, log_scale=False):
    x_konum = np.arange(len(matris_listesi))
    for v in VARYANT_SIRA:
        xs, ys = [], []
        for i, m in enumerate(matris_listesi):
            kayit = veri.get((m, v))
            if kayit is None:
                continue
            y = y_alici(kayit)
            if y is None:
                continue
            xs.append(i)
            ys.append(y)
        ax.plot(xs, ys, linestyle="none", marker=MARKER[v], markersize=8 if v != "md_light" else 11,
                markerfacecolor=DOLGU[v], markeredgecolor=RENK[v], markeredgewidth=1.6,
                label=VARYANT_AD[v], alpha=0.9)
    for i in range(len(matris_listesi)):
        ax.axvline(i, color="#dddddd", linewidth=0.6, zorder=0)
    ax.set_xticks(x_konum)
    ax.set_xticklabels(matris_listesi, rotation=60, ha="right", fontsize=7)
    ax.set_ylabel(y_etiket)
    if log_scale:
        ax.set_yscale("log")
    eksenleri_boya(ax)

def kisa_ad(m):
    m = m.replace("block_diagonal_", "blok_").replace("generated_", "gen_")
    m = m.replace("erdos_renyi_", "er_")
    return m

def sekil_tamamlanma(kayitlar):
    ana = {}
    er = {}
    for s in kayitlar:
        anahtar = (s["matris"], s["varyant"])
        if s["kume"] == "K" or s["kume"] == "G":
            ana[anahtar] = s
        elif s["kume"] == "E":
            er[anahtar] = s
    bilinen_sira = [m for m in
        ["band_n100_k2","band_n100_k3","band_n500_k3","band_n1000_k2","band_n1000_k3",
         "block_diagonal_n500_p0p1","block_diagonal_n500_p0p3","block_diagonal_n500_p0p5","block_diagonal_n500_p0p9",
         "block_diagonal_n1000_p0p1","block_diagonal_n1000_p0p3","block_diagonal_n1000_p0p5","block_diagonal_n1000_p0p9",
         "generated_500_494_1","generated_500_499_4","generated_1000_989_4","generated_1000_995_5"]
        if any(k[0] == m for k in ana)]
    gercek_sira = sorted({k[0] for k in ana} - set(bilinen_sira),
                         key=lambda m: (int(ana[(m, "ras")]["n"]) if (m, "ras") in ana else 0, m))
    fig, axs = plt.subplots(2, 1, figsize=(13, 9))
    d1 = {(kisa_ad(m), v): ana[(m, v)] for (m, v) in ana if m in bilinen_sira}
    dikey_hizali_ciz(axs[0], [kisa_ad(m) for m in bilinen_sira], d1,
                     lambda s: deger(s, "pm_orani"), "tam eslesme orani")
    axs[0].set_title("Tamamlanma orani (perfect matching) — bilinen-permanentli siniflar, 10.000 ornek", fontsize=10)
    axs[0].set_ylim(-0.03, 1.05)
    axs[0].legend(ncol=6, fontsize=8, loc="lower left")
    d2 = {(kisa_ad(m), v): ana[(m, v)] for (m, v) in ana if m in gercek_sira}
    dikey_hizali_ciz(axs[1], [kisa_ad(m) for m in gercek_sira], d2,
                     lambda s: deger(s, "pm_orani"), "tam eslesme orani")
    axs[1].set_title("Tamamlanma orani — SuiteSparse GERCEK matrisler (structural rank = FULL), 10.000 ornek", fontsize=10)
    axs[1].set_ylim(-0.03, 1.05)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig1_tamamlanma_ana.png", dpi=150)
    plt.close(fig)

    gruplar = [(n, k) for n in (500, 1000) for k in (3, 4, 5)]
    fig, ax = plt.subplots(figsize=(12, 5))
    etiketler = []
    x = 0
    x_konumlar = {}
    for n, k in gruplar:
        for kopya in range(1, 6):
            m = f"er_n{n}_k{k}_{kopya}"
            x_konumlar[m] = x
            etiketler.append(f"n{n} K{k} #{kopya}")
            x += 1
        x += 1
    for v in VARYANT_SIRA:
        xs, ys = [], []
        for m, konum in x_konumlar.items():
            s = er.get((m, v))
            if s is None:
                continue
            y = deger(s, "pm_orani")
            xs.append(konum)
            ys.append(y)
        ax.plot(xs, ys, linestyle="none", marker=MARKER[v], markersize=7 if v != "md_light" else 10,
                markerfacecolor=DOLGU[v], markeredgecolor=RENK[v], markeredgewidth=1.4,
                label=VARYANT_AD[v], alpha=0.9)
    for m, konum in x_konumlar.items():
        ax.axvline(konum, color="#eeeeee", linewidth=0.5, zorder=0)
    ax.set_xticks(list(x_konumlar.values()))
    ax.set_xticklabels(etiketler, rotation=75, ha="right", fontsize=6)
    ax.set_ylabel("tam eslesme orani")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("Tamamlanma orani — Erdos-Renyi (K=3,4,5; her ayar 5 farkli rastgele matris), 10.000 ornek", fontsize=10)
    ax.legend(ncol=6, fontsize=8, loc="upper right")
    eksenleri_boya(ax)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig2_tamamlanma_erdos.png", dpi=150)
    plt.close(fig)

def sekil_mutlak_hata(kayitlar):
    veri = {}
    for s in kayitlar:
        if s["kume"] == "K" and s["mutlak_bagil_hata"] != "":
            veri[(kisa_ad(s["matris"]), s["varyant"])] = s
    sira = sorted({k[0] for k in veri}, key=lambda m: (("band" not in m) * 1 + ("blok" in m) * 0 + ("gen" in m) * 2, m))
    fig, ax = plt.subplots(figsize=(13, 5.5))
    def hata_al(s):
        h = deger(s, "mutlak_bagil_hata")
        if h is None:
            return None
        return max(h, 1e-6)
    dikey_hizali_ciz(ax, sira, veri, hata_al, "|tahmin - gercek| / gercek", log_scale=True)
    ax.axhline(1.0, color="#999999", linestyle="--", linewidth=1)
    ax.text(0.1, 1.15, "hata = 1  (cokmus: tahmin 0)", fontsize=8, color="#555555")
    ax.set_title("Mutlak bagil hata |Est - Gercek| / Gercek — bilinen permanentli 17 matris, 10.000 ornek, seed=42", fontsize=10)
    ax.legend(ncol=6, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig3_mutlak_hata.png", dpi=150)
    plt.close(fig)

def sekil_canli_nnz():
    hedefler = ["generated_1000_995_5", "er_n1000_k3_1", "block_diagonal_n1000_p0p3", "dwt_503"]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8))
    for eksen, m in zip(axs.flat, hedefler):
        for v in VARYANT_SIRA:
            yol = f"results/D__{m}__{v}.egri"
            if not os.path.exists(yol):
                continue
            ts, ys = [], []
            for satir in open(yol):
                a, b = satir.split()
                ts.append(int(a))
                ys.append(float(b))
            n = ts[-1]
            eksen.plot([t / n for t in ts], ys, linewidth=1.4, color=RENK[v],
                       marker=MARKER[v], markevery=max(1, n // 12), markersize=6,
                       markerfacecolor=DOLGU[v], markeredgecolor=RENK[v], label=VARYANT_AD[v])
        eksen.set_title(kisa_ad(m), fontsize=10)
        eksen.set_xlabel("iterasyon / n")
        eksen.set_ylabel("canli nnz (0-dolgulu ortalama)")
        eksenleri_boya(eksen)
    axs[0, 0].legend(fontsize=8)
    fig.suptitle("Canli NNZ vs iterasyon — 2.000 ornek ortalamasi (cokusten sonra 0 ile doldurulur)", fontsize=11)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig4_canli_nnz.png", dpi=150)
    plt.close(fig)

def dogum_oku(m, v):
    yol = f"results/D__{m}__{v}.dogum"
    if not os.path.exists(yol):
        return []
    kayitlar = []
    for satir in open(yol):
        a, b, c, d = satir.split()
        kayitlar.append((int(a), int(b), int(c), float(d)))
    return kayitlar

def sekil_dogum_scaling():
    matrisler = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    for eksen, v in zip(axs, ["once", "light"]):
        for m in matrisler:
            degerler = [d for (_, _, _, d) in dogum_oku(m, v) if d >= 0]
            if not degerler:
                continue
            ust = max(0.2, min(1.0, max(degerler) * 1.05))
            kovalar = np.arange(0.0, ust + 0.01, 0.01)
            eksen.hist(degerler, bins=kovalar, histtype="step", linewidth=1.5, label=kisa_ad(m), density=True)
        eksen.set_title(f"Cokusu doguran kenarin scaling degeri — {VARYANT_AD[v]} (kova genisligi 0.01)", fontsize=9)
        eksen.set_xlabel("s = rv * a * cv (secim anindaki deger)")
        eksen.set_ylabel("yogunluk")
        eksen.legend(fontsize=8)
        eksenleri_boya(eksen)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig5_dogum_scaling_kova001.png", dpi=150)
    plt.close(fig)

def sekil_dogum_gecikme():
    matrisler = ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1",
                 "er_n1000_k3_1", "block_diagonal_n1000_p0p3", "dwt_503", "olm500", "band_n1000_k3"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    dogum_ort = {v: [] for v in VARYANT_SIRA}
    gecikme_ort = {v: [] for v in VARYANT_SIRA}
    for m in matrisler:
        for v in VARYANT_SIRA:
            kayitlar = dogum_oku(m, v)
            if not kayitlar:
                continue
            n = kayitlar[0][2]
            dogumlar = [a / n for (a, _, _, _) in kayitlar]
            gecikmeler = [(b - a) / n for (a, b, _, _) in kayitlar]
            dogum_ort[v].append(float(np.mean(dogumlar)))
            gecikme_ort[v].append(float(np.mean(gecikmeler)))
    x = np.arange(len(VARYANT_SIRA))
    for eksen, kaynak, baslik, etiket in [
        (axs[0], dogum_ort, "Cokusun dogum yeri (iterasyon / n) — matris ortalamalarinin dagilimi",
         "dogum konumu / n"),
        (axs[1], gecikme_ort, "Fark edilme gecikmesi: (tespit - dogum) / n — 'bosa harcanan is'",
         "gecikme / n")]:
        for i, v in enumerate(VARYANT_SIRA):
            noktalar = kaynak[v]
            if not noktalar:
                continue
            eksen.plot([i] * len(noktalar), noktalar, linestyle="none", marker=MARKER[v],
                       markersize=8, markerfacecolor=DOLGU[v], markeredgecolor=RENK[v],
                       markeredgewidth=1.5, alpha=0.85)
            eksen.plot([i - 0.25, i + 0.25], [np.mean(noktalar)] * 2, color=RENK[v], linewidth=2.5)
        eksen.set_xticks(x)
        eksen.set_xticklabels([VARYANT_AD[v] for v in VARYANT_SIRA], rotation=25, ha="right", fontsize=8)
        eksen.set_title(baslik, fontsize=9)
        eksen.set_ylabel(etiket)
        eksenleri_boya(eksen)
    fig.suptitle("Cokus analizi — cizgi: 8 matrisin ortalamasi, nokta: tek matris ortalamasi (2.000 ornek/koşu)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig6_dogum_gecikme.png", dpi=150)
    plt.close(fig)

def sekil_adaptif_tassa(kayitlar):
    modlar = ["never", "adaptive", "always"]
    mod_ad = {"never": "TASSA yok", "adaptive": "Adaptif TASSA", "always": "Her adimda TASSA"}
    mod_renk = {"never": "#888888", "adaptive": "#d62728", "always": "#1f77b4"}
    veri = {}
    for s in kayitlar:
        if s["kume"] == "A":
            veri[(s["matris"], s["tassa_modu"])] = s
    matrisler = sorted({k[0] for k in veri})
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    x = np.arange(len(matrisler))
    genislik = 0.26
    for j, mod in enumerate(modlar):
        pm = [deger(veri[(m, mod)], "pm_orani") if (m, mod) in veri else 0 for m in matrisler]
        axs[0].bar(x + (j - 1) * genislik, pm, genislik, label=mod_ad[mod], color=mod_renk[mod])
        ms = [deger(veri[(m, mod)], "ornek_ms") if (m, mod) in veri else 0 for m in matrisler]
        axs[1].bar(x + (j - 1) * genislik, ms, genislik, label=mod_ad[mod], color=mod_renk[mod])
    for j, mod in enumerate(["adaptive", "always"]):
        pay = []
        for m in matrisler:
            s = veri.get((m, mod))
            if s is None:
                pay.append(0)
                continue
            t = deger(s, "sure_toplam") or 1
            pay.append(100.0 * (deger(s, "sure_tassa") or 0) / t)
        axs[2].bar(x + (j - 0.5) * genislik, pay, genislik, label=mod_ad[mod], color=mod_renk[mod])
    basliklar = ["Tam eslesme orani", "Ornek basina sure (ms)", "TASSA'nin toplam suredeki payi (%)"]
    for eksen, b in zip(axs, basliklar):
        eksen.set_xticks(x)
        eksen.set_xticklabels([kisa_ad(m) for m in matrisler], rotation=60, ha="right", fontsize=7)
        eksen.set_title(b + (" [log]" if "sure (ms)" in b else ""), fontsize=9)
        eksen.legend(fontsize=7)
        eksenleri_boya(eksen)
    axs[1].set_yscale("log")
    fig.suptitle("Adaptif (kosullu) TASSA — kural: derece=1 ya da en buyuk olasilikli aday secildiyse TASSA yapma; kucugu sectiysen yap (taban: MD+Light)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig7_adaptif_tassa.png", dpi=150)
    plt.close(fig)

def sekil_scaling_maliyeti(kayitlar):
    veri = {}
    for s in kayitlar:
        if s["kume"] in ("K", "C"):
            veri[(s["matris"], s["varyant"])] = s
    matrisler = ["band_n1000_k3", "block_diagonal_n1000_p0p3", "generated_500_494_1", "dwt_503"]
    varyantlar = ["md", "md_once", "md_light", "md_every"]
    adlar = {"md": "MD", "md_once": "MD+Once", "md_light": "MD+Light", "md_every": "MD+HerAdimTamScaling"}
    fig, ax = plt.subplots(figsize=(10, 4.6))
    x = np.arange(len(matrisler))
    genislik = 0.2
    renkler = ["#ff7f0e", "#9467bd", "#d62728", "#1f77b4"]
    for j, v in enumerate(varyantlar):
        ys = []
        for m in matrisler:
            s = veri.get((m, v))
            ys.append(deger(s, "ornek_ms") if s else 0)
        cubuklar = ax.bar(x + (j - 1.5) * genislik, ys, genislik, label=adlar[v], color=renkler[j])
        for c, y in zip(cubuklar, ys):
            if y:
                ax.text(c.get_x() + c.get_width() / 2, y * 1.05, f"{y:.2f}", ha="center", fontsize=7)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([kisa_ad(m) for m in matrisler], fontsize=8)
    ax.set_ylabel("ornek basina sure (ms, log)")
    ax.set_title("Scaling maliyeti: her adimda TAM Sinkhorn (5 ic iterasyon) vs Light vs Once — ayni tohum", fontsize=10)
    ax.legend(fontsize=8)
    eksenleri_boya(ax)
    fig.tight_layout()
    fig.savefig(f"{CIKTI}/fig8_scaling_maliyeti.png", dpi=150)
    plt.close(fig)

def tablo_t4(kayitlar):
    gruplar = {"band": [], "blok": [], "generated": [], "erdos": [], "gercek": []}
    for s in kayitlar:
        if s["kume"] not in ("K", "E", "G"):
            continue
        m = s["matris"]
        if m.startswith("band"): g = "band"
        elif m.startswith("block"): g = "blok"
        elif m.startswith("generated"): g = "generated"
        elif m.startswith("er_"): g = "erdos"
        else: g = "gercek"
        gruplar[g].append(s)
    satirlar = []
    for v in VARYANT_SIRA:
        satir = [VARYANT_AD[v]]
        for g in ["band", "blok", "generated", "erdos", "gercek"]:
            pmler = [deger(s, "pm_orani") for s in gruplar[g] if s["varyant"] == v]
            pmler = [p for p in pmler if p is not None]
            if not pmler:
                satir += ["-", "-", "-"]
                continue
            satir += [f"{np.median(pmler):.3f}", f"{np.mean(pmler):.3f}",
                      str(sum(1 for p in pmler if p >= 0.99))]
        satirlar.append(satir)
    with open(f"{CIKTI}/tablo_t4.tsv", "w") as f:
        f.write("varyant\t" + "\t".join(f"{g}_{c}" for g in ["band", "blok", "gen", "erdos", "gercek"]
                                        for c in ["medyan", "ort", "n99"]) + "\n")
        for s in satirlar:
            f.write("\t".join(s) + "\n")
    print("tablo_t4.tsv yazildi")

def main():
    os.makedirs(CIKTI, exist_ok=True)
    kayitlar = tabloyu_oku()
    sekil_tamamlanma(kayitlar)
    sekil_mutlak_hata(kayitlar)
    sekil_canli_nnz()
    sekil_dogum_scaling()
    sekil_dogum_gecikme()
    sekil_adaptif_tassa(kayitlar)
    sekil_scaling_maliyeti(kayitlar)
    tablo_t4(kayitlar)
    print("tum grafikler uretildi")

if __name__ == "__main__":
    main()
