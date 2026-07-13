import csv

def oku():
    with open("results/toplu.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def bul(kayitlar, kume, matris, mod_veya_varyant):
    for s in kayitlar:
        if s["kume"] == kume and s["matris"] == matris:
            if kume == "A" and s["tassa_modu"] == mod_veya_varyant:
                return s
            if kume != "A" and s["varyant"] == mod_veya_varyant:
                return s
    return None

def yuzde(s):
    return f"{100.0 * float(s['pm_orani']):.1f}"

def ms(s):
    return f"{float(s['ornek_ms']):.2f}"

def hata(s):
    return f"{float(s['mutlak_bagil_hata']):.3f}"

def main():
    k = oku()
    m = {}
    m["AGenBinNever"] = yuzde(bul(k, "A", "generated_1000_995_5", "never"))
    m["AGenBinAdaptif"] = yuzde(bul(k, "A", "generated_1000_995_5", "adaptive"))
    m["AGenBinAlways"] = yuzde(bul(k, "A", "generated_1000_995_5", "always"))
    m["AGenBesNever"] = yuzde(bul(k, "A", "generated_500_494_1", "never"))
    m["AGenBesAdaptif"] = yuzde(bul(k, "A", "generated_500_494_1", "adaptive"))
    m["AGenBesAlways"] = yuzde(bul(k, "A", "generated_500_494_1", "always"))
    m["AErBinNever"] = yuzde(bul(k, "A", "er_n1000_k3_1", "never"))
    m["AErBinAdaptif"] = yuzde(bul(k, "A", "er_n1000_k3_1", "adaptive"))
    m["AErBinAlways"] = yuzde(bul(k, "A", "er_n1000_k3_1", "always"))
    m["AGenBesMsNever"] = ms(bul(k, "A", "generated_500_494_1", "never"))
    m["AGenBesMsAdaptif"] = ms(bul(k, "A", "generated_500_494_1", "adaptive"))
    m["AGenBesMsAlways"] = ms(bul(k, "A", "generated_500_494_1", "always"))
    m["AGenBinMsNever"] = ms(bul(k, "A", "generated_1000_995_5", "never"))
    m["AGenBinMsAdaptif"] = ms(bul(k, "A", "generated_1000_995_5", "adaptive"))
    m["AGenBinMsAlways"] = ms(bul(k, "A", "generated_1000_995_5", "always"))
    m["AGenBesHataNever"] = hata(bul(k, "A", "generated_500_494_1", "never"))
    m["AGenBesHataAdaptif"] = hata(bul(k, "A", "generated_500_494_1", "adaptive"))
    m["AGenBesHataAlways"] = hata(bul(k, "A", "generated_500_494_1", "always"))
    s = bul(k, "A", "generated_500_494_1", "adaptive")
    m["AGenBesTassaPay"] = f"{100.0 * float(s['sure_tassa']) / float(s['sure_toplam']):.0f}"
    s = bul(k, "A", "generated_1000_995_5", "adaptive")
    m["AGenBinTassaPay"] = f"{100.0 * float(s['sure_tassa']) / float(s['sure_toplam']):.0f}"

    m["KBandBinIkiHata"] = hata(bul(k, "K", "band_n1000_k2", "md_light"))
    m["KBandBinUcHata"] = hata(bul(k, "K", "band_n1000_k3", "md_light"))
    m["KBlokBinBirHata"] = hata(bul(k, "K", "block_diagonal_n1000_p0p1", "md_light"))
    m["KGenBinHata"] = hata(bul(k, "K", "generated_1000_995_5", "md_light"))
    m["KGenBesHata"] = hata(bul(k, "K", "generated_500_494_1", "md_light"))

    every = bul(k, "C", "band_n1000_k3", "md_every")
    m["CBandEveryMs"] = ms(every)
    m["CBandLightMs"] = ms(bul(k, "K", "band_n1000_k3", "md_light"))
    m["CBandOnceMs"] = ms(bul(k, "K", "band_n1000_k3", "md_once"))
    m["CBandMdMs"] = ms(bul(k, "K", "band_n1000_k3", "md"))
    m["CGenEveryMs"] = ms(bul(k, "C", "generated_500_494_1", "md_every"))
    m["CGenLightMs"] = ms(bul(k, "K", "generated_500_494_1", "md_light"))

    grup_pm = {}
    for s in k:
        if s["kume"] not in ("K", "E", "G"):
            continue
        ad = s["matris"]
        if ad.startswith("er_"):
            g = "erdos"
        elif ad.startswith("band"):
            g = "band"
        elif ad.startswith("block"):
            g = "blok"
        elif ad.startswith("generated"):
            g = "gen"
        else:
            g = "gercek"
        grup_pm.setdefault((g, s["varyant"]), []).append(float(s["pm_orani"]))
    for (g, v), lst in grup_pm.items():
        anahtar = f"PM{g}{v}".replace("_", "")
        m[anahtar] = f"{100.0 * sum(lst) / len(lst):.1f}"

    dogum = {}
    for mat in ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]:
        for v in ["ras", "md_light"]:
            yol = f"results/D__{mat}__{v}.dogum"
            try:
                satirlar = [l.split() for l in open(yol)]
            except FileNotFoundError:
                continue
            n = int(satirlar[0][2])
            gec = sum((int(b) - int(a)) / n for a, b, _, _ in satirlar) / len(satirlar)
            dogum.setdefault(v, []).append(gec)
    m["GecikmeRas"] = f"{sum(dogum['ras']) / len(dogum['ras']):.2f}"
    m["GecikmeMdLight"] = f"{sum(dogum['md_light']) / len(dogum['md_light']):.3f}"

    toplam_sure = 0.0
    for satir in open("results/kosu_sureleri.tsv"):
        toplam_sure += float(satir.split("\t")[1])
    m["ToplamKosuDakika"] = f"{toplam_sure / 60.0:.0f}"
    m["ToplamKosuSayisi"] = str(sum(1 for _ in open("results/kosu_sureleri.tsv")))

    def tr(x):
        return x.replace(".", "{,}")
    with open("rapor/rakamlar.tex", "w") as f:
        for ad, deger in sorted(m.items()):
            f.write(f"\\newcommand{{\\{ad}}}{{{tr(deger)}}}\n")
    print("rakamlar.tex:", len(m), "makro")

if __name__ == "__main__":
    main()
