import glob
import math
import os

def known_log10_tablosu():
    tablo = {}
    for satir in open("known_log10.tsv"):
        parca = satir.rstrip("\n").split("\t")
        tablo[parca[0]] = float(parca[3])
    return tablo

def ozet_oku(yol):
    kayit = {}
    for satir in open(yol):
        parca = satir.rstrip("\n").split("\t")
        if len(parca) >= 2:
            kayit[parca[0]] = parca[1]
    return kayit

def mutlak_hata(tahmin_log10, bilinen_log10):
    if tahmin_log10 is None:
        return 1.0
    fark = tahmin_log10 - bilinen_log10
    if fark > 30:
        return 10.0 ** 30
    return abs(10.0 ** fark - 1.0)

def main():
    bilinen = known_log10_tablosu()
    cikti = open("results/toplu.tsv", "w")
    basliklar = ["kume", "matris", "varyant", "tassa_modu", "n", "nnz_tassa", "ornek",
                 "tam_eslesme", "pm_orani", "ort_kard", "tahmin_log10", "bilinen_log10",
                 "mutlak_bagil_hata", "sure_toplam", "sure_scaling", "sure_tassa",
                 "tassa_cagri", "tassa_silinen", "tassa_firsat", "ornek_ms"]
    cikti.write("\t".join(basliklar) + "\n")
    for yol in sorted(glob.glob("results/*.ozet")):
        ad = os.path.basename(yol)[:-5]
        parcalar = ad.split("__")
        kume = parcalar[0]
        matris = parcalar[1]
        etiket = parcalar[2] if len(parcalar) > 2 else ""
        k = ozet_oku(yol)
        if "sonuc" in k:
            continue
        varyant = k.get("varyant", etiket)
        tassa_modu = k.get("tassa_modu", "never")
        if kume == "A":
            tassa_modu = etiket
        tahmin = k.get("tahmin_log10", "YOK")
        tahmin_deger = None if tahmin == "YOK" else float(tahmin)
        bilinen_deger = bilinen.get(matris)
        hata = ""
        if bilinen_deger is not None:
            hata = f"{mutlak_hata(tahmin_deger, bilinen_deger):.6e}"
        satir = [kume, matris, varyant, tassa_modu,
                 k.get("n", ""), k.get("nnz_tassa_sonrasi", ""), k.get("ornek", ""),
                 k.get("tam_eslesme", ""), k.get("tam_eslesme_orani", ""),
                 k.get("ortalama_kardinalite", ""),
                 "" if tahmin_deger is None else f"{tahmin_deger:.6f}",
                 "" if bilinen_deger is None else f"{bilinen_deger:.6f}",
                 hata,
                 k.get("sure_toplam_sn", ""), k.get("sure_scaling_sn", ""), k.get("sure_tassa_sn", ""),
                 k.get("tassa_cagri", ""), k.get("tassa_silinen", ""), k.get("tassa_firsat", ""),
                 k.get("ornek_basina_ms", "")]
        cikti.write("\t".join(satir) + "\n")
    cikti.close()
    print("toplu.tsv yazildi")

if __name__ == "__main__":
    main()
