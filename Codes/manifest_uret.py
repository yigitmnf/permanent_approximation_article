import os

satirlar = []

BILINEN = ["band_n100_k2","band_n100_k3","band_n500_k3","band_n1000_k2","band_n1000_k3",
 "block_diagonal_n500_p0p1","block_diagonal_n500_p0p3","block_diagonal_n500_p0p5","block_diagonal_n500_p0p9",
 "block_diagonal_n1000_p0p1","block_diagonal_n1000_p0p3","block_diagonal_n1000_p0p5","block_diagonal_n1000_p0p9",
 "generated_500_494_1","generated_500_499_4","generated_1000_989_4","generated_1000_995_5"]
V = ["ras","once","light","md","md_once","md_light"]
for m in BILINEN:
    for v in V:
        satirlar.append((f"results/K__{m}__{v}", f"mats/{m}.mtx {v} 10000 120 42 results/K__{m}__{v}"))

ADAPTIF = ["generated_500_494_1","generated_1000_995_5","er_n500_k3_1","er_n500_k5_1",
           "er_n1000_k3_1","block_diagonal_n1000_p0p1","dwt_503","olm500"]
for m in ADAPTIF:
    p = " --pattern" if m in ("dwt_503","olm500") else ""
    satirlar.append((f"results/A__{m}__never",    f"mats/{m}.mtx md_light 10000 150 42 results/A__{m}__never{p} --tassa never"))
    satirlar.append((f"results/A__{m}__adaptive", f"mats/{m}.mtx md_light 10000 150 42 results/A__{m}__adaptive{p} --tassa adaptive"))
    satirlar.append((f"results/A__{m}__always",   f"mats/{m}.mtx md_light 2000 150 42 results/A__{m}__always{p} --tassa always"))

DOGUM = ["generated_500_494_1","generated_1000_995_5","er_n500_k3_1","er_n1000_k3_1",
         "block_diagonal_n1000_p0p3","band_n1000_k3","dwt_503","olm500"]
EGRI = {"generated_1000_995_5","er_n1000_k3_1","block_diagonal_n1000_p0p3","dwt_503"}
for m in DOGUM:
    p = " --pattern" if m in ("dwt_503","olm500") else ""
    e = " --egri" if m in EGRI else ""
    for v in V:
        satirlar.append((f"results/D__{m}__{v}", f"mats/{m}.mtx {v} 2000 90 42 results/D__{m}__{v}{p} --dogum{e}"))

for n in (500,1000):
    for k in (3,4,5):
        for c in range(1,6):
            m = f"er_n{n}_k{k}_{c}"
            for v in V:
                satirlar.append((f"results/E__{m}__{v}", f"mats/{m}.mtx {v} 10000 60 42 results/E__{m}__{v}"))

GERCEK = ["ibm32","bcspwr01","bcspwr02","curtis54","dwt_59","chesapeake","mycielskian6","polbooks",
 "nos4","olm100","tub100","ck104","breasttissue_10NN","pivtol","grid1","can_256","sphere3",
 "dwt_503","Trefethen_500","olm500","tomography","nos3","nos5","nos6","nos7","gr_30_30",
 "jagmesh1","dwt_918","dwt_992","lshp_577","lshp_778","685_bus","662_bus"]
for m in GERCEK:
    for v in V:
        satirlar.append((f"results/G__{m}__{v}", f"mats/{m}.mtx {v} 10000 60 42 results/G__{m}__{v} --pattern"))

MALIYET = ["band_n1000_k3","block_diagonal_n1000_p0p3","generated_500_494_1","dwt_503"]
for m in MALIYET:
    p = " --pattern" if m == "dwt_503" else ""
    satirlar.append((f"results/C__{m}__every",    f"mats/{m}.mtx every 500 150 42 results/C__{m}__every{p}"))
    satirlar.append((f"results/C__{m}__md_every", f"mats/{m}.mtx md_every 500 150 42 results/C__{m}__md_every{p}"))

with open("manifest.tsv","w") as f:
    for pref, cmd in satirlar:
        f.write(pref + "\t" + cmd + "\n")
print("toplam kosu:", len(satirlar))
