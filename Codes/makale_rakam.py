import csv
import glob
import math
import os
import statistics

def known():
    t = {}
    for line in open("known_log10.tsv"):
        p = line.rstrip("\n").split("\t")
        t[p[0]] = float(p[3])
    return t

def read_ozet(path):
    d = {}
    for line in open(path):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            d[p[0]] = p[1]
    return d

def relerr(est_log10, true_log10):
    if est_log10 is None:
        return 1.0
    diff = est_log10 - true_log10
    if diff > 30:
        return 10.0 ** 30
    return abs(10.0 ** diff - 1.0)

def rows():
    with open("results/toplu.tsv") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def A(rs, mat, mode, field):
    for r in rs:
        if r["kume"] == "A" and r["matris"] == mat and r["tassa_modu"] == mode:
            return float(r[field])
    return None

def K(rs, mat, var, field):
    for r in rs:
        if r["kume"] == "K" and r["matris"] == mat and r["varyant"] == var:
            return float(r[field])
    return None

def C(rs, mat, var, field):
    for r in rs:
        if r["kume"] == "C" and r["matris"] == mat and r["varyant"] == var:
            return float(r[field])
    return None

def group_stats(rs):
    g = {}
    for r in rs:
        if r["kume"] not in ("K", "E", "G"):
            continue
        m = r["matris"]
        if m.startswith("er_"):
            grp = "Er"
        elif m.startswith("band"):
            grp = "Band"
        elif m.startswith("block"):
            grp = "Blk"
        elif m.startswith("generated"):
            grp = "Gen"
        else:
            grp = "Real"
        g.setdefault((grp, r["varyant"]), []).append(float(r["pm_orani"]))
    return g

def rep_stats(prefix, mode):
    kt = known()
    pms, errs = [], []
    for path in sorted(glob.glob(f"results/P__{prefix}__rep_{mode}_s*.ozet")) + \
                [f"results/A__{prefix}__{mode}.ozet"]:
        if not os.path.exists(path):
            continue
        d = read_ozet(path)
        pms.append(float(d["tam_eslesme_orani"]))
        mat = prefix
        if mat in kt:
            est = d.get("tahmin_log10", "YOK")
            errs.append(relerr(None if est == "YOK" else float(est), kt[mat]))
    return pms, errs

def eq_stats(mat, budget_tag):
    kt = known()
    out = {}
    for mode in ("never", "adaptive", "always"):
        recs = []
        for path in sorted(glob.glob(f"results/P__{mat}__{budget_tag}_{mode}_s*.ozet")):
            d = read_ozet(path)
            est = d.get("tahmin_log10", "YOK")
            recs.append({"samples": int(d["ornek"]),
                         "pm": float(d["tam_eslesme_orani"]),
                         "err": relerr(None if est == "YOK" else float(est), kt[mat]),
                         "t": float(d["sure_toplam_sn"])})
        out[mode] = recs
    return out

def fmt(x, nd=1):
    return f"{x:.{nd}f}"

def main():
    rs = rows()
    kt = known()
    m = {}

    hard = [("generated_500_494_1", "GenA"), ("generated_1000_995_5", "GenB"),
            ("er_n500_k3_1", "ErA"), ("er_n500_k5_1", "ErC"), ("er_n1000_k3_1", "ErB"),
            ("block_diagonal_n1000_p0p1", "BlkA"), ("dwt_503", "DwtA"), ("olm500", "OlmA")]
    for mat, tag in hard:
        for mode, mtag in (("never", "Nev"), ("adaptive", "Ada"), ("always", "Alw")):
            pm = A(rs, mat, mode, "pm_orani")
            ms = A(rs, mat, mode, "ornek_ms")
            if pm is not None:
                m[f"pm{tag}{mtag}"] = fmt(100 * pm, 1)
                m[f"ms{tag}{mtag}"] = fmt(ms, 2)
    for mat, tag in (("generated_500_494_1", "GenA"), ("generated_1000_995_5", "GenB")):
        for mode, mtag in (("never", "Nev"), ("adaptive", "Ada"), ("always", "Alw")):
            r = [x for x in rs if x["kume"] == "A" and x["matris"] == mat and x["tassa_modu"] == mode]
            if r and r[0]["tahmin_log10"]:
                m[f"err{tag}{mtag}"] = fmt(relerr(float(r[0]["tahmin_log10"]), kt[mat]), 3)
        r = [x for x in rs if x["kume"] == "A" and x["matris"] == mat and x["tassa_modu"] == "adaptive"][0]
        m[f"share{tag}Ada"] = fmt(100 * float(r["sure_tassa"]) / float(r["sure_toplam"]), 0)

    m["ratioAdaAlwGenA"] = fmt(float(m["msGenAAda"]) / float(m["msGenAAlw"]) * 100, 0)
    m["ratioAdaAlwGenB"] = fmt(float(m["msGenBAda"]) / float(m["msGenBAlw"]) * 100, 0)

    thr = {}
    for mat, tag in (("generated_500_494_1", "GenA"), ("generated_1000_995_5", "GenB"),
                     ("er_n1000_k3_1", "ErB")):
        for mode, mtag in (("never", "Nev"), ("adaptive", "Ada"), ("always", "Alw")):
            pm = A(rs, mat, mode, "pm_orani")
            ms = A(rs, mat, mode, "ornek_ms")
            thr[f"{tag}{mtag}"] = 1000.0 * pm / ms
            m[f"thr{tag}{mtag}"] = fmt(1000.0 * pm / ms, 0)

    g = group_stats(rs)
    for (grp, var), lst in g.items():
        vtag = {"ras": "Ras", "once": "Once", "light": "Light", "md": "Md",
                "md_once": "MdOnce", "md_light": "MdLight"}[var]
        m[f"pm{grp}{vtag}Mean"] = fmt(100 * sum(lst) / len(lst), 1)
        m[f"pm{grp}{vtag}Med"] = fmt(100 * statistics.median(lst), 1)
        m[f"n{grp}{vtag}Hi"] = str(sum(1 for p in lst if p >= 0.99))

    for mat, var, tag in (("band_n1000_k2", "md_light", "BandBkTwo"),
                          ("band_n1000_k3", "md_light", "BandBkThree"),
                          ("block_diagonal_n1000_p0p1", "md_light", "BlkBpOne"),
                          ("block_diagonal_n500_p0p9", "md_light", "BlkFpNine"),
                          ("generated_500_494_1", "md_light", "GenFive"),
                          ("generated_1000_995_5", "md_light", "GenTen")):
        m[f"err{tag}"] = fmt(K(rs, mat, var, "mutlak_bagil_hata"), 3)

    m["msEveryBand"] = fmt(C(rs, "band_n1000_k3", "md_every", "ornek_ms"), 2)
    m["msLightBand"] = fmt(K(rs, "band_n1000_k3", "md_light", "ornek_ms"), 2)
    m["msOnceBand"] = fmt(K(rs, "band_n1000_k3", "md_once", "ornek_ms"), 2)
    m["msMdBand"] = fmt(K(rs, "band_n1000_k3", "md", "ornek_ms"), 2)
    m["slowEvery"] = fmt(float(m["msEveryBand"]) / float(m["msLightBand"]), 0)

    lagsum = {}
    for mat in ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1",
                "block_diagonal_n1000_p0p3", "dwt_503", "olm500", "band_n1000_k3"]:
        for var in ("ras", "once", "light", "md", "md_once", "md_light"):
            path = f"results/D__{mat}__{var}.dogum"
            if not os.path.exists(path):
                continue
            recs = [l.split() for l in open(path)]
            if not recs:
                continue
            n = int(recs[0][2])
            lag = sum((int(b) - int(a)) / n for a, b, _, _ in recs) / len(recs)
            birth = sum(int(a) / n for a, b, _, _ in recs) / len(recs)
            lagsum.setdefault(var, {"lag": [], "birth": []})
            lagsum[var]["lag"].append(lag)
            lagsum[var]["birth"].append(birth)
    for var, tag in (("ras", "Ras"), ("light", "Light"), ("md", "Md"), ("md_light", "MdLight")):
        m[f"lag{tag}"] = fmt(sum(lagsum[var]["lag"]) / len(lagsum[var]["lag"]), 2)
        m[f"birth{tag}"] = fmt(sum(lagsum[var]["birth"]) / len(lagsum[var]["birth"]), 2)

    frac = {}
    for var in ("once", "light"):
        vals = []
        for mat in ["generated_500_494_1", "generated_1000_995_5", "er_n500_k3_1", "er_n1000_k3_1"]:
            path = f"results/D__{mat}__{var}.dogum"
            if not os.path.exists(path):
                continue
            vals += [float(l.split()[3]) for l in open(path) if float(l.split()[3]) >= 0]
        frac[var] = sum(1 for v in vals if v <= 0.25) / len(vals)
    m["fracLowOnce"] = fmt(100 * frac["once"], 0)
    m["fracLowLight"] = fmt(100 * frac["light"], 0)

    for mat, tag in (("generated_500_494_1", "GenA"), ("er_n500_k3_1", "ErA")):
        for mode, mtag in (("never", "Nev"), ("adaptive", "Ada")):
            pms, errs = rep_stats(mat, mode)
            m[f"repN{tag}{mtag}"] = str(len(pms))
            m[f"repPm{tag}{mtag}"] = fmt(100 * statistics.mean(pms), 1)
            m[f"repPmSd{tag}{mtag}"] = fmt(100 * statistics.stdev(pms), 1)
            if errs:
                m[f"repErr{tag}{mtag}"] = fmt(statistics.mean(errs), 3)
                m[f"repErrSd{tag}{mtag}"] = fmt(statistics.stdev(errs), 3)

    eqA = eq_stats("generated_500_494_1", "eq30")
    for mode, mtag in (("never", "Nev"), ("adaptive", "Ada"), ("always", "Alw")):
        recs = eqA[mode]
        m[f"eqA{mtag}Samp"] = f"{int(statistics.mean([r['samples'] for r in recs])):,}".replace(",", "\\,")
        m[f"eqA{mtag}Pm"] = fmt(100 * statistics.mean([r["pm"] for r in recs]), 1)
        m[f"eqA{mtag}Err"] = fmt(statistics.mean([r["err"] for r in recs]), 3)
        m[f"eqA{mtag}ErrSd"] = fmt(statistics.stdev([r["err"] for r in recs]), 3)
    eqB = eq_stats("generated_1000_995_5", "eq90")
    for mode, mtag in (("never", "Nev"), ("adaptive", "Ada"), ("always", "Alw")):
        recs = eqB[mode]
        m[f"eqB{mtag}N"] = str(len(recs))
        m[f"eqB{mtag}Samp"] = f"{int(statistics.mean([r['samples'] for r in recs])):,}".replace(",", "\\,")
        m[f"eqB{mtag}Pm"] = fmt(100 * statistics.mean([r["pm"] for r in recs]), 1)
        m[f"eqB{mtag}Err"] = fmt(statistics.mean([r["err"] for r in recs]), 3)
        if len(recs) > 1:
            m[f"eqB{mtag}ErrSd"] = fmt(statistics.stdev([r["err"] for r in recs]), 3)
            m[f"eqB{mtag}ErrMax"] = fmt(max(r["err"] for r in recs), 3)
            m[f"eqB{mtag}ErrMin"] = fmt(min(r["err"] for r in recs), 3)

    for mat, tag in (("band_n1000_k3_weighted", "WBandBk"), ("band_n500_k3_weighted", "WBandFh"),
                     ("block_diagonal_n1000_p0p3_weighted", "WBlkBk"),
                     ("block_diagonal_n500_p0p3_weighted", "WBlkFh")):
        d = read_ozet(f"results/P__{mat}__w_md_light.ozet")
        est = d.get("tahmin_log10", "YOK")
        if mat in kt:
            m[f"err{tag}"] = fmt(relerr(None if est == "YOK" else float(est), kt[mat]), 3)
        else:
            m[f"err{tag}"] = "---"
        m[f"pm{tag}"] = fmt(100 * float(d["tam_eslesme_orani"]), 1)

    tassa_cut = []
    for r in rs:
        if r["kume"] == "G" and r["varyant"] == "ras":
            pass
    for path in glob.glob("results/G__*__ras.ozet"):
        d = read_ozet(path)
        h, t = d.get("nnz_ham"), d.get("nnz_tassa_sonrasi")
        if h and t and int(h) > 0:
            tassa_cut.append(1.0 - int(t) / int(h))
    m["staticCutMax"] = fmt(100 * max(tassa_cut), 0)
    m["staticCutMed"] = fmt(100 * statistics.median(tassa_cut), 1)

    total = 0.0
    nrun = 0
    for line in open("results/kosu_sureleri.tsv"):
        total += float(line.split("\t")[1])
        nrun += 1
    m["totalRuns"] = str(nrun)
    m["totalMinutes"] = fmt(total / 60.0, 0)

    with open("makale/numbers.tex", "w") as f:
        for k2, v2 in sorted(m.items()):
            f.write(f"\\newcommand{{\\{k2}}}{{{v2}}}\n")
    print("numbers.tex:", len(m), "macros")

if __name__ == "__main__":
    main()
