import random
import sys

def uret(n, k, tohum, yol):
    rng = random.Random(tohum)
    girdiler = set((i, i) for i in range(1, n + 1))
    for _ in range(n * k):
        i = rng.randint(1, n)
        j = rng.randint(1, n)
        girdiler.add((i, j))
    sirali = sorted(girdiler)
    with open(yol, "w") as f:
        f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
        f.write(f"{n} {n} {len(sirali)}\n")
        for i, j in sirali:
            f.write(f"{i} {j} 1\n")
    return len(sirali)

if __name__ == "__main__":
    n = int(sys.argv[1])
    k = int(sys.argv[2])
    adet = int(sys.argv[3])
    hedef_klasor = sys.argv[4]
    for kopya in range(1, adet + 1):
        tohum = 1000 * n + 100 * k + kopya
        yol = f"{hedef_klasor}/er_n{n}_k{k}_{kopya}.mtx"
        nnz = uret(n, k, tohum, yol)
        print(yol, "nnz =", nnz)
