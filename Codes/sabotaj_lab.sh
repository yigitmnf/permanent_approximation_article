set -e
KAYNAK_DIZIN=${1:-.}
CALISMA=sabotaj_calisma
mkdir -p $CALISMA
cp "$KAYNAK_DIZIN/rasmussen_generic_weighted.cpp" "$KAYNAK_DIZIN/dm.cpp" "$KAYNAK_DIZIN/dm.h" $CALISMA/
cd $CALISMA

python3 - << 'PY'
import random, itertools
random.seed(5)
n = 9
while True:
    a = [[1 if random.random() < 0.45 else 0 for _ in range(n)] for _ in range(n)]
    toplam = 0
    for p in itertools.permutations(range(n)):
        v = 1
        for i in range(n):
            v *= a[i][p[i]]
            if v == 0:
                break
        toplam += v
    if toplam > 0:
        break
girdiler = [(i + 1, j + 1) for i in range(n) for j in range(n) if a[i][j]]
with open("t9w.mtx", "w") as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n")
    f.write(f"{n} {n} {len(girdiler)}\n")
    for i, j in girdiler:
        f.write(f"{i} {j} 1\n")
print("gercek permanent =", toplam)
PY

g++ -O2 -std=c++17 -o taban rasmussen_generic_weighted.cpp dm.cpp

cp rasmussen_generic_weighted.cpp s1.cpp
sed -i 's|prod = prod \* (sum_unmatched / cv_sample\[chosen_column\]);|prod = prod * sum_unmatched;|' s1.cpp
g++ -O2 -std=c++17 -o sabotaj1 s1.cpp dm.cpp

cp rasmussen_generic_weighted.cpp s2.cpp
python3 - << 'PY'
s = open("s2.cpp").read()
eski = """\t\tfor(int ptr = A.row_ptr[current_vertex]; ptr < A.row_ptr[current_vertex + 1]; ptr++) {
\t\t\tint c = A.col_idx[ptr];
\t\t\tif(c_is_matched[c] == false) {"""
yeni = eski.replace("if(c_is_matched[c] == false)", "if(true)")
assert s.count(eski) == 1
open("s2.cpp", "w").write(s.replace(eski, yeni, 1))
PY
g++ -O2 -std=c++17 -o sabotaj2 s2.cpp dm.cpp

cp rasmussen_generic_weighted.cpp s3.cpp
python3 - << 'PY'
s = open("s3.cpp").read()
eski = """\t\t\tif(no_unmatched == 0) {
\t\t\t\tprod = 0;
\t\t\t\tbreak; //if everything is matched sample contribution is 0
\t\t\t}

\t\t\tint chosen_column = -1;
\t\t\tRT rand_val = sum_unmatched * dis(gen);
\t\t\tchosen_column = unmatched_c[0]; """
yeni = eski.replace("prod = 0;\n\t\t\t\tbreak; //if everything is matched sample contribution is 0", "break;")
assert s.count(eski) == 1
open("s3.cpp", "w").write(s.replace(eski, yeni, 1))
PY
g++ -O2 -std=c++17 -o sabotaj3 s3.cpp dm.cpp

for ad in taban sabotaj1 sabotaj2 sabotaj3; do
  echo "================ $ad  (gercek per = 115; cikti sirasi: olceksiz, once, her-adim, light)"
  ./$ad t9w.mtx 0.35 1 42 2>/dev/null | grep -m4 "Estimation:" | cut -f1
done
