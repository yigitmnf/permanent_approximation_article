BASLANGIC=$(date +%s)
BUTCE=${1:-215}
TAMAM=0; ATLANDI=0
while IFS=$'\t' read -r pref cmd; do
  if [ -f "$pref.ozet" ]; then ATLANDI=$((ATLANDI+1)); continue; fi
  T0=$(date +%s.%N)
  ./deney_kosucu $cmd > /dev/null 2>&1
  T1=$(date +%s.%N)
  echo -e "$pref\t$(echo "$T1 - $T0" | bc)" >> results/kosu_sureleri.tsv
  TAMAM=$((TAMAM+1))
  GECEN=$(( $(date +%s) - BASLANGIC ))
  if [ $GECEN -gt $BUTCE ]; then echo "BUTCE_DOLDU tamam=$TAMAM atlandi=$ATLANDI"; exit 0; fi
done < makale_manifest.tsv
echo "MANIFEST_BITTI tamam=$TAMAM atlandi=$ATLANDI"
