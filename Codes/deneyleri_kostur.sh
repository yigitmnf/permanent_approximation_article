mkdir -p results
R=results
K=./deney_kosucu
SEED=42

BILINEN="band_n100_k2 band_n100_k3 band_n500_k3 band_n1000_k2 band_n1000_k3 block_diagonal_n500_p0p1 block_diagonal_n500_p0p3 block_diagonal_n500_p0p5 block_diagonal_n500_p0p9 block_diagonal_n1000_p0p1 block_diagonal_n1000_p0p3 block_diagonal_n1000_p0p5 block_diagonal_n1000_p0p9 generated_500_494_1 generated_500_499_4 generated_1000_989_4 generated_1000_995_5"
VARYANTLAR="ras once light md md_once md_light"

for m in $BILINEN; do
  for v in $VARYANTLAR; do
    $K mats/$m.mtx $v 10000 120 $SEED $R/K__${m}__${v} > /dev/null 2>&1
  done
  echo "K $m bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
done
echo "FAZ1_BITTI" >> $R/ilerleme.log

ADAPTIF="generated_500_494_1 generated_1000_995_5 er_n500_k3_1 er_n500_k5_1 er_n1000_k3_1 block_diagonal_n1000_p0p1 dwt_503 olm500"
for m in $ADAPTIF; do
  P=""
  case $m in dwt_503|olm500) P="--pattern";; esac
  $K mats/$m.mtx md_light 10000 150 $SEED $R/A__${m}__never $P --tassa never > /dev/null 2>&1
  $K mats/$m.mtx md_light 10000 150 $SEED $R/A__${m}__adaptive $P --tassa adaptive > /dev/null 2>&1
  $K mats/$m.mtx md_light 2000 150 $SEED $R/A__${m}__always $P --tassa always > /dev/null 2>&1
  echo "A $m bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
done
echo "FAZ2_BITTI" >> $R/ilerleme.log

DOGUM="generated_500_494_1 generated_1000_995_5 er_n500_k3_1 er_n1000_k3_1 block_diagonal_n1000_p0p3 band_n1000_k3 dwt_503 olm500"
EGRISET="generated_1000_995_5 er_n1000_k3_1 block_diagonal_n1000_p0p3 dwt_503"
for m in $DOGUM; do
  P=""
  case $m in dwt_503|olm500) P="--pattern";; esac
  E=""
  for em in $EGRISET; do if [ "$m" = "$em" ]; then E="--egri"; fi; done
  for v in $VARYANTLAR; do
    $K mats/$m.mtx $v 2000 90 $SEED $R/D__${m}__${v} $P --dogum $E > /dev/null 2>&1
  done
  echo "D $m bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
done
echo "FAZ3_BITTI" >> $R/ilerleme.log

for n in 500 1000; do
  for k in 3 4 5; do
    for c in 1 2 3 4 5; do
      m=er_n${n}_k${k}_${c}
      for v in $VARYANTLAR; do
        $K mats/$m.mtx $v 10000 60 $SEED $R/E__${m}__${v} > /dev/null 2>&1
      done
    done
    echo "E n$n k$k bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
  done
done
echo "FAZ4_BITTI" >> $R/ilerleme.log

GERCEK="ibm32 bcspwr01 bcspwr02 curtis54 dwt_59 chesapeake mycielskian6 polbooks nos4 olm100 tub100 ck104 breasttissue_10NN pivtol grid1 can_256 sphere3 dwt_503 Trefethen_500 olm500 tomography nos3 nos5 nos6 nos7 gr_30_30 jagmesh1 dwt_918 dwt_992 lshp_577 lshp_778 685_bus 662_bus"
for m in $GERCEK; do
  for v in $VARYANTLAR; do
    $K mats/$m.mtx $v 10000 60 $SEED $R/G__${m}__${v} --pattern > /dev/null 2>&1
  done
  echo "G $m bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
done
echo "FAZ5_BITTI" >> $R/ilerleme.log

MALIYET="band_n1000_k3 block_diagonal_n1000_p0p3 generated_500_494_1 dwt_503"
for m in $MALIYET; do
  P=""
  case $m in dwt_503) P="--pattern";; esac
  $K mats/$m.mtx every 500 150 $SEED $R/C__${m}__every $P > /dev/null 2>&1
  $K mats/$m.mtx md_every 500 150 $SEED $R/C__${m}__md_every $P > /dev/null 2>&1
  echo "C $m bitti $(date +%H:%M:%S)" >> $R/ilerleme.log
done
echo "HEPSI_BITTI" >> $R/ilerleme.log
