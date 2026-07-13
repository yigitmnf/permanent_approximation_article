# permanent_approximation_article

# Seyrek Permanent Tahmini — Deney Araçları

MD sıralı + Sinkhorn ölçekli Rasmussen SIS örnekleyicisi, dinamik/adaptif TASSA budaması ve makaledeki tüm deney/analiz boru hattı. Tek çekirdek, C++17, harici bağımlılık yok (Python tarafı: matplotlib, numpy).

## Envanter

| Dosya | Görev |
|---|---|
| `deney_kosucu.cpp` | Çalıştırıcı v1: 7 varyant (ras/once/light/every/md/md_once/md_light), `--tassa never\|adaptive\|always`, `--pattern`, `--egri`, `--dogum` (doom enstrümantasyonu) |
| `deney_kosucu_v2.cpp` | v1 + `--esik <tau>` (eşikli adaptif tetik) ve `--dsizle` (çifte-stokastik sapma izleme); bayraksız kullanımda v1 ile alan-alan aynı çıktı |
| `erdos_renyi_uretici.py` | ER matris üretici; tohum = `1000n + 100k + kopya` (deterministik) |
| `manifest_uret.py` | 560 satırlık ana deney manifestini üretir (`manifest.tsv`) |
| `surdur.sh` / `msurdur.sh` / `ysurdur.sh` | Kaldığı yerden sürdürülebilir koşucular: `manifest.tsv` / `makale_manifest.tsv` / `yeni_manifest.tsv`; `.ozet` varsa atlar; argüman: saniye bütçesi |
| `deneyleri_kostur.sh` | Tarihsel: manifest sistemi öncesi ilk parti koşu dizisi (yeniden üretim için `surdur.sh` yeterlidir) |
| `fetch.sh` | SuiteSparse gerçek matris indirici (`dl/` altına) |
| `toplu_ozet.py` | Tüm `results/*.ozet` dosyalarını `results/toplu.tsv`'ye birleştirir |
| `grafik_uret.py` | Ana rapor figürleri (8 adet) |
| `rakam_cikar.py` | Rapor makroları → `rapor/rakamlar.tex` |
| `makale_figur.py` | Makale figürleri (vektörel PDF, `figs/`) |
| `makale_rakam.py` | Makale makroları → `numbers.tex` (236 makro) |
| `yeni_grafikler.py` | Ek deney figürleri (τ taraması, eşit-süre, yakınsama, hassasiyet, DS-sapma, ağırlıklı ER) + `rakamlar2.tex` |
| `sabotaj_lab.sh` | Öğretici: temel kodu 3 şekilde bozup yanlılık/yansızlık farkını ölçer |

## Derleme

```
g++ -O2 -std=c++17 -o deney_kosucu    deney_kosucu.cpp
g++ -O2 -std=c++17 -o deney_kosucu_v2 deney_kosucu_v2.cpp
```

## Çalıştırıcı kullanımı

```
./deney_kosucu_v2 <matris.mtx> <varyant> <ornek> <sure_sn> <tohum> <cikti_oneki>
                  [--pattern] [--egri] [--dogum]
                  [--tassa never|adaptive|always] [--esik tau] [--dsizle]
```

Varyantlar: `ras once light every md md_once md_light`.
Çıktılar: `<onek>.ozet` (TSV özet), isteğe bağlı `.egri` (yakınsama), `.dogum` (doom kayıtları), `.dsat` (DS sapması).

Örnek — adaptif budamalı MD+Light, 10.000 örnek:

```
./deney_kosucu_v2 mats/generated_500_494_1.mtx md_light 10000 300 42 sonuc --tassa adaptive
```

## Boru hattı (sıfırdan yeniden üretim)

```
python3 erdos_renyi_uretici.py 500 3 5 mats
python3 erdos_renyi_uretici.py 500 4 5 mats
python3 erdos_renyi_uretici.py 500 5 5 mats
python3 erdos_renyi_uretici.py 1000 3 5 mats
python3 erdos_renyi_uretici.py 1000 4 5 mats
python3 erdos_renyi_uretici.py 1000 5 5 mats
bash fetch.sh
python3 manifest_uret.py
bash surdur.sh 3000
bash msurdur.sh 3000
bash ysurdur.sh 3000
python3 toplu_ozet.py
python3 rakam_cikar.py
python3 grafik_uret.py
python3 makale_rakam.py
python3 makale_figur.py
python3 yeni_grafikler.py
```

Sürdürücüler bütçe dolunca `BUTCE_DOLDU`, iş bitince `MANIFEST_BITTI` basar; aynı komut tekrar verilerek kalan koşular tamamlanır.

## Manifest formatı

Her satır: `cikti_oneki<TAB>komut_argumanlari`. `.ozet` mevcutsa satır atlanır; bu sayede koşular kesintiye dayanıklıdır.

## Veri

Matrisler (`mats/`), bilinen permanentler (`known_log10.tsv`, 60 basamak hassasiyetle türetilmiş log10 değerleri) ve üç manifest, makale artifact paketinde yer alır; ER matrisleri yukarıdaki üreticiyle bire bir yeniden üretilebilir.
