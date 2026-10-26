# Adsorpsiyon kolonu: parametre kestirimi ve ölçüm tasarımı

Tek boyutlu bir adsorpsiyon kolonu üzerine sentetik bir çalışma. Depoda sonlu hacim referans çözücüsü, taşınım parametrelerinin ters kestirimi, ölçüm tasarımı deneyleri ve henüz deneme aşamasında olan bir PINN var.

## Çalıştırma

```bash
python -m pip install -r requirements.txt
python run_demo.py
python run_demo.py --mode results
```

Spyder'da `run_demo.py` dosyasını açıp çalıştırmak da yeterli. Varsayılan mod ileri çözümü çalıştırır ve çıktıları `generated/` klasörüne yazar. `python run_demo.py --mode inverse` küçük ve kaba bir ters çözüm örneğidir; arşivdeki sonuçları yeniden üretmez.

## Kısaca

- Tüm veriler sentetik. Modelde radyoaktif bozunma, pH ve iyon rekabeti yok.
- Sıvı faz verisi yığın yoğunluğu ile doygunluk kapasitesini ayrı ayrı belirleyemez; kestirilen büyüklük P = yoğunluk × kapasite.
- İki giriş derişimi ortalamada daha iyi sonuç verdi, ama her gürültü örneğinde değil.
- Yerel standart sapmalar güven aralığı olarak okunmamalı.
- PINN kısmı tamamlanmadı.
- Prajitno vd. (2021) makalesindeki parametrelerle çizilen Thomas eğrisi ölçülmüş veri değildir; makalenin verileri depoda yer almıyor.

Ayrıntılar için [README.md](README.md) ve [deney açıklamaları](experiments/README.md).
