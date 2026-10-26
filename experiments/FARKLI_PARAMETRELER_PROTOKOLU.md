# Sonraki sınama: farklı doğru parametreler

22 Ekim 2026'da, bu yeni vakaların sonuçları hesaplanmadan belirlenen geliştirme protokolü. Nihai kör test veya istatistiksel güç analizi değildir. Önceki tek-nokta çapraz pilot ayrı tutulacak.

## Soru ve vakalar

İki giriş derişiminin eşit gözlem sayısındaki kazanımı, kinetik veya afinite değiştiğinde korunuyor mu? Önceki doğru k=.001, P=48, KL=.04, eps=.4 noktasıyla karşılaştırılmak üzere dört yeni vaka seçildi:

| Vaka | k (min^-1) | P | KL (L/mg) | eps |
|---|---:|---:|---:|---:|
| Yavaş kinetik | .0005 | 48 | .04 | .4 |
| Hızlı kinetik | .002 | 48 | .04 | .4 |
| Düşük afinite | .001 | 48 | .01 | .4 |
| Yüksek afinite | .001 | 48 | .16 | .4 |

Her vakada yalnız bir fiziksel parametre değişir. P, porozite, dispersiyon veya debinin genel etkisi bu taramada sınanmış olmayacak. Bunlar literatürden belirli bir izotopa kalibre edilmiş değerler değildir.

## Sabit tutulacak kararlar

- Tek C0=100 için 40; C0=100/200 için 20+20 gözlem. theta=.005..1, Q=.157 mL/min, T=24000 dakika. Gürültü sigma=.02 C/C0. Önceki model ve sınır koşulları korunur.
- İki yeni derişim gürültüsü seed 20261121/20261122; iki bağımsız porozite ölçümü seed 20261131/20261132, log-sigma=.05. Seed'ler tarih değildir. Standart hatalar vakalar ve tasarımlar arasında eşlenir; fiziksel ölçümlerin aynı olduğu iddia edilmez. Sonuçlara veya hata işaretine bakarak seed değiştirilmez.
- Her vaka için 2×2 çapraz gürültü/yardımcı ölçüm, iki tasarım ve önceki iki sabit başlangıç: 16 fit/vaka, toplam 64 yeni fit. İki faktör düzeyi çok sınırlı bir dayanıklılık taramasıdır; güven aralığı veya kapsama çalışması değildir.
- Başlangıçlar [.0007,38.4,.024,.32] ve [.0015,57.6,.072,.5] olarak kalır; gizli gerçek değerlere göre yakınlaştırılmaz. Parametre sınırları, log-parametreleme, ftol=xtol=1e-8, gtol=1e-6 ve en fazla 30 amaç değerlendirmesi korunur. Başlangıç seçimi yalnız gözlem amacıyla yapılır.
- Sabit yoğunluk bilgisi .8254263138594846; P ana kapasite çıktısıdır. Yoğunluk gürültüsü yeni bir faktör olarak eklenmez.
- Her yeni vakada 100/150/200 mg/L için 1024 hücre referans ve 2048 kontrolü yapılır. Önceki örnek-zaman birleşiminde outlet RMS<.001, maksimum fark<.002 ve bağıl kütle hatası<1e-6 şartları sağlanmadan fit başlamaz. Kontrol başarısızsa vaka çıkarılmaz; sayısal çözünürlük kararı ve gerekçesi sonuçlardan önce kaydedilir. Fit 256 hücrede.
- Her fit hemen kaydedilir; başarısız, sınıra dayanan veya değerlendirme bütçesi dolan koşular saklanır. Bunlar sonradan daha kolay vakalarla değiştirilmez.

## Raporlama ve karar

Her vaka ayrı: işaretli/mutlak k,P,KL,eps hatası; gözlem maliyeti; iki başlangıç farkı; durum/optimality ve süre. Ortak 150 mg/L tahmini seçilen parametrelerle 1024 hücrede hesaplanır. Bu, daha önce araştırma sırasında kullanılan geliştirme koşuludur, nihai kör test değildir. İki tasarımın eşlenik farkları raporlanır; dört hücre bağımsız dört veri tekrarı sayılmaz. Vakalar tek medyanda birleştirilerek başarısız rejimler gizlenmez.

Kazanım kaybolursa veya tersine dönerse araştırma sonucu hangi rejimde işe yaradığını açıklamaya daraltılır; yalnız iyi vakalar seçilmez. Başlangıçlar anlaşmazsa amaç ve parametre farkları ayrı gösterilir; iki başlangıcın anlaşması küresel optimum kanıtı değildir. Bu tarama bittikten sonra nihai tekrar/belirsizlik bütçesi belirlenir. Makale özgünlüğü en yakın çalışmalarla ayrıca doğrulanmalıdır.

Durum: protokol hazır; 64 yeni fit henüz başlatılmadı. Önce bu dört vaka için referans kontrollerinin uygulanması gerekiyor.
