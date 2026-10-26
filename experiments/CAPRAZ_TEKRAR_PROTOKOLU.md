# Derişim ve yardımcı ölçüm hataları: sonraki çapraz pilot

Amaç: tek gürültü örneğinde görülen iki-C0 kazanımını, derişim ve yardımcı porozite gözleminin rastgele hatalarını birbirine karıştırmadan sınamak. Bu plan sonuçlar görülmeden yazıldı; nihai makale deneyi veya kör test kaydı değildir.

- Tasarımlar: 40 nokta C0=100; 20+20 nokta C0=100/200 mg/L. Önceki sabit zaman aralığı, debi ve model varsayımları korunur. Eşit gözlem sayısı, eşit toplam deney veya CPU maliyeti değildir.
- Derişim için üç yeni standart Gaussian vektör: seed 20261101/20261102/20261103. Her tasarımda sigma=.02 C/C0; tasarımlar aynı standart vektörü kullanır, fiziksel zaman/koşul eşleşmeleri aynı değildir.
- Porozite için üç ayrı gözlem: seed 20261111/20261112/20261113; log(eps_obs)=log(.4)+.05*z. Üç gözlem her derişim veri örneğiyle çaprazlanır. Ana pilotta sistematik yanlılık sıfırdır; önceki yanlılık taraması ayrı deney olarak kalır. Seed sayıları yalnız rastgele üretim anahtarıdır, veri toplama tarihleri değildir.
- 3×3=9 kombinasyon/tasarım; iki sabit optimizasyon başlangıcıyla toplam 36 fit planlanır. Dokuz kombinasyon dokuz bağımsız derişim verisi değildir: yalnız üç derişim ve üç yardımcı gözlem gerçekleşmesi vardır. Hata kaynaklarının etkisi faktörlerine göre ayrı raporlanır.
- İlk aşamada yoğunluk gözlemi önceki sabit değerde tutulur; P=rho*qmax ana kapasite çıktısıdır. Ayrı qmax sonucu bu sabit yoğunluk bilgisine koşulludur. Yoğunluk ölçümünün rastgele hatası bu pilotta üçüncü bir faktör olarak karıştırılmaz.
- Parametre sınırları ve gözlem modeli önceki fitlerle aynı; en çok 30 amaç değerlendirmesi/başlangıç. İki başlangıç arasında seçim yalnız gözlem amacına göre; doğru parametreye göre seçim yapılmaz. Bütün başarısız/durmamış ve sınıra dayanan fitler raporlanır.
- Referans verisi için 1024 hücre ve seçilmiş 2048 hücre kontrolü, fit için 256 hücre hedeflenir. Bu ön referans kontrolü maliyetli veya yetersiz çıkarsa deney çalıştırılmadan protokol güncellenip gerekçesi kaydedilir. Önceki 512 hücreyle aynı-ızgara fit yaklaşımı nihai bağımsız doğrulama diye kullanılmaz.
- Ana çıktılar: P/KL/k/eps işaretli ve mutlak parametre hataları, ortak 150 mg/L koşulunda çıkış tahmini, aktif sınırlar, durma/optimality, başlangıç farkı ve hesap maliyeti. Üç faktör düzeyinden kalibre edilmiş güven aralığı çıkarılmaz.

Tek doğru parametre noktası pilotudur. Nihai çalışma için daha fazla bağımsız veri/başlangıç, farklı doğru parametreler ve yakın literatürle açık katkı ayrımı gerekir. PINN karşılaştırması ancak kendi ileri çözüm doğrulaması sağlanırsa güvenilir biçimde eklenebilir. Yeni bir PINN kabul edilmiş değildir.
