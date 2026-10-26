"""Local linear uncertainty diagnostic at four saved 512-cell fits; not coverage."""
import hashlib,json,time
import numpy as np
from refine_selected_regimes import ROOT,OUT as SOURCE,NAMES,DESIGNS
from fast_column_sensitivity import sensitive
from parameter_regimes import save

OUT=ROOT/'runs_local_regime_uncertainty'/'pilot_v1'

def geometry(J):
    _,s,Vt=np.linalg.svd(J,full_matrices=False)
    if s[-1]<=np.finfo(float).eps*max(J.shape)*s[0]:return dict(singular_values=s.tolist(),numerically_full_rank=False)
    cov=(Vt.T*(1/s**2))@Vt;std=np.sqrt(np.diag(cov));corr=cov/np.outer(std,std)
    return dict(singular_values=s.tolist(),numerically_full_rank=True,condition=float(s[0]/s[-1]),log_parameter_std=std.tolist(),correlation=corr.tolist(),covariance=cov.tolist())

def main():
    OUT.mkdir(parents=True,exist_ok=True);data=json.loads((SOURCE/'result.json').read_text())
    provenance={n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['local_regime_uncertainty.py','fast_column_sensitivity.py','parameter_regimes.py','refine_selected_regimes.py','physics_model.py']};provenance['source_result']=hashlib.sha256((SOURCE/'result.json').read_bytes()).hexdigest()
    if (OUT/'provenance.json').exists():assert json.loads((OUT/'provenance.json').read_text())==provenance
    else:save(OUT/'provenance.json',provenance)
    save(OUT/'protocol.json',dict(source='runs_regime_grid/pilot_v1',fit_cells=512,sigma_c=.02,sigma_log_porosity=.05,parameter_order=NAMES,scope='Observed-design local Gauss-Newton geometry at saved fitted parameters; fixed known noise scales. Not new fits, nonlinear confidence regions, coverage, blind validation or global identifiability. Without-auxiliary geometry uses SAME fitted state, not reestimated parameters. Density fixed; no separate qmax uncertainty.'))
    rows=json.loads((OUT/'partial.json').read_text()) if (OUT/'partial.json').exists() else []
    for r in data['results']:
        if any(a['key']==r['key'] for a in rows):continue
        tic=time.perf_counter();x=np.log([r['estimate'][n] for n in NAMES]);J=np.vstack([sensitive(x,feed,np.linspace(.005,1,count),512)[1] for feed,count in DESIGNS[r['design']]])/.02
        row=dict(key=r['key'],case=r['case'],design=r['design'],estimate=r['estimate'],concentration_only_same_state=geometry(J),with_porosity=geometry(np.vstack([J,[0,0,0,20]])),seconds=time.perf_counter()-tic);rows.append(row);save(OUT/'partial.json',rows);print('SAVED',r['key'],row['with_porosity']['log_parameter_std'],flush=True)
    lines=[]
    for r in rows:
        g=r['with_porosity'];lines.append('| '+r['case']+' | '+r['design']+' | '+' | '.join(f'{v:.5f}' for v in g.get('log_parameter_std',[float('nan')]*4))+' |')
    report='''# Seçilmiş 512 hücre fitlerinde yerel belirsizlik tanısı

Dört kayıtlı fitte normalize çıkış duyarlılıkları hesaplandı. x=log(k,P,KL,eps), J gözlem standart sapmasına bölünmüş Jacobian; yardımcı porozite gözlemi satırı [0,0,0,20]. Tam rank durumunda yerel kovaryans (J^T J)^-1 olarak SVD ile hesaplandı. Gürültü ölçekleri bilinen sentetik değerlerde tutuldu; artık maliyetinden yeniden ölçeklenmedi. Bu Gauss–Newton yaklaşımıdır, doğrusal olmayan kesin Hessian değildir.

| Vaka | Tasarım | log(k) std | log(P) std | log(KL) std | log(eps) std |
|---|---|---:|---:|---:|---:|
'''+ '\n'.join(lines)+'''

Değerler yerel lineer modelde log-parametre standart sapmasıdır; gerçek parametre hatası, kalibre edilmiş güven aralığı veya başarısızlık olasılığı değildir. Küçük değerlerde yaklaşık bağıl belirsizlik olarak okunabilir; büyük değerler için bu yaklaşık yorum güvenilmez. Sayısal rank tek başına pratik tanımlanabilirlik kanıtı değildir. Sınırlar, doğrusal olmayanlık, model/ızgara hatası ve yardımcı ölçüm yanlılığı bu kovaryans formülünde kapsanmaz.

result.json tekil değerler, koşul sayıları ve korelasyon matrislerini içerir. Yardımcı satır çıkarıldığında hesaplanan geometri aynı kestirim noktasındadır; yardımcı ölçüm olmadan yeniden fit edilmiş sonuç değildir. Tek bir gürültü/yardımcı gözlem eşleşmesi kullanıldı; gürültü tekrarları veya kapsama deneyi yapılmadı. İki tasarım farklı kestirim noktalarında değerlendirildiğinden bu tablo tek başına saf tasarım etkisini ayırmaz. Tasarımın gerçek parametredeki yerel bilgi pilotu ayrıca mevcuttur.

Sonraki kullanım: korelasyonları ve büyük yerel belirsizlikleri aday profil yönlerini seçmek için kullanmak. Güven aralığı ilan etmeden önce doğrusal olmayan profil ve yeterli bağımsız tekrar gerekir. Yeni PINN modeli kabul edilmedi.
'''
    save(OUT/'result.json',dict(results=rows));(OUT/'RAPOR.md').write_text(report,encoding='utf-8');(ROOT/'YEREL_REJIM_BELIRSIZLIGI.md').write_text(report,encoding='utf-8');print('DONE',OUT,flush=True)

if __name__=='__main__':main()
