"""Rank-one local information update at fixed fitted states; no refits."""
import hashlib,json
import numpy as np
from pathlib import Path
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

ROOT=Path(__file__).resolve().parent
SOURCE=ROOT/'runs_local_regime_uncertainty'/'pilot_v1'/'result.json'
OUT=ROOT/'runs_auxiliary_precision_geometry'/'pilot_v1'

def main():
    OUT.mkdir(parents=True,exist_ok=True);data=json.loads(SOURCE.read_text());rows=[];maxcheck=0.
    for r in data['results']:
        C=np.array(r['with_porosity']['covariance']);I=np.linalg.inv(C)
        for sigma in [.01,.02,.05,.10,.20]:
            dp=1/sigma**2-1/.05**2;den=1+dp*C[3,3];assert den>0
            updated=C-dp/den*np.outer(C[:,3],C[3,:]);info=I.copy();info[3,3]+=dp
            direct=np.linalg.solve(info,np.eye(4));check=float(np.linalg.norm(updated-direct)/np.linalg.norm(direct));maxcheck=max(maxcheck,check);assert check<1e-7
            assert np.linalg.eigvalsh(updated).min()>0
            rows.append(dict(key=r['key'],case=r['case'],design=r['design'],sigma_log_porosity=sigma,log_parameter_std=np.sqrt(np.diag(updated)).tolist(),direct_solve_relative_check=check))
    # Stronger independent information cannot increase local marginal variances.
    for key in {r['key'] for r in rows}:
        group=sorted([r for r in rows if r['key']==key],key=lambda r:r['sigma_log_porosity']);assert np.min(np.diff(np.array([r['log_parameter_std'] for r in group]),axis=0))>-1e-9
    (OUT/'result.json').write_text(json.dumps(dict(results=rows,max_direct_solve_relative_check=maxcheck),indent=2))
    (OUT/'provenance.json').write_text(json.dumps(dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
    fig=Figure(figsize=(10,4),layout='constrained');FigureCanvasAgg(fig);axes=fig.subplots(1,2);lines=[]
    for ax,case in zip(axes,['fast_kinetics','high_affinity']):
        for design in ['one_condition','two_concentrations']:
            group=[r for r in rows if r['case']==case and r['design']==design]
            ax.plot([r['sigma_log_porosity'] for r in group],[r['log_parameter_std'][2] for r in group],'o-',label=design)
            for r in group:lines.append('| '+case+' | '+design+f" | {r['sigma_log_porosity']:.2f} | "+' | '.join(f'{v:.5f}' for v in r['log_parameter_std'])+' |')
        ax.set(title=case,xlabel='Assumed log-porosity measurement sigma',ylabel='Local log(KL) std');ax.legend(fontsize=8)
    fig.savefig(OUT/'precision_effect.png',dpi=170)
    report='''# Yardımcı porozite ölçüm hassasiyeti: sabit noktada yerel bilgi hesabı

Önceki dört 512 hücre kestirimindeki yerel kovaryanslar kullanıldı. Her kestirim noktası ve derişim Jacobian'ı sabit tutulup bağımsız porozite gözleminin bilgi katkısı değiştirildi: I_new=I_base+(1/sigma_new^2-1/.05^2)*e_eps*e_eps^T. Kovaryans rank-bir güncellemesiyle hesaplandı ve doğrudan matris çözümüyle kontrol edildi. Daha hassas yardımcı bilgi için yerel marjinal varyansların artmadığı doğrulandı. Yeni ODE çözümü, gözlem üretimi veya parametre optimizasyonu yapılmadı.

| Vaka | Tasarım | Varsayılan log-eps sigma | log(k) std | log(P) std | log(KL) std | log(eps) std |
|---|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(lines)+'''

Bu tablo aynı kestirim noktasındaki bilgi ağırlığının etkisini ayırır. Ölçüm hassasiyeti değiştiğinde gözlem değerinin ve yeniden kestirilen parametrelerin nasıl değişeceğini simüle etmez; sistematik yanlılığı kapsamaz. Tasarımlar farklı kestirim noktalarında olduğundan tasarımlar arası fark yalnız bilgi bütçesine atfedilemez. Log-parametre std değerleri gerçek yüzde hata veya kalibre edilmiş güven aralığı değildir. Özellikle büyük değerlerde yaklaşık yüzde yorumundan kaçınılmalı. Önceki doğrusal olmayan KL profili yerel kuadratik yaklaşımdan sapma gösterdiği için bu tablodan nihai güven aralığı türetilmemeli.

Bu ucuz tanı, hangi hassasiyet aralığının sonraki doğrusal olmayan deney için değer taşıdığını seçmeye yardımcı olur; tek başına genel optimum ölçüm hassasiyeti belirlemez. Yoğunluk sabit tutuldu, qmax ayrı belirsizliği hesaplanmadı. Veriler sentetik; yeni PINN kabulü ve deneysel doğrulama yoktur.
'''
    (OUT/'RAPOR.md').write_text(report,encoding='utf-8');(ROOT/'POROZITE_HASSASIYETI_YEREL_BILGI.md').write_text(report,encoding='utf-8');print(report)

if __name__=='__main__':main()
