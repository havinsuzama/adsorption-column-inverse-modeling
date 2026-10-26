"""Post-pilot grid diagnostic: four warm-start refits, not new replicates."""
import hashlib,json,time
import numpy as np
from scipy.optimize import least_squares
from parameter_regimes import ROOT,NAMES,DESIGNS,column,save
from fast_column_sensitivity import sensitive
from physics_model import solve

SOURCE=ROOT/'runs_parameter_regimes'/'pilot_v1'
OUT=ROOT/'runs_regime_grid'/'pilot_v1'

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    files=['refine_selected_regimes.py','parameter_regimes.py','fast_column_sensitivity.py','physics_model.py']
    provenance={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in files}
    for case in ['fast_kinetics','high_affinity']:
        for name in ['result.json','observations.npz','reference.npz']:provenance[f'{case}/{name}']=hashlib.sha256((SOURCE/case/name).read_bytes()).hexdigest()
    if (OUT/'provenance.json').exists():assert json.loads((OUT/'provenance.json').read_text())==provenance
    else:save(OUT/'provenance.json',provenance)
    protocol=dict(cases=['fast_kinetics','high_affinity'],concentration_index=0,porosity_index=0,source_cells=256,fit_cells=512,data_cells=1024,max_nfev=15,initialization='Observed-objective selected 256-cell fit, warm start',scope='Post-pilot selected-case discretization diagnostic, not independent noise/start replication. Selection stated before these refits. Same bounds, objective and stopping tolerances; no acceptance threshold changed.')
    save(OUT/'protocol.json',protocol)
    rows=json.loads((OUT/'fits.json').read_text()) if (OUT/'fits.json').exists() else []
    for case in protocol['cases']:
        old=json.loads((SOURCE/case/'result.json').read_text());obs=np.load(SOURCE/case/'observations.npz');ref=np.load(SOURCE/case/'reference.npz')
        for design,conditions in DESIGNS.items():
            key=case+'_'+design
            if any(r['key']==key for r in rows):continue
            before=next(r for r in old['selected'] if r['design']==design and r['concentration_index']==0 and r['porosity_index']==0);start=np.log([before['estimate'][n] for n in NAMES]);lastx=None;cache=None;tic=time.perf_counter()
            def evaluate(x):
                nonlocal lastx,cache
                if lastx is None or not np.array_equal(lastx,x):
                    pairs=[sensitive(x,feed,np.linspace(.005,1,count),512) for feed,count in conditions];pred=np.concatenate([p[0] for p in pairs]);J=np.vstack([p[1] for p in pairs])/.02
                    cache=(np.r_[(pred-obs[design+'_0'])/.02,(x[3]-np.log(obs['porosity_observed'][0]))/.05],np.vstack([J,[0,0,0,20]]));lastx=x.copy()
                return cache
            initial=float(.5*np.sum(evaluate(start)[0]**2));print('REFIT',key,'initial cost',initial,flush=True)
            fit=least_squares(lambda x:evaluate(x)[0],start,jac=lambda x:evaluate(x)[1],bounds=np.log([[1e-4,5,.005,.15],[.01,200,.4,.75]]),max_nfev=15,ftol=1e-8,xtol=1e-8,gtol=1e-6);v=np.exp(fit.x);times=np.linspace(0,1,101)
            curve=solve(column(v,150),1024,times=times,rtol=2e-9,atol=1e-11)['c'][:,-1];target=np.interp(times,ref['times'],ref['150']);rms=float(np.sqrt(np.mean((curve-target)**2)))
            row=dict(key=key,case=case,design=design,estimate=dict(zip(NAMES,v.tolist())),grid256=before['estimate'],relative_parameter_change=(v/np.exp(start)-1).tolist(),relative_errors=(v/np.array(list(old['truth'].values()))-1).tolist(),initial_cost=initial,cost=float(fit.cost),success=bool(fit.success),status=int(fit.status),message=fit.message,nfev=fit.nfev,optimality=float(fit.optimality),active_bounds=fit.active_mask.tolist(),common_150_rmse=rms,grid256_common_150_rmse=before['common_150_rmse'],seconds=time.perf_counter()-tic)
            assert fit.cost<=initial+1e-8;rows.append(row);save(OUT/'fits.json',rows);print('SAVED',key,json.dumps(row),flush=True)
    save(OUT/'result.json',dict(protocol=protocol,results=rows))
    lines=['| '+r['case']+' | '+r['design']+' | '+' | '.join(f'{100*x:+.3f}%' for x in r['relative_parameter_change'])+f" | {r['grid256_common_150_rmse']:.7f} → {r['common_150_rmse']:.7f} |" for r in rows]
    report='''# Seçilmiş rejimlerde 256→512 hücre yeniden kestirimi

İlk derişim gürültüsü ve ilk porozite gözlemi eşleşmesinde, hızlı kinetik ve yüksek afinite vakalarının iki tasarımı 512 hücrede yeniden optimize edildi. Seçim önceki pilot görüldükten sonra tanı amacıyla yapıldı; bu bağımsız tekrar veya bütün rejimlerin yakınsama çalışması değildir. Başlangıçlar önceki 256 hücre kestirimleridir. Veri 1024 hücrede üretilmişti; fit aynı veri üretim ızgarasına taşınmadı. Ortak 150 mg/L tahminleri her iki parametre kümesinde de 1024 hücrede değerlendirildi.

| Vaka | Tasarım | k kayması | P kayması | KL kayması | eps kayması | Ortak RMS: önce → sonra |
|---|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(lines)+'''

Tablodaki yüzdeler gerçek parametre hatası değil, 256→512 parametre değişimidir. Gerçek hatalar, gözlem amaçları, durma/optimality ve süreler result.json içinde saklanır. ftol/xtol ile durma küresel optimum kanıtı değildir. En fazla 15 amaç değerlendirmesi kullanıldı; başarısız durmalar çıkarılmadı. Küçük kaymalar tek başına bütün gürültü örneklerinde tasarım sıralamasının korunmasını veya güven aralığı kalibrasyonunu göstermez. Referans hâlâ aynı fizik/ayrıklaştırma ailesine aittir; deneysel doğrulama yoktur. Yeni PINN modeli kabul edilmedi.
'''
    (OUT/'RAPOR.md').write_text(report,encoding='utf-8');(ROOT/'SECILMIS_REJIM_IZGARA_KONTROLU.md').write_text(report,encoding='utf-8');print('DONE',OUT,flush=True)

if __name__=='__main__':main()
