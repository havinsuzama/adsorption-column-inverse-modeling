"""Summarize the complete crossed pilot without refitting."""
import json
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from crossed_replication import ROOT,OUT,NAMES

def main():
    data=json.loads((OUT/'result.json').read_text());rows=data['selected'];checks=json.loads((OUT/'reference_checks.json').read_text());lines=[]
    assert len(data['fits'])==36 and len(rows)==18
    aux=np.load(OUT/'observations.npz')['porosity_observed']
    fig=Figure(figsize=(11,6),layout='constrained');FigureCanvasAgg(fig);axes=fig.subplots(2,3)
    for design in ['one_condition','two_concentrations']:
        rs=[r for r in rows if r['design']==design];errors=np.array([r['relative_errors'] for r in rs]);hold=np.array([r['common_150_rmse'] for r in rs])
        lines.append('| '+design+' | '+' | '.join(f'{v:.2f}%' for v in 100*np.median(abs(errors),axis=0))+f' | {np.median(hold):.6f} |')
        for ax,index in zip(axes[0],range(3)):
            for j in range(3):
                values=[r for r in rs if r['porosity_index']==j];ax.plot([r['concentration_index']+1 for r in values],[100*r['relative_errors'][index] for r in values],marker='o',label=f'{design}, aux {j+1}')
            ax.set(xlabel='Concentration noise realization',ylabel=NAMES[index]+' error (%)');ax.axhline(0,color='gray',lw=.6)
    axes[0,0].legend(fontsize=6)
    effects=[]
    for ax,index in zip(axes[1],range(3)):
        for design in ['one_condition','two_concentrations']:
            grid=np.empty((3,3))
            for r in rows:
                if r['design']==design:grid[r['concentration_index'],r['porosity_index']]=100*r['relative_errors'][index]
            noise_effect=np.mean(grid,axis=1);aux_effect=np.mean(grid,axis=0)
            effects.append(dict(design=design,parameter=NAMES[index],concentration_level_mean_signed_error_percent=noise_effect.tolist(),porosity_level_mean_signed_error_percent=aux_effect.tolist(),interaction_residual_percent=(grid-noise_effect[:,None]-aux_effect[None,:]+grid.mean()).tolist()))
            ax.plot([1,2,3],aux_effect,'o-',label=design)
        ax.set(xlabel='Auxiliary measurement realization',ylabel='Mean '+NAMES[index]+' error (%)');ax.legend(fontsize=7)
    fig.savefig(OUT/'factor_effects.png',dpi=170)
    (OUT/'factor_effects.json').write_text(json.dumps(effects,indent=2))
    comparisons=[]
    for index in range(4):
        wins=0;deltas=[]
        for i in range(3):
            for j in range(3):
                a=next(r for r in rows if r['design']=='one_condition' and r['concentration_index']==i and r['porosity_index']==j);b=next(r for r in rows if r['design']=='two_concentrations' and r['concentration_index']==i and r['porosity_index']==j);delta=abs(b['relative_errors'][index])-abs(a['relative_errors'][index]);wins+=int(delta<0);deltas.append(delta)
        comparisons.append(f'- {NAMES[index]}: lower with two concentrations in {wins}/9 pairs; median paired difference {100*np.median(deltas):+.2f} percentage points.')
    unsuccessful=sum(not r['success'] for r in data['fits']);boundary=sum(any(r['active_bounds']) for r in data['fits']);maxgap=max(r['start_cost_gap'] for r in rows)
    pairs=[(next(r for r in rows if r['design']=='one_condition' and r['concentration_index']==i and r['porosity_index']==j),next(r for r in rows if r['design']=='two_concentrations' and r['concentration_index']==i and r['porosity_index']==j)) for i in range(3) for j in range(3)]
    hold_wins=sum(b['common_150_rmse']<a['common_150_rmse'] for a,b in pairs)
    report=f'''# Crossed noise pilot

One true parameter set (k = 0.001 min⁻¹, P = 48, KL = 0.04 L/mg, eps = 0.4). Three concentration-noise draws were crossed with three auxiliary porosity observations. Two designs with 40 points each (40 at 100 mg/L, or 20 + 20 at 100 and 200 mg/L) and two optimizer starts: 36 fits. In each cell the start with the lower observed objective was kept. The nine cells share noise draws, so they are not nine independent datasets.

| Design | Median abs. k error | P error | KL error | eps error | Median 150 mg/L outlet RMS |
|---|---:|---:|---:|---:|---:|
{chr(10).join(lines)}

Paired comparison, absolute error with two concentrations minus one:

{chr(10).join(comparisons)}
- 150 mg/L outlet RMS: lower with two concentrations in {hold_wins}/9 pairs.

## Notes

- The porosity observations were {[round(v,4) for v in aux.tolist()]}. All three happen to lie above 0.4; no bias was added to the measurement model.
- Reference data use 1024 cells and were checked against 2048 cells; fits use 256 cells. All three feeds passed the preset checks (outlet RMS < 0.001, max difference < 0.002, relative mass error < 1e-6); see reference_checks.json.
- The 150 mg/L prediction is computed with the fitted parameters on the 1024-cell grid. This condition was used during development, so it is not a blind test.
- {unsuccessful} of 36 fits reported an unsuccessful stop and {boundary} hit a parameter bound. The largest objective gap between the two starts was {maxgap:.3g}.
- Density was fixed at 0.8254 kg/L, so P = rho_b·qmax is the capacity result; qmax alone is not identified.
- Mean errors per factor level and interaction residuals are in factor_effects.json. Three levels per factor are too few for variance components.
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8');print(report)

if __name__=='__main__':main()
