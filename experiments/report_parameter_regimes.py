"""Report all four regimes after completion; never fit or discard failed starts."""
import json
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from parameter_regimes import ROOT,OUT,CASES,NAMES

def main():
    assert (OUT/'completion.json').exists(),'Experiment not complete; inspect per-case fits.json for progress.'
    lines=[];details=[];summaries=[];fig=Figure(figsize=(11,7),layout='constrained');FigureCanvasAgg(fig);axes=fig.subplots(2,2)
    for ax,(case,truth) in zip(axes.ravel(),CASES.items()):
        data=json.loads((OUT/case/'result.json').read_text());rows=data['selected'];assert len(data['fits'])==16
        for design in ['one_condition','two_concentrations']:
            selected=[r for r in rows if r['design']==design and 'relative_errors' in r]
            if not selected:
                lines.append(f'| {case} | {design} | no finite fit | — | — | — |');continue
            errors=np.array([r['relative_errors'] for r in selected]);hold=[r['common_150_rmse'] for r in selected]
            lines.append('| '+case+' | '+design+' | '+' | '.join(f'{x:.2f}%' for x in 100*np.median(abs(errors[:,:3]),axis=0))+f' | {np.median(hold):.6f} |')
            ax.plot(range(1,len(hold)+1),hold,'o-',label=design)
        comparisons={}
        for name,index in zip(NAMES,range(4)):
            pairs=[]
            for i in range(2):
                for j in range(2):
                    a=next(r for r in rows if r['design']=='one_condition' and r['concentration_index']==i and r['porosity_index']==j);b=next(r for r in rows if r['design']=='two_concentrations' and r['concentration_index']==i and r['porosity_index']==j)
                    if 'relative_errors' in a and 'relative_errors' in b:pairs.append(abs(b['relative_errors'][index])-abs(a['relative_errors'][index]))
            comparisons[name]=dict(available_pairs=len(pairs),two_concentration_lower_error=sum(v<0 for v in pairs),paired_absolute_error_differences=pairs)
        hold_deltas=[]
        for i in range(2):
            for j in range(2):
                pair=[next(r for r in rows if r['design']==d and r['concentration_index']==i and r['porosity_index']==j) for d in ['one_condition','two_concentrations']]
                if all('common_150_rmse' in r for r in pair):hold_deltas.append(pair[1]['common_150_rmse']-pair[0]['common_150_rmse'])
        failures=sum(not r['success'] for r in data['fits']);bounds=sum(any(r.get('active_bounds',[])) for r in data['fits']);gap=max((r.get('start_cost_gap',0) for r in rows),default=0)
        details.append(dict(case=case,truth=truth,comparisons=comparisons,common_150_paired_differences=hold_deltas,unsuccessful=failures,exception_count=sum('error' in r for r in data['fits']),active_bound_count=bounds,max_start_cost_gap=gap,fit_seconds=sum(r['seconds'] for r in data['fits'])))
        summaries.append(f"- {case}: two concentrations better in {comparisons[NAMES[1]]['two_concentration_lower_error']}/4 pairs for P, {comparisons[NAMES[2]]['two_concentration_lower_error']}/4 for KL and {sum(x<0 for x in hold_deltas)}/{len(hold_deltas)} for the 150 mg/L prediction. Unsuccessful stops {failures}/16, active bounds {bounds}/16, largest start gap {gap:.3g}.")
        ax.set(title=case,xlabel='Concentration / auxiliary realization',ylabel='150 mg/L outlet RMS');ax.set_xticks([1,2,3,4],['C1/A1','C1/A2','C2/A1','C2/A2']);ax.legend(fontsize=8)
    fig.savefig(OUT/'regime_predictions.png',dpi=170)
    (OUT/'comparison.json').write_text(json.dumps(details,indent=2))
    report='''# Four parameter regimes

Each regime changes one parameter from the base case (k = 0.001, P = 48, KL = 0.04, eps = 0.4): slow kinetics k = 0.0005, fast kinetics k = 0.002, low affinity KL = 0.01, high affinity KL = 0.16. Per regime: 2 concentration-noise draws × 2 porosity observations × 2 designs × 2 starts = 16 fits; 64 in total. The same standard noise draws are used across regimes and designs, so the four cells per regime are not independent replicates. Both designs have 40 observations.

| Regime | Design | Median abs. k error | P error | KL error | Median 150 mg/L outlet RMS |
|---|---|---:|---:|---:|---:|
'''+ '\n'.join(lines)+'\n\n'+ '\n'.join(summaries)+'''

## Notes

- In each cell the start with the lower observed objective was kept, never the one closer to the true parameters. All starts, including unsuccessful ones, are in each regime's result.json; paired differences are in comparison.json.
- Reference data use 1024 cells and were checked against 2048 cells for 100, 150 and 200 mg/L before any fit started (reference_checks.json in each regime folder). Fits use 256 cells.
- The 150 mg/L prediction is computed on the 1024-cell grid. This condition was used during development, so it is not a blind test.
- P and porosity are the same in every regime, so their effect is not tested. qmax is not identified separately from bulk density.
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8');print(report)

if __name__=='__main__':main()
