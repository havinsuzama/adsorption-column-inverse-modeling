"""Report saved sparse profile diagnostic, including unsuccessful attempts."""
import json
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from high_affinity_profile import ROOT,OUT

def main():
    data=json.loads((OUT/'result.json').read_text());assert len(data['fits'])==14
    fig=Figure(figsize=(10,4),layout='constrained');FigureCanvasAgg(fig);axes=fig.subplots(1,2);lines=[]
    for ax,design in zip(axes,['one_condition','two_concentrations']):
        rows=sorted([r for r in data['points'] if r['design']==design],key=lambda r:r['actual_local_z'])
        ax.plot([r['actual_local_z'] for r in rows],[r['twice_delta_cost'] for r in rows],'o-',label='Reoptimized nuisance parameters')
        ax.plot([r['actual_local_z'] for r in rows],[r['local_quadratic_prediction'] for r in rows],'s--',label='Local quadratic prediction')
        ax.set(title=design,xlabel='Log KL shift / local std',ylabel='2 × observed objective increase');ax.legend(fontsize=8)
        for r in rows:
            lines.append(f"| {design} | {r['z_requested']:+d} | {r['target_KL']:.6f} | {r['actual_local_z']:+.4f} | {r['twice_delta_cost']:.5f} | {r['local_quadratic_prediction']:.5f} | {r['available_starts']} | {r['start_cost_gap']:.3g} |")
    fig.savefig(OUT/'profile_comparison.png',dpi=170)
    failed=sum(not r['success'] for r in data['fits']);exceptions=sum('error' in r for r in data['fits']);negative=sum(r['twice_delta_cost'] < -1e-6 for r in data['points'])
    report=f'''# Sparse KL profile at high affinity

High-affinity regime, first concentration-noise draw and first porosity observation, 512-cell fits. KL was fixed at the fitted value plus 0, ±1 and ±2 local standard deviations (clipped to the original bounds), and k, P and eps were re-optimized at each point. Five points per design, with a second start at the two ends: 14 fits. Points and starts were written to protocol.json before fitting.

With δ = log(KL) − log(KL_hat) and s the local log(KL) standard deviation, a quadratic objective would rise by (δ/s)². The table compares this with the observed rise 2ΔJ relative to the re-optimized centre.

| Design | Requested z | Fixed KL | Actual z | 2ΔJ observed | Quadratic prediction | Starts | Start gap |
|---|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(lines)}

## Notes

- In both designs the objective rises faster than the quadratic prediction toward low KL and more slowly toward high KL, so the local standard deviation is not a symmetric interval here.
- When a point hits the KL upper bound (0.4), its actual z is below the requested value and the quadratic prediction uses the actual z.
- Points with an objective below the centre (beyond 1e-6): {negative}. Unsuccessful stops: {failed}/14; numerical exceptions: {exceptions}.
- Five points per design show the shape of the profile; they are not a confidence interval.
'''
    (OUT/'REPORT.md').write_text(report,encoding='utf-8');print(report)

if __name__=='__main__':main()
