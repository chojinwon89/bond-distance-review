"""Refresh results and figures for the matched CPU initialization experiment."""
from pathlib import Path
import csv,json,datetime
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'benchmarks/ag_glycerol_cpu_contact_20261002'
rows=[json.loads(p.read_text()) for p in sorted((OUT/'runs').glob('*/summary.json'))]
status=dict(updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),completed=len(rows),planned=20,
            valid=sum(r['valid'] for r in rows),finished=len(rows)==20)
print(json.dumps(status),flush=True)
(OUT/'status.json').write_text(json.dumps(status,indent=2))
if not rows:raise SystemExit(0)
with (OUT/'all_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
pairs=[]
for seed in sorted(set(r['seed'] for r in rows)):
    rr={r['mode']:r for r in rows if r['seed']==seed}
    if len(rr)<2:continue
    a,b=rr['baseline'],rr['contact25']
    pairs.append(dict(seed=seed,baseline_E_ads_eV=a['final_E_ads_eV'],contact25_E_ads_eV=b['final_E_ads_eV'],
        contact25_minus_baseline_eV=b['final_E_ads_eV']-a['final_E_ads_eV'],baseline_O_Ag_A=a['nearest_O_Ag_A'],
        contact25_O_Ag_A=b['nearest_O_Ag_A'],both_valid=a['valid'] and b['valid']))
if pairs:
    with (OUT/'paired_comparison.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]));w.writeheader();w.writerows(pairs)
fig,axs=plt.subplots(1,2,figsize=(11,5.3))
for ax,key,title in zip(axs,['final_E_ads_eV','nearest_O_Ag_A'],['MLIP adsorption energy (eV)','Final nearest O–Ag distance (Å)']):
    for p in pairs:
        rr={r['mode']:r for r in rows if r['seed']==p['seed']}
        if p['both_valid']:
            ax.plot([0,1],[rr['baseline'][key],rr['contact25'][key]],color='gray',alpha=.4,lw=.8)
    for mode,x,col in [('baseline',0,'#2374AB'),('contact25',1,'#DC6B26')]:
        group=sorted([r for r in rows if r['mode']==mode],key=lambda r:r[key])
        span=max(r[key] for r in rows)-min(r[key] for r in rows)
        gap=max(span*.037,1e-5)
        labels=[]
        for r in group:
            labels.append(max(r[key],labels[-1]+gap) if labels else r[key])
        side=-1 if mode=='baseline' else 1
        for r,label_y in zip(group,labels):
            ax.scatter(x,r[key],c=col,marker='o' if r['valid'] else 'x',s=40,zorder=3)
            ax.annotate(str(r['seed']),(x,r[key]),xytext=(x+side*.13,label_y),
                        textcoords='data',fontsize=7,ha='right' if side<0 else 'left',va='center',
                        arrowprops=dict(arrowstyle='-',color=col,lw=.5,alpha=.6))
    ax.set_xticks([0,1],['Original GOAD','Initial O–Ag = 2.5 Å']);ax.set_ylabel(title)
    ax.set_xlim(-.3,1.4);ax.grid(axis='y',alpha=.2)
fig.suptitle('Ag(111)–glycerol · GOAD + SevenNet-MF · CPU',weight='bold')
fig.text(.1,.025,f"{len(rows)}/20 runs complete. Lines connect matched valid seeds; × marks invalid/unconverged endpoints.\n2.5 Å applies only at initialization; all final O–Ag distances are unconstrained. These are MLIP results, not new DFT results.",fontsize=9)
fig.tight_layout(rect=[0,.11,1,.94])
for ext in ['png','pdf']:fig.savefig(OUT/f'comparison.{ext}',dpi=200)
best={mode:min([r for r in rows if r['mode']==mode and r['valid']],key=lambda r:r['final_E_ads_eV'],default=None) for mode in ['baseline','contact25']}
(OUT/'best_endpoints.json').write_text(json.dumps(best,indent=2))
