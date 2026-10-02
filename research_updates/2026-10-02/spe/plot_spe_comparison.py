"""Reproduce the live website's MLIP versus DFT SPE pairs without audit filtering."""
from pathlib import Path
import csv
import json
import math
import re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from bs4 import BeautifulSoup
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = Path(__file__).resolve().parent
FUNCS = [('pbe', 'PBE'), ('pbe_d3', 'PBE+D3'), ('r2scan', 'r²SCAN'), ('beef_vdw', 'BEEF-vdW')]
COLORS = {'Ag':'#2374AB', 'Au':'#B68B12', 'Cu':'#DC6B26', 'Ir':'#8062AB',
          'Pd':'#299477', 'Pt':'#C14F79', 'Rh':'#657584'}
rows = list(csv.DictReader((ROOT/'sources/dft_comparison_singlepoint.csv').open()))
soup = BeautifulSoup((ROOT/'sources/dft_comparison.html').read_text(), 'html.parser')
cards = {(g['data-surf'],g['data-mol']):g for g in soup.select('div.g[data-surf][data-mol]')}
pairs, missing = [], []
checked = 0
for r in rows:
    if not r['E_ads_ML_SP'] or not r['E_ads_DFT_SP']:
        missing.append(r)
        continue
    x,y = float(r['E_ads_DFT_SP']),float(r['E_ads_ML_SP'])
    assert math.isfinite(x) and math.isfinite(y)
    assert abs((y-x)-float(r['delta_SP'])) < 2e-7
    card = cards[r['surface'],r['molecule']]
    name = dict(FUNCS)[r['functional']]
    tr = next(tr for tr in card.select('table.mini tr') if tr.select_one('td.l') and tr.select_one('td.l').get_text(strip=True)==name)
    cells = tr.select('td')
    assert abs(float(cells[4].get_text(strip=True))-x) <= 0.000501
    assert abs(float(cells[2].get_text(strip=True))-y) <= 0.000501
    checked += 1
    pairs.append(dict(r,metal=re.match('[A-Za-z]+',r['surface'])[0],
                      molecule_label='DME (CH3OCH3)' if r['molecule']=='DME' else r['molecule'].replace('_',' '),
                      x_dft_spe_eV=x,y_mlip_eV=y))
assert len(pairs)==len(set((r['surface'],r['molecule'],r['functional']) for r in pairs))
for filename, data in [('plotted_pairs.csv',pairs),('missing_ml_pairs.csv',missing)]:
    with (ROOT/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)

stats=[]
for f,label in FUNCS:
    rr=[r for r in pairs if r['functional']==f]
    diff=np.array([r['y_mlip_eV']-r['x_dft_spe_eV'] for r in rr])
    stats.append(dict(functional=f,n=len(rr),SPE_review=sum(r['SPE_review']=='true' for r in rr),
                      outside_zoom=sum(not(-5<=r['x_dft_spe_eV']<=1 and -5<=r['y_mlip_eV']<=1) for r in rr),
                      all_published_MAE_eV=float(np.mean(abs(diff))),
                      all_published_RMSE_eV=float(np.sqrt(np.mean(diff**2)))))
(ROOT/'validation.json').write_text(json.dumps(dict(html_pairs_checked=checked,missing_ml=len(missing),panels=stats),indent=2))

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':16,
                     'axes.labelsize':12,'svg.fonttype':'none','pdf.fonttype':42})
for zoom in [False,True]:
    fig,axs=plt.subplots(2,2,figsize=(13.5,11.8))
    fig.subplots_adjust(left=.09,right=.975,top=.83,bottom=.19,wspace=.23,hspace=.36)
    fig.suptitle('GOAD + SevenNet-Omni vs DFT single-point energies',x=.09,y=.974,ha='left',fontsize=20,weight='bold')
    fig.text(.09,.937,'Adsorption energies from the published comparison table • 2 October 2026',fontsize=12,color='#4C5969')
    mode='Zoom: −5 to +1 eV; counts identify points outside the view' if zoom else 'Full range: every available pair is shown; x-axis ranges differ by panel'
    fig.text(.09,.910,mode,fontsize=11,color='#4C5969')
    legend=[Line2D([],[],marker='o',ls='',color=c,label=m,markersize=7) for m,c in COLORS.items()]
    legend += [Line2D([],[],marker='o',ls='',color='#555555',mfc='none',label='SPE review flag',markersize=7),Line2D([],[],ls='--',color='#4E5966',label='ML = DFT')]
    fig.legend(handles=legend,loc='upper left',bbox_to_anchor=(.08,.89),ncol=9,frameon=False,columnspacing=1.15,handlelength=1.5)
    for idx,((f,label),ax,st) in enumerate(zip(FUNCS,axs.flat,stats)):
        rr=[r for r in pairs if r['functional']==f]
        for metal,color in COLORS.items():
            for flagged in [False,True]:
                sub=[r for r in rr if r['metal']==metal and (r['SPE_review']=='true')==flagged]
                ax.scatter([r['x_dft_spe_eV'] for r in sub],[r['y_mlip_eV'] for r in sub],s=30 if flagged else 26,
                           facecolors='none' if flagged else color,edgecolors=color if flagged else 'white',
                           linewidths=.85 if flagged else .3,alpha=.8 if flagged else .72,zorder=3)
        xmin=-5 if zoom else min(-5,math.floor(min(r['x_dft_spe_eV'] for r in rr)-.3))
        ax.set(xlim=(xmin,1),ylim=(-5,1),xlabel=r'DFT SPE adsorption energy, $E_{\mathrm{ads}}^{\mathrm{DFT,SP}}$ (eV)',
               ylabel=r'MLIP adsorption energy, $E_{\mathrm{ads}}^{\mathrm{ML}}$ (eV)')
        ax.plot([-5,1],[-5,1],ls='--',color='#4E5966',lw=1.2,zorder=2)
        ax.axhline(0,color='#BFC7D0',lw=.7);ax.axvline(0,color='#BFC7D0',lw=.7)
        ax.grid(color='#E6EBF0',lw=.65,zorder=0)
        for spine in ax.spines.values():spine.set_color('#BCC5CF')
        ax.set_title(f'{chr(97+idx)}   {label}',loc='left',weight='bold',pad=13)
        detail=f"n = {st['n']}  |  SPE flagged: {st['SPE_review']}"
        if zoom:detail+=f"  |  outside view: {st['outside_zoom']}"
        ax.text(.025,.025,detail,transform=ax.transAxes,fontsize=9.2,color='#3E4B58',
                bbox=dict(fc='white',ec='none',alpha=.9,pad=4))
        # Highlight the same requested molecular example in each panel.
        target=next((r for r in rr if r['surface']=='Cu111' and r['molecule']=='DME'),None)
        if target:
            ax.annotate('Cu(111) · DME',xy=(target['x_dft_spe_eV'],target['y_mlip_eV']),xytext=(-95,53),textcoords='offset points',
                        fontsize=9,color='#A64E1C',arrowprops=dict(arrowstyle='-',color='#A64E1C',lw=.85),
                        bbox=dict(boxstyle='round,pad=.25',fc='white',ec='none',alpha=.9),zorder=5)
        if not zoom and f=='beef_vdw':
            extreme=min(rr,key=lambda r:r['x_dft_spe_eV'])
            ax.annotate('Pt(100) · NO\n−42.800 eV',xy=(extreme['x_dft_spe_eV'],extreme['y_mlip_eV']),
                        xytext=(18,-65),textcoords='offset points',fontsize=9,color='#657080',
                        arrowprops=dict(arrowstyle='-',color='#657080'),bbox=dict(fc='white',ec='none',alpha=.9))
    fig.text(.09,.04,'Filled: no SPE review flag in the source table. Hollow: source SPE review flag retained.\n'
             'All finite pairs are retained; 74 entries lack an MLIP value. Flags are not a new validation of the calculations.\n'
             'SPE uses the MLIP geometry with slab/molecule references; it is not the DFT-relaxed adsorption energy.',
             fontsize=10,color='#4C5969',linespacing=1.65)
    stem='mlip_vs_dft_spe_zoom' if zoom else 'mlip_vs_dft_spe_full_range'
    for ext in ['png','pdf','svg']:fig.savefig(ROOT/f'{stem}.{ext}',dpi=220,facecolor='white')
    plt.close(fig)

# Self-contained interactive companion; hover identifies every molecule and surface.
fig=make_subplots(rows=2,cols=2,subplot_titles=[label for _,label in FUNCS],vertical_spacing=.16,horizontal_spacing=.12)
for i,(f,label) in enumerate(FUNCS):
    row,col=i//2+1,i%2+1
    for metal,color in COLORS.items():
        rr=[r for r in pairs if r['functional']==f and r['metal']==metal]
        fig.add_trace(go.Scatter(x=[r['x_dft_spe_eV'] for r in rr],y=[r['y_mlip_eV'] for r in rr],mode='markers',
            name=metal,legendgroup=metal,showlegend=i==0,
            marker=dict(color=color,size=7,opacity=.8,symbol=['circle-open' if r['SPE_review']=='true' else 'circle' for r in rr]),
            customdata=[[r['surface'],r['molecule_label'],r['SPE_review'],r['provenance'],r['complex_directory']] for r in rr],
            hovertemplate='%{customdata[0]} · %{customdata[1]}<br>DFT SPE: %{x:.6f} eV<br>MLIP: %{y:.6f} eV<br>SPE review: %{customdata[2]}<br>%{customdata[3]}<br>%{customdata[4]}<extra></extra>'),row=row,col=col)
    fig.add_trace(go.Scatter(x=[-5,1],y=[-5,1],mode='lines',line=dict(color='gray',dash='dash'),showlegend=False,hoverinfo='skip'),row=row,col=col)
    fig.update_xaxes(title_text='DFT SPE adsorption energy (eV)',range=[-5,1],row=row,col=col)
    fig.update_yaxes(title_text='MLIP adsorption energy (eV)',range=[-5,1],row=row,col=col)
full={};zoom={}
for i,(f,_) in enumerate(FUNCS):
    k='xaxis'+(str(i+1) if i else '')
    full[k+'.range']=[min(-5,math.floor(min(r['x_dft_spe_eV'] for r in pairs if r['functional']==f)-.3)),1]
    zoom[k+'.range']=[-5,1]
fig.update_layout(title='GOAD + SevenNet-Omni vs DFT SPE — hover for molecule and surface',height=950,width=1300,
                  template='plotly_white',margin=dict(t=140,b=65),legend=dict(orientation='h',y=1.12),
                  updatemenus=[dict(type='buttons',direction='right',x=.55,y=1.16,buttons=[
                      dict(label='Zoom −5 to +1 eV',method='relayout',args=[zoom]),dict(label='Full published range',method='relayout',args=[full])])])
html=fig.to_html(include_plotlyjs=True,full_html=True)
html=html.replace('<body>','<body><p style="font:15px sans-serif;margin:20px">Source: published comparison table, 2 October 2026. Hollow points retain SPE review flags. All 1,450 finite pairs are included; 108 lie outside the initial zoom. DME = CH3OCH3. Missing MLIP values are not filled from a different dataset.</p>')
(ROOT/'mlip_vs_dft_spe_interactive.html').write_text(html)
print(json.dumps(stats,indent=2))
print('Validated all',checked,'plotted pairs against the live HTML table.')
