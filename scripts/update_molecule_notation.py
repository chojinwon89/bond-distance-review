#!/usr/bin/env python3
"""Expose exact formula/name aliases in the comparison cards, selector and search."""
import argparse
import csv
import html
import json
from pathlib import Path
import re
from bs4 import BeautifulSoup
from molecule_names import ALIASES, GROUPS, FORMULAS, canonical

ROOT=Path(__file__).resolve().parents[1]
CARD=re.compile(r'(<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">)(.*?)(?=<div class="g" data-surf=|<script>)',re.S)


def update(root=ROOT,project=Path('/pscratch/sd/j/jcho5/VASP')):
    root=Path(root);page=(root/'dft_comparison.html').read_text()
    tables=re.findall(r'<table\b.*?</table>',page,re.S)
    molecules={c['data-mol'] for c in BeautifulSoup(page,'html.parser').select('.g')}
    def label(name):
        formula=FORMULAS.get(name,name)
        return name.replace('_',' ') + (' ('+formula+')' if formula!=name else '')
    def card(match):
        prefix,s,m,body=match.groups();s,m=html.unescape(s),html.unescape(m)
        body=re.sub(r'<h3>.*?(?=<span class="tag")', '<h3>'+html.escape(s)+' &middot; '+html.escape(label(m))+' ',body,count=1)
        return prefix+body
    page=CARD.sub(card,page)
    def options(match):
        return re.sub(r'(<option value="([^"]+)">).*?(</option>)',
                      lambda m:m[1]+html.escape(label(html.unescape(m[2])))+m[3],match[0])
    page=re.sub(r'<select id="fmol">.*?</select>',options,page,flags=re.S)
    search={m:' '.join(dict.fromkeys([m,FORMULAS.get(m,m),*GROUPS.get(m,())])) for m in sorted(molecules)}
    page=re.sub(r'  // molecule-aliases-start.*?  // molecule-aliases-end\n','',page,flags=re.S)
    block='  // molecule-aliases-start\n  var moleculeAliases='+json.dumps(search,ensure_ascii=True)+';\n'
    targets={a.replace('_',' ').lower():c for a,c in ALIASES.items()}
    block+='  var moleculeQueryTargets='+json.dumps(targets,ensure_ascii=True)+';\n'
    block+='  function normalizeMoleculeQuery(value){return value.replace(/[₀-₉]/g,function(c){return String("₀₁₂₃₄₅₆₇₈₉".indexOf(c));}).replace(/_/g," ").trim().toLowerCase();}\n  // molecule-aliases-end\n'
    page=page.replace('  function apply(){',block+'  function apply(){')
    page=page.replace('q=fq.value.trim().toLowerCase()', 'q=normalizeMoleculeQuery(fq.value)')
    page=page.replace("(cs+' '+cm).toLowerCase().indexOf(q)","normalizeMoleculeQuery(cs+' '+cm+' '+(moleculeAliases[cm]||'')).indexOf(q)")
    fallback="normalizeMoleculeQuery(cs+' '+cm+' '+(moleculeAliases[cm]||'')).indexOf(q)>-1"
    page=page.replace('(!q||'+fallback+')','(!q||(moleculeQueryTargets[q]?cm===moleculeQueryTargets[q]:'+fallback+'))')
    page=re.sub(r'<!-- molecule-notation-note -->.*?<!-- /molecule-notation-note -->\n?','',page,flags=re.S)
    note='''<!-- molecule-notation-note -->
<div class="note info"><b>Molecule notation:</b> DME = CH3OCH3 (dimethyl ether); ethanol = CH3CH2OH = C2H5OH;
methanol = CH3OH; acetaldehyde = CH3CHO; acetic acid = CH3COOH; formaldehyde = H2CO; formic acid = HCOOH.
Names, condensed formulas and subscript digits work in the search box. Ethanol and DME remain distinct despite both having formula C2H6O;
formate and formic acid, and propanol and isopropanol, also remain distinct.
<br>Converged binding energies are displayed even outside the ±5 eV review range, with amber review markers.
Missing or failed references and slab-size mismatches are separate from unusual energy magnitudes.
<a href="molecule_aliases.csv" download>Notation equivalences</a> &middot;
<a href="molecule_name_audit.csv" download>Directory-to-page name audit</a> &middot;
<a href="dme_energy_status.csv" download>DME energy and reference audit</a>.</div>
<!-- /molecule-notation-note -->
'''
    page=page.replace('<h2>Per-system structure',note+'<h2>Per-system structure')
    assert re.findall(r'<table\b.*?</table>',page,re.S)==tables
    (root/'dft_comparison.html').write_text(page)
    with (root/'molecule_aliases.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['alias','canonical_molecule','display_formula']);w.writeheader()
        w.writerows(dict(alias=a,canonical_molecule=c,display_formula=FORMULAS.get(c,c)) for a,c in sorted(ALIASES.items()))
    audit=[]
    for tree in ['dft_jobs','vasp_mol']:
        for d in sorted((project/tree).iterdir()):
            if not d.is_dir():continue
            token=d.name
            if tree=='dft_jobs':
                match=re.fullmatch(r'(.+)_([A-Z][a-z]?\d+)(?:_.+)?',token)
                if not match:continue
                token=match[1]
            name=canonical(token)
            audit.append(dict(tree=tree,directory=str(d),molecule_token=token,canonical_molecule=name,appears_on_page=name in molecules))
    with (root/'molecule_name_audit.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    dme=[]
    for kind,file in [('relaxed','dft_perlmutter_energy_audit.csv'),('SPE','dft_perlmutter_singlepoint_audit.csv')]:
        with (root/file).open(newline='') as f: records=list(csv.DictReader(f))
        for r in records:
            if canonical(r['molecule'])!='DME':continue
            dme.append(dict(calculation=kind,surface=r['surface'],functional=r['functional'],system=r['system'],status=r['status'],
                raw_binding_energy_eV=r.get('E_ads_raw',r.get('E_ads_SPE','')),convergence=r.get('convergence',r.get('complex_status','')),
                complex_source=r.get('complex_outcar',r.get('complex_directory','')),slab_source=r.get('slab_outcar',r.get('slab_directory','')),
                molecule_source=r.get('molecule_outcar',r.get('molecule_directory','')),note=r['note']))
    with (root/'dme_energy_status.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dme[0]));w.writeheader();w.writerows(dme)
    print('Notation:',len(molecules),'page molecules;',len(ALIASES),'explicit aliases;',len(audit),'directory mappings')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--project-root',type=Path,default=Path('/pscratch/sd/j/jcho5/VASP'))
    args=p.parse_args();update(project=args.project_root)
