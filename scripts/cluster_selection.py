#!/usr/bin/env python3
"""Perlmutter-first candidate selection, with geometry-verified Kestrel fallback."""
from collections import Counter,defaultdict
import csv
from decimal import Decimal
import html
import json
from pathlib import Path
import re

from component_store import latest, read_store, reference_candidates
from molecule_names import canonical
from shared_reference_energies import derive, FIELDS
from structure_identity import same_structure

ROOT=Path(__file__).resolve().parents[1]


def choose(perlmutter,kestrel,mode):
    """Never prefer an energy because it is lower or nearer the ML prediction."""
    if perlmutter:
        if perlmutter['status']=='ok':return perlmutter,'Perlmutter first',False
        if kestrel and kestrel['status']=='ok':
            same,proof=same_structure(perlmutter['_complex'],kestrel['_complex'],mode)
            if same:return kestrel,'Kestrel fallback: Perlmutter energy review; same structure verified by '+proof,True
            return perlmutter,'Perlmutter energy review retained; Kestrel structure equivalence unverified',False
        return perlmutter,'Perlmutter energy review retained; no unflagged Kestrel alternative',False
    if kestrel:return kestrel,'Kestrel: no usable Perlmutter binding-energy candidate',False
    return None,'No complete matched component set',False


def write_csv(path,rows,fields):
    with path.open('w',newline='') as out:
        w=csv.DictWriter(out,fieldnames=fields);w.writeheader();w.writerows(rows)


def read_csv(path):
    with path.open(newline='') as source:return list(csv.DictReader(source))


def build(root=ROOT):
    records=latest(read_store(root/'dft_component_store.jsonl.gz'))
    byid={r['snapshot_id']:r for r in records}; grouped=defaultdict(list);observed=defaultdict(list)
    coverage=read_csv(root/'dft_completion_coverage.csv')
    targets={(r['surface'],canonical(r['molecule']),r['functional']) for r in coverage}
    spe_allowed=set()
    for filename in ['dft_perlmutter_singlepoint_audit.csv','dft_kestrel_singlepoint_audit.csv']:
        for r in read_csv(root/filename):
            if r['status'] in ('ok','energy_review','reference_unavailable','potential_mismatch'):spe_allowed.add(r['complex_directory'])
    for comp in records:
        c=comp['calculation'];key=(comp['surface'],comp['molecule'],c['functional'])
        if comp['role']!='complex' or key not in targets:continue
        mode='SPE' if c.get('nsw')==0 else 'relaxed' if (c.get('nsw') or 0)>0 else 'unknown'
        observed[key+(mode,)].append(comp)
        if c['status']!='converged' or mode=='unknown':continue
        if mode=='SPE' and c['directory'] not in spe_allowed:continue
        refs=reference_candidates(comp,records)['candidate_reference_ids']
        def rank(sid):
            r=byid[sid]
            return (r['cluster']!=comp['cluster'],r['provenance'].get('kind')!='live-files',r['calculation']['directory'])
        candidate=None
        for sid in sorted(refs['slab'],key=rank):
            for gid in sorted(refs['molecule'],key=rank):
                try:candidate=derive(comp,byid[sid],byid[gid],mode)
                except ValueError:continue
                candidate['_complex']=comp;break
            if candidate:break
        if candidate:grouped[key+(mode,)].append(candidate)
    selections=[];duplicates=[];candidate_rows=[]
    for entry in coverage:
        for mode in ['relaxed','SPE']:
            key=(entry['surface'],canonical(entry['molecule']),entry['functional'],mode)
            candidates=sorted(grouped.get(key,[]),key=lambda r:(len(Path(r['complex_directory']).parts),r['complex_directory']))
            primary=next((r for r in candidates if r['complex_cluster']=='perlmutter'),None)
            backup=next((r for r in candidates if r['complex_cluster']=='kestrel'),None)
            # Among Kestrel duplicate directories, a verified same-structure match
            # outranks an unverified candidate, without using the energy as a rank.
            if primary:
                match=next((r for r in candidates if r['complex_cluster']=='kestrel' and same_structure(primary['_complex'],r['_complex'],mode)[0]),None)
                if match:backup=match
            selected,reason,fallback=choose(primary,backup,mode)
            if primary and backup:
                same,proof=same_structure(primary['_complex'],backup['_complex'],mode)
                duplicates.append(dict(surface=key[0],molecule=key[1],functional=key[2],mode=mode,
                    E_ads_perlmutter=primary['E_ads'],E_ads_kestrel=backup['E_ads'],
                    complex_total_difference_eV=str(Decimal(str(primary['E_complex']))-Decimal(str(backup['E_complex']))),
                    perlmutter_directory=primary['complex_directory'],kestrel_directory=backup['complex_directory'],
                    same_structure=str(same).lower(),geometry_evidence=proof,selection_reason=reason))
            if selected and not primary:
                local=[r for r in observed.get(key,[]) if r['cluster']=='perlmutter']
                if local and not any(same_structure(r,selected['_complex'],mode)[0] for r in local):
                    selected=None;reason='Kestrel candidate held: matching Perlmutter structure not verified'
            row={field:'' for field in FIELDS}
            row.update(surface=key[0],molecule=key[1],functional=key[2],mode=mode)
            if selected:row.update({k:v for k,v in selected.items() if not k.startswith('_')})
            row.update(selection_reason=reason,kestrel_fallback=str(fallback).lower(),
                previous_E_ads=entry['E_ads_'+mode],candidate_count=len(candidates))
            selections.append(row)
            candidate_rows.extend({k:v for k,v in c.items() if not k.startswith('_')} for c in candidates)
    write_csv(root/'dft_cluster_candidates.csv',candidate_rows,FIELDS)
    write_csv(root/'dft_cluster_selection.csv',selections,FIELDS+['selection_reason','kestrel_fallback','previous_E_ads','candidate_count'])
    duplicate_fields=['surface','molecule','functional','mode','E_ads_perlmutter','E_ads_kestrel','complex_total_difference_eV','perlmutter_directory','kestrel_directory','same_structure','geometry_evidence','selection_reason']
    write_csv(root/'dft_cluster_duplicates.csv',duplicates,duplicate_fields)
    summary=dict(policy='Require matching PAW TITEL identities for all component species. Withhold historical numbers lacking validated references. Perlmutter first. If |E_ads| > 5 eV, prefer an unflagged Kestrel result only with verified matching structure. Directory names, composition, similar energies and matching cells alone do not prove equal coordinates. Preserve nonselected history; do not substitute gas energies across functionals.',
        selections=dict(Counter(r['selection_reason'] for r in selections)),candidate_results=len(candidate_rows),duplicate_pairs=len(duplicates),
        geometry_verified_pairs=sum(r['same_structure']=='true' for r in duplicates),
        normal_existing_to_review=sum(bool(r['previous_E_ads']) and abs(Decimal(r['previous_E_ads']))<=5 and r['status']=='energy_review' for r in selections),
        kestrel_record_provenance=dict(Counter(r['provenance']['kind'] for r in records if r['cluster']=='kestrel')),
        access_note='This build reads the local component store. Archived exports do not constitute a fresh Kestrel filesystem search.')
    (root/'dft_cluster_selection_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))


def apply(root=ROOT):
    from update_kestrel_singlepoint import CARD,FUNCTIONALS,number,combined
    selection_path=root/'dft_cluster_selection.csv'
    if not selection_path.exists():return
    selection_rows=read_csv(selection_path)
    selection={(r['surface'],r['molecule'],r['functional'],r['mode']):r for r in selection_rows}
    potential_audit=root/'dft_potential_mismatches.csv'
    mismatched={(r['surface'],r['molecule'],r['functional'],r['mode']) for r in read_csv(potential_audit)} if potential_audit.exists() else set()
    fresh_audit=root/'dft_perlmutter_singlepoint_audit.csv'
    if potential_audit.exists() and fresh_audit.exists():
        mismatched.update((r['surface'],canonical(r['molecule']),r['functional'],'SPE') for r in read_csv(fresh_audit) if r['status']=='potential_mismatch')
    recovery_path=root/'dft_gas_reference_recovery.json'
    recovery={r['functional']:r for r in json.loads(recovery_path.read_text())['jobs']} if recovery_path.exists() else {}
    path=root/'dft_comparison.html';page=path.read_text();sp_path=root/'dft_comparison_singlepoint.csv'
    previous=read_csv(sp_path);sp_fields=list(previous[0]);sp={(r['surface'],r['molecule'],r['functional']):r for r in previous}
    snapshot=root/'dft_cluster_selection_before.csv'
    if not snapshot.exists():
        from bs4 import BeautifulSoup
        before=[]
        for card_node in BeautifulSoup(page,'html.parser').select('.g'):
            s,m=card_node['data-surf'],canonical(card_node['data-mol'])
            for tr in card_node.select('table tr'):
                cells=tr.select('td')
                if len(cells)!=7:continue
                f=FUNCTIONALS[cells[0].get_text()];saved=sp.get((s,m,f),{})
                for mode,index,delta_index in [('relaxed',1,3),('SPE',4,6)]:
                    before.append(dict(surface=s,molecule=m,functional=f,mode=mode,E_ads=cells[index].get_text(),E_ads_ML=cells[2].get_text(),delta=cells[delta_index].get_text(),
                        provenance=saved.get('provenance','') if mode=='SPE' else cells[index].get('title',''),
                        complex_directory=saved.get('complex_directory','') if mode=='SPE' else ''))
        write_csv(snapshot,before,list(before[0]))
    changes=[]
    def card(match):
        prefix,surface,molecule,body=match.groups();surface,molecule=html.unescape(surface),canonical(html.unescape(molecule))
        def row(match):
            cells=re.findall(r'<td\b[^>]*>.*?</td>',match[1],re.S)
            if len(cells)!=7:return match[0]
            functional=FUNCTIONALS[html.unescape(re.sub(r'<[^>]+>','',cells[0]))];key=(surface,molecule,functional)
            ml=number(cells[2]);selected_spe=None
            for mode,index,delta_index in [('relaxed',1,3),('SPE',4,6)]:
                selected=selection.get(key+(mode,))
                if not selected:continue
                old=number(cells[index])
                selected['applied_E_ads']=str(old) if old is not None else ''
                selected['application_status']='existing value retained' if old is not None else 'missing'
                if not selected['E_ads']:
                    if key+(mode,) in mismatched or (potential_audit.exists() and (old is not None or 'refs unverified' in cells[index])):
                        cls='sp-energy energy-review' if mode=='SPE' else 'energy-review'
                        label='potential mismatch' if key+(mode,) in mismatched else 'refs unverified'
                        cells[index]=f'<td class="{cls}" title="No fully validated replacement with matching PAW potentials is selected; previous numbers remain in the history downloads">{label}</td>'
                        cells[delta_index]='<td class="sp-energy">&mdash;</td>' if mode=='SPE' else '<td>&mdash;</td>'
                        selected['applied_E_ads']='';selected['application_status']='withheld: '+label
                        if mode=='SPE':sp.pop(key,None);cells[5]='<td class="sp-energy sp-combined">&mdash;</td>'
                        continue
                    retry=recovery.get(functional) if molecule=='DME' else None
                    if old is None and retry and retry['scheduler_state'] in ('PENDING','RUNNING'):
                        label='gas ref queued' if retry['scheduler_state']=='PENDING' else 'gas ref running'
                        cls=' class="sp-energy"' if mode=='SPE' else ''
                        title=html.escape('Gas reference recovery '+retry['job_id']+' '+retry['scheduler_state']+' at snapshot; other component statuses remain in the audit. '+retry['directory'],quote=True)
                        cells[index]=f'<td{cls} title="{title}">{label}</td>'
                    elif old is None and retry and 'gas ref' in cells[index]:
                        cls=' class="sp-energy"' if mode=='SPE' else ''
                        cells[index]=f'<td{cls} title="Gas retry status changed; inspect the component audit for remaining blockers">see audit</td>'
                    continue
                value=Decimal(selected['E_ads']);old=number(cells[index]);review=selected['status']=='energy_review'
                if not potential_audit.exists() and old is not None and selected['complex_cluster']=='kestrel' and selected['kestrel_fallback']!='true':
                    selected['application_status']='existing value retained; no verified replacement required'
                    continue
                # Existing unflagged values are held while an unusual replacement
                # lacks a geometry-verified Kestrel comparison. Raw candidates stay downloadable.
                if not potential_audit.exists() and old is not None and abs(old)<=5 and 'energy-review' not in cells[index] and review:
                    selected['application_status']='existing unflagged value retained pending geometry verification'
                    continue
                selected['applied_E_ads']=str(value);selected['application_status']='selected candidate applied'
                cls=('sp-energy ' if mode=='SPE' else '')+('energy-review' if review else '')
                attr=f' class="{cls.strip()}"' if cls.strip() else ''
                title=html.escape(selected['selection_reason']+'; '+selected['complex_directory']+'; slab '+selected['slab_directory']+'; molecule '+selected['molecule_directory'],quote=True)
                cells[index]=f'<td{attr} title="{title}">{value:.3f}</td>'
                cells[delta_index]=f'<td{attr} title="Displayed ML minus {mode} DFT">{ml-value:+.3f}</td>' if ml is not None else ('<td class="sp-energy">&mdash;</td>' if mode=='SPE' else '<td>&mdash;</td>')
                if old is None or abs(value-old)>Decimal('0.0005'):
                    changes.append(dict(surface=surface,molecule=molecule,functional=functional,mode=mode,old_E_ads=str(old) if old is not None else '',new_E_ads=str(value),reason=selected['selection_reason']))
                if mode=='SPE':selected_spe=selected
            relaxed_value=number(cells[1]);spe_value=number(cells[4]);saved=sp.get(key)
            if spe_value is not None:
                saved=dict(saved or {field:'' for field in sp_fields})
                saved.update(surface=surface,molecule=molecule,functional=functional)
                if selected_spe:
                    spe_value=Decimal(selected_spe['E_ads']);saved.update(E_ads_DFT_SP=str(spe_value),E_ads_ML_SP=str(ml) if ml is not None else '',delta_SP=str(ml-spe_value) if ml is not None else '',
                        provenance='Cluster policy: '+selected_spe['selection_reason'],complex_directory=selected_spe['complex_directory'])
                elif saved.get('E_ads_DFT_SP'):spe_value=Decimal(saved['E_ads_DFT_SP'])
                total=combined(ml,spe_value,relaxed_value)
                cells[5]=f'<td class="sp-energy sp-combined">{total:+.3f}</td>' if total is not None else '<td class="sp-energy sp-combined">&mdash;</td>'
                saved.update(ML_minus_SPE_minus_relaxed_DFT=str(total) if total is not None else '',E_ads_ML_relaxed_displayed=str(ml) if ml is not None else '',
                    E_ads_DFT_relaxed_displayed=str(relaxed_value) if relaxed_value is not None else '',relaxed_review=str('energy-review' in cells[1]).lower(),
                    SPE_review=str('energy-review' in cells[4]).lower(),analysis_eligible=str(ml is not None and relaxed_value is not None and 'energy-review' not in cells[1] and 'energy-review' not in cells[4]).lower())
                sp[key]=saved
            return '<tr>'+''.join(cells)+'</tr>'
        return prefix+re.sub(r'<tr>(.*?)</tr>',row,body,flags=re.S)
    page=CARD.sub(card,page)
    summary_text=f'<b>Single-point adsorption energies:</b> {len(sp)} functional results across {len({(k[0],k[1]) for k in sp})} of the 415 systems below; source choices follow the Perlmutter-first policy and retain explicit provenance.'
    page=re.sub(r'<b>Single-point adsorption energies:</b>.*?(?=\n  <br><b>Perlmutter)',summary_text,page,flags=re.S)
    page=re.sub(r'<!-- cluster-policy -->.*?<!-- /cluster-policy -->\n?','',page,flags=re.S)
    note='''<!-- cluster-policy -->
<div class="note info"><b>Source preference:</b> Perlmutter first. For an unusual result (|E<sub>ads</sub>| &gt; 5 eV), a Kestrel alternative must pass completion/reference checks and match the structure fingerprint.
<br><b>Potential compatibility:</b> Complex and references must use matching PAW potential identities, including variant and dataset date.
Known mixed-potential values are withheld until compatible references are available. Historical values remain in the
<a href="dft_potential_mismatches.csv">potential mismatch audit</a>; <a href="dft_potential_audit_summary.json">counts by surface</a>.
<a href="dft_potential_slab_recovery.json">Matching-potential slab recovery inputs and submission</a>.
Matching names or near-identical energies alone are insufficient. Under the matching-potential policy, old values lacking a validated component set are withheld rather than assumed compatible.
Archived Kestrel sources remain labeled; no fresh Kestrel login is claimed.
<br><b>DME reference recovery:</b> PBE and PBE+D3 gas-reference retries are tracked separately; their totals are used only after convergence.
<a href="dft_gas_reference_recovery.json">Gas-reference job status</a>.
<br><a href="dft_cluster_selection.csv" download>Source selection for every cell</a> &middot;
<a href="dft_cluster_duplicates.csv" download>Perlmutter/Kestrel duplicate comparisons</a> &middot;
<a href="dft_cluster_candidates.csv" download>All matched energy candidates</a> &middot;
<a href="dft_cluster_selection_summary.json">Selection policy and counts</a>.</div>
<!-- /cluster-policy -->
'''
    page=page.replace('<h2>Per-system structure',note+'<h2>Per-system structure');path.write_text(page)
    write_csv(sp_path,[sp[k] for k in sorted(sp)],sp_fields)
    selection_fields=[k for k in selection_rows[0] if k not in ('applied_E_ads','application_status')]+['applied_E_ads','application_status']
    for r in selection_rows:
        r.setdefault('applied_E_ads','');r.setdefault('application_status','not on page')
    write_csv(selection_path,selection_rows,selection_fields)
    history=root/'dft_cluster_selection_changes.csv';old=read_csv(history) if history.exists() else []
    fields=['surface','molecule','functional','mode','old_E_ads','new_E_ads','reason']
    seen={tuple(r[k] for k in fields) for r in old}
    for r in changes:
        if tuple(r[k] for k in fields) not in seen:old.append(r);seen.add(tuple(r[k] for k in fields))
    write_csv(history,old,fields)
    summary_path=root/'dft_cluster_selection_summary.json'
    summary=json.loads(summary_path.read_text())
    summary['application_status']=dict(Counter(r['application_status'] for r in selection_rows))
    summary['applied_candidate_cluster']=dict(Counter(r['complex_cluster'] for r in selection_rows if r['application_status']=='selected candidate applied'))
    summary['numeric_changes_recorded']=len(old)
    summary_path.write_text(json.dumps(summary,indent=2)+'\n')
    print('Applied cluster policy:',len(changes),'numeric additions/changes')


if __name__=='__main__':build()
