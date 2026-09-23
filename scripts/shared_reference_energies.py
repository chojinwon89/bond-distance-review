#!/usr/bin/env python3
"""Fill missing adsorption results using audited, size-matched stored references."""
import csv
from decimal import Decimal
import html
import json
from pathlib import Path
import re

from component_store import read_store, latest, reference_candidates, relaxed
from molecule_names import canonical

ROOT = Path(__file__).resolve().parents[1]
FIELDS = ['surface','molecule','functional','mode','status','E_ads','E_complex','E_slab','E_molecule',
          'complex_cluster','slab_cluster','molecule_cluster','complex_directory','slab_directory','molecule_directory',
          'complex_snapshot_id','slab_snapshot_id','molecule_snapshot_id','provenance','reference_status','validation_notes']


def read_rows(root=ROOT):
    path=root/'dft_shared_reference_energies.csv'
    return list(csv.DictReader(path.open())) if path.exists() else []


def compatible_settings(comp, slab, gas):
    """Preserve project reference conventions while checking stored core settings."""
    for reference in [slab,gas]:
        if Decimal(comp.get('ENCUT','0')) != Decimal(reference.get('ENCUT','-1')):
            return False
    for key,default in [('ISPIN','1'),('ISMEAR','1'),('SIGMA','0.2')]:
        if Decimal(comp.get(key,default)) != Decimal(slab.get(key,default)):
            return False
    for reference in [slab,gas]:
        if comp.get('LASPH','.FALSE.').strip('.').upper() != reference.get('LASPH','.FALSE.').strip('.').upper():
            return False
    return True


def derive(comp, slab, gas, mode, allow_unverified=False):
    c,s,g=[r['calculation'] for r in [comp,slab,gas]]
    if c['status']!='converged' or (c.get('nsw')!=0 if mode=='SPE' else not relaxed(c)):
        raise ValueError('Complex does not pass the requested calculation mode and completion checks')
    settings_ok=compatible_settings(c.get('settings',{}),s.get('settings',{}),g.get('settings',{}))
    if not allow_unverified and not settings_ok:
        raise ValueError('Stored cutoff, slab spin/smearing or LASPH settings differ')
    refs=reference_candidates(comp,[slab,gas],allow_unverified=allow_unverified)['candidate_reference_ids']
    if slab['snapshot_id'] not in refs['slab'] or gas['snapshot_id'] not in refs['molecule']:
        raise ValueError('References fail convergence, functional, identity, atom-count or full-cell matching')
    energy=Decimal(str(c['energy']))-Decimal(str(s['energy']))-Decimal(str(g['energy']))
    result=dict(surface=comp['surface'],molecule=comp['molecule'],functional=c['functional'],mode=mode,
        status='energy_review' if abs(energy)>5 else 'ok',E_ads=str(energy),
        E_complex=c['energy'],E_slab=s['energy'],E_molecule=g['energy'],
        provenance='Converged components; matching PAW TITEL identities, functional, composition, full slab cell and stored core settings. Archived Kestrel provenance retained. No energy scaling or offsets.')
    from potential_matching import potential_match
    from component_store import cell_matches
    notes=[]
    if not cell_matches(c.get('cell'),s.get('cell')):notes.append('slab cell differs or is unverified')
    if not potential_match(c,s) or not potential_match(c,g):notes.append('potential identities differ or are unverified')
    if not settings_ok:notes.append('core settings differ or are unverified')
    result['reference_status']='unverified' if notes else 'matched'
    result['validation_notes']='; '.join(notes)
    if notes:result['provenance']='Converged components; same functional and composition. Displayed under user-requested permissive reference policy. '+result['validation_notes']
    for role,row in [('complex',comp),('slab',slab),('molecule',gas)]:
        result[role+'_cluster']=row['cluster'];result[role+'_directory']=row['calculation']['directory']
        result[role+'_snapshot_id']=row['snapshot_id']
    return result


def build(root=ROOT):
    records=read_store(root/'dft_component_store.jsonl.gz'); current=latest(records)
    byid={r['snapshot_id']:r for r in records}
    # Retain already-derived rows when later coverage lists no longer call them gaps.
    selected={(r['surface'],r['molecule'],r['functional'],r['mode']):r for r in read_rows(root)}
    requests=json.loads((root/'dft_kestrel_search_requests.json').read_text())['requests']
    spe_allowed=set()
    for filename in ['dft_perlmutter_singlepoint_audit.csv','dft_kestrel_singlepoint_audit.csv']:
        for r in csv.DictReader((root/filename).open()):
            if r['status'] in ('ok','energy_review','reference_unavailable','potential_mismatch'):
                spe_allowed.add(r['complex_directory'])
    audit=[]
    for request in requests:
        key=(request['surface'],canonical(request['molecule']),request['functional'],request['mode'])
        candidates=[byid[c['complex_snapshot_id']] for c in request['completed_complexes']]
        # Prefer the current Perlmutter complex, then the base layout; never select
        # a candidate because its adsorption energy is particularly low.
        candidates.sort(key=lambda r:(r['cluster']!='perlmutter',len(Path(r['calculation']['directory']).parts),r['calculation']['directory']))
        for comp in candidates:
            if request['mode']=='SPE' and comp['calculation']['directory'] not in spe_allowed:
                continue
            references=reference_candidates(comp,current)['candidate_reference_ids']
            def order(snapshot):
                row=byid[snapshot]
                return (row['cluster']!=comp['cluster'],row['calculation']['directory'])
            found=None
            for sid in sorted(references['slab'],key=order):
                for gid in sorted(references['molecule'],key=order):
                    try:found=derive(comp,byid[sid],byid[gid],request['mode'])
                    except ValueError as error:
                        audit.append(dict(surface=key[0],molecule=key[1],functional=key[2],mode=key[3],reason=str(error)))
                    if found:break
                if found:break
            if found:
                selected.setdefault(key,found)
                break
    # Revalidate stored selections by immutable IDs on every regeneration.
    for key,row in list(selected.items()):
        try:selected[key]=derive(*(byid[row[role+'_snapshot_id']] for role in ['complex','slab','molecule']),row['mode'])
        except ValueError as error:
            audit.append(dict(surface=key[0],molecule=key[1],functional=key[2],mode=key[3],reason=str(error)))
            del selected[key]
    with (root/'dft_shared_reference_energies.csv').open('w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=FIELDS);writer.writeheader();writer.writerows(selected[k] for k in sorted(selected))
    (root/'dft_shared_reference_audit.json').write_text(json.dumps(dict(
        active_kestrel_project_root='/scratch/jcho5/goad-global-optimization',
        historical_kestrel_root='/kfs3/scratch/jcho5/goad-global-optimization',
        note='User supplied active Kestrel path. Archived source paths remain unchanged; root spellings are not silently equated or merged. This run uses saved records, not a live Kestrel filesystem read.',
        policy='Only fill missing cells. Require full slab cell within 1e-5 Angstrom and matching atom counts, functional, completion and stored core settings. SPE complexes must also pass the original cluster audit through the complex checks. KPOINTS and selective-dynamics constraints are not contained in the archived component export; no new numerical-convergence claim is made.',
        selected=len(selected),rejected_candidates=audit),indent=2)+'\n')
    print('Shared-reference results:',len(selected),'relaxed:',sum(r['mode']=='relaxed' for r in selected.values()),
          'SPE:',sum(r['mode']=='SPE' for r in selected.values()),'energy review:',sum(r['status']=='energy_review' for r in selected.values()))


def fill_relaxed(page, rows):
    from update_kestrel_singlepoint import CARD, FUNCTIONALS, number
    lookup={(r['surface'],canonical(r['molecule']),r['functional']):r for r in rows if r['mode']=='relaxed'}
    def card(match):
        prefix,surface,molecule,body=match.groups()
        def row(match):
            cells=re.findall(r'<td\b[^>]*>.*?</td>',match[1],re.S)
            if len(cells)<4:return match[0]
            label=html.unescape(re.sub(r'<[^>]+>','',cells[0]));functional=FUNCTIONALS.get(label)
            result=lookup.get((html.unescape(surface),canonical(html.unescape(molecule)),functional))
            if not result or number(cells[1]) is not None:return match[0]
            value=Decimal(result['E_ads']);ml=number(cells[2]);review=result['status']=='energy_review'
            marker=' class="energy-review"' if review else ''
            title=html.escape(result['provenance']+' Slab: '+result['slab_directory']+(' ENERGY REVIEW: |E_ads| > 5 eV.' if review else ''),quote=True)
            cells[1]=f'<td{marker} title="{title}">{value:.3f}</td>'
            cells[3]=f'<td{marker} title="ML minus relaxed DFT">{ml-value:+.3f}</td>' if ml is not None else '<td>&mdash;</td>'
            return '<tr>'+''.join(cells)+'</tr>'
        return prefix+re.sub(r'<tr>(.*?)</tr>',row,body,flags=re.S)
    return CARD.sub(card,page)


if __name__=='__main__':build()
