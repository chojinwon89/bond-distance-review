"""Portable read-only HPC export and append-only local import, including coordinates."""
import argparse,gzip,hashlib,json
from pathlib import Path
from component_store import collect,merge,validate
ROOT=Path(__file__).resolve().parents[1]


def export(cluster,project,output):
    records=collect(cluster,project);coordinates=[]
    for r in records:
        c=r['calculation']
        if r['role']!='complex' or not (c.get('nsw') or 0)>0:continue
        path=Path(c['directory'])/'CONTCAR'
        if not path.exists() or not path.stat().st_size:continue
        data=path.read_bytes()
        coordinates.append(dict(cluster=cluster,directory=c['directory'],snapshot_id=r['snapshot_id'],sha256=hashlib.sha256(data).hexdigest(),vasp=data.decode()))
    payload=dict(schema_version=1,cluster=cluster,components=records,coordinates=coordinates)
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_bytes(gzip.compress(json.dumps(payload).encode(),mtime=0))
    print('Exported',len(records),'components and',len(coordinates),'coordinate files to',output)


def import_bundle(path,root=ROOT):
    payload=json.loads(gzip.decompress(path.read_bytes()))
    if payload['schema_version']!=1:raise ValueError('Unknown bundle schema')
    for r in payload['components']:validate(r)
    byid={r['snapshot_id']:r for r in payload['components']}
    for c in payload['coordinates']:
        r=byid[c['snapshot_id']]
        if (c['cluster'],c['directory'])!=(r['cluster'],r['calculation']['directory']):raise ValueError('Coordinate source mismatch')
        if hashlib.sha256(c['vasp'].encode()).hexdigest()!=c['sha256']:raise ValueError('Coordinate checksum mismatch')
    folder=root/'data/adsorption_results';folder.mkdir(exist_ok=True,parents=True)
    index=folder/'coordinate_sources.json';sources=json.loads(index.read_text()) if index.exists() else {}
    for c in payload['coordinates']:
        destination=folder/'cluster_coordinates'/c['sha256']/'CONTCAR';destination.parent.mkdir(exist_ok=True,parents=True)
        if not destination.exists():destination.write_text(c['vasp'])
        sources[c['cluster']+'|'+c['directory']]=dict(path=str(destination.relative_to(root)),sha256=c['sha256'],snapshot_id=c['snapshot_id'])
    # All input data were checked before the first write; old snapshots are never removed.
    merge(root/'dft_component_store.jsonl.gz',payload['components'])
    index.write_text(json.dumps(sources,indent=2)+'\n')
    print('Merged',len(payload['components']),'component observations;',len(sources),'cached coordinate sources')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    e=sub.add_parser('export');e.add_argument('--cluster',choices=['kestrel','perlmutter'],required=True);e.add_argument('--project-root',type=Path,required=True);e.add_argument('--output',type=Path,required=True)
    i=sub.add_parser('import');i.add_argument('--input',type=Path,required=True)
    a=p.parse_args()
    if a.command=='export':export(a.cluster,a.project_root,a.output)
    else:import_bundle(a.input)
