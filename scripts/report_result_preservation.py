"""Regenerate the local archive manifest and its publication recovery summary."""
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from bs4 import BeautifulSoup
from publication_store import read, REL
from plot_energy_decomposition import load_rows

ROOT=Path(__file__).resolve().parents[1]


def main(root=ROOT):
    records=read(root)
    with (root/'dft_cluster_selection.csv').open() as f:selection=list(csv.DictReader(f))
    historic=[r for r in selection if r['application_status'].startswith('historical retained:')]
    selected=[r for r in selection if r['application_status']=='selected candidate applied']
    components=json.loads((root/'dft_component_store_summary.json').read_text())
    files=[REL,Path('dft_component_store.jsonl.gz')]
    geometry_path=root/'functional_geometry_audit.json'
    geometry=json.loads(geometry_path.read_text()) if geometry_path.exists() else None
    if geometry:files.extend([Path('functional_geometry.json'),Path(geometry['coordinate_snapshot'])])
    manifest=dict(schema_version=1,description='Local immutable observations; active validation remains separate',
        historical_baseline_commit='23048650b883c5118871d1829a317187a65f0989',
        refreshed_baseline_commit='9d88d3812666d5d0822f160999b9bdf4f0925a02',
        files=[dict(path=str(p),sha256=hashlib.sha256((root/p).read_bytes()).hexdigest(),bytes=(root/p).stat().st_size) for p in files],
        publication_observations=len(records),component_snapshots=components['component_snapshots'],
        current_component_records_by_cluster=components['by_cluster'],
        historical_cells_by_mode=dict(Counter(r['mode'] for r in historic)),
        historical_cells_by_current_audit=dict(Counter(r['application_status'].split(': ',1)[1] for r in historic)),
        selected_component_derived_cells_by_mode=dict(Counter(r['mode'] for r in selected)),
        selected_reference_status=dict(Counter(r.get('reference_status','matched') for r in selected)),
        selected_energy_review_cells=sum(r['status']=='energy_review' for r in selected),
        validated_paired_plot_points=len(load_rows(root/'dft_comparison_singlepoint.csv')),
        caveats=['Current audit labels describe selection blockers, not proof of the exact references used in every old publication.',
                 'Historical values are excluded from validated figures.', 'Kestrel records are archived imports; this recovery did not log in to Kestrel.',
                 'The ledger contains extracted results and provenance, not complete raw VASP output files.'])
    if geometry:manifest['geometry']=dict(structures=geometry['structures'],coordinate_snapshot=geometry['coordinate_snapshot'])
    (root/'data/adsorption_results/manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    n=len(historic)
    note=f'''<!-- result-preservation -->
<div class="note info"><b>Results preserved · 22 September 2026:</b> {len(selected)} component-derived values are displayed under the current reference policy. Another {n} previously published values are retained: {sum(r['mode']=='relaxed' for r in historic)} relaxed and {sum(r['mode']=='SPE' for r in historic)} SPE. Unverified references and historical values are excluded from the validated plots; the current paired analysis contains {manifest['validated_paired_plot_points']} points.
<br>“Refs unverified” means this audit cannot verify a compatible component set; “potential mismatch” identifies incompatible potential candidates in the current audit. Neither label deletes the stored observation. <a href="data/adsorption_results/README.md">Archive and explanation</a> · <a href="data/adsorption_results/manifest.json">Recovery counts and file checksums</a>.</div>
<!-- /result-preservation -->
'''
    for name in ['index.html','dft_comparison.html']:
        p=root/name;s=p.read_text()
        s=re.sub(r'<!-- result-preservation -->.*?<!-- /result-preservation -->\n?','',s,flags=re.S)
        s=s.replace('Current validated tables contain 772 relaxed and 1010 SPE energies, including review-marked values.', 'That refresh selected 772 relaxed and 1010 SPE energies with matched components, including magnitude-review values; historical results have since been restored below.')
        s=s.replace('<!-- /refresh-20260922 -->','<!-- /refresh-20260922 -->\n'+note)
        p.write_text(s)
    p=root/'dft_data_gaps.html';s=p.read_text()
    if '<!-- preservation-caveat -->' not in s:
        s=s.replace('<body>','<body><!-- preservation-caveat --><p><strong>Dated validation-gap snapshot:</strong> Previously blank cells may now show historical numbers. Those numbers still need compatible verified references; these search requests remain relevant. See <a href="dft_comparison.html">the updated tables</a> and <a href="dft_completion_coverage.csv">current coverage</a>.</p>',1)
        p.write_text(s)
    p=root/'dft_cluster_selection_summary.json';summary=json.loads(p.read_text())
    summary['policy']=summary['policy'].replace('Withhold historical numbers lacking validated references.', 'Preserve historical numbers with explicit audit labels when validated references are unavailable; exclude them from validated statistics.')
    p.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
