#!/usr/bin/env python3
"""Use exact, already-published base GOAD images and contacts where Omni images are absent."""
import csv
import hashlib
import html
import json
from pathlib import Path
import re
from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[1]
CARD=re.compile(r'(<div class="g" data-surf="([^"]+)" data-mol="([^"]+)">)(.*?)(?=<div class="g" data-surf=|<script>)',re.S)


def update(root=ROOT):
    root=Path(root)
    page=(root/'dft_comparison.html').read_text()
    tables=re.findall(r'<table\b.*?</table>',page,re.S)
    records={(r['surface'],r['molecule']):r for r in csv.DictReader((root/'bond_distances.csv').open())}
    audit_path=root/'mlip_base_structure_sources.json'
    audit=json.loads(audit_path.read_text()) if audit_path.exists() else {'policy':'Use exact published base GOAD PNG and its bond_distances.csv record; never substitute a 5m or other relaxation variant. Energy tables unchanged. Coordination sites and source coordinates were not archived with these base records.','images':[]}
    def card(match):
        prefix,s,m,body=match.groups();s,m=html.unescape(s),html.unescape(m)
        if '<div class="ph">no MLIP image</div>' not in body:return match.group(0)
        row=records[(s,m)];image=row['png']
        if image!=s+'_'+m+'.png' or not (root/image).is_file():raise ValueError('Missing exact base image')
        source=html.escape(image,quote=True)
        old='<div><div class="ph">no MLIP image</div><div class="cap">GOAD+SevenNet-OMNI</div></div>'
        replacement=f'<div><img src="{source}" alt="{html.escape(s+" "+m,quote=True)} published base GOAD"><div class="cap">GOAD (published base) &middot; <a href="mlip_base_structure_sources.json">source</a></div></div>'
        if old not in body:raise ValueError('Unexpected placeholder structure')
        body=body.replace(old,replacement)
        contact=f'MLIP contact (published base): <b>{html.escape(row["pair"])} {float(row["min_dist"]):.3f} &#8491;</b> &middot; site not archived'
        body=re.sub(r'MLIP contact:.*?(?=<br>|</div>)',contact,body,count=1,flags=re.S)
        audit['images'].append(dict(surface=s,molecule=m,image=image,image_sha256=hashlib.sha256((root/image).read_bytes()).hexdigest(),
                                   source_csv='bond_distances.csv',source_record=row,geometry_variant='published base GOAD'))
        return prefix+body
    page=CARD.sub(card,page)
    page=re.sub(r'<!-- base-geometry-note -->.*?<!-- /base-geometry-note -->\n?', '',page,flags=re.S)
    note=f'''<!-- base-geometry-note -->
<div class="note info"><b>Geometry coverage:</b> all 415 systems now have images on both sides.
{len(audit['images'])} MLIP slots use the exact <b>published base GOAD</b> image and its archived contact distance because an Omni image is unavailable;
the captions identify this variant. Coordination sites and original coordinates for these base records are not archived here.
These additions provide structural context; they do not change the geometries or values associated with the energy tables.
<a href="mlip_base_structure_sources.json">Base-image and contact provenance</a>.</div>
<!-- /base-geometry-note -->
'''
    page=page.replace('<h2>Per-system structure',note+'<h2>Per-system structure')
    status='''<div class="status">
<span class="live">MLIP IMAGES 415/415</span> 348 SevenNet-OMNI images and 67 published base GOAD images.
<span class="live">MLIP CONTACTS 415/415</span>
<span class="live">DFT IMAGES 415/415</span>
<span class="live pend">DFT CONTACTS 400/415</span>
The remaining 15 DFT contact measurements require the original Kestrel CONTCARs used to render their images.
</div>'''
    page=re.sub(r'<div class="status">.*?</div>',status,page,count=1,flags=re.S)
    assert re.findall(r'<table\b.*?</table>',page,re.S)==tables
    soup=BeautifulSoup(page,'html.parser')
    assert len(soup.select('.g .pair img'))==830
    assert not soup.select('.ph')
    audit_path.write_text(json.dumps(audit,indent=2)+'\n')
    (root/'dft_comparison.html').write_text(page)
    print('Filled base-image coverage:',len(audit['images']),'records; all 830 structure images present')


if __name__=='__main__':update()
