"""Render archived coordinates with the original ASE figure style and cache surface atlases."""
import io
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
from ase import Atoms
from ase.io import write
from PIL import Image, ImageOps

ROOT=Path(__file__).resolve().parents[1]
METHODS=['pbe','pbe_d3','r2scan','beef_vdw']
ROTATION='-70x,20y,10z'
WIDTH,HEIGHT=240,640


def build(root=ROOT):
    systems=json.loads((root/'functional_geometry.json').read_text())['systems']
    folder=root/'dftcmp/functional';folder.mkdir(parents=True,exist_ok=True)
    manifest=dict(schema_version=1,style='Original ASE PNG renderer; unmodified archived coordinates',rotation=ROTATION,show_unit_cell=2,cell_width=WIDTH,cell_height=HEIGHT,systems={})
    surfaces=sorted({s['surface'] for s in systems.values()})
    for surface in surfaces:
        keys=sorted(k for k,s in systems.items() if s['surface']==surface)
        atlas=Image.new('RGB',(WIDTH*len(METHODS),HEIGHT*len(keys)),'white')
        count=0
        for row,key in enumerate(keys):
            rendered={}
            for col,method in enumerate(METHODS):
                s=systems[key]['structures'].get(method)
                if not s:continue
                atoms=Atoms(symbols=s['symbols'],positions=s['positions'],cell=s['cell'],pbc=s['pbc'])
                buffer=io.BytesIO();write(buffer,atoms,format='png',rotation=ROTATION,show_unit_cell=2)
                buffer.seek(0);image=Image.open(buffer).convert('RGB')
                image=ImageOps.contain(image,(WIDTH,HEIGHT),Image.Resampling.LANCZOS)
                atlas.paste(image,(col*WIDTH+(WIDTH-image.width)//2,row*HEIGHT+(HEIGHT-image.height)//2))
                rendered[method]=dict(atlas=f'dftcmp/functional/{surface}.png',column=col,row=row,columns=len(METHODS),rows=len(keys),source_sha256=s['source_sha256'])
                count+=1
            manifest['systems'][key]=rendered
        path=folder/(surface+'.png');atlas.quantize(colors=256,method=Image.Quantize.FASTOCTREE).save(path,optimize=True)
        print(surface,count,'static figures',flush=True)
    (root/'functional_figures.json').write_text(json.dumps(manifest,separators=(',',':'))+'\n')


if __name__=='__main__':build()
