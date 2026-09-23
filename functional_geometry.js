/* Coordinates and measurements come from saved VASP/MLIP files, not illustrative models. */
'use strict';
(() => {
 const methods=['mlip','pbe','pbe_d3','r2scan','beef_vdw'];
 const names={mlip:'MLIP',pbe:'PBE',pbe_d3:'PBE+D3',r2scan:'r²SCAN',beef_vdw:'BEEF-vdW'};
 const colors={H:'#eeeeee',C:'#50565f',N:'#5288eb',O:'#ec5559',S:'#efcf51',Ag:'#aab6c6',Au:'#d9bc59',Cu:'#c58b62',Pt:'#b9c0c9',Pd:'#9bc8c8',Rh:'#aaaedb',Ir:'#9baac1'};
 const el=(tag,cls,text)=>{const n=document.createElement(tag);if(cls)n.className=cls;if(text!==undefined)n.textContent=text;return n;};
 const num=(n,d=3)=>Number.isFinite(n)?n.toFixed(d):'—';
 function table(caption,headers,rows){const wrap=el('div','geometry-scroll'),t=el('table','geometry-table');t.append(el('caption','',caption));const h=el('tr');headers.forEach(v=>h.append(el('th','',v)));t.append(h);rows.forEach(values=>{const tr=el('tr');values.forEach(v=>tr.append(el('td','',String(v))));t.append(tr);});wrap.append(t);return wrap;}
 function download(s,label){const symbols=[],counts=[];s.symbols.forEach(x=>{if(symbols.at(-1)!==x){symbols.push(x);counts.push(1);}else counts[counts.length-1]++;});
   const text=[label,'1.0',...s.cell.map(v=>v.join(' ')),symbols.join(' '),counts.join(' '),'Cartesian',...s.positions.map(v=>v.join(' '))].join('\n')+'\n';
   const a=el('a');a.href=URL.createObjectURL(new Blob([text],{type:'text/plain'}));a.download=label+'.vasp';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);
 }
 function viewer(canvas,s,metal){let yaw=.35,pitch=1.08;const ads=new Set(s.adsorbate_indices),pos=s.display_positions;
   function draw(){const box=canvas.getBoundingClientRect(),w=Math.max(100,box.width),h=Math.max(100,box.height),dpr=Math.min(devicePixelRatio||1,2);canvas.width=w*dpr;canvas.height=h*dpr;const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);ctx.fillStyle='#11151e';ctx.fillRect(0,0,w,h);
     const points=pos.map(([x,y,z],i)=>{const xx=x*Math.cos(yaw)-y*Math.sin(yaw),yy=x*Math.sin(yaw)+y*Math.cos(yaw);return {i,x:xx,y:yy*Math.cos(pitch)-z*Math.sin(pitch),depth:yy*Math.sin(pitch)+z*Math.cos(pitch)};});
     const xmin=Math.min(...points.map(p=>p.x)),xmax=Math.max(...points.map(p=>p.x)),ymin=Math.min(...points.map(p=>p.y)),ymax=Math.max(...points.map(p=>p.y));const scale=Math.min((w-45)/Math.max(xmax-xmin,1),(h-50)/Math.max(ymax-ymin,1));points.forEach(p=>{p.px=w/2+(p.x-(xmin+xmax)/2)*scale;p.py=h/2+(p.y-(ymin+ymax)/2)*scale;});
     ctx.lineWidth=2;ctx.strokeStyle='#bcc6d2';(s.molecular_bonds||[]).forEach(([i,j])=>{const a=points[i],b=points[j];ctx.beginPath();ctx.moveTo(a.px,a.py);ctx.lineTo(b.px,b.py);ctx.stroke();});
     const focus=s.contacts.filter(c=>c.element!=='H').slice(0,3);if(!focus.length)focus.push(...s.contacts.slice(0,2));ctx.setLineDash([3,3]);ctx.strokeStyle='#6dccdb';ctx.lineWidth=1;
     focus.forEach(c=>{const a=points[c.atom],b=points[c.metal_atom];const displayed=Math.hypot(...pos[c.atom].map((v,k)=>v-pos[c.metal_atom][k]));if(Math.abs(displayed-c.distance_A)<.01){ctx.beginPath();ctx.moveTo(a.px,a.py);ctx.lineTo(b.px,b.py);ctx.stroke();}});ctx.setLineDash([]);
     points.sort((a,b)=>a.depth-b.depth).forEach(p=>{const symbol=s.symbols[p.i],r=Math.max(3,Math.min(11,scale*(symbol==='H'?.22:ads.has(p.i)?.4:.64)));ctx.beginPath();ctx.arc(p.px,p.py,r,0,2*Math.PI);ctx.fillStyle=colors[symbol]||'#ab9fcc';ctx.fill();ctx.strokeStyle=ads.has(p.i)?'#e2e8f0':'#68768c';ctx.lineWidth=.6;ctx.stroke();if(focus.some(c=>c.atom===p.i)){ctx.fillStyle='#fff';ctx.font='bold 10px sans-serif';ctx.fillText(s.labels[p.i],p.px+r+2,p.py-4);}});
   }
   let dragging=false,last;canvas.addEventListener('pointerdown',e=>{dragging=true;last=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);});canvas.addEventListener('pointermove',e=>{if(!dragging)return;yaw+=(e.clientX-last[0])*.015;pitch+=(e.clientY-last[1])*.015;last=[e.clientX,e.clientY];draw();});canvas.addEventListener('pointerup',()=>dragging=false);canvas.addEventListener('pointercancel',()=>dragging=false);
   new ResizeObserver(draw).observe(canvas);draw();return view=>{pitch=view==='top'?0:1.08;yaw=.35;draw();};
 }
 function populate(card,system){if(card.classList.contains('geometry-ready'))return;card.classList.add('geometry-ready');const gallery=el('div','functional-gallery'),controls=el('div','geometry-tools'),views=[];
   ['side','top'].forEach(view=>{const b=el('button','',view==='side'?'Side / tilted view':'Top view');b.type='button';b.onclick=()=>views.forEach(f=>f(view));controls.append(b);});
   const message=el('p','geometry-note','Drag a structure to rotate. Atom labels identify the molecular contacts. Distances use periodic boundaries; dashed lines are geometric contacts, not bond-order assignments.');
   methods.forEach(method=>{const s=system.structures[method],fig=el('figure'),cap=el('figcaption','',names[method]);
     if(s){const canvas=el('canvas');canvas.setAttribute('aria-label',system.surface+' '+system.molecule+' '+names[method]+' atomic structure');canvas.setAttribute('role','img');fig.append(canvas);cap.append(el('div',s.status==='converged'||method==='mlip'?'geometry-note':'geometry-status',method==='mlip'?s.model:(s.cluster+' · '+(s.status==='converged'?'relaxed; converged':'last saved geometry; '+s.status))));
       if(method!=='mlip'&&!s.matches_energy_complex)cap.append(el('div','geometry-note','Geometry source differs from the selected energy source.'));
       cap.title=s.source_path;const b=el('button','geometry-download','Download coordinates');b.type='button';b.onclick=()=>download(s,system.surface+'_'+system.molecule+'_'+method);cap.append(b);requestAnimationFrame(()=>views.push(viewer(canvas,s,system.metal)));
     }else fig.append(el('div','missing-geometry','No readable final coordinates available for '+names[method]+'. A stored energy does not supply atomic coordinates.'));
     fig.append(cap);gallery.append(fig);
   });
   const detail=el('div','geometry-data'),summaryRows=methods.map(method=>{const s=system.structures[method];return [names[method],s?num(s.nearest_contact_A):'—',s?num(s.nearest_heavy_contact_A):'—',s?s.contacts[0].metal_label+'–'+s.contacts[0].label:'—'];});
   detail.append(table('Surface–molecule geometry (Å)',['Method','Nearest contact','Nearest heavy atom','Closest atom pair'],summaryRows));
   const contactDetails=el('details');contactDetails.append(el('summary','','Compare specific contacting atoms · DFT versus MLIP'));
   methods.slice(1).forEach(method=>{const s=system.structures[method];if(!s)return;const pairs=s.paired_contacts||[];
     contactDetails.append(el('p','geometry-note',names[method]+': '+(s.atom_mapping_note||'No MLIP coordinates available')));
     if(pairs.length)contactDetails.append(table(names[method]+' contacts',['ML atom ↔ DFT atom','MLIP (Å)','DFT (Å)','DFT − MLIP','Metal sites (ML / DFT)'],pairs.map(p=>[p.ML_atom+' ↔ '+p.DFT_atom,num(p.ML_distance_A),num(p.DFT_distance_A),num(p.delta_A),p.ML_metal+' / '+p.DFT_metal])));
   });
   detail.append(contactDetails);
   const angles=el('details');angles.append(el('summary','','Bond angles and torsions · matched molecular atoms'));
   let angularCount=0;
   methods.slice(1).forEach(method=>{const s=system.structures[method];if(!s)return;['surface_angles','angles','torsions'].forEach(kind=>{const values=s['paired_'+kind]||[];angularCount+=values.length;if(values.length)angles.append(table(names[method]+' · '+(kind==='surface_angles'?'metal–anchor–neighbor angles':kind==='angles'?'bond angles':'torsions')+' (°)',['ML atom path ↔ DFT atom path','MLIP','DFT','DFT − MLIP'],values.map(p=>[p.label+' ↔ '+p.DFT_label,num(p.ML_degrees,1),num(p.DFT_degrees,1),num(p.delta_degrees,1)])));});});
   const rawAngles=el('details');rawAngles.append(el('summary','','Available molecular angles and torsions by method'));methods.forEach(method=>{const s=system.structures[method];if(!s)return;['angles','torsions'].forEach(kind=>{const values=s[kind]||[];if(values.length)rawAngles.append(table(names[method]+' '+kind,['Atom path','Angle (°)'],values.map(v=>[v.label,num(v.degrees,1)])));});});angles.append(rawAngles);
   angles.append(el('p','geometry-note',angularCount?'Angles follow molecular bonds; terminal hydrogen atoms provide torsions where a heavy-atom path is unavailable. Torsion differences wrap to −180°…180°. Symmetry-equivalent hydrogens can give equivalent torsion choices.':'No comparable molecular angle or torsion is defined for the available structures.'));
   detail.append(angles);
   const all=el('details');all.append(el('summary','','All atom-to-surface distances, including hydrogen'));
   methods.forEach(method=>{const s=system.structures[method];if(s)all.append(table(names[method],['Molecular atom','Nearest metal','Distance (Å)','Contact candidate'],s.contacts.map(c=>[c.label,c.metal_label,num(c.distance_A),c.contact_candidate?'yes':'no'])));});detail.append(all);
   const oldDFT=card.querySelector('.pair > div:nth-child(2)');if(oldDFT){const archive=el('details');archive.append(el('summary','','Previously published DFT image'));const oldMeta=card.querySelector('.m');if(oldMeta)archive.append(oldMeta.cloneNode(true));archive.append(oldDFT.cloneNode(true));detail.append(archive);}
   const first=card.querySelector('.energy-scroll');card.insertBefore(controls,first);card.insertBefore(message,first);card.insertBefore(gallery,first);card.append(detail);
 }
 const intro=el('div','geometry-preamble','Four functional-specific DFT structures and atom-resolved geometry comparisons. Expand the contact and angle tables for details.');
 const heading=Array.from(document.querySelectorAll('h2')).find(n=>n.textContent.startsWith('Per-system structure'));if(heading)heading.after(intro);
 fetch('functional_geometry.json').then(r=>{if(!r.ok)throw Error('Geometry data unavailable');return r.json();}).then(data=>{
   intro.append(document.createTextNode(' '));[['Geometry audit','functional_geometry_audit.json'],['Contact distances CSV','functional_geometry_contacts.csv']].forEach(([label,url])=>{const a=el('a','',label);a.href=url;intro.append(a,document.createTextNode(' · '));});
   const observer=new IntersectionObserver(entries=>entries.forEach(entry=>{if(!entry.isIntersecting)return;const card=entry.target,key=card.dataset.surf+'|'+card.dataset.mol;const system=data.systems[key];if(system)populate(card,system);observer.unobserve(card);}),{rootMargin:'500px'});
   document.querySelectorAll('.g').forEach(card=>observer.observe(card));
 }).catch(error=>{intro.textContent='Functional geometry could not load: '+error.message+'. The saved energy table remains available.';});
})();
