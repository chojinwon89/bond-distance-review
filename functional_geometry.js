/* Coordinates and measurements come from saved VASP/MLIP files, not illustrative models. */
'use strict';
(() => {
 const methods=['mlip','pbe','pbe_d3','r2scan','beef_vdw'];
 const names={mlip:'MLIP',pbe:'PBE',pbe_d3:'PBE+D3',r2scan:'r²SCAN',beef_vdw:'BEEF-vdW'};
 const el=(tag,cls,text)=>{const n=document.createElement(tag);if(cls)n.className=cls;if(text!==undefined)n.textContent=text;return n;};
 const num=(n,d=3)=>Number.isFinite(n)?n.toFixed(d):'—';
 function table(caption,headers,rows){const wrap=el('div','geometry-scroll'),t=el('table','geometry-table');t.append(el('caption','',caption));const h=el('tr');headers.forEach(v=>h.append(el('th','',v)));t.append(h);rows.forEach(values=>{const tr=el('tr');values.forEach(v=>tr.append(el('td','',String(v))));t.append(tr);});wrap.append(t);return wrap;}
 function populate(card,system,figures,restarts,completion){if(card.classList.contains('geometry-ready'))return;card.classList.add('geometry-ready');const gallery=el('div','functional-gallery');
   methods.forEach(method=>{const s=system.structures[method],fig=el('figure'),cap=el('figcaption','',names[method]);
     const saved=figures[method];const oldML=card.querySelector('.pair > div:first-child img');
     if(method==='mlip'&&oldML){const image=oldML.cloneNode(true);image.className='original-mlip-figure';image.loading='lazy';fig.append(image);}
     else if(s&&saved&&s.source_sha256===saved.source_sha256){const frame=el('div','static-figure-frame'),image=el('div','static-structure');
       image.setAttribute('role','img');image.setAttribute('aria-label',system.surface+' '+system.molecule+' DFT '+names[method]+' structure');
       image.style.backgroundImage='url("'+saved.atlas+'")';
       image.style.backgroundSize=(saved.columns*100)+'% '+(saved.rows*100)+'%';
       image.style.backgroundPosition=(saved.column/(saved.columns-1)*100)+'% '+(saved.rows>1?saved.row/(saved.rows-1)*100:0)+'%';
       frame.append(image);fig.append(frame);
     }else fig.append(el('div','missing-geometry',s?'Figure awaiting regeneration for '+names[method]+'.':'No readable final coordinates available for '+names[method]+'.'));
     if(s){cap.append(el('div',s.status==='converged'||method==='mlip'?'geometry-note':'geometry-status',method==='mlip'?s.model:(s.cluster+' · '+(s.status==='converged'?'relaxed; converged':'last saved geometry; '+s.status))));
       const restart=restarts.jobs&&restarts.jobs[s.source_path];
       if(method!=='mlip'&&s.status!=='converged'&&restart&&restart.source_sha256===s.source_sha256){
         cap.append(el('div','geometry-note','Previous run: '+restart.previous_cause+'.'));
         const note=el('div','geometry-note','Restart '+restart.restart_status+' · '+restart.job_id);
         note.append(el('div','', 'Checked '+restarts.checked_at.slice(0,16).replace('T',' ')+' UTC'));
         note.title=restart.restart_directory;cap.append(note);
         if(restart.restart_status==='converged')cap.append(el('div','geometry-note','New result verified; this image still shows the previous run until the next geometry refresh.'));
       }
       if(method!=='mlip'&&!s.matches_energy_complex)cap.append(el('div','geometry-note','Geometry source differs from the selected energy source.'));
       cap.title=s.source_path;
     }
     const pending=(completion.jobs||[]).find(j=>j.role==='relaxed'&&j.surface===system.surface&&j.molecule===system.molecule&&j.functional===method);
     if(pending&&(!s||s.status!=='converged')){cap.append(el('div','geometry-note','DFT relaxation '+pending.status+' · '+pending.job_id));cap.append(el('div','geometry-note','Checked '+completion.checked_at.slice(0,16).replace('T',' ')+' UTC'));}
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
   const first=card.querySelector('.energy-scroll');card.insertBefore(gallery,first);card.append(detail);
 }
 const intro=el('div','geometry-preamble','Original structure figure style, with separate images for all four DFT functionals. Bond distances, angles and torsions are reported in the comparison tables below.');
 const heading=Array.from(document.querySelectorAll('h2')).find(n=>n.textContent.startsWith('Per-system structure'));if(heading)heading.after(intro);
 const restartData=fetch('relaxation_restart_status.json',{cache:'no-cache'}).then(r=>{if(!r.ok)throw Error('Restart audit unavailable');return r.json();}).catch(()=>({jobs:{}}));
 const completionData=fetch('missing_component_status.json',{cache:'no-cache'}).then(r=>{if(!r.ok)throw Error('Completion status unavailable');return r.json();}).catch(()=>({jobs:[]}));
 Promise.all([...['functional_geometry.json','functional_figures.json'].map(url=>fetch(url).then(r=>{if(!r.ok)throw Error('Geometry data unavailable');return r.json();})),restartData,completionData]).then(([data,figures,restarts,completion])=>{
   intro.append(document.createTextNode(' '));[['Geometry audit','functional_geometry_audit.json'],['Contact distances CSV','functional_geometry_contacts.csv']].forEach(([label,url])=>{const a=el('a','',label);a.href=url;intro.append(a,document.createTextNode(' · '));});
   if(restarts.summary){const a=el('a','','Restart audit ('+restarts.summary.total+' jobs)');a.href='relaxation_restart_status.json';intro.append(a);intro.append(el('div','geometry-note','Restart status is a saved snapshot from '+restarts.checked_at.slice(0,16).replace('T',' ')+' UTC. Kestrel comparison covers archived results; live verification is pending SSH access.'));}
   if(completion.summary){const line=el('div','geometry-note'),a=el('a','','Missing-component batch '+completion.array_job_id);a.href='missing_component_status.json';line.append(a,document.createTextNode(': '+completion.summary.roles.molecule+' gas references, '+completion.summary.roles.slab+' clean slabs, '+completion.summary.roles.spe+' SPE jobs and '+completion.summary.roles.relaxed+' relaxations. Status checked '+completion.checked_at.slice(0,16).replace('T',' ')+' UTC.'));intro.append(line);}
   const observer=new IntersectionObserver(entries=>entries.forEach(entry=>{if(!entry.isIntersecting)return;const card=entry.target,key=card.dataset.surf+'|'+card.dataset.mol;const system=data.systems[key];if(system)populate(card,system,figures.systems[key]||{},restarts,completion);observer.unobserve(card);}),{rootMargin:'500px'});
   document.querySelectorAll('.g').forEach(card=>observer.observe(card));
 }).catch(error=>{intro.textContent='Functional geometry could not load: '+error.message+'. The saved energy table remains available.';});
})();
