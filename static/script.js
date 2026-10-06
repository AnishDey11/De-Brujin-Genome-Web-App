const $ = id => document.getElementById(id);
let dataset = {};
let currentData = null;

function setStatus(message, ok=true){
  const el = $('status');
  el.textContent = message;
  el.classList.toggle('error', !ok);
}

function cleanDNA(text){ return String(text || '').toUpperCase().replace(/[^ACGT]/g, ''); }
function formatDNA(seq, width=80){
  const clean = cleanDNA(seq);
  const lines = [];
  for(let i=0;i<clean.length;i+=width) lines.push(clean.slice(i,i+width));
  return lines.join('\n');
}

async function loadSequences(){
  try{
    const r = await fetch('/api/sequences');
    const data = await r.json();
    dataset = data.sequences || {};
    const select = $('sequence');
    select.innerHTML = '';
    Object.entries(dataset).forEach(([name,len])=>{
      const option=document.createElement('option');
      option.value=name;
      option.textContent=`${name}  •  ${len} bp`;
      select.appendChild(option);
    });
    if(!select.options.length) throw new Error('No sequencing reads found.');
    await refreshSelectedSequence();
    setStatus(`${select.options.length} sequencing read(s) loaded.`, true);
    await analyzeSelected();
  }catch(err){ setStatus(err.message, false); }
}

async function refreshSelectedSequence(){
  const id=$('sequence').value;
  if(!id) return;
  const r=await fetch(`/api/sequence-text?sequence_id=${encodeURIComponent(id)}`);
  const d=await r.json();
  if(!r.ok) throw new Error(d.error || 'Could not load sequence.');
  $('selectedSequence').textContent=formatDNA(d.sequence);
  $('selectedMeta').textContent=`${id} • ${d.sequence.length} bases`;
}

async function analyzeSelected(){
  const id=$('sequence').value;
  const k=parseInt($('k').value,10);
  if(!id) return;
  if(!Number.isInteger(k) || k<2) return setStatus('k must be at least 2.', false);
  setStatus(`Analyzing ${id} using only ${id}…`, true);
  try{
    const r=await fetch(`/api/analyze?sequence_id=${encodeURIComponent(id)}&k=${k}`);
    const d=await r.json();
    if(!r.ok) throw new Error(d.error || 'Analysis failed.');
    currentData=d;
    render(d);
    setStatus(`Analysis complete: ${id}. Weights use only this input.`, true);
  }catch(err){ setStatus(err.message,false); }
}

async function analyzeCustom(){
  const sequence=cleanDNA($('customSequence').value);
  const k=parseInt($('k').value,10);
  if(!sequence) return setStatus('Enter a DNA sequence first.',false);
  if(!Number.isInteger(k) || k<2) return setStatus('k must be at least 2.',false);
  $('selectedSequence').textContent=formatDNA(sequence);
  $('selectedMeta').textContent=`CUSTOM • ${sequence.length} bases`;
  setStatus('Analyzing your typed DNA only…',true);
  try{
    const r=await fetch('/api/analyze-custom',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({sequence,k})});
    const d=await r.json();
    if(!r.ok) throw new Error(d.error || 'Custom analysis failed.');
    currentData=d;
    render(d);
    setStatus('Custom sequence analyzed. All weights use only this typed input.',true);
  }catch(err){ setStatus(err.message,false); }
}

async function uploadDataset(){
  const file=$('file').files[0];
  if(!file) return setStatus('Choose a TXT or FASTA file first.',false);
  const fd=new FormData(); fd.append('file',file);
  setStatus('Uploading dataset…',true);
  try{
    const r=await fetch('/api/upload',{method:'POST',body:fd});
    const d=await r.json();
    if(!r.ok) throw new Error(d.error || 'Upload failed.');
    dataset=d.sequences || {};
    const select=$('sequence'); select.innerHTML='';
    Object.entries(dataset).forEach(([name,len])=>{
      const o=document.createElement('option'); o.value=name; o.textContent=`${name}  •  ${len} bp`; select.appendChild(o);
    });
    await refreshSelectedSequence();
    setStatus(`Uploaded ${Object.keys(dataset).length} read(s). Choose a read to analyze.`,true);
    await analyzeSelected();
  }catch(err){ setStatus(err.message,false); }
}

function render(d){
  $('length').textContent=d.sequence_length;
  $('edges').textContent=d.edges_count;
  $('unique').textContent=d.unique_kmers;
  $('nodes').textContent=d.nodes_count;
  $('maxWeight').textContent=d.max_weight;
  $('eulerian').textContent=d.eulerian?'Yes':'No';
  $('bestScore').textContent=d.best_score;
  $('bestEdges').textContent=d.best_edges_used;
  $('bestMethod').textContent=d.best_method;
  $('bestPath').textContent=(d.best_path_nodes||[]).join('  →  ') || '-';
  $('bestReconstruction').textContent=d.best_reconstruction || '-';
  $('eulerScore').textContent=d.euler_score;
  $('eulerEdges').textContent=d.euler_edges_used;
  $('eulerReconstruction').textContent=d.euler_reconstruction || '-';
  renderEdgeTable(d.all_edges||[]);
  renderFrequencyTable(d.common_kmers||[]);
  renderGraph(d.nodes||[],d.edges||[]);
}

function renderEdgeTable(edges){
  const body=$('edgeTable'); body.innerHTML='';
  edges.forEach(e=>{
    const tr=document.createElement('tr');
    tr.innerHTML=`<td>${e.id}</td><td class="mono">${escapeHtml(e.kmer)}</td><td class="mono">${escapeHtml(e.from)}</td><td class="mono">${escapeHtml(e.to)}</td><td><strong>${e.weight}</strong></td><td>${e.normalized_weight}</td>`;
    body.appendChild(tr);
  });
}
function renderFrequencyTable(items){
  const body=$('frequencyTable'); body.innerHTML='';
  items.forEach((x,i)=>{
    const tr=document.createElement('tr');
    tr.innerHTML=`<td>${i+1}</td><td class="mono">${escapeHtml(x.kmer)}</td><td><strong>${x.frequency}</strong></td>`;
    body.appendChild(tr);
  });
}

function renderGraph(nodes,edges){
  const box=$('graph'); box.innerHTML='';
  if(!nodes.length){box.textContent='No graph to display.';return;}

  // Clean, readable layered layout. Nodes are arranged in rows and parallel
  // edges are curved separately so that repeated k-mers remain visible.
  const width=1500;
  const nodeR=30;
  const gapX=150;
  const gapY=175;
  const marginX=90;
  const marginY=85;
  const perRow=Math.max(5, Math.min(9, Math.floor((width-2*marginX)/gapX)+1));
  const rows=Math.ceil(nodes.length/perRow);
  const height=Math.max(620, rows*gapY+150);
  const positions={};

  nodes.forEach((n,i)=>{
    const row=Math.floor(i/perRow);
    const col=i%perRow;
    const countInRow=Math.min(perRow,nodes.length-row*perRow);
    const usable=width-2*marginX;
    const step=countInRow>1 ? usable/(countInRow-1) : 0;
    positions[n.id]={
      x: countInRow===1 ? width/2 : marginX+col*step,
      y: marginY+row*gapY
    };
  });

  const pairCounts={};
  edges.forEach(e=>{
    const key=`${e.from}|||${e.to}`;
    pairCounts[key]=(pairCounts[key]||0)+1;
  });
  const pairSeen={};

  let svg=`<svg class="graph-svg" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" aria-label="Weighted De Bruijn graph">`;
  svg+=`<defs>
    <marker id="smallArrow" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="5" markerHeight="5" orient="auto">
      <path class="arrowhead" d="M0,0 L8,4 L0,8 Z"></path>
    </marker>
    <filter id="labelShadow" x="-20%" y="-20%" width="140%" height="140%">
      <feDropShadow dx="0" dy="1" stdDeviation="1.2" flood-opacity=".55"/>
    </filter>
  </defs>`;

  // Draw edges first so nodes stay visually on top of connections.
  edges.forEach((e,idx)=>{
    const a=positions[e.from], b=positions[e.to];
    if(!a||!b) return;

    const key=`${e.from}|||${e.to}`;
    const total=pairCounts[key];
    const seen=pairSeen[key]||0;
    pairSeen[key]=seen+1;
    const center=(total-1)/2;
    const offset=seen-center;

    const weight=Math.max(0,e.normalized_weight||0);
    const strokeWidth=1.8+3.4*weight;
    const high=weight>=0.65 ? ' high' : '';
    let path='';
    let labelX, labelY;

    if(e.from===e.to){
      const loopR=48+Math.abs(offset)*18;
      const lx=a.x+loopR;
      const ly=a.y-loopR;
      path=`M ${a.x+nodeR*0.72} ${a.y-nodeR*0.72}
            C ${lx} ${a.y-loopR*1.35}, ${lx} ${ly}, ${a.x+nodeR*0.72} ${a.y+nodeR*0.72}`;
      labelX=lx+5;
      labelY=ly-7;
    }else{
      const dx=b.x-a.x;
      const dy=b.y-a.y;
      const dist=Math.max(1,Math.hypot(dx,dy));
      const ux=dx/dist, uy=dy/dist;
      const sx=a.x+ux*nodeR;
      const sy=a.y+uy*nodeR;
      const ex=b.x-ux*(nodeR+7);
      const ey=b.y-uy*(nodeR+7);

      // Perpendicular bend separates parallel/nearby edges.
      const px=-uy, py=ux;
      const bend=offset*34;
      const mx=(sx+ex)/2+px*bend;
      const my=(sy+ey)/2+py*bend;
      path=`M ${sx} ${sy} Q ${mx} ${my} ${ex} ${ey}`;
      labelX=mx+px*8;
      labelY=my+py*8;
    }

    svg+=`<path class="edge-line${high}" d="${path}" stroke-width="${strokeWidth.toFixed(2)}" fill="none" marker-end="url(#smallArrow)" opacity="0.94"></path>`;

    // Compact label sitting on a dark pill keeps the k-mer/weight readable.
    const label=`${e.kmer} · w=${e.weight}`;
    const labelW=Math.max(54,label.length*7.1+18);
    svg+=`<g class="edge-label-group" transform="translate(${labelX.toFixed(1)},${labelY.toFixed(1)})" filter="url(#labelShadow)">
      <rect class="edge-label-bg" x="${(-labelW/2).toFixed(1)}" y="-12" width="${labelW.toFixed(1)}" height="22" rx="8"></rect>
      <text class="weight-label" x="0" y="4" text-anchor="middle">${escapeHtml(label)}</text>
    </g>`;
  });

  nodes.forEach(n=>{
    const p=positions[n.id];
    if(!p) return;
    svg+=`<g class="graph-node">
      <circle class="node-circle" cx="${p.x}" cy="${p.y}" r="${nodeR}"></circle>
      <text class="node-label" x="${p.x}" y="${p.y+5}">${escapeHtml(n.id.length>14?n.id.slice(0,14)+'…':n.id)}</text>
      <text class="node-degree" x="${p.x}" y="${p.y+50}">in ${n.in} · out ${n.out}</text>
    </g>`;
  });

  svg+='</svg>';
  box.innerHTML=svg;
}
function escapeHtml(s){return String(s).replace(/[&<>'"]/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;',"'":'&#39;','"':'&quot;'}[c]));}

$('sequence').addEventListener('change', async()=>{await refreshSelectedSequence(); await analyzeSelected();});
$('buildBtn').addEventListener('click', analyzeSelected);
$('customBtn').addEventListener('click', analyzeCustom);
$('uploadBtn').addEventListener('click', uploadDataset);
$('loadSelectedBtn').addEventListener('click', async()=>{
  const id=$('sequence').value;
  const r=await fetch(`/api/sequence-text?sequence_id=${encodeURIComponent(id)}`);
  const d=await r.json();
  $('customSequence').value=d.sequence || '';
  $('customSequence').focus();
});
$('clearCustomBtn').addEventListener('click',()=>{$('customSequence').value='';$('customSequence').focus();});
$('k').addEventListener('keydown',e=>{if(e.key==='Enter') analyzeSelected();});

loadSequences();
