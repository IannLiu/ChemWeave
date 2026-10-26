'use strict';

const $ = id => document.getElementById(id);
const escapeHTML = value => String(value ?? '').replace(/[&<>"']/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
const moleculeURL = (smiles, mapped = false) => `/api/molecule?smiles=${encodeURIComponent(smiles)}&mapped=${mapped}`;
const ruleNames = {alkane_dehydrogenation:'Alkane dehydrogenation',olefin_double_bond_shift:'Olefin double-bond shift',terminal_olefin_coupling:'Terminal-olefin coupling',olefin_scission:'Olefin scission'};
const reasons = {exhausted:'Network exhausted',max_depth:'Depth limit reached',targets_reached:'All targets found',max_species:'Species limit reached',max_reactions:'Reaction limit reached',user_stopped:'Stopped by user'};
let builtinRules = [], availableRules = [], network = null, speciesById = new Map(), reactionsById = new Map();
let selectedId = null, focusId = null, loadedRevision = -1, busy = false, polling = false;
let analysis = null, pathwayReactionIds = null;
const emptyInspector = $('details').innerHTML;

const cy = cytoscape({
  container: $('graph'), elements: [], minZoom: .12, maxZoom: 2.5, wheelSensitivity: .2,
  style: [
    {selector:'node',style:{'label':'data(label)','color':'#536a56','font-family':'-apple-system, sans-serif','font-size':11,'text-wrap':'wrap','text-max-width':120,'text-valign':'bottom','text-margin-y':5,'border-width':1,'border-color':'#cfdbc8','background-color':'#ffffff','overlay-opacity':0}},
    {selector:'node[kind="species"]',style:{'shape':'round-rectangle','width':105,'height':72,'background-image':'data(image)','background-fit':'contain','background-width':'94%','background-height':'92%','background-clip':'none'}},
    {selector:'node[kind="reaction"]',style:{'shape':'diamond','width':18,'height':18,'background-color':'#8fa583','border-color':'#8fa583','font-size':8,'color':'#8a9d7d','text-margin-y':6}},
    {selector:'.seed',style:{'border-width':2,'border-color':'#5a947c','background-color':'#f4f9f4'}},
    {selector:'.target',style:{'border-width':2,'border-color':'#c6a352','background-color':'#fffdf3'}},
    {selector:'edge',style:{'width':1.15,'line-color':'#bacbb0','target-arrow-color':'#a9bd9c','target-arrow-shape':'triangle','arrow-scale':.65,'curve-style':'bezier','label':'data(label)','font-size':9,'color':'#839577','text-background-color':'#f7f9f6','text-background-opacity':1,'text-background-padding':2}},
    {selector:':selected',style:{'border-color':'#235f56','border-width':3,'background-color':'#f0f5e9'}},
    {selector:'.dim',style:{'opacity':.16}},
    {selector:'.search-hit',style:{'border-width':3,'border-color':'#ce9a36'}}
  ]
});
new ResizeObserver(() => {
  cy.resize();
  if (network) {
    cy.fit(cy.elements(),55);
    if (cy.zoom() > 1.25) { cy.zoom(1.25); cy.center(); }
  }
}).observe($('graph'));
cy.on('tap','node', event => selectNode(event.target.id()));

async function api(path, options = {}) {
  const response = await fetch(path, {...options,headers:{'Content-Type':'application/json',...options.headers}});
  if (!response.ok) {
    const error = await response.json().catch(() => ({}));
    const detail = Array.isArray(error.detail) ? error.detail.map(e => `${e.loc.slice(1).join('.')}: ${e.msg}`).join('\n') : error.detail;
    throw new Error(detail || `Request failed (${response.status})`);
  }
  return response.json();
}

function message(text = '', success = false) {
  $('message').textContent = text;
  $('message').hidden = !text;
  $('message').classList.toggle('success', success);
}

function addSpecies(containerId, record = {}) {
  const row = document.createElement('div');
  row.className = 'species-row';
  row.innerHTML = `<div class="species-main"><img alt="Molecular structure" hidden><div class="species-fields"><input data-field="smiles" aria-label="${containerId === 'seeds' ? 'Seed' : 'Target'} species SMILES" placeholder="SMILES" required spellcheck="false"><input data-field="name" aria-label="Species name" placeholder="Name (optional)"></div><button class="remove-species" type="button" aria-label="Remove species">×</button></div><details class="row-state"><summary>Phase and site</summary><div><select data-field="phase" aria-label="Phase"><option value="gas">Gas</option><option value="surface">Surface</option><option value="liquid">Liquid</option></select><input data-field="site" aria-label="Site" placeholder="Site (optional)"></div></details>`;
  for (const key of ['smiles','name','phase','site']) {
    const input = row.querySelector(`[data-field="${key}"]`);
    if (record[key] != null) {
      if (key === 'phase' && !Array.from(input.options).some(o => o.value === record[key])) input.add(new Option(record[key],record[key]));
      input.value = record[key];
    }
  }
  row.querySelector('.remove-species').onclick = () => row.remove();
  const smilesInput = row.querySelector('[data-field="smiles"]');
  smilesInput.addEventListener('input', () => { row.classList.remove('invalid'); smilesInput.setCustomValidity(''); row.querySelector('img').hidden = true; });
  smilesInput.addEventListener('blur', async () => {
    const smiles = smilesInput.value.trim();
    if (!smiles) return;
    try {
      const result = await api('/api/preview',{method:'POST',body:JSON.stringify({smiles})});
      if (smilesInput.value.trim() !== smiles) return;
      row.querySelector('img').src = moleculeURL(result.smiles);
      row.querySelector('img').hidden = false;
      row.classList.remove('invalid');
      smilesInput.setCustomValidity('');
    } catch (error) {
      if (smilesInput.value.trim() !== smiles) return;
      row.classList.add('invalid');
      smilesInput.setCustomValidity(error.message);
    }
  });
  $(containerId).append(row);
  if (record.smiles) smilesInput.dispatchEvent(new Event('blur'));
  else smilesInput.focus();
}

function renderRules(rules, label = 'Built-in demo rules') {
  availableRules = rules;
  $('rule-set-name').textContent = label;
  $('reset-rules').hidden = label === 'Built-in demo rules';
  $('rules').innerHTML = rules.map((rule,index) => {
    const sides = rule.smarts.split('>>');
    return `<div class="rule-item"><label><input type="checkbox" data-rule-index="${index}" checked><span>${escapeHTML(ruleNames[rule.id] || rule.id)}</span><span class="arity">${sides[0].split('.').length} → ${sides[1]?.split('.').length ?? '?'}</span></label><details><summary>View rule</summary><pre>${escapeHTML(rule.smarts)}</pre><p>${escapeHTML(rule.notes || rule.family || '')}</p></details></div>`;
  }).join('');
}

function fillSettings(request) {
  $('seeds').replaceChildren();
  $('targets').replaceChildren();
  for (const record of request.seeds) addSpecies('seeds',record);
  for (const record of request.targets || []) addSpecies('targets',record);
  const config = request.config || {};
  for (const [id,key,fallback] of [['max-carbon','max_carbon',6],['max-depth','max_depth',null],['max-species','max_species',200],['max-reactions','max_reactions',1000]]) $(id).value = config[key] ?? fallback ?? '';
  $('stop-on-targets').checked = config.stop_when_targets_found ?? true;
  renderRules(request.rules || builtinRules,request.rules ? 'Project rules' : 'Built-in demo rules');
}

function setBusy(value) {
  busy = value;
  $('configuration-fields').disabled = value;
  $('pathway-fields').disabled = value;
  $('pathway-view').disabled = value;
  $('export-pathways').disabled = value;
  $('run').hidden = value;
  $('stop').hidden = !value;
  $('stop').disabled = false;
  for (const id of ['new-project','open-project','run-example']) $(id).disabled = value;
  $('graph-busy').hidden = !value;
  for (const id of ['save-project','export','export-format']) $(id).disabled = value || !network;
}

async function startGeneration(event) {
  event?.preventDefault();
  if (busy || !$('settings').reportValidity()) return;
  const records = id => Array.from($(id).children).map(row => Object.fromEntries(Array.from(row.querySelectorAll('[data-field]')).map(input => [input.dataset.field,input.value.trim() || (input.dataset.field === 'site' ? null : '')])));
  const seeds = records('seeds');
  const rules = Array.from($('rules').querySelectorAll('input:checked')).map(input => availableRules[Number(input.dataset.ruleIndex)]);
  if (!seeds.length) return message('Add at least one seed species.');
  if (!rules.length) return message('Select at least one reaction rule.');
  const request = {seeds,targets:records('targets'),rules,config:{stop_when_targets_found:$('stop-on-targets').checked,max_carbon:Number($('max-carbon').value),max_depth:$('max-depth').value === '' ? null : Number($('max-depth').value),max_species:Number($('max-species').value),max_reactions:Number($('max-reactions').value)}};
  message();
  setBusy(true);
  try {
    await api('/api/run',{method:'POST',body:JSON.stringify(request)});
    network = null;
    analysis = pathwayReactionIds = null;
    $('pathway-panel').hidden = true;
    selectedId = focusId = null;
    cy.elements().remove();
    $('details').innerHTML = emptyInspector;
    $('empty-state').hidden = true;
    $('network-count').textContent = 'Generating';
    loadedRevision = -1;
  } catch (error) { setBusy(false); message(error.message); }
}

function loadNetwork(payload, restoreSettings = false) {
  network = payload.network;
  speciesById = new Map(network.species.map(s => [s.id,s]));
  reactionsById = new Map(network.reactions.map(r => [r.id,r]));
  selectedId = focusId = null;
  analysis = payload.analysis || null;
  pathwayReactionIds = null;
  $('search').value = '';
  $('details').innerHTML = emptyInspector;
  const maxDepth = Math.max(0,...network.species.map(s => s.depth ?? 0));
  $('depth-filter').innerHTML = '<option value="all">All depths</option>' + Array.from({length:maxDepth+1},(_,i) => `<option value="${i}">Depth ≤ ${i}</option>`).join('');
  $('depth-filter').disabled = network.species.every(s => s.depth == null);
  $('depth-filter').value = network.species.length > 80 ? String(Math.min(2,maxDepth)) : 'all';
  $('network-count').textContent = `${network.species.length} species · ${network.reactions.length} reactions`;
  $('empty-state').hidden = true;
  if (restoreSettings && payload.request) fillSettings(payload.request);
  renderPathways();
  if (!analysis) renderGraph();
  if (network.targets.length) {
    const found = network.targets.filter(t => t.species_id).length;
    $('details').innerHTML = `<div class="detail-content"><div class="target-progress">Targets found: ${found} / ${network.targets.length}<ul>${network.targets.map(t => `<li>${t.species_id ? '✓' : '○'} ${escapeHTML(t.name || t.smiles)}${t.phase !== 'gas' ? ' · '+escapeHTML(t.phase) : ''}${t.site ? ' · '+escapeHTML(t.site) : ''}</li>`).join('')}</ul></div>${emptyInspector}</div>`;
  }
  setBusy(false);
}

function renderGraph() {
  if (!network) return;
  const bound = $('depth-filter').value === 'all' ? Infinity : Number($('depth-filter').value);
  let species = network.species.filter(s => s.depth == null || s.depth <= bound);
  const allowed = new Set(species.map(s => s.id));
  let reactions = network.reactions.filter(r => (r.depth == null || r.depth <= bound) && [...r.reactants,...r.products].every(([sid]) => allowed.has(sid)));
  if (pathwayReactionIds !== null) {
    reactions = network.reactions.filter(r => pathwayReactionIds.has(r.id));
    const included = new Set([analysis.target, ...reactions.flatMap(r => [...r.reactants,...r.products].map(([sid]) => sid))]);
    species = network.species.filter(s => included.has(s.id));
  }
  if (focusId) {
    reactions = reactions.filter(r => r.id === focusId || [...r.reactants,...r.products].some(([sid]) => sid === focusId));
    const nearby = new Set(reactions.flatMap(r => [...r.reactants,...r.products].map(([sid]) => sid)));
    nearby.add(focusId);
    species = species.filter(s => nearby.has(s.id));
  }
  const targets = new Set(network.targets.map(t => t.species_id).filter(Boolean));
  if (pathwayReactionIds !== null) targets.add(analysis.target);
  const nodes = species.map(s => ({data:{id:s.id,kind:'species',label:`${s.name || s.formula}\n${s.id}`,image:moleculeURL(s.smiles),column:(s.depth ?? 0)*2},classes:[network.seeds.includes(s.id) ? 'seed' : '',targets.has(s.id) ? 'target' : ''].join(' ')}));
  const edges = [];
  for (const reaction of reactions) {
    nodes.push({data:{id:reaction.id,kind:'reaction',label:reaction.id,column:Math.max(1,(reaction.depth ?? 1)*2-1)}});
    for (const [side,terms] of [['left',reaction.reactants],['right',reaction.products]]) for (const [sid,n] of terms) edges.push({data:{id:`${reaction.id}-${side}-${sid}`,source:side === 'left' ? sid : reaction.id,target:side === 'left' ? reaction.id : sid,label:n > 1 ? String(n) : ''}});
  }
  const columns = new Map();
  for (const node of nodes) {
    const col = node.data.column;
    if (!columns.has(col)) columns.set(col,[]);
    columns.get(col).push(node);
  }
  for (const [column,items] of columns) items.forEach((node,index) => { node.position = {x:column*105,y:(index-(items.length-1)/2)*125}; });
  cy.elements().remove();
  cy.add([...nodes,...edges]);
  cy.layout({name:network.species.every(s => s.depth == null) ? 'cose' : 'preset',animate:false,padding:65}).run();
  cy.fit(cy.elements(),65);
  if (cy.zoom() > 1.25) { cy.zoom(1.25); cy.center(); }
  if (selectedId && cy.getElementById(selectedId).length) cy.getElementById(selectedId).select();
  $('visible-count').textContent = `Showing ${species.length} / ${network.species.length} species · ${reactions.length} / ${network.reactions.length} reactions`;
  $('depth-filter').disabled = pathwayReactionIds !== null || network.species.every(s => s.depth == null);
  $('view-note').textContent = pathwayReactionIds !== null ? 'Formation pathways · Complete reactions' : network.species.every(s => s.depth == null) ? 'Depth is not recorded in this legacy file' : (focusId ? 'Neighborhood view · Complete reactions' : 'Arranged by minimum depth');
  $('neighbors').classList.toggle('active',Boolean(focusId));
  searchGraph();
}

function searchGraph() {
  const query = $('search').value.trim().toLowerCase();
  cy.elements().removeClass('dim search-hit');
  if (!query || !network) return;
  const hits = cy.nodes().filter(node => {
    const s = speciesById.get(node.id());
    return s && [s.id,s.name,s.smiles,s.formula].some(value => value.toLowerCase().includes(query));
  });
  cy.elements().addClass('dim');
  hits.add(hits.neighborhood()).removeClass('dim');
  hits.addClass('search-hit');
}

function selectNode(id) {
  selectedId = id;
  cy.nodes().unselect();
  cy.getElementById(id).select();
  $('neighbors').disabled = false;
  if (speciesById.has(id)) {
    const s = speciesById.get(id);
    const related = network.reactions.filter(r => [...r.reactants,...r.products].some(([sid]) => sid === id));
    $('details').innerHTML = `<div class="detail-content"><span class="detail-kicker">SPECIES · ${escapeHTML(s.id)}</span><h3>${escapeHTML(s.name || s.formula)}</h3>${network.seeds.includes(id) ? '<span class="detail-tag">Seed species</span>' : ''}${network.targets.some(t => t.species_id === id) ? '<span class="detail-tag target">Target species</span>' : ''}<img class="structure-image" src="${escapeHTML(moleculeURL(s.smiles))}" alt="${escapeHTML(s.smiles)}"><label class="map-toggle"><input id="species-map" type="checkbox">Show species atom-map numbers</label><dl><dt>SMILES</dt><dd><code>${escapeHTML(s.smiles)}</code></dd><dt>Formula / minimum depth</dt><dd>${escapeHTML(s.formula)} · ${s.depth == null ? 'Unknown' : 'Depth '+s.depth}</dd><dt>Phase / site</dt><dd>${escapeHTML(s.phase)} / ${escapeHTML(s.site || '—')}</dd></dl><button id="analyze-species" type="button" class="secondary">Analyze formation</button><div class="detail-label"><h4>Related reactions</h4><span class="optional">${related.length} reactions</span></div><div class="related-list">${related.map(r => `<button type="button" data-related="${escapeHTML(r.id)}">${escapeHTML(r.id)} · ${escapeHTML(ruleNames[r.rule_id] || r.rule_id)}</button>`).join('') || '<p class="detail-note">No related reactions in this network.</p>'}</div><p class="detail-note">These atom-map numbers belong to the species representative. Inspect a reaction to see its single-step mapping.</p></div>`;
    $('analyze-species').disabled = busy;
    $('analyze-species').onclick = () => { $('pathway-target').value = id; $('pathway-settings').requestSubmit(); };
    $('species-map').onchange = event => { $('details').querySelector('.structure-image').src = moleculeURL(event.target.checked ? s.mapped_smiles : s.smiles,event.target.checked); };
    $('details').querySelectorAll('[data-related]').forEach(button => { button.onclick = () => selectNode(button.dataset.related); });
  } else if (reactionsById.has(id)) showReaction(reactionsById.get(id));
}

function showReaction(reaction, channelIndex = 0, mapped = false) {
  const channels = reaction.channels || [];
  const channel = channels[channelIndex];
  const sides = channel ? [channel.reactants,channel.products] : [reaction.reactants,reaction.products].map(terms => terms.flatMap(([sid,n]) => Array(n).fill(sid)));
  const cards = (ids,side) => ids.map((sid,index) => {
    const s = speciesById.get(sid);
    const smiles = mapped && channel ? channel[side === 0 ? 'mapped_reactant_smiles' : 'mapped_product_smiles'][index] : s.smiles;
    return `<div class="reaction-molecule"><img src="${escapeHTML(moleculeURL(smiles,mapped))}" alt="${escapeHTML(smiles)}"><small>${escapeHTML(s.name || s.formula)} · ${escapeHTML(sid)}</small></div>`;
  }).join('');
  $('details').innerHTML = `<div class="detail-content"><span class="detail-kicker">REACTION · ${escapeHTML(reaction.id)}</span><h3>${escapeHTML(ruleNames[reaction.rule_id] || reaction.rule_id)}</h3><span class="detail-tag">Depth ${reaction.depth ?? 'Unknown'}</span><span class="detail-tag">${reaction.kind === 'overall' ? 'Overall transformation' : escapeHTML(reaction.kind)}</span>${channels.length ? `<select id="channel" class="channel-select" aria-label="Mapping channel">${channels.map((c,i) => `<option value="${i}" ${i === channelIndex ? 'selected' : ''}>Channel ${i+1} · ${c.match_count} template matches</option>`).join('')}</select><label class="map-toggle"><input id="reaction-map" type="checkbox" ${mapped ? 'checked' : ''}>Show reaction atom-map numbers</label>` : '<p class="detail-note">This legacy record has no complete mapping channels.</p>'}<div class="reaction-side">${cards(sides[0],0)}</div><div class="reaction-arrow">↓</div><div class="reaction-side">${cards(sides[1],1)}</div><dl><dt>Rule ID / reaction family</dt><dd>${escapeHTML(reaction.rule_id)}<br>${escapeHTML(reaction.family)}</dd></dl>${mapped && channel ? `<details><summary>Mapped reaction SMILES</summary><code>${escapeHTML(channel.mapped_reaction_smiles)}</code></details>` : ''}<p class="detail-note">Each channel stores one representative of symmetry-equivalent mappings. Template-match counts are not kinetic degeneracy factors and may be incomplete after an early stop.</p></div>`;
  if (channels.length) {
    $('channel').onchange = event => showReaction(reaction,Number(event.target.value),$('reaction-map').checked);
    $('reaction-map').onchange = event => showReaction(reaction,Number($('channel').value),event.target.checked);
  }
}

async function pollStatus() {
  if (polling) return;
  polling = true;
  try {
    const state = await api('/api/status');
    const running = ['running','analyzing','stopping'].includes(state.status);
    if (running !== busy) setBusy(running);
    $('stop').disabled = state.status === 'stopping';
    $('stop').textContent = state.status === 'stopping' ? 'Stopping…' : 'Stop';
    $('graph-busy').querySelector('strong').textContent = state.status === 'analyzing' ? 'Analyzing formation pathways' : 'Generating your reaction network';
    $('graph-busy').querySelector('p').textContent = state.status === 'analyzing' ? 'You can stop the search and keep the pathways already found.' : 'Statistics update during the run. The graph appears when it finishes.';
    const s = state.summary;
    $('status-text').textContent = state.status === 'analyzing' ? 'Analyzing pathways' : state.status === 'running' ? 'Generating' : state.status === 'stopping' ? 'Waiting for the current calculation' : state.status === 'error' ? 'Calculation failed' : state.status === 'finished' ? (state.analysis_termination ? ({exhausted:'Pathway analysis complete',max_pathways:'Pathway limit reached',user_stopped:'Pathway analysis stopped'}[state.analysis_termination]) : (reasons[s.termination] || s.termination)) : 'Ready';
    $('status-dot').className = running ? 'running' : state.status === 'error' ? 'error' : '';
    for (const [id,key] of [['stat-species','species'],['stat-reactions','reactions'],['stat-depth','depth_reached']]) $(id).textContent = s[key] ?? '—';
    $('stat-targets').textContent = s.targets_total ? `${s.targets_found}/${s.targets_total}` : '—';
    $('elapsed').textContent = state.status === 'idle' ? 'All calculations run locally' : `Elapsed ${state.elapsed.toFixed(1)} s · Local computation`;
    if (state.status === 'finished' && state.revision !== loadedRevision) {
      const payload = await api('/api/network');
      loadNetwork(payload,true);
      loadedRevision = state.revision;
      if (s.missing_targets?.length) message(`Run finished with ${s.missing_targets.length} target(s) still missing. Adjust the depth, search limits or rules and run again.`);
    } else if (state.status === 'error' && state.revision !== loadedRevision) {
      loadedRevision = state.revision;
      message(state.error || 'Calculation failed');
      $('empty-state').hidden = Boolean(network);
      $('graph-busy').hidden = true;
    }
  } catch (error) { message(`Cannot connect to the local service: ${error.message}`); }
  finally { polling = false; }
}

async function startPathwayAnalysis(event) {
  event.preventDefault();
  if (busy || !$('pathway-settings').reportValidity()) return;
  const request = {target:$('pathway-target').value,max_reactions:Number($('pathway-max-reactions').value),max_pathways:Number($('pathway-max-count').value)};
  message();
  setBusy(true);
  try {
    await api('/api/pathways',{method:'POST',body:JSON.stringify(request)});
    loadedRevision = -1;
    await pollStatus();
  } catch (error) { setBusy(false); message(error.message); }
}

function renderPathways() {
  $('pathway-panel').hidden = false;
  const preferred = analysis?.target || network.targets.find(t => t.species_id)?.species_id;
  $('pathway-target').replaceChildren(...network.species.map(s => new Option(`${s.name || s.smiles} · ${s.id} · ${s.phase}${s.site ? ' / '+s.site : ''}`,s.id)));
  $('pathway-target').value = preferred || network.species.find(s => !network.seeds.includes(s.id))?.id || network.seeds[0];
  $('pathway-results').hidden = !analysis;
  $('pathway-steps').replaceChildren();
  if (!analysis) return;
  $('pathway-max-reactions').value = analysis.max_reactions;
  $('pathway-max-count').value = analysis.max_pathways;
  const target = speciesById.get(analysis.target);
  $('pathway-summary').textContent = `${analysis.pathways.length} pathway(s) to ${target.name || target.smiles} · up to ${analysis.max_reactions} reactions per route.`;
  const notices = [];
  if (!analysis.network_search_complete) notices.push(`Network search ended early (${reasons[analysis.network_termination] || analysis.network_termination}); alternative routes may be missing.`);
  if (analysis.termination === 'max_pathways') notices.push('Pathway limit reached. Enumeration is not confirmed complete; results are not ranked.');
  else if (analysis.termination === 'user_stopped') notices.push('Analysis stopped. The pathways found so far are retained.');
  else notices.push('Enumeration complete within this stored network and reaction-count limit.');
  if (!analysis.pathways.length) notices.push('No route found within these limits.');
  $('pathway-warning').textContent = notices.join(' ');
  $('pathway-view').replaceChildren(new Option('Full network','network'),new Option('Union of found pathways','all'),...analysis.pathways.map((p,i) => new Option(`Pathway ${i+1} · ${p.reactions.length} reactions`,String(i))));
  $('pathway-view').value = 'all';
  showPathway();
}

function showPathway() {
  const view = $('pathway-view').value;
  focusId = null;
  $('search').value = '';
  $('pathway-steps').replaceChildren();
  const chosen = view === 'network' ? [] : view === 'all' ? analysis.pathways : [analysis.pathways[Number(view)]];
  pathwayReactionIds = view === 'network' ? null : new Set(chosen.flatMap(p => p.reactions));
  if (view !== 'network' && view !== 'all') {
    const path = chosen[0];
    const side = terms => terms.map(([sid,n]) => `${n > 1 ? n+' ' : ''}${speciesById.get(sid).name || speciesById.get(sid).smiles}`).join(' + ');
    for (const rid of path.reactions) {
      const r = reactionsById.get(rid);
      const item = document.createElement('li');
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = `${rid}: ${side(r.reactants)} → ${side(r.products)}`;
      button.onclick = () => selectNode(rid);
      item.append(button);
      $('pathway-steps').append(item);
    }
    if (!path.reactions.length) $('pathway-steps').textContent = 'The target is already a seed species; no reaction is required.';
  }
  renderGraph();
}

$('pathway-settings').addEventListener('submit',startPathwayAnalysis);
$('pathway-view').onchange = showPathway;
$('export-pathways').onclick = () => { window.location.href = '/api/export/pathways'; };

$('settings').addEventListener('submit',startGeneration);
$('run-example').onclick = () => $('settings').requestSubmit();
$('add-seed').onclick = () => addSpecies('seeds');
$('add-target').onclick = () => addSpecies('targets');
$('add-bulk').onclick = () => { for (const smiles of $('bulk-smiles').value.split(/\r?\n/).map(s => s.trim()).filter(Boolean)) addSpecies('seeds',{smiles}); $('bulk-smiles').value = ''; };
$('stop').onclick = async () => { try { await api('/api/stop',{method:'POST'}); await pollStatus(); } catch (error) { message(error.message); } };
$('search').addEventListener('input',searchGraph);
$('depth-filter').onchange = renderGraph;
$('neighbors').onclick = () => { focusId = focusId ? null : selectedId; renderGraph(); };
$('show-all').onclick = () => { pathwayReactionIds = null; $('pathway-view').value = 'network'; $('pathway-steps').replaceChildren(); focusId = null; $('depth-filter').value = 'all'; $('search').value = ''; renderGraph(); };
$('zoom-in').onclick = () => cy.zoom({level:Math.min(cy.zoom()*1.25,2.5),renderedPosition:{x:cy.width()/2,y:cy.height()/2}});
$('zoom-out').onclick = () => cy.zoom({level:Math.max(cy.zoom()/1.25,.12),renderedPosition:{x:cy.width()/2,y:cy.height()/2}});
$('fit').onclick = () => cy.fit(cy.elements(),65);
$('save-project').onclick = () => { window.location.href = '/api/export/project'; };
$('export').onclick = () => {
  const format = $('export-format').value;
  if (format === 'png') { const a = document.createElement('a'); a.href = cy.png({full:true,bg:'#f7f9f6',scale:2}); a.download = 'chemweave-network.png'; a.click(); }
  else window.location.href = `/api/export/${format}`;
};
$('open-project').onclick = () => $('project-file').click();
$('project-file').onchange = async event => {
  const file = event.target.files[0];
  if (!file) return;
  try {
    const payload = await api('/api/open',{method:'POST',body:await file.text()});
    const state = await api('/api/status');
    loadedRevision = state.revision;
    loadNetwork(payload,true);
    message(`Opened ${file.name}`,true);
    await pollStatus();
  } catch (error) { message(`Could not open file: ${error.message}`); }
  event.target.value = '';
};
$('import-rules').onclick = () => $('rules-file').click();
$('rules-file').onchange = async event => {
  const file = event.target.files[0];
  if (!file) return;
  try {
    const rules = JSON.parse(await file.text());
    if (!Array.isArray(rules) || !rules.length || rules.some(r => !r.id || !r.smarts)) throw new Error('Select a JSON array of rules containing id and smarts fields.');
    renderRules(rules,file.name);
    message();
  } catch (error) { message(`Could not import rules: ${error.message}`); }
  event.target.value = '';
};
$('reset-rules').onclick = () => renderRules(builtinRules);
$('new-project').onclick = () => $('new-dialog').showModal();
$('cancel-new').onclick = () => $('new-dialog').close();
$('confirm-new').onclick = async () => {
  try { await api('/api/new',{method:'POST'}); window.location.reload(); }
  catch (error) { $('new-dialog').close(); message(error.message); }
};

async function initialize() {
  try {
    const result = await api('/api/rules');
    builtinRules = result.rules;
    fillSettings({seeds:[{smiles:'CCC',name:'Propane'}],config:{max_carbon:6,max_depth:4,max_species:200,max_reactions:1000}});
    await pollStatus();
    window.setInterval(pollStatus,700);
  } catch (error) { message(error.message); }
}
initialize();
