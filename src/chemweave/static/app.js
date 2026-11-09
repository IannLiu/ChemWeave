'use strict';

const $ = id => document.getElementById(id);
const escapeHTML = value => String(value ?? '').replace(/[&<>"']/g, char => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[char]));
const moleculeURL = (smiles, mapped = false) => `/api/molecule?smiles=${encodeURIComponent(smiles)}&mapped=${mapped}`;
const ruleNames = {user_metal_H2_adsorption:'H₂ adsorption on metal',methane_h_abstraction:'Methane H/D atom abstraction',alkane_dehydrogenation:'Alkane dehydrogenation',olefin_double_bond_shift:'Olefin double-bond shift',terminal_olefin_coupling:'Terminal-olefin coupling',olefin_scission:'Olefin scission'};
const reasons = {exhausted:'Network exhausted',max_depth:'Depth limit reached',targets_reached:'All targets found',max_species:'Species limit reached',max_reactions:'Reaction limit reached',user_stopped:'Stopped by user'};
let builtinExample = null, builtinRules = [], availableRules = [], network = null, speciesById = new Map(), reactionsById = new Map();
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
  row.innerHTML = `<div class="species-main"><img alt="Molecular structure" hidden><div class="species-fields"><input data-field="smiles" aria-label="${containerId === 'seeds' ? 'Seed' : 'Target'} species SMILES" placeholder="SMILES" required spellcheck="false"><input data-field="name" aria-label="Species name" placeholder="Name (optional)"></div><button class="remove-species" type="button" aria-label="Remove species">×</button></div><details class="row-state"><summary>Phase and catalyst sites</summary><div><select data-field="phase" aria-label="Phase"><option value="gas">Gas</option><option value="surface">Surface</option><option value="liquid">Liquid</option></select><input data-field="site" aria-label="Site" placeholder="Legacy site label (optional)"></div><label>Explicit site atoms (JSON by map number)<textarea data-field="site_atoms" rows="3" aria-label="Explicit site atoms" spellcheck="false" placeholder='{"1":{"kind":"metal","type":"Fe_model","capacity":1}}'></textarea></label><p class="field-hint">Use [*:1] in species SMILES, [#0:1] in reaction SMARTS. Brønsted model: [*:1]O with kind bronsted.</p></details>`;
  for (const key of ['smiles','name','phase','site','site_atoms']) {
    const input = row.querySelector(`[data-field="${key}"]`);
    if (record[key] != null) {
      if (key === 'phase' && !Array.from(input.options).some(o => o.value === record[key])) input.add(new Option(record[key],record[key]));
      input.value = key === 'site_atoms' ? JSON.stringify(record[key],null,2) : record[key];
    }
  }
  row.querySelector('.remove-species').onclick = () => row.remove();
  const smilesInput = row.querySelector('[data-field="smiles"]');
  smilesInput.addEventListener('input', () => { row.classList.remove('invalid'); smilesInput.setCustomValidity(''); row.querySelector('img').hidden = true; });
  smilesInput.addEventListener('blur', async () => {
    const smiles = smilesInput.value.trim();
    if (!smiles) return;
    try {
      const result = await api('/api/preview',{method:'POST',body:JSON.stringify({smiles,
        phase:row.querySelector('[data-field="phase"]').value,
        site:row.querySelector('[data-field="site"]').value.trim() || null,
        site_atoms:JSON.parse(row.querySelector('[data-field="site_atoms"]').value || '{}')})});
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
  for (const field of ['phase','site','site_atoms']) row.querySelector(`[data-field="${field}"]`).addEventListener('change', () => {
    smilesInput.setCustomValidity('');
    smilesInput.dispatchEvent(new Event('blur'));
  });
  $(containerId).append(row);
  if (record.smiles) smilesInput.dispatchEvent(new Event('blur'));
  else smilesInput.focus();
}

function renderRules(rules, label = 'Built-in demo rules') {
  availableRules = rules;
  $('custom-rules').value = JSON.stringify(rules,null,2);
  $('rule-set-name').textContent = rules.length ? label : 'No templates loaded';
  $('reset-rules').hidden = !rules.length || label === 'Built-in demo rules';
  $('rules').innerHTML = rules.map((rule,index) => {
    const sides = rule.smarts.split('>>');
    const arities = sides.map(side => {
      let depth = 0, count = 1;
      for (const char of side) {
        if (char === '(') depth++;
        else if (char === ')') depth--;
        else if (char === '.' && depth === 0) count++;
      }
      return count;
    });
    return `<div class="rule-item"><label><input type="checkbox" data-rule-index="${index}" checked><span>${escapeHTML(ruleNames[rule.id] || rule.id)}</span><span class="arity">${arities[0]} → ${arities[1] ?? '?'}</span></label><details><summary>View rule</summary><pre>${escapeHTML(rule.smarts)}</pre><p>${escapeHTML(rule.notes || rule.family || '')}</p></details></div>`;
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
  $('single-event-enabled').checked = request.single_event != null;
  $('single-event-options').value = JSON.stringify(request.single_event || {calculator:'mapped_event_v1'},null,2);
  renderRules(request.rules || [], 'Project rules');
}

function setBusy(value) {
  busy = value;
  $('configuration-fields').disabled = value;
  $('pathway-fields').disabled = value;
  $('apply-single-event').disabled = value || !network;
  $('pathway-view').disabled = value;
  $('export-pathways').disabled = value;
  $('run').hidden = value;
  $('stop').hidden = !value;
  $('stop').disabled = false;
  for (const id of ['new-project','load-example','open-project','run-example']) $(id).disabled = value;
  $('graph-busy').hidden = !value;
  for (const id of ['save-project','export','export-format']) $(id).disabled = value || !network;
}

async function startGeneration(event) {
  event?.preventDefault();
  if (busy || !$('settings').reportValidity()) return;
  const records = id => Array.from($(id).children).map(row => Object.fromEntries(Array.from(row.querySelectorAll('[data-field]')).map(input => [input.dataset.field,input.dataset.field === 'site_atoms' ? JSON.parse(input.value || '{}') : input.value.trim() || (input.dataset.field === 'site' ? null : '')])));
  let seeds, targets;
  try { seeds = records('seeds'); targets = records('targets'); }
  catch (error) { return message(`Invalid site attributes JSON: ${error.message}`); }
  const rules = Array.from($('rules').querySelectorAll('input:checked')).map(input => availableRules[Number(input.dataset.ruleIndex)]);
  if (!seeds.length) return message('Add at least one seed species.');
  if (!rules.length) return message('Select at least one reaction rule.');
  const request = {seeds,targets,rules,config:{stop_when_targets_found:$('stop-on-targets').checked,max_carbon:Number($('max-carbon').value),max_depth:$('max-depth').value === '' ? null : Number($('max-depth').value),max_species:Number($('max-species').value),max_reactions:Number($('max-reactions').value)}};
  try { request.single_event = $('single-event-enabled').checked ? JSON.parse($('single-event-options').value) : null; }
  catch (error) { return message(`Invalid single-event JSON: ${error.message}`); }
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
  $('apply-single-event').disabled = busy;
  const seChannels = network.reactions.flatMap(r => r.channels.map(c => c.single_event)).filter(Boolean);
  $('single-event-summary').textContent = network.single_event ? `${seChannels.filter(c => ['calculated','imported'].includes(c.status)).length} annotated / ${seChannels.length} mapped channels · ${Object.keys(network.single_event.parameters).length} shared parameter classes · ${network.single_event.unmapped_reactions.length} unmapped reactions` : '';

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
    $('details').innerHTML = `<div class="detail-content"><span class="detail-kicker">SPECIES · ${escapeHTML(s.id)}</span><h3>${escapeHTML(s.name || s.formula)}</h3>${network.seeds.includes(id) ? '<span class="detail-tag">Seed species</span>' : ''}${network.targets.some(t => t.species_id === id) ? '<span class="detail-tag target">Target species</span>' : ''}<img class="structure-image" src="${escapeHTML(moleculeURL(s.smiles))}" alt="${escapeHTML(s.smiles)}"><label class="map-toggle"><input id="species-map" type="checkbox">Show species atom-map numbers</label><dl><dt>SMILES</dt><dd><code>${escapeHTML(s.smiles)}</code></dd><dt>Formula / minimum depth</dt><dd>${escapeHTML(s.formula)} · ${s.depth == null ? 'Unknown' : 'Depth '+s.depth}</dd><dt>Phase / site</dt><dd>${escapeHTML(s.phase)} / ${escapeHTML(s.site || "—")}</dd><dt>Explicit catalyst nodes</dt><dd><pre>${escapeHTML(JSON.stringify(s.site_atoms || {},null,2))}</pre></dd></dl><button id="analyze-species" type="button" class="secondary">Analyze formation</button><div class="detail-label"><h4>Related reactions</h4><span class="optional">${related.length} reactions</span></div><div class="related-list">${related.map(r => `<button type="button" data-related="${escapeHTML(r.id)}">${escapeHTML(r.id)} · ${escapeHTML(ruleNames[r.rule_id] || r.rule_id)}</button>`).join('') || '<p class="detail-note">No related reactions in this network.</p>'}</div><p class="detail-note">These atom-map numbers belong to the species representative. Inspect a reaction to see its single-step mapping.</p></div>`;
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
  $('details').innerHTML = `<div class="detail-content"><span class="detail-kicker">REACTION · ${escapeHTML(reaction.id)}</span><h3>${escapeHTML(ruleNames[reaction.rule_id] || reaction.rule_id)}</h3><span class="detail-tag">Depth ${reaction.depth ?? 'Unknown'}</span><span class="detail-tag">Elementary input</span>${channels.length ? `<select id="channel" class="channel-select" aria-label="Mapping channel">${channels.map((c,i) => `<option value="${i}" ${i === channelIndex ? 'selected' : ''}>Channel ${i+1} · ${c.match_count} template matches</option>`).join('')}</select><label class="map-toggle"><input id="reaction-map" type="checkbox" ${mapped ? 'checked' : ''}>Show reaction atom-map numbers</label>` : '<p class="detail-note">This legacy record has no complete mapping channels.</p>'}<div class="reaction-side">${cards(sides[0],0)}</div><div class="reaction-arrow">↓</div><div class="reaction-side">${cards(sides[1],1)}</div><dl><dt>Rule ID / reaction family</dt><dd>${escapeHTML(reaction.rule_id)}<br>${escapeHTML(reaction.family)}</dd></dl>${mapped && channel ? `<details><summary>Mapped reaction SMILES</summary><code>${escapeHTML(channel.mapped_reaction_smiles)}</code></details>` : ''}${singleEventDetails(channel)}<p class="detail-note">Each channel stores one representative of symmetry-equivalent mappings. Template-match counts are not kinetic degeneracy factors and may be incomplete after an early stop.</p></div>`;
  if (channels.length) {
    $('channel').onchange = event => showReaction(reaction,Number(event.target.value),$('reaction-map').checked);
    $('reaction-map').onchange = event => showReaction(reaction,Number($('channel').value),event.target.checked);
  }
}

function singleEventDetails(channel) {
  const se = channel?.single_event;
  if (!se) return '<p class="detail-note">Path factors have not been calculated. Enable Single-event annotation or apply it to this network.</p>';
  const center = se.reaction_center;
  const parameter = network.single_event?.parameters?.[se.parameter_class];
  const ok = ['calculated','imported'].includes(se.status);
  const counts = se.counting;
  const atoms = center.atom_changes.map(c => `Atom ${c.atom}: ` + Object.keys(c.before).filter(k => c.before[k] !== c.after[k]).map(k => `${k === 'radicals' ? 'radical electrons' : k} ${c.before[k] === '' ? 'unspecified' : c.before[k]} → ${c.after[k] === '' ? 'unspecified' : c.after[k]}`).join(', ')).join('; ');
  const bonds = center.bond_changes.map(c => `${c.atoms.join('–')}: ${c.before?.order ?? 'none'} → ${c.after?.order ?? 'none'}`).join('; ');
  return `<div class="single-event-details"><h4>Single-event annotation</h4><dl>
    <dt>Status</dt><dd>${escapeHTML(se.status)}</dd>
    <dt>Method / algorithm</dt><dd>${escapeHTML(se.method || 'Legacy annotation')} / ${escapeHTML(se.algorithm || '—')}</dd>
    ${ok ? `<dt>Independent single-event factor (n_e)</dt><dd>${escapeHTML(se.factor)}</dd>` : `<dt>Reason</dt><dd>${escapeHTML(se.reason)}${se.reason_code ? '<br>'+escapeHTML(se.reason_code) : ''}</dd>`}
    ${ok && counts ? `<dt>Path counting</dt><dd>${escapeHTML(counts.unique_labeled_paths)} distinct labelled events ÷ ${escapeHTML(counts.identical_reactant_divisor)} identical-reactant divisor = ${escapeHTML(se.factor)}<br>Statistical correction applied once</dd>` : ''}
    <dt>Reaction center (atom maps)</dt><dd>${escapeHTML(center.atom_maps.join(', '))}</dd>
    <dt>Bond changes</dt><dd>${escapeHTML(bonds || 'None')}</dd><dt>Atom-state changes</dt><dd><code>${escapeHTML(atoms || 'None')}</code></dd>
    ${se.family ? `<dt>Template family (optional label)</dt><dd>${escapeHTML(se.family)}</dd>` : ''}
    ${counts || se.event ? `<dt>Calculation evidence</dt><dd><details><summary>Show event and counting evidence</summary><pre>${escapeHTML(JSON.stringify({event:se.event,counting:counts},null,2))}</pre></details></dd>` : ''}
    ${se.sigma_reactants ? `<dt>Legacy reactant / TS symmetry numbers</dt><dd>${escapeHTML(se.sigma_reactants.join(' × '))} / ${escapeHTML(se.sigma_ts)}</dd>` : ''}
    ${se.kinetic_channel ? `<dt>Unique event channel</dt><dd><code>${escapeHTML(se.kinetic_channel)}</code></dd>` : ''}
    ${parameter ? `<dt>Optional shared parameter class</dt><dd>${escapeHTML(se.parameter_class)}<br>A = ${escapeHTML(parameter.A ?? 'Unfitted')} ${escapeHTML(parameter.A_units || '')}<br>n = ${escapeHTML(parameter.n ?? 'Unfitted')}<br>Ea = ${escapeHTML(parameter.Ea ?? 'Unfitted')} ${escapeHTML(parameter.Ea_units)}<br>${escapeHTML(parameter.status)} · ${escapeHTML(parameter.source)}</dd>` : ''}
    ${se.convention ? `<dt>Counting convention</dt><dd>${escapeHTML(se.convention)}</dd>` : ''}
    ${se.source ? `<dt>Source</dt><dd>${escapeHTML(se.source)}</dd>` : ''}
  </dl></div>`;
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

document.querySelectorAll('[data-path-example]').forEach(button => {
  button.onclick = async () => {
    try { fillSettings(await api(`/api/single-event/path-example/${button.dataset.pathExample}`)); message('Path-counting example loaded. Generate to enumerate labelled paths; no lattice or orientation multiplicity is included.',true); }
    catch (error) { message(error.message); }
  };
});

$('single-event-example').onclick = async () => {
  try { fillSettings(await api('/api/single-event/example')); message('Methane + H atom example loaded. Generate to validate sigma_R / sigma_TS = 12 / 3 = 4.',true); }
  catch (error) { message(error.message); }
};
$('apply-single-event').onclick = async () => {
  if (busy || !network) return;
  try {
    const options = JSON.parse($('single-event-options').value);
    const payload = await api('/api/single-event',{method:'POST',body:JSON.stringify(options)});
    const state = await api('/api/status');
    loadedRevision = state.revision;
    loadNetwork(payload,true);
    message('Reaction-path factors updated.',true);
  } catch (error) { message(error.message); }
};
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
$('apply-custom-rules').onclick = async () => {
  try {
    const rules = await api('/api/rules/validate', {method:'POST',body:$('custom-rules').value});
    renderRules(rules.rules, 'Custom elementary templates');
    message('Custom templates loaded. Set seed species, enable path factors if wanted, and Generate.',true);
  } catch (error) { message(`Could not load templates: ${error.message}`); }
};
$('import-rules').onclick = () => $('rules-file').click();
$('rules-file').onchange = async event => {
  const file = event.target.files[0];
  if (!file) return;
  try {
    const rules = JSON.parse(await file.text());
    if (!Array.isArray(rules) || !rules.length || rules.some(r => !r.id || !r.smarts)) throw new Error('Select a JSON array of rules containing id and smarts fields.');
    const checked = await api('/api/rules/validate',{method:'POST',body:JSON.stringify(rules)});
    renderRules(checked.rules,file.name);
    message();
  } catch (error) { message(`Could not import rules: ${error.message}`); }
  event.target.value = '';
};
$('reset-rules').onclick = () => renderRules(builtinRules);
async function resetWorkspace(request) {
  await api('/api/new',{method:'POST'});
  network = analysis = pathwayReactionIds = selectedId = focusId = null;
  speciesById.clear(); reactionsById.clear(); loadedRevision = -1;
  cy.elements().remove();
  $('details').innerHTML = emptyInspector;
  $('pathway-panel').hidden = true;
  $('single-event-summary').textContent = '';
  $('network-count').textContent = 'Not generated';
  $('visible-count').textContent = 'Drag to pan · Scroll to zoom · Click for details';
  $('view-note').textContent = '';
  $('depth-filter').disabled = false;
  $('search').value = '';
  $('depth-filter').innerHTML = '<option value="all">All depths</option>';
  $('neighbors').disabled = true;
  $('empty-state').hidden = false;
  $('bulk-smiles').value = '';
  fillSettings(request);
  setBusy(false);
  message();
  await pollStatus();
}
$('load-example').onclick = async () => {
  try {
    await resetWorkspace(builtinExample);
    $('input-hint').textContent = 'Propane + metal sites loaded · Five templates · Click Generate when ready.';
    message('Gas and catalytic example loaded. Edit any inputs, then click Generate.',true);
  } catch (error) { message(error.message); }
};
$('new-project').onclick = () => $('new-dialog').showModal();
$('cancel-new').onclick = () => $('new-dialog').close();
$('confirm-new').onclick = async () => {
  try {
    await resetWorkspace({seeds:[],rules:[],config:{max_depth:4}});
    $('input-hint').textContent = 'Start with your own species and templates, or click Load example.';
    $('new-dialog').close();
  }
  catch (error) { $('new-dialog').close(); message(error.message); }
};

async function initialize() {
  try {
    const result = await api('/api/rules');
    builtinRules = result.rules;
    builtinExample = {seeds:result.seeds,rules:builtinRules,config:{max_carbon:6,max_depth:4,max_species:200,max_reactions:1000},single_event:{calculator:'mapped_event_v1'}};
    fillSettings({seeds:[],rules:[],config:{max_depth:4}});
    await pollStatus();
    window.setInterval(pollStatus,700);
  } catch (error) { message(error.message); }
}
initialize();
