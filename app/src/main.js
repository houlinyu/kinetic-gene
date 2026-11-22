import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'three/addons/renderers/CSS2DRenderer.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/addons/postprocessing/UnrealBloomPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import GUI from 'lil-gui';
import './style.css';
import { Lane, CDS_STUB } from './sim.js';
import { LaneView, COLORS } from './world.js';
import { PARAMS, SOURCES, values } from './params.js';
import * as pace from './pace.js';
import { expected } from './analytic.js';
import { Scorecard } from './scorecard.js';
import { FoldView } from './folded.js';
import { Struct3DView, ASSETS, structNoteHTML } from './struct3d.js';
import { chipMotifs, chipText, strengthScale, motifStyle, summaryLine, stripHTML, tipHTML, legendHTML } from './motifs.js';
import { whyHTML, fitCardHTML } from './why.js';
import * as targetTab from './target/panel.js'; // therapeutic levers: mounted on BOTH pages, filtered

const STRESS_TC = 0.15; // matches findings.js
const CALM_TC = 1;
const $ = (id) => document.getElementById(id);
const fmt = (x, d = 1) => (x == null || !Number.isFinite(x) ? '—' : Number(x).toFixed(d));
const fmtR = (x) => (!Number.isFinite(x) ? '—' : x >= 10 ? x.toFixed(1) : x >= 1 ? x.toFixed(2) : x.toPrecision(2));
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// ---------- honesty: ratios are rank order, not fold change ----------
// Measured on this data set: predicted output spans 4.3e6-fold across leaders and within-gene
// sibling ratios reach 3.7e6x. No cell does that, so every printed ratio is clamped to
// >100x / <0.01x and labelled as a rank comparison rather than a fold change.
const RCAP_HI = 100, RCAP_LO = 0.01;
const RANK_TIP = "the model's dynamic range is unphysically wide; treat these as rank order, not fold-change";
// "×8.3" inside the bounds, ">100×" / "<0.01×" outside them
const fmtRC = (x) => (!Number.isFinite(x) ? '—' : x > RCAP_HI ? '&gt;100×' : x < RCAP_LO ? '&lt;0.01×' : '×' + fmtR(x));
const rankWord = (r) => (!Number.isFinite(r) ? 'not comparable'
  : r > RCAP_HI ? 'B ranks far above A' : r > 1.05 ? 'B ranks above A'
  : r < RCAP_LO ? 'B ranks far below A' : r < 0.95 ? 'B ranks below A' : 'A and B rank together');

// ---------- renderer / scene ----------
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.outputColorSpace = THREE.SRGBColorSpace;
renderer.toneMapping = THREE.ACESFilmicToneMapping;
document.body.appendChild(renderer.domElement);
const labels = new CSS2DRenderer();
labels.setSize(innerWidth, innerHeight);
Object.assign(labels.domElement.style, { position: 'fixed', top: 0, left: 0, pointerEvents: 'none' });
document.body.appendChild(labels.domElement);

const scene = new THREE.Scene();
const BG_CALM = new THREE.Color(COLORS.bg), BG_STRESS = new THREE.Color(0x1a0b0e);
scene.background = BG_CALM.clone();
scene.fog = new THREE.Fog(COLORS.bg, 30, 70);
const camera = new THREE.PerspectiveCamera(40, innerWidth / innerHeight, 0.1, 200);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true; controls.enablePan = false;
controls.minDistance = 4; controls.maxDistance = 150;
scene.add(new THREE.HemisphereLight(0xdde6ff, 0x20242e, 1.0));
const dl = new THREE.DirectionalLight(0xffffff, 1.0); dl.position.set(4, 8, 6); scene.add(dl);

// signature effect: bloom on the only HDR things (ribosomes, CDS stub, protein bursts)
const composer = new EffectComposer(renderer);
composer.addPass(new RenderPass(scene, camera));
const bloom = new UnrealBloomPass(new THREE.Vector2(innerWidth / 2, innerHeight / 2), 0.45, 0.25, 0.95);
composer.addPass(bloom);
composer.addPass(new OutputPass());

// ---------- state ----------
const P = values();
let data, genes = [], gene, laneA, laneB, viewA, viewB, paused = false;
// ANIMATION PACE. `sim.speed` is now RELATIVE: 1 = the watchable pace solved by pace.js, which
// maps the model's real nt/s onto a wall-clock tempo that is neither frozen nor a blur at any
// setting. The physical parameters keep their units and still produce every printed number;
// only the number of model seconds per second of your time changes, and that factor is printed
// on the page (#paceNote).
const sim = { speed: 1 };
const ANIM = { min: 0.25, max: 4, def: 1 };
let paceT = { ts: 1, crossSim: 1, crossWall: 1 };
const LANE_Y = [1.0, -2.6];
let idxA = 0, idxB = 1;
const target = { A: 0, B: 0, ratio: 1, stress: 0 };
const shown = { A: 0, B: 0, ratio: 1, stress: 0 };
const stressed = () => P.tc < 0.5;
let viewMode = 'linear', stacked = false;
let simple = true, allMotifs = false;
try { simple = localStorage.getItem('utr-mode') !== 'detail'; } catch {}
const FOLD_S = 6;
const foldBox = () => (stacked ? [9, 6] : [6, 8]);
try { const v = localStorage.getItem('utr-view'); if (['folded', 'struct'].includes(v)) viewMode = v; } catch {}

// ---------- layout: fit the strands into the free area (right of the panel, left of the findings) ----------
function layout() {
  const W = innerWidth, H = innerHeight, mobile = W <= 760;
  renderer.setSize(W, H); labels.setSize(W, H); composer.setSize(W, H);
  let x0 = 0, x1 = W, y0 = 0, y1 = H;
  const pr = $('panel').getBoundingClientRect(), sr = $('side').getBoundingClientRect();
  // the page switcher + the view toggle sit above the stage (see #tabs in style.css)
  if (mobile) { y0 = 118; y1 = pr.top - 8; }
  else { x0 = pr.right + 8; x1 = (getComputedStyle($('side')).display === 'none' ? W : sr.left) - 8; y0 = 92; y1 = H - 44; }
  const rw = Math.max(100, x1 - x0), rh = Math.max(100, y1 - y0);
  const wantStack = rw / rh < 1.25;
  if (wantStack !== stacked) { stacked = wantStack; if (viewMode !== 'linear' && laneA) setup(); }
  const rcx = (x0 + x1) / 2, rcy = (y0 + y1) / 2;
  const dx = W / 2 - rcx, dy = H / 2 - rcy;
  const fullW = W + 2 * Math.abs(dx), fullH = H + 2 * Math.abs(dy);
  camera.aspect = fullW / fullH;
  camera.setViewOffset(fullW, fullH, Math.abs(dx) + dx, Math.abs(dy) + dy, W, H);
  // scene bounds: strands x ∈ [-7, 7] plus "cap ▶" label; lanes + titles in y
  let sceneW = mobile ? 17.2 : 15.6, sceneH = (LANE_Y[0] + 1.9) - (LANE_Y[1] - 1.5);
  let cx = 0, cy = (LANE_Y[0] + 1.9 + LANE_Y[1] - 1.5) / 2;
  if (viewMode !== 'linear') {
    const [bw, bh] = foldBox();
    if (stacked) { sceneW = bw + 1.4; sceneH = 2 * bh + 2.6; cy = 0.3; }
    else { sceneW = 2 * bw + 2.2; sceneH = bh + 1.6; cy = 0.3; }
  }
  const t = Math.tan(THREE.MathUtils.degToRad(camera.fov / 2));
  const z = Math.max(sceneW * fullH / (2 * t * rw), sceneH * fullH / (2 * t * rh));
  controls.target.set(cx, cy, 0); camera.position.set(cx, cy, z);
  scene.fog.near = z * 1.2; scene.fog.far = z * 3; camera.far = z * 4;
  camera.updateProjectionMatrix(); controls.update();
  const hint = $('hint'); hint.style.left = rcx + 'px'; hint.style.maxWidth = rw - 16 + 'px';
  $('stressPill').style.left = mobile ? '' : rcx + 'px';
}
addEventListener('resize', layout);

// ---------- UI ----------
const isMane = (iso) => iso.tags?.includes('MANE_Select');
function isoOptions(current) {
  return gene.isoforms.map((iso, i) =>
    `<option value="${i}" ${i === current ? 'selected' : ''}>${iso.name} · ${iso.len} nt${isMane(iso) ? ' · MANE' : ''}</option>`).join('');
}
// Optimus 5′ only reads the last 100 nt of a leader, so two different isoforms can get byte-identical
// inputs — and then the same number must not be printed twice as if it were two predictions.
const optimusSame = (a, b) => a && b && a !== b
  && a.pred?.optimus_input === 'last100' && b.pred?.optimus_input === 'last100'
  && typeof a.utr5 === 'string' && typeof b.utr5 === 'string'
  && a.utr5.slice(-100) === b.utr5.slice(-100);
function extraLines(iso, other) {
  const out = [];
  const m = iso.measured;
  if (m && m.values && Object.keys(m.values).length) {
    const cell = 'HEK293T' in m.values ? 'HEK293T' : Object.keys(m.values)[0];
    const v = m.values[cell];
    // an ambiguous TE is gene-level: never show it as an isoform number
    if (v != null && !m.ambiguous) out.push(`Measured TE (Cao 2021, ${esc(cell)}): <b>${fmt(v, 2)}</b>`);
  }
  const mrl = iso.pred?.optimus_mrl;
  if (mrl != null) out.push(optimusSame(iso, other)
    ? `Optimus 5′ MRL: <b>${fmt(mrl, 2)}</b> <span class="cav">— identical Optimus input: the model reads only the last 100 nt, which are the same in both lanes, so this number cannot separate A from B</span>`
    : `Optimus 5′ MRL: <b>${fmt(mrl, 2)}</b>`);
  return out.length ? `<div class="meas">${out.join('<br>')}</div>` : '';
}
// ---------- experimental evidence (data/evidence.json audit) ----------
let evidence = null;
const evOf = (g) => evidence?.genes?.[g?.gene] ?? null;
// The curated evidence entry decides first; an ambiguous measured TE is an extra caveat on top of
// it, not a reason to throw the tier-A credit away (DDIT3-201 / ATF5-201 are both tested leaders).
function honestyLine(iso) {
  const m = iso.measured;
  const te = m && m.ambiguous ? '<div class="honest">measured TE is gene-level (isoforms not separable)</div>' : '';
  const ev = evOf(gene);
  if (!ev || ev.isoform_resolved === false || !isMane(iso))
    return '<div class="honest">model prediction — no experimental data for this isoform</div>' + te;
  return `<div class="honest ok">leader tested experimentally (${citeHTML(ev.cite, ev.url)})</div>` + te;
}
// cite keys the params.js SOURCES list does not carry (huez1998, lee2009, mihailovich2007,
// nanbru1997, ray2006) must still link out and still appear in the Sources drawer.
function citeHTML(key, url) {
  if (!key) return 'see evidence';
  const src = SOURCES[key];
  const href = url ?? src?.url;
  const ttl = src?.text ?? 'citation from data/evidence.json';
  return href ? `<a href="${esc(href)}" target="_blank" rel="noopener" title="${esc(ttl)}">${esc(key)}</a>` : `<span title="${esc(ttl)}">${esc(key)}</span>`;
}
function evDetailHTML(g, ev) {
  if (!ev) return '';
  const link = ev.url ? `<a href="${esc(ev.url)}" target="_blank" rel="noopener">${esc(ev.cite ?? 'paper')}</a>` : esc(ev.cite ?? '');
  return `<div class="tt-h">${esc(g)} · evidence tier ${esc(ev.tier)}</div>
    ${ev.quote ? `<div class="tt-note">“${esc(ev.quote)}”</div><div class="tt-pos">${esc(ev.quote_source ?? '')}</div>` : ''}
    ${ev.assays?.length ? `<div class="tt-note"><b>Assays:</b> ${ev.assays.map(esc).join('; ')}</div>` : ''}
    ${ev.cells ? `<div class="tt-note"><b>Cells:</b> ${esc(ev.cells)}</div>` : ''}
    ${ev.caveat ? `<div class="tt-note cav"><b>Caveat:</b> ${esc(ev.caveat)}</div>` : ''}
    <div class="tt-src">Source: ${link}</div>`;
}
function renderEvidence() {
  const ev = evOf(gene), b = $('evBadge');
  if (!ev) {   // absence must look like absence, not like a broken badge
    b.innerHTML = `<span class="evb tnone" title="data/evidence.json has no entry for this gene, so nothing here has been audited against a paper">○ no curated evidence entry</span>`;
    $('geneEv').innerHTML = `<div class="ftitle">Evidence for ${esc(gene?.gene ?? '—')}</div>
      <div class="tt-note">No entry in <code>data/evidence.json</code>: this gene was added after the evidence audit, so there is no tier, no quote and no assay list for it. Everything shown for it is model prediction only.</div>`;
    return;
  }
  b.innerHTML = `<span class="evb t${esc(ev.tier)}">${ev.tier === 'A' ? '● experimentally characterised leader' : '◐ indirect evidence (IRES/other)'}${ev.show_as_example ? '' : ' · no isoform-level evidence'} <u>why?</u></span>`;
  const el = b.firstElementChild;
  el.onpointerenter = (e) => { if (!tip.classList.contains('pinned')) showRich(e, evDetailHTML(gene.gene, ev) + '<div class="tt-hint">click to pin</div>'); };
  el.onpointerleave = () => { if (!tip.classList.contains('pinned')) tip.classList.remove('on', 'rich'); };
  el.onclick = (e) => { e.stopPropagation(); showRich(e, evDetailHTML(gene.gene, ev), true); };
  $('geneEv').innerHTML = `<div class="ftitle">Evidence for ${esc(gene.gene)}</div>${evDetailHTML(gene.gene, ev)}`;
}
function fillGenes() {
  const all = $('allGenes').checked;
  const shown = (g) => {
    const ev = evOf(g);
    if (evidence?.genes) return all || ev?.show_as_example === true;
    // no evidence file: old rule (the model can't start at a non-AUG MANE start codon)
    return all || !(g.isoforms.find(isMane) ?? g.isoforms[0])?.tags?.some(t => t.startsWith('non_AUG_start'));
  };
  const cur = gene;
  genes = data.genes.filter(shown);
  $('gene').innerHTML = genes.map((g, i) => {
    const ev = evOf(g), extra = ev && !ev.show_as_example ? ` · tier ${ev.tier}, no isoform-level evidence` : '';
    return `<option value="${i}" title="${esc(ev?.caveat ?? '')}">${g.gene} · ${g.isoforms.length} isoforms${extra}</option>`;
  }).join('');
  const keep = genes.indexOf(cur);
  if (keep >= 0) $('gene').value = keep;
  return keep >= 0;
}
function cardHTML(tag, iso, idx, other) {
  return `<div class="lhead"><span class="tag">${tag}</span><select data-lane="${tag}">${isoOptions(idx)}</select></div>
    <div class="rate"><span id="rate${tag}">—</span> <small>proteins / min per mRNA</small></div>${honestyLine(iso)}
    <div class="cs" id="cs${tag}"></div>${summaryLine(iso)}${extraLines(iso, other)}${stripHTML(iso)}`;
}
function statsHTML(tag, iso) {
  const nU = iso.uorfs?.length || 0, nOver = (iso.uorfs || []).filter(u => u.type === 'overlapping').length;
  return `<div><span class="tag">${tag}</span> ${esc(iso.name)}: ${iso.len} nt · ${nU} uORF${nU === 1 ? '' : 's'}${nOver ? ` (${nOver} overlapping)` : ''} · Kozak ${iso.main_kozak?.class ?? '—'} ·
    MFE ${fmt(iso.fold?.mfe)} kcal/mol · cap-50nt ΔG ${fmt(iso.fold?.cap_dG)} · GC ${fmt((iso.gc ?? NaN) * 100, 0)}%${Number.isFinite(iso.fold?.acc_start) ? ` · start accessibility ${fmt(iso.fold.acc_start, 2)}` : ''}${iso.fold?.method ? ` · fold: ${esc(iso.fold.method)}` : ''}</div>`;
}

function pickDefaultPair() {
  const n = gene.isoforms.length;
  idxA = Math.max(0, gene.isoforms.findIndex(isMane));
  const A = gene.isoforms[idxA];
  let cand = [...gene.isoforms.keys()].filter(i => i !== idxA && gene.isoforms[i].len <= 2.5 * A.len && gene.isoforms[i].len >= A.len / 2.5);
  if (!cand.length) cand = [...gene.isoforms.keys()].filter(i => i !== idxA);
  idxB = n > 1 ? cand.sort((a, b) => Math.abs(gene.isoforms[b].len - A.len) - Math.abs(gene.isoforms[a].len - A.len))[0] : 0;
}

function label(view, cls, text, x, y, cx, cy) {
  const el = document.createElement('div'); el.className = cls; el.textContent = text;
  const o = new CSS2DObject(el); o.position.set(x, y, 0); o.center.set(cx, cy);
  view.group.add(o); return o;
}
function dropView(view) {
  if (!view) return;
  view.group.traverse(o => { if (o.isCSS2DObject) o.element.remove(); });
  view.dispose();
}
function decorate(view, tag, iso) {
  const fold = viewMode !== 'linear';
  if (view.attachStallTags) view.attachStallTags((text, pos) => {
    const el = document.createElement('div'); el.className = 'stall-tag'; el.textContent = text;
    const o = new CSS2DObject(el); o.position.copy(pos); o.position.z += view.R * 3; o.center.set(0.5, 1.4);
    view.content.add(o); return o;
  });
  if (view.scaleLabel) {
    const el = document.createElement('div'); el.className = 'scale-tag'; el.textContent = view.scaleLabel.text;
    const o = new CSS2DObject(el); o.position.copy(view.scaleLabel.pos); o.center.set(0.5, -0.6); view.content.add(o);
  }
  const [tx, ty] = view.titlePos, [cx, cy, cz] = view.capPos, [dx, dy, dz] = view.cdsPos;
  label(view, 'lane-title', `${tag} · ${iso.name} (${isMane(iso) ? 'MANE, standard' : 'alternative'})${fold ? ` · ${iso.len} nt` : ''}`, tx, ty, fold ? 0.5 : 0, 1);
  label(view, 'lane-cap', fold ? '5′ cap' : 'cap ▶', cx, cy, fold ? 0.5 : 1, fold ? -0.4 : 0.5).position.z = cz ?? 0;
  label(view, 'lane-cds', fold ? 'AUG → protein' : 'protein-coding ▶', dx, dy, 0.5, fold ? -0.3 : 1).position.z = dz ?? 0;
  // motif chips
  // simple mode: only decisive motifs (strong uAUG, overlapping-uORF starts, main Kozak)
  const overStarts = new Set((iso.uorfs || []).filter(u => u.type === 'overlapping' || u.end > iso.len).map(u => u.start));
  const decisive = (m) => m.type === 'kozak_main' || (m.type === 'uAUG' && (m.strength === 'strong' || overStarts.has(m.start)));
  const chips = chipMotifs(iso, innerWidth <= 760 ? 6 : fold ? 10 : 12, simple && !allMotifs ? decisive : null);
  for (const m of chips) {
    const a = view.motifAnchor(m);
    const el = document.createElement('div');
    el.className = 'mchip' + (m.strength === 'weak' ? ' weak' : '');
    const col = motifStyle(m.type).color;
    el.style.setProperty('--mc', col);
    el.style.fontSize = (10.5 * strengthScale(m)).toFixed(1) + 'px';
    el.textContent = chipText(m);
    el.onpointerenter = (e) => { if (!tip.classList.contains('pinned')) showRich(e, tipHTML(m, iso, data)); };
    el.onpointerleave = () => { if (!tip.classList.contains('pinned')) tip.classList.remove('on', 'rich'); };
    el.onclick = (e) => { e.stopPropagation(); showRich(e, tipHTML(m, iso, data), true); };
    const o = new CSS2DObject(el);
    const p = Array.isArray(a.pos) ? a.pos : a.pos.toArray();
    o.position.set(p[0], p[1], p[2] ?? 0); o.center.set(a.center[0], a.center[1]);
    view.group.add(o);
  }
}
function positionFold() {
  if (viewMode === 'linear' || !(viewA instanceof FoldView)) return;
  const [bw, bh] = foldBox(), off = (stacked ? bh : bw) / 2 + 0.6;
  viewA.group.position.set(stacked ? 0 : -off, stacked ? off + 0.2 : 0, 0);
  viewB.group.position.set(stacked ? 0 : off, stacked ? -off + 0.2 : 0, 0);
  viewA.group.rotation.set(-0.3, stacked ? 0 : 0.2, 0);
  viewB.group.rotation.set(-0.3, stacked ? 0 : -0.2, 0);
}
function setView(mode) {
  viewMode = mode;
  try { localStorage.setItem('utr-view', mode); } catch {}
  document.querySelectorAll('#viewToggle button').forEach(b => b.classList.toggle('on', b.dataset.v === mode));
  setHint();
  if (0) $('hint').textContent = mode === 'folded'
    ? 'folded 2D structure · click amber uORF = mutate AUG · click a stem = melt it · hover chips = motif · drag = orbit · 1/2/3 = view · S = stress'
    : mode === 'struct'
      ? 'coarse 3D model · click amber backbone = mutate uORF AUG · click a helix = melt it · hover chips = motif · drag = orbit/zoom · 1/2/3 = view · S = stress'
      : 'click orange uORF = mutate AUG · click blue arcs = melt hairpin · hover chips = motif · drag = orbit · 1/2/3 = view · S = stress';
  $('viewNote').innerHTML = mode === 'struct' ? structNoteHTML() : '';
  if (0) $('viewNote').innerHTML = mode === 'struct'
    ? `<b>3D model built from the predicted secondary structure with ideal A-form helices — not a predicted tertiary structure.</b><br>
       White blob = scanning 43S; helices open ahead of it and re-close behind; <span class="hot">red</span> = blocked by a stable hairpin (green = eIF4A helicase working).
       <span class="pk"><i style="background:#4fc3f7"></i>G·C <i style="background:#9ad18b"></i>A·U <i style="background:#c39bff"></i>G·U pair</span>`
    : '';
  $('viewNote').style.display = mode === 'struct' ? '' : 'none';
  if (gene) { setup(); }
  layout();
}
function setHint() {
  if (simple) { $('hint').textContent = 'click a uORF, or the rate-limiting hairpin (hover to find it) · S = stress · D = detail'; return; }
  $('hint').textContent = viewMode === 'folded'
    ? 'folded 2D structure · click amber uORF = mutate AUG · click a stem = melt it (only the rate-limiting one changes the number — hover to see which) · drag = orbit · 1/2/3 = view · S = stress · D = simple'
    : viewMode === 'struct'
      ? 'coarse 3D model · click amber backbone = mutate uORF AUG · click a helix = melt it · hover chips = motif · drag = orbit/zoom · 1/2/3 = view · D = simple'
      : 'click orange uORF = mutate AUG · click blue arcs = melt hairpin (only the rate-limiting one moves the number — hover to see which) · drag = orbit · 1/2/3 = view · S = stress · D = simple';
}
function setMode(isSimple) {
  simple = isSimple;
  try { localStorage.setItem('utr-mode', simple ? 'simple' : 'detail'); } catch {}
  document.body.classList.toggle('simple', simple);
  document.body.classList.remove('side-open');
  document.querySelectorAll('#modeToggle button').forEach(b => b.classList.toggle('on', (b.dataset.m === 'simple') === simple));
  setHint();
  if (gene) setup();
  layout();
  if (data) onParams('tc');
  spliceTab?.setSimple(simple);
}
document.querySelectorAll('#modeToggle button').forEach(b => b.onclick = () => setMode(b.dataset.m === 'simple'));
$('allMotifs').onchange = (e) => { allMotifs = e.target.checked; if (gene) setup(); };
function showRich(e, html, pin) {
  tip.innerHTML = html; tip.classList.add('on', 'rich'); tip.classList.remove('strong');
  tip.classList.toggle('pinned', !!pin);
  clearTimeout(tipTimer);
  const r = tip.getBoundingClientRect();
  let x = e.clientX + 14, y = e.clientY + 12;
  if (x + r.width > innerWidth - 8) x = Math.max(8, e.clientX - r.width - 14);
  if (y + r.height > innerHeight - 8) y = Math.max(8, e.clientY - r.height - 12);
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
}
document.addEventListener('pointerdown', (e) => {
  if (tip.classList.contains('pinned') && !tip.contains(e.target) && !e.target.classList?.contains('mchip')) tip.classList.remove('on', 'rich', 'pinned');
});
document.querySelectorAll('#viewToggle button').forEach(b => b.onclick = () => setView(b.dataset.v));

function prewarm(lane, secs = 45) {
  for (let t = 0; t < secs; t += 0.05) lane.step(0.05, P);
}

function setup() {
  dropView(viewA); dropView(viewB);
  const isoA = gene.isoforms[idxA], isoB = gene.isoforms[idxB] ?? isoA;
  const maxLen = Math.max(isoA.len, isoB.len);
  laneA = new Lane(isoA); laneB = new Lane(isoB);
  prewarm(laneA); prewarm(laneB);
  if (viewMode !== 'linear') {
    const V = viewMode === 'struct' ? Struct3DView : FoldView;
    const sz = (L) => foldBox().map(v => v * THREE.MathUtils.clamp(Math.sqrt(L / maxLen), 0.6, 1));
    viewA = new V(scene, sz(isoA.len)); viewB = new V(scene, sz(isoB.len));
    positionFold();
  } else {
    viewA = new LaneView(scene, LANE_Y[0], maxLen); viewB = new LaneView(scene, LANE_Y[1], maxLen);
  }
  viewA.build(laneA); decorate(viewA, 'A', isoA);
  viewB.build(laneB);
  if (viewA.u && viewB.u) { // true scale: both lanes share one Å -> world factor
    const um = Math.min(viewA.u, viewB.u);
    viewA.group.scale.setScalar(um / viewA.u); viewB.group.scale.setScalar(um / viewB.u);
  }
  decorate(viewB, 'B', isoB);
  if (gene.isoforms.length < 2) {            // nothing to race: do not draw a second copy of A
    viewB.group.visible = false;
    if (Array.isArray(viewB.pickables)) viewB.pickables.length = 0;
  }
  laneA.events.length = 0; laneB.events.length = 0;
  $('cardA').innerHTML = cardHTML('A', isoA, idxA, isoB);
  // One isoform in the data, or the same one picked twice, is not a comparison: say so instead of
  // racing a leader against itself and printing ×1.00.
  const same = gene.isoforms.length < 2 || idxB === idxA;
  $('cardB').innerHTML = gene.isoforms.length < 2
    ? `<div class="lhead"><span class="tag">B</span> <span class="mute">— none</span></div>
       <div class="sameiso"><b>${esc(gene.gene)} has only one 5′ leader in this data set</b>, so there is no second isoform to race it against. Pick another gene to compare.</div>`
    : cardHTML('B', isoB, idxB, isoA);
  document.body.classList.toggle('one-iso', same);
  $('stats').innerHTML = statsHTML('A', isoA) + statsHTML('B', isoB);
  renderEvidence();
  renderXLink();              // does this gene exist on the ① Transcription page?
  $('why').innerHTML = whyHTML(isoA, isoB, isoA.name, isoB.name);
  document.querySelectorAll('.lane-card select').forEach(s => s.onchange = () => {
    if (s.dataset.lane === 'A') idxA = +s.value; else idxB = +s.value; setup();
  });
  recalc(true);
}

// ---------- exact numbers (analytic model) ----------
const optsOf = (lane) => ({
  disabledU: new Set(lane.uorfs.filter(u => !u.enabled).map(u => u.idx)),
  melted: new Set(lane.hairpins.filter(h => h.melted).map(h => h.idx)),
});
const laneOut = (lane) => expected(lane.iso, P, optsOf(lane)).perMin;
// Which hairpin currently gates this lane: melt each one alone and see which melt moves the output.
// (In this model only the single longest stall limits throughput, so usually exactly one does.)
function limitingHP(lane) {
  const base = laneOut(lane);
  let best = -1, bestGain = 0;
  lane.hairpins.forEach((h, i) => {
    if (h.melted) return;
    const o = optsOf(lane);
    o.melted = new Set([...o.melted, h.idx]);
    const gain = expected(lane.iso, P, o).perMin - base;
    if (gain > bestGain * 1.000001 && gain > base * 1e-3) { bestGain = gain; best = i; }
  });
  return best;
}
function hairpinTip(lane, idx) {
  const h = lane.hairpins[idx];
  if (h.melted) return 'melted hairpin · click to restore';
  const dg = `hairpin ΔG ${fmt(h.dG)} kcal/mol`;
  if (!(h.dG < 0)) return `${dg} — not a favourable hairpin (ΔG ≥ 0); the model still stalls on it · click to melt`;
  return limitingHP(lane) === idx
    ? `${dg} · RATE-LIMITING right now — melting this one raises the output · click to melt`
    : `${dg} · not rate-limiting: only the longest stall gates throughput, so melting it changes the number by <0.1% · click to melt`;
}
let fatesCache = null;
function recalc(snap = false) {
  if (!laneA) return;
  const res = {};
  for (const [tag, lane] of [['A', laneA], ['B', laneB]]) {
    const o = optsOf(lane);
    const now = expected(lane.iso, P, o);
    const c = expected(lane.iso, { ...P, tc: CALM_TC }, o).perMin;
    const s = expected(lane.iso, { ...P, tc: STRESS_TC }, o).perMin;
    res[tag] = now; target[tag] = now.perMin;
    const r = s / Math.max(c, 1e-12), arrow = r > 1.05 ? '<b class="up">▲</b>' : r < 0.95 ? '<b class="down">▼</b>' : '<b>=</b>';
    const cs = $('cs' + tag); if (!cs) continue;
    cs.innerHTML = `<span class="${stressed() ? '' : 'on'}">calm ${fmtR(c)}</span> → <span class="${stressed() ? 'on' : ''}">stress ${fmtR(s)}</span> ${arrow} <span class="mute" title="${esc(RANK_TIP)}">${fmtRC(r)}</span>`;
  }
  target.ratio = res.B.perMin / Math.max(res.A.perMin, 1e-12);
  fatesCache = res;
  renderFates();
  if (snap) { shown.A = target.A; shown.B = target.B; shown.ratio = target.ratio; }
  renderNumbers();
}

const FCOL = { main: '#e85d9c', uorfDrop: '#ffb547', overlapLoss: '#ff6b4a', leaky: '#7c8494' };
const FNAME = { main: 'made protein', uorfDrop: 'lost after uORF', overlapLoss: 'stuck in overlapping uORF', leaky: 'scanned past start' };
function renderFates() {
  const bar = (r, tag) => {
    const f = { main: r.pMain, uorfDrop: r.lost.uorfDrop, overlapLoss: r.lost.overlapLoss };
    f.leaky = Math.max(0, 1 - f.main - f.uorfDrop - f.overlapLoss);
    return `<div class="row"><span class="tag sm">${tag}</span><div class="bar">${Object.keys(FCOL).map(k => `<span title="${FNAME[k]}: ${(100 * f[k]).toFixed(0)}%" style="width:${100 * f[k]}%;background:${FCOL[k]}"></span>`).join('')}</div><span class="pct">${(100 * f.main).toFixed(0)}%</span></div>`;
  };
  const one = gene?.isoforms?.length < 2;
  // This bar is PER RIBOSOME THAT LOADS (an efficiency); the verdict above is per mRNA per minute
  // (a flux). They can point opposite ways, so each says which quantity it is.
  $('fates').innerHTML = `<div class="ftitle">What happens to each ribosome that loads <span class="mute">— % per loaded ribosome, not per mRNA per minute</span></div>${bar(fatesCache.A, 'A')}${one ? '' : bar(fatesCache.B, 'B')}
    <div class="fkey">${Object.keys(FCOL).map(k => `<span><i style="background:${FCOL[k]}"></i>${FNAME[k]}</span>`).join('')}</div>
    <div class="fnote">A lane can start a higher <i>fraction</i> of its ribosomes and still make less protein per minute: the verdict divides by how fast ribosomes load and how long they stall.</div>`;
}

function renderNumbers() {
  const ra = $('rateA'), rb = $('rateB');       // lane B's card is absent for single-isoform genes
  if (ra) ra.textContent = fmtR(shown.A);
  if (rb) rb.textContent = fmtR(shown.B);
  const r = shown.ratio;
  const rs = fmtRC(r);
  if (gene && (gene.isoforms.length < 2 || idxA === idxB)) {
    $('verdict').innerHTML = `<div class="vsame">${gene.isoforms.length < 2
      ? `<b>Only one 5′ leader for ${esc(gene.gene)} in this data set</b> — there is no second isoform to compare, so no ratio is shown.`
      : '<b>The same isoform is selected in both lanes</b> — pick a different one in lane B to get a comparison.'}</div>`;
    return;
  }
  // RANK comparison, not a fold change: the printed number is capped (see RANK_TIP)
  $('verdict').innerHTML = `<div class="vbig" title="${esc(RANK_TIP)}"><b${rs.length > 5 ? ' class="long"' : ''}>${rs}</b><span>${rankWord(r)}<br><span class="mute">rank comparison of protein per mRNA per minute</span></span></div>
    <div class="vsub">exact model · ${stressed() ? '<span class="hot">stressed cell (eIF2·TC low)</span>' : 'calm cell'}</div>
    <div class="vcap" title="${esc(RANK_TIP)}">⚠ The model's dynamic range is unphysically wide — predicted output spans ~4×10⁶-fold across these leaders, and sibling isoforms of one gene differ by up to ~4×10⁶. Printed ratios are capped at ×100 / ×0.01 and mean <b>rank order, not fold-change</b>.</div>`;
}

// ---------- stress ----------
function setStress(on) {
  P.tc = on ? STRESS_TC : CALM_TC;
  gui.controllersRecursive().forEach(c => c.updateDisplay());
  onParams('tc');
}
function syncStressUI() {
  const on = stressed();
  document.body.classList.toggle('stress', on);
  target.stress = on ? 1 : 0;
  $('stressBtn').textContent = on ? '☀ Remove stress' : '⚡ Apply stress';
}
$('stressBtn').onclick = () => setStress(!stressed());

// ---------- scorecard / feasibility ----------
const card = new Scorecard($('score'), {
  onPick(lam, beta) {
    P.tcRate = +(lam * P.scanSpeed).toPrecision(3); P.loadTCexp = +beta.toFixed(2);
    gui.controllersRecursive().forEach(c => c.updateDisplay());
    onParams('tcRate');
  },
});
let mapTimer = null;
function onParams(k) {
  syncStressUI();
  recalc();
  if (!data) return;
  const s = card.update(data, P);
  const open = document.body.classList.contains('side-open');
  // same denominator as the scorecard itself: every direction check that ran, failures included
  const tested = s.res.length, ok = s.res.filter(f => f.result.pass).length;
  $('sideToggle').innerHTML = simple && innerWidth > 760
    ? `<b>${ok} of ${tested}</b> direction checks at <i>our chosen</i> parameter point
       <span class="mute">— all five together in 0.70% of re-sampled parameter sets</span> <span class="det">${open ? 'hide ▴' : 'details ▸'}</span>`
    : `${ok}/${tested} checks ${open ? '▾' : '▸'}`;
  $('sideToggle').classList.toggle('ok', ok === tested);
  if (!['tc', 'tcRate', 'loadTCexp'].includes(k)) { clearTimeout(mapTimer); mapTimer = setTimeout(() => card.recompute(data, P), 250); }
}
$('sideToggle').onclick = () => { document.body.classList.toggle('side-open'); onParams('tc'); };

// ---------- parameter GUI (Advanced, closed) ----------
const gui = new GUI({ title: 'Advanced: all parameters', container: $('guiHost') });
gui.close();
// Relative to the watchable pace, not to real time: 1 = the pace pace.js solved for this leader
// and these rates. The absolute compression it works out to is printed live in #paceNote, and
// pace.js clamps the result so neither end of this slider can freeze or blur the picture.
gui.add(sim, 'speed', ANIM.min, ANIM.max, 0.25).name('animation speed × watchable pace');
for (const k of Object.keys(PARAMS)) {
  const d = PARAMS[k];
  gui.add(P, k, d.min, d.max, d.step).name(`${d.label}${d.cite ? '' : ' (assumed)'}`).onChange(() => onParams(k));
}
function renderSources() {
  // evidence.json cites papers that params.js SOURCES does not list; add them from the audit itself
  // rather than letting a bare key like "lee2009" appear with no reference anywhere.
  const extra = [];
  for (const [g, e] of Object.entries(evidence?.genes ?? {})) {
    if (!e.cite || SOURCES[e.cite] || extra.some(x => x.k === e.cite)) continue;
    extra.push({ k: e.cite, url: e.url, text: `cited in data/evidence.json for ${g}${e.quote ? ` — “${e.quote}”` : ''}` });
  }
  $('sources').innerHTML = `<div class="ftitle">Sources</div><ul>${Object.entries(SOURCES).map(([k, s]) =>
    `<li><b>${esc(k)}</b>: ${s.url ? `<a href="${esc(s.url)}" target="_blank" rel="noopener">${esc(s.text)}</a>` : esc(s.text)}</li>`).join('')}
    ${extra.map(s => `<li><b>${esc(s.k)}</b>: ${s.url ? `<a href="${esc(s.url)}" target="_blank" rel="noopener">${esc(s.text)}</a>` : esc(s.text)}</li>`).join('')}</ul>
    <div class="fine">Assumed (no citation yet): ${Object.entries(PARAMS).filter(([, d]) => !d.cite).map(([, d]) => esc(d.label)).join(' · ')}</div>`;
}
renderSources();

// ---------- picking: hover cursor + floating tag ----------
const ray = new THREE.Raycaster(), mouse = new THREE.Vector2();
let downAt = null, tipTimer = null;
const tip = $('tip');
function pick(e) {
  if (!viewA) return null;
  mouse.set((e.clientX / innerWidth) * 2 - 1, -(e.clientY / innerHeight) * 2 + 1);
  ray.setFromCamera(mouse, camera);
  const hits = ray.intersectObjects([...viewA.pickables, ...viewB.pickables], false);
  // uORF sleeves sit inside hairpin pick volumes in the folded view: prefer them
  return (hits.find(h => h.object.userData.kind === 'uorf') ?? hits[0])?.object.userData ?? null;
}
function showTip(e, text, hold) {
  tip.classList.remove('rich', 'pinned');
  tip.textContent = text; tip.style.left = e.clientX + 14 + 'px'; tip.style.top = e.clientY - 10 + 'px';
  tip.classList.add('on'); tip.classList.toggle('strong', !!hold);
  clearTimeout(tipTimer);
  if (hold) tipTimer = setTimeout(() => tip.classList.remove('on', 'strong'), 1800);
}
renderer.domElement.addEventListener('pointermove', e => {
  if (tip.classList.contains('strong') || tip.classList.contains('rich')) return;
  const u = pick(e);
  renderer.domElement.style.cursor = u ? 'pointer' : '';
  if (!u) { tip.classList.remove('on'); return; }
  const lane = u.lane.lane;
  showTip(e, u.kind === 'uorf'
    ? (lane.uorfs[u.idx].enabled ? `uORF (${lane.uorfs[u.idx].type}) · click to mutate its AUG` : 'mutated uORF · click to restore')
    : hairpinTip(lane, u.idx));
});
renderer.domElement.addEventListener('pointerdown', e => downAt = [e.clientX, e.clientY]);
renderer.domElement.addEventListener('pointerup', e => {
  if (!downAt || Math.hypot(e.clientX - downAt[0], e.clientY - downAt[1]) > 5) return;
  const u = pick(e); if (!u) return;
  const { kind, idx, lane: view } = u, lane = view.lane;
  let msg;
  if (kind === 'uorf') { const on = !lane.uorfs[idx].enabled; lane.setUorf(idx, on); msg = on ? 'uORF restored' : 'AUG mutated: uORF removed'; }
  else {
    const m = !lane.hairpins[idx].melted;
    const before = laneOut(lane);
    lane.setHairpin(idx, m);
    const after = laneOut(lane);
    const moved = Math.abs(after - before) > Math.abs(before) * 1e-3;
    // only the single largest hairpin stall gates throughput in this model, so most melts are
    // no-ops. Say so instead of leaving the user guessing.
    msg = m ? (moved ? 'hairpin melted — output rises' : 'melted — not rate-limiting, output unchanged')
      : (moved ? 'hairpin restored — output falls' : 'restored — not rate-limiting, output unchanged');
  }
  prewarm(lane, 30); lane.events.length = 0;
  view.refreshToggles();
  showTip(e, msg, true);
  recalc();
});
addEventListener('keydown', e => {
  if (['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
  if (e.code === 'KeyT') { setPage(page === 'translation' ? 'transcription' : 'translation'); return; }
  // the levers panel is an overlay on whichever page opened it, and it owns the keyboard while
  // it is up (Escape closes it) — the stage underneath must not react to view/stress keys
  if (document.body.classList.contains('tgt-open')) return;
  if (page === 'transcription') { // page 1 owns its own shortcuts
    if (e.code === 'Digit1') spliceTab?.setView('diagram');
    if (e.code === 'Digit2') spliceTab?.setView('struct3d');
    if (e.code === 'KeyF') spliceTab?.setView(spliceTab.view === 'diagram' ? 'struct3d' : 'diagram');
    if (e.code === 'KeyD') setMode(!simple);
    return;
  }
  // never let page 2's shortcuts fire anywhere else
  if (page !== 'translation') return;
  if (e.code === 'Space') { paused = !paused; e.preventDefault(); }
  if (e.code === 'KeyS') setStress(!stressed());
  const VM = ['linear', 'folded', 'struct'];
  if (e.code === 'KeyD') setMode(!simple);
  if (e.code === 'KeyF') setView(VM[(VM.indexOf(viewMode) + 1) % 3]);
  if (['Digit1', 'Digit2', 'Digit3'].includes(e.code)) setView(VM[+e.code.slice(5) - 1]);
});

// ---------- animation pace (view only: it never touches a printed number) ----------
// The lanes are drawn to one shared scale (world.js: WIDTH / (maxLen + CDS_STUB)), so both
// ribosomes move at the same px/s. One traversal of that drawn width is the thing we keep inside
// the watchable window.
function paceSolve(mul = sim.speed, rate = P.scanSpeed) {
  const nt = laneA ? Math.max(laneA.iso.len, laneB ? laneB.iso.len : 0) + CDS_STUB : 1;
  // how long the lane stays empty between two 43S loadings (hairpin blocking excluded: it is
  // lane-specific, and we must not speed the whole clock up because one lane is gated)
  const load = Math.max(1e-6, P.loadRate * Math.pow(Math.max(P.tc, 1e-6), P.loadTCexp ?? 1));
  return pace.solve({
    nt, rate, mul,
    refNt: 340,                        // a ~250 nt leader + the 90 nt CDS stub: a typical case
    refRate: PARAMS.scanSpeed.value,   // the default scanning speed lands on PACE.MID
    gapSim: 1 / load,
    fastRate: P.elongSpeed,            // elongating ribosomes (uORF / CDS) must not blur either
  });
}
let paceEl = null;
function renderPaceNote() {
  if (!paceEl) {
    paceEl = document.createElement('div');
    paceEl.id = 'paceNote'; paceEl.className = 'pacenote';
    $('cards').before(paceEl);
  }
  const p = paceT;
  if (!laneA || !(p.ts > 0)) { paceEl.textContent = ''; return; }
  const lo = paceSolve(sim.speed, PARAMS.scanSpeed.min), hi = paceSolve(sim.speed, PARAMS.scanSpeed.max);
  const mins = Math.floor(p.crossSim / 60), secs = Math.round(p.crossSim % 60);
  paceEl.innerHTML = `<b>⏱ ${esc(pace.label(p.ts))}</b>
    <span>A ribosome crosses this leader in <b>${fmt(p.crossWall)} s</b> of your time —
    ${p.crossSim >= 60 ? `${mins} min ${secs} s` : `${fmt(p.crossSim)} s`} in the cell
    (${Math.round(p.crossSim * P.scanSpeed)} nt at ${fmt(P.scanSpeed, 0)} nt/s).</span>
    <span>${esc(pace.ratioNote('Scanning speed', lo, hi))}${p.limit ? ` · at this setting the picture is held at the <b>${p.limit}</b> watchable pace, so the on-screen difference here is smaller than the real one` : ''}</span>
    <span class="fine">Only the wall-clock tempo is compressed. Every number on this page is the exact model at the real rates.</span>`;
}

// ---------- loop ----------
let last = performance.now(), uiClock = 0;
const fogCol = new THREE.Color();
function tick(now) {
  const dt = Math.min((now - last) / 1000, 0.05); last = now;
  if (laneA) paceT = paceSolve();
  if (laneA && !paused) {
    // sub-stepping keeps the stochastic model accurate however much wall time is compressed;
    // capped so a 100x+ compression cannot cost a frame (0.02 s -> at most 0.06 s per sub-step)
    const simDt = dt * paceT.ts, sub = Math.min(240, Math.max(1, Math.ceil(simDt / 0.02)));
    for (let i = 0; i < sub; i++) { laneA.step(simDt / sub, P); laneB.step(simDt / sub, P); }
    viewA.update(dt); viewB.update(dt);
  }
  // tweens
  const k = 1 - Math.exp(-dt * 5);
  // snap on settle: easing alone leaves a visible 0.1-0.5% residue (26.5 -> 26.3) after un-stressing
  for (const key of ['A', 'B', 'ratio']) {
    shown[key] += (target[key] - shown[key]) * k;
    if (Math.abs(target[key] - shown[key]) <= Math.abs(target[key]) * 1e-4 + 1e-9) shown[key] = target[key];
  }
  shown.stress += (target.stress - shown.stress) * (1 - Math.exp(-dt * 3));
  scene.background.copy(BG_CALM).lerp(BG_STRESS, shown.stress);
  scene.fog.color.copy(fogCol.copy(BG_CALM).lerp(BG_STRESS, shown.stress));
  uiClock += dt;
  if (uiClock > 0.06 && laneA) { uiClock = 0; renderNumbers(); renderPaceNote(); }
  controls.update();
  composer.render(); labels.render(scene, camera);
}
renderer.setAnimationLoop(tick);

// ---------- pacing diagnostics (read-only) ----------
// Used by scratchpad/pace/measure.mjs to MEASURE on-screen px/s and traversal time instead of
// judging it by eye. Reports where the markers are in CSS pixels; changes nothing.
const _pv = new THREE.Vector3();
const toPx = (grp, x) => {
  _pv.set(x, 0, 0).applyMatrix4(grp.matrixWorld).project(camera);
  return (_pv.x + 1) / 2 * renderer.domElement.clientWidth;
};
window.__kgPace = () => {
  const o = { page, canvasW: renderer.domElement.clientWidth, timeScale: paceT?.ts, mul: sim.speed, anim: ANIM, paused };
  if (laneA && viewA?.X && viewMode === 'linear') {
    const L = laneA.iso.len;
    o.lane = { iso: laneA.iso.name, utrLen: L, spanNt: L + CDS_STUB, x0: toPx(viewA.group, viewA.X(0)), x1: toPx(viewA.group, viewA.X(L + CDS_STUB)) };
    o.marks = laneA.ribos.map(r => ({ id: r.id, px: toPx(viewA.group, viewA.X(r.pos)), nt: r.pos, mode: r.mode }));
    o.pace = paceT;
  }
  if (spliceTab?.paceInfo) o.splice = spliceTab.paceInfo();
  return o;
};

// ---------- data ----------
// A data failure must not blank the page: it renders a visible error and still mounts the tabs.
async function loadData() {
  const tryJSON = async (u) => {
    const r = await fetch(u);
    if (!r.ok) throw new Error(`${u}: HTTP ${r.status}`);
    return await r.json();            // throws if Vite served index.html instead of JSON
  };
  try { return await tryJSON('./utr_data.json'); }
  catch (e1) {
    console.warn('main: utr_data.json unusable', e1);
    try { return await tryJSON('./mock_data.json'); }
    catch (e2) { throw new Error(`${e1.message}; fallback ${e2.message}`); }
  }
}
function dataError(msg) {
  $('verdict').innerHTML = `<div class="vcap"><b>No transcript data.</b> <code>app/public/utr_data.json</code> could not be loaded or parsed (${esc(msg)}),
    so the ② Translation page has nothing to run. The ① Transcription page still works.</div>`;
  $('cards').innerHTML = ''; $('fates').innerHTML = ''; $('why').innerHTML = '';
  $('stressBtn').disabled = true;
  renderer.domElement.style.display = 'none';
}
async function load() {
  data = await loadData();
  try { const r = await fetch('./evidence.json'); if (r.ok) evidence = await r.json(); } catch { evidence = null; }
  card.evidence = evidence;
  renderSources();          // now that the audit's own citations are available
  const getJ = async (u) => { try { const r = await fetch(u); return r.ok ? await r.json() : null; } catch { return null; } };
  [ASSETS.ribo, ASSETS.constants] = await Promise.all([getJ('./struct/43s.json'), getJ('./struct/constants.json')]);
  if (viewMode === 'struct') $('viewNote').innerHTML = structNoteHTML();
  // which genes the ① Transcription page can model: needed for the cross-link, before that page
  // (and its 3D view) is ever loaded
  const sj = await getJ('./splice_data.json');
  spliceGenes = new Set((sj?.genes ?? []).filter(g => (g.exons ?? []).length >= 2 && (g.isoforms ?? []).length).map(g => g.gene));
  const sel = $('gene');
  fillGenes();
  sel.onchange = () => { geneTouched = true; xMsg = null; gene = genes[+sel.value]; pickDefaultPair(); setup(); };
  $('allGenes').onchange = () => { if (!fillGenes()) { gene = genes[0]; pickDefaultPair(); setup(); } };
  $('fit').innerHTML = fitCardHTML(data.meta);
  $('mlegend').innerHTML = data.genes.some(g => g.isoforms.some(i => Array.isArray(i.motifs))) ? legendHTML() : '';
  gene = genes[0]; pickDefaultPair(); setup();
  $('prov').textContent = `Structures: ViennaRNA MFE. Transcripts: ${data.meta?.gencode_or_ensembl_release ?? 'Ensembl'}. Numbers come from the exact (expected-value) scanning model; ribosomes are a stochastic animation of the same rules. Model, not measurement.`;
  onParams('init');
  targetTab.update({ data, evidence });   // both levers panels read the same transcripts + audit
  layout();
}
// ================= pages =================
// Two pages, split by biology, each with its own URL so it can be opened directly and shared:
//   #/transcription — Pol II + co-transcriptional splicing (src/splice/view.js): which isoforms exist
//   #/translation   — ribosome scanning (this file): protein per mRNA, per isoform
// The three-tab version's URLs (#/make, #/use) and its stored localStorage values ('make', 'use',
// 'target') still resolve, so existing links and restored state do not break.
const PAGE_ALIAS = {
  transcription: 'transcription', translation: 'translation',
  make: 'transcription', use: 'translation', target: 'translation',
};
const normPage = (s) => PAGE_ALIAS[String(s ?? '').toLowerCase().replace(/^#/, '').replace(/^\//, '').replace(/\/+$/, '')] ?? null;
const hashPage = () => normPage(location.hash);
let page = 'translation', spliceTab = null, spliceLoading = null;
// Page 1 opens on its own clean default gene; it only follows page 2 once the user has picked a gene
// there (so the demo path does not land on a gene where the splicing model is incoherent).
let geneTouched = false;
try { page = normPage(localStorage.getItem('utr-page')) ?? normPage(localStorage.getItem('utr-tab')) ?? page; } catch {}
page = hashPage() ?? page;

async function setPage(p, opts = {}) {
  p = normPage(p) ?? 'translation';
  page = p;
  try { localStorage.setItem('utr-page', p); localStorage.setItem('utr-tab', p === 'transcription' ? 'make' : 'use'); } catch {}
  const want = '#/' + p;
  // keep the URL in step without pushing a history entry per switch
  if (location.hash !== want) history.replaceState(null, '', location.pathname + location.search + want);
  document.body.classList.toggle('tab-make', p === 'transcription');
  document.querySelectorAll('#tabs .pagebtn').forEach(b => {
    const on = b.dataset.page === p;
    b.classList.toggle('on', on); b.setAttribute('aria-selected', on ? 'true' : 'false');
  });
  const on1 = p === 'translation';
  renderer.domElement.style.display = on1 ? '' : 'none';
  labels.domElement.style.display = on1 ? '' : 'none';
  renderer.setAnimationLoop(on1 ? tick : null);   // only one stage animates
  // the levers panel is mounted on both pages; only this page's instance is live
  targetTab.setPage(p, opts.levers ? { levers: true } : {});
  if (on1) { spliceTab?.setActive(false); renderXLink(); layout(); return; }
  if (!spliceTab) {
    // one loader, however many times the user switches while it is in flight
    spliceLoading ??= (async () => { const m = await import('./splice/view.js'); const t = new m.SpliceTab(); await t.init(); return t; })();
    const t = await spliceLoading;
    if (!spliceTab) { spliceTab = t; wireSpliceLinks(); }
  }
  if (page !== 'transcription') { spliceTab.setActive(false); return; }  // switched away while loading
  spliceTab.setSimple(simple);
  if (gene && geneTouched) spliceTab.matchGene(gene.gene);
  spliceTab.setActive(true);
}
addEventListener('hashchange', () => { const p = hashPage(); if (p && p !== page) setPage(p); });
addEventListener('resize', () => spliceTab?.resize());

// ---------- the two pages are linked: the same isoform, seen at both stages ----------
let spliceGenes = null;   // gene symbols that have a transcription model (splice_data.json)
let xMsg = null;          // "this isoform is not on the other page" note, if any
// page 2 → page 1: "← where does this isoform come from?"
function renderXLink() {
  const el = $('xlink');
  if (!el) return;
  if (!gene || !spliceGenes) { el.innerHTML = xMsg ? `<span class="xnone">${xMsg}</span>` : ''; return; }
  const iso = gene.isoforms[idxA];
  el.innerHTML = (spliceGenes.has(gene.gene)
    ? `<button type="button" class="xgo" title="show ${esc(gene.gene)} on the ① Transcription page: Pol II + splicing, where this isoform is made">← where does this isoform come from?</button>`
    : `<span class="xnone"><b>${esc(gene.gene)} is on this page only.</b> <code>splice_data.json</code> carries no exon/intron model for it,
         so the ① Transcription page cannot show where ${esc(iso?.name ?? 'this isoform')} comes from.</span>`)
    + (xMsg ? `<span class="xnone">${xMsg}</span>` : '');
  el.querySelector('.xgo')?.addEventListener('click', () => { geneTouched = true; setPage('transcription'); });
}
// page 1 → page 2: "→ translate this isoform". Returns null on success, else the reason, which
// page 1 prints next to its isoform rows (a gene on one page and not the other must say so).
function openIsoform(geneSym, isoName) {
  if (!data) return 'the translation data (<code>utr_data.json</code>) failed to load, so this isoform cannot be translated here';
  const g = data.genes.find(x => x.gene === geneSym);
  if (!g) return `<b>${esc(geneSym)} is on the ① Transcription page only.</b> <code>utr_data.json</code> has no 5′ leader for it, so the ② Translation page cannot run the scanning model on it.`;
  let i = genes.indexOf(g);
  if (i < 0) { $('allGenes').checked = true; fillGenes(); i = genes.indexOf(g); }
  if (i < 0) return `<b>${esc(geneSym)} cannot be selected on the ② Translation page.</b>`;
  geneTouched = true;
  $('gene').value = i; gene = genes[i];
  pickDefaultPair();
  const j = gene.isoforms.findIndex(iso => iso.name === isoName);
  if (j >= 0 && j !== idxA) { if (j === idxB) idxB = idxA; idxA = j; }
  xMsg = j >= 0 ? null
    : `<b>${esc(isoName)} has no 5′ leader in <code>utr_data.json</code></b>, so lane A shows ${esc(gene.isoforms[idxA]?.name ?? '—')} instead — this is <b>not</b> the isoform selected on the ① Transcription page.`;
  setup();
  setPage('translation');
  return null;
}
function wireSpliceLinks() {
  if (!spliceTab) return;
  spliceTab.onTranslate = (geneSym, isoName) => openIsoform(geneSym, isoName);
}

// ---------- therapeutic levers: ONE module, mounted on both pages (src/target/panel.js) ----------
targetTab.mount(document.body, {
  setPage,
  openGene(name, fromPage) {            // "open this gene on this page"
    if (fromPage === 'transcription') { setPage('transcription').then(() => spliceTab?.matchGene(name)); return true; }
    let i = genes.findIndex(g => g.gene === name);
    if (i < 0) { $('allGenes').checked = true; fillGenes(); i = genes.findIndex(g => g.gene === name); }
    if (i < 0) return false;
    geneTouched = true; xMsg = null;
    $('gene').value = i; gene = genes[i]; pickDefaultPair(); setup();
    setPage('translation');
    return true;
  },
});

syncStressUI();
document.body.classList.toggle('simple', simple);
document.querySelectorAll('#modeToggle button').forEach(b => b.classList.toggle('on', (b.dataset.m === 'simple') === simple));
document.body.classList.toggle('tab-make', page === 'transcription');
document.querySelectorAll('#tabs .pagebtn').forEach(b => b.classList.toggle('on', b.dataset.page === page));
setView(viewMode);
load().catch(e => { console.error('main: data load failed', e); dataError(e.message); })
  .finally(() => setPage(page));   // both pages mount even when the translation data is missing
