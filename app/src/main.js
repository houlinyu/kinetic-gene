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
import { expected } from './analytic.js';
import { Scorecard } from './scorecard.js';

const STRESS_TC = 0.15; // matches findings.js
const CALM_TC = 1;
const $ = (id) => document.getElementById(id);
const fmt = (x, d = 1) => (x == null || !Number.isFinite(x) ? '—' : Number(x).toFixed(d));
const fmtR = (x) => (!Number.isFinite(x) ? '—' : x >= 10 ? x.toFixed(1) : x >= 1 ? x.toFixed(2) : x.toPrecision(2));
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

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
const sim = { speed: 3 };
const LANE_Y = [1.0, -2.6];
let idxA = 0, idxB = 1;
const target = { A: 0, B: 0, ratio: 1, stress: 0 };
const shown = { A: 0, B: 0, ratio: 1, stress: 0 };
const stressed = () => P.tc < 0.5;

// ---------- layout: fit the strands into the free area (right of the panel, left of the findings) ----------
function layout() {
  const W = innerWidth, H = innerHeight, mobile = W <= 760;
  renderer.setSize(W, H); labels.setSize(W, H); composer.setSize(W, H);
  let x0 = 0, x1 = W, y0 = 0, y1 = H;
  const pr = $('panel').getBoundingClientRect(), sr = $('side').getBoundingClientRect();
  if (mobile) { y0 = 48; y1 = pr.top - 8; }
  else { x0 = pr.right + 8; x1 = (getComputedStyle($('side')).display === 'none' ? W : sr.left) - 8; y0 = 40; y1 = H - 44; }
  const rw = Math.max(100, x1 - x0), rh = Math.max(100, y1 - y0);
  const rcx = (x0 + x1) / 2, rcy = (y0 + y1) / 2;
  const dx = W / 2 - rcx, dy = H / 2 - rcy;
  const fullW = W + 2 * Math.abs(dx), fullH = H + 2 * Math.abs(dy);
  camera.aspect = fullW / fullH;
  camera.setViewOffset(fullW, fullH, Math.abs(dx) + dx, Math.abs(dy) + dy, W, H);
  // scene bounds: strands x ∈ [-7, 7] plus "cap ▶" label; lanes + titles in y
  const sceneW = mobile ? 17.2 : 15.6, sceneH = (LANE_Y[0] + 1.9) - (LANE_Y[1] - 0.9);
  const cx = 0, cy = (LANE_Y[0] + 1.9 + LANE_Y[1] - 0.9) / 2;
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
function extraLines(iso) {
  const out = [];
  const m = iso.measured;
  if (m && m.values && Object.keys(m.values).length) {
    const cell = 'HEK293T' in m.values ? 'HEK293T' : Object.keys(m.values)[0];
    const v = m.values[cell];
    if (v != null) out.push(`Measured TE (Cao 2021, ${esc(cell)}): <b>${fmt(v, 2)}</b>${m.ambiguous ? ' <span class="amb" title="UTR shared with other isoforms in the measurement">(ambiguous)</span>' : ''}`);
  }
  const mrl = iso.pred?.optimus_mrl;
  if (mrl != null) out.push(`Optimus 5′ MRL: <b>${fmt(mrl, 2)}</b>`);
  return out.length ? `<div class="meas">${out.join('<br>')}</div>` : '';
}
function cardHTML(tag, iso, idx) {
  return `<div class="lhead"><span class="tag">${tag}</span><select data-lane="${tag}">${isoOptions(idx)}</select></div>
    <div class="rate"><span id="rate${tag}">—</span> <small>proteins / min per mRNA</small></div>
    <div class="cs" id="cs${tag}"></div>${extraLines(iso)}`;
}
function statsHTML(tag, iso) {
  const nU = iso.uorfs?.length || 0, nOver = (iso.uorfs || []).filter(u => u.type === 'overlapping').length;
  return `<div><span class="tag">${tag}</span> ${esc(iso.name)}: ${iso.len} nt · ${nU} uORF${nU === 1 ? '' : 's'}${nOver ? ` (${nOver} overlapping)` : ''} · Kozak ${iso.main_kozak?.class ?? '—'} ·
    MFE ${fmt(iso.fold?.mfe)} kcal/mol · cap-50nt ΔG ${fmt(iso.fold?.cap_dG)} · GC ${fmt((iso.gc ?? NaN) * 100, 0)}%</div>`;
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
  label(view, 'lane-title', `${tag} · ${iso.name} (${isMane(iso) ? 'MANE, standard' : 'alternative'})`, view.x0, 1.55, 0, 1);
  label(view, 'lane-cap', 'cap ▶', view.x0 - 0.08, 0, 1, 0.5);
  label(view, 'lane-cds', 'protein-coding ▶', view.X(iso.len + CDS_STUB / 2), 0.5, 0.5, 1);
}

function prewarm(lane, secs = 45) {
  for (let t = 0; t < secs; t += 0.05) lane.step(0.05, P);
}

function setup() {
  dropView(viewA); dropView(viewB);
  const isoA = gene.isoforms[idxA], isoB = gene.isoforms[idxB] ?? isoA;
  const maxLen = Math.max(isoA.len, isoB.len);
  laneA = new Lane(isoA); laneB = new Lane(isoB);
  prewarm(laneA); prewarm(laneB);
  viewA = new LaneView(scene, LANE_Y[0], maxLen); viewA.build(laneA); decorate(viewA, 'A', isoA);
  viewB = new LaneView(scene, LANE_Y[1], maxLen); viewB.build(laneB); decorate(viewB, 'B', isoB);
  laneA.events.length = 0; laneB.events.length = 0;
  $('cardA').innerHTML = cardHTML('A', isoA, idxA);
  $('cardB').innerHTML = cardHTML('B', isoB, idxB);
  $('stats').innerHTML = statsHTML('A', isoA) + statsHTML('B', isoB);
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
    $('cs' + tag).innerHTML = `<span class="${stressed() ? '' : 'on'}">calm ${fmtR(c)}</span> → <span class="${stressed() ? 'on' : ''}">stress ${fmtR(s)}</span> ${arrow} <span class="mute">×${fmtR(r)}</span>`;
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
  $('fates').innerHTML = `<div class="ftitle">What happens to each ribosome that loads</div>${bar(fatesCache.A, 'A')}${bar(fatesCache.B, 'B')}
    <div class="fkey">${Object.keys(FCOL).map(k => `<span><i style="background:${FCOL[k]}"></i>${FNAME[k]}</span>`).join('')}</div>`;
}

function renderNumbers() {
  $('rateA').textContent = fmtR(shown.A);
  $('rateB').textContent = fmtR(shown.B);
  const r = shown.ratio;
  $('verdict').innerHTML = `<div class="vbig"><b>${fmtR(r)}×</b><span>as much protein<br>from <span class="tag">B</span> as from <span class="tag">A</span></span></div>
    <div class="vsub">exact model · ${stressed() ? '<span class="hot">stressed cell (eIF2·TC low)</span>' : 'calm cell'}</div>`;
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
  $('sideToggle').textContent = `${s.passed}/${s.total} findings ${document.body.classList.contains('side-open') ? '▾' : '▸'}`;
  if (!['tc', 'tcRate', 'loadTCexp'].includes(k)) { clearTimeout(mapTimer); mapTimer = setTimeout(() => card.recompute(data, P), 250); }
}
$('sideToggle').onclick = () => { document.body.classList.toggle('side-open'); onParams('tc'); };

// ---------- parameter GUI (Advanced, closed) ----------
const gui = new GUI({ title: 'Advanced: all parameters', container: $('guiHost') });
gui.close();
gui.add(sim, 'speed', 0.25, 10, 0.25).name('animation speed ×');
for (const k of Object.keys(PARAMS)) {
  const d = PARAMS[k];
  gui.add(P, k, d.min, d.max, d.step).name(`${d.label}${d.cite ? '' : ' (assumed)'}`).onChange(() => onParams(k));
}
$('sources').innerHTML = `<div class="ftitle">Sources</div><ul>${Object.entries(SOURCES).map(([k, s]) =>
  `<li><b>${k}</b>: ${s.url ? `<a href="${s.url}" target="_blank" rel="noopener">${esc(s.text)}</a>` : esc(s.text)}</li>`).join('')}</ul>
  <div class="fine">Assumed (no citation yet): ${Object.entries(PARAMS).filter(([, d]) => !d.cite).map(([, d]) => esc(d.label)).join(' · ')}</div>`;

// ---------- picking: hover cursor + floating tag ----------
const ray = new THREE.Raycaster(), mouse = new THREE.Vector2();
let downAt = null, tipTimer = null;
const tip = $('tip');
function pick(e) {
  if (!viewA) return null;
  mouse.set((e.clientX / innerWidth) * 2 - 1, -(e.clientY / innerHeight) * 2 + 1);
  ray.setFromCamera(mouse, camera);
  const hits = ray.intersectObjects([...viewA.pickables, ...viewB.pickables], false);
  return hits[0]?.object.userData ?? null;
}
function showTip(e, text, hold) {
  tip.textContent = text; tip.style.left = e.clientX + 14 + 'px'; tip.style.top = e.clientY - 10 + 'px';
  tip.classList.add('on'); tip.classList.toggle('strong', !!hold);
  clearTimeout(tipTimer);
  if (hold) tipTimer = setTimeout(() => tip.classList.remove('on', 'strong'), 1800);
}
renderer.domElement.addEventListener('pointermove', e => {
  if (tip.classList.contains('strong')) return;
  const u = pick(e);
  renderer.domElement.style.cursor = u ? 'pointer' : '';
  if (!u) { tip.classList.remove('on'); return; }
  const lane = u.lane.lane;
  showTip(e, u.kind === 'uorf'
    ? (lane.uorfs[u.idx].enabled ? `uORF (${lane.uorfs[u.idx].type}) · click to mutate its AUG` : 'mutated uORF · click to restore')
    : (lane.hairpins[u.idx].melted ? 'melted hairpin · click to restore' : `hairpin ΔG ${fmt(lane.hairpins[u.idx].dG)} kcal/mol · click to melt`));
});
renderer.domElement.addEventListener('pointerdown', e => downAt = [e.clientX, e.clientY]);
renderer.domElement.addEventListener('pointerup', e => {
  if (!downAt || Math.hypot(e.clientX - downAt[0], e.clientY - downAt[1]) > 5) return;
  const u = pick(e); if (!u) return;
  const { kind, idx, lane: view } = u, lane = view.lane;
  let msg;
  if (kind === 'uorf') { const on = !lane.uorfs[idx].enabled; lane.setUorf(idx, on); msg = on ? 'uORF restored' : 'AUG mutated: uORF removed'; }
  else { const m = !lane.hairpins[idx].melted; lane.setHairpin(idx, m); msg = m ? 'hairpin melted' : 'hairpin restored'; }
  prewarm(lane, 30); lane.events.length = 0;
  view.refreshToggles();
  showTip(e, msg, true);
  recalc();
});
addEventListener('keydown', e => {
  if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
  if (e.code === 'Space') { paused = !paused; e.preventDefault(); }
  if (e.code === 'KeyS') setStress(!stressed());
});

// ---------- loop ----------
let last = performance.now(), uiClock = 0;
const fogCol = new THREE.Color();
renderer.setAnimationLoop((now) => {
  const dt = Math.min((now - last) / 1000, 0.05); last = now;
  if (laneA && !paused) {
    const simDt = dt * sim.speed, sub = Math.ceil(simDt / 0.02);
    for (let i = 0; i < sub; i++) { laneA.step(simDt / sub, P); laneB.step(simDt / sub, P); }
    viewA.update(dt); viewB.update(dt);
  }
  // tweens
  const k = 1 - Math.exp(-dt * 5);
  for (const key of ['A', 'B', 'ratio']) shown[key] += (target[key] - shown[key]) * k;
  shown.stress += (target.stress - shown.stress) * (1 - Math.exp(-dt * 3));
  scene.background.copy(BG_CALM).lerp(BG_STRESS, shown.stress);
  scene.fog.color.copy(fogCol.copy(BG_CALM).lerp(BG_STRESS, shown.stress));
  uiClock += dt;
  if (uiClock > 0.06 && laneA) { uiClock = 0; renderNumbers(); }
  controls.update();
  composer.render(); labels.render(scene, camera);
});

// ---------- data ----------
async function load() {
  let res = await fetch('./utr_data.json');
  if (!res.ok) res = await fetch('./mock_data.json');
  data = await res.json();
  // the model can't start at a non-AUG MANE start codon (MYC, VEGFA use CUG): hide those genes
  genes = data.genes.filter(g => !(g.isoforms.find(isMane) ?? g.isoforms[0])?.tags?.some(t => t.startsWith('non_AUG_start')));
  const sel = $('gene');
  sel.innerHTML = genes.map((g, i) => `<option value="${i}">${g.gene} · ${g.isoforms.length} isoforms</option>`).join('');
  sel.onchange = () => { gene = genes[+sel.value]; pickDefaultPair(); setup(); };
  gene = genes[0]; pickDefaultPair(); setup();
  $('prov').textContent = `Structures: ViennaRNA MFE. Transcripts: ${data.meta?.gencode_or_ensembl_release ?? 'Ensembl'}. Numbers come from the exact (expected-value) scanning model; ribosomes are a stochastic animation of the same rules. Model, not measurement.`;
  onParams('init');
  layout();
}
syncStressUI();
layout();
load();
