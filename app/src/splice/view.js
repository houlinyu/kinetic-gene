// Tab 2 "Make it": transcription + co-transcriptional splicing.
// SIMPLE view = an SVG gene diagram (promoter + TSS, exon boxes, introns, splice sites with their
// scores, junction arcs whose thickness is the predicted usage), a Pol II bead per nascent
// transcript with the RNA trailing behind it, introns looping out and being cut, finished
// transcripts dropping into per-isoform bins.
// 3D view = app/src/splice/view3d.js (true molecular scale).
//
// Every DISPLAYED number comes from model_exact.expected(); the beads are model.js's stochastic
// simulation of the same rules. Data: app/public/splice_data.json (data/SPEC_SPLICE.md).
// model.js / model_exact.js / params.js are owned by the SPLICE-MODEL agent and are picked up at
// build time; if any of them is missing this file falls back to its own clearly-labelled
// competition model so the page still renders.

const MODS = import.meta.glob('./*.js');
// ISOFORM DIVERSITY: the combinatorial space (promoter × splice-site choices), the live
// per-transcript identity chips, the divergence watch and the ② Translation join.
// Pure module owned alongside this file; every probability it prints comes from expected().
import * as DIV from './diversity.js';
// Wall-clock pacing (view only). The model keeps nt/s and seconds; this decides how many model
// seconds pass per second of your time, and the factor is printed live in #spaceNote.
import * as pace from '../pace.js';
const $ = (id) => document.getElementById(id);
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmt = (x, d = 1) => (x == null || !Number.isFinite(x) ? '—' : Number(x).toFixed(d));
const pct = (x) => (Number.isFinite(x) ? (100 * x).toFixed(0) + '%' : '—');
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const SVGNS = 'http://www.w3.org/2000/svg';
const ntSort = (a, b) => a.start - b.start;
// seconds -> "5 h 19 min" / "4 min 34 s" / "42 s", for the time-compression line
const dur = (s) => (s >= 3600 ? `${Math.floor(s / 3600)} h ${Math.round((s % 3600) / 60)} min`
  : s >= 60 ? `${Math.floor(s / 60)} min ${Math.round(s % 60)} s` : `${s.toFixed(s < 10 ? 1 : 0)} s`);

// ---------------------------------------------------------------- fallback parameters
// Only used when splice/params.js is missing. Same contract as app/src/params.js:
// every value is sourced (cite) or labelled "assumed".
const FALLBACK_PARAMS = {
  elongRate:  { value: 25,   min: 8,   max: 67,  step: 1,    unit: 'nt/s', label: 'Pol II elongation rate', cite: 'jonkers2014' },
  pausePos:   { value: 50,   min: 0,   max: 200, step: 5,    unit: 'nt',   label: 'Promoter-proximal pause site', cite: 'core2019' },
  pauseEscape:{ value: 0.002, min: 0.0002, max: 0.2, step: 0.0002, unit: '/s', label: 'Pause escape rate', cite: null },
  initRate:   { value: 0.02, min: 0.002, max: 0.5, step: 0.002, unit: '/s', label: 'Initiation rate', cite: null },
  scoreSlope: { value: 1.0,  min: 0,   max: 4,   step: 0.05, unit: '',     label: 'Splice-site strength scale', cite: null },
  commitScale:{ value: 0.02, min: 0.002, max: 0.3, step: 0.002, unit: '/s', label: 'Commitment rate of a median-strength site', cite: 'singh2009' },
  removalTime:{ value: 480,  min: 30,  max: 1200, step: 10,  unit: 's',    label: 'Intron removal time after commitment', cite: 'singh2009' },
  postTx:     { value: 0.2,  min: 0,   max: 1,   step: 0.01, unit: '',     label: 'Fraction spliced post-transcriptionally (1 = no kinetic coupling)', cite: null },
};
const FALLBACK_SOURCES = {
  jonkers2014: { text: 'Jonkers, Kwak & Lis 2014 eLife: Pol II elongation rates of ~0.5–4 kb/min across genes', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3977343/' },
  core2019: { text: 'Core & Adelman 2019 Genes Dev: promoter-proximal Pol II pauses ~20–60 nt downstream of the TSS', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6690705/' },
  singh2009: { text: 'Singh & Padgett 2009 Nat Struct Mol Biol: introns are removed co-transcriptionally, within minutes of synthesis', url: 'https://pubmed.ncbi.nlm.nih.gov/19377482/' },
  delamata2003: { text: 'de la Mata et al. 2003 Mol Cell: a slow Pol II mutant increases inclusion of the fibronectin EDI exon — the kinetic-coupling experiment this tab animates', url: 'https://pubmed.ncbi.nlm.nih.gov/14536091/' },
};
// which parameter plays which role in the four headline controls
const ROLES = {
  elong:    ['elongRate', 'elong', 'elongation', 'vPol'],
  pause:    ['pauseEscape', 'pauseDwell', 'pauseProb'],
  strength: ['scoreSlope', 'ssScale', 'commitScale'],
  cotxn:    ['postTx', 'coTxn'],
};
const kbmin = (v, unit) => (String(unit).includes('nt/s') ? v * 60 / 1000 : v);

// ---------------------------------------------------------------- fallback model (only if the real one is missing)
const kOf = (j, P) => {
  const sc = [j.donor?.score, j.acceptor?.score].filter(Number.isFinite);
  const mean = sc.length ? sc.reduce((a, b) => a + b, 0) / sc.length : 8;
  return Math.max(1e-9, (P.commitScale ?? 0.02) * Math.exp((P.scoreSlope ?? 1) * (mean - 8) / 3));
};
const conflict = (a, b) => a !== b && a.start < b.end && b.start < a.end;
function fbGroups(J) {
  const seen = new Set(), out = [];
  for (const i of J) {
    if (seen.has(i)) continue;
    const grp = [i]; seen.add(i);
    for (const j of J) if (!seen.has(j) && grp.some(k => conflict(k, j))) { grp.push(j); seen.add(j); }
    out.push(grp);
  }
  return out;
}
// closed-form race: an intron can only commit once Pol II has passed its 3'ss, commitment is
// exponential, the first one to commit kills its competitors. Slow Pol II stretches the head
// start a proximal site gets -> kinetic coupling.
function fbRace(ks, ts, horizon) {
  const m = ks.length, ord = [...ks.keys()].sort((a, b) => ts[a] - ts[b]);
  const p = new Array(m).fill(0);
  let C = 1, K = 0; const act = [];
  for (let s = 0; s < m; s++) {
    const i = ord[s];
    act.push(i); K += ks[i];
    const T2 = s + 1 < m ? ts[ord[s + 1]] : horizon;
    const dT = Math.max(0, Math.min(T2, horizon) - ts[i]);
    if (K > 0 && dT > 0) {
      const f = (1 - Math.exp(-K * dT)) / K;
      for (const j of act) p[j] += ks[j] * C * f;
      C *= Math.exp(-K * dT);
    }
    if (T2 >= horizon) break;
  }
  return { p, none: Math.max(0, C) };
}
function fbExpected(g, J, P) {
  const v = P.elongRate ?? 25, span = Math.max(...J.map(j => j.end), ...(g.exons ?? []).map(e => e.end), 1);
  const tEnd = span / Math.max(v, 1e-6), horizon = tEnd + (P.removalTime ?? 480) * 2;
  const win = new Map(), none = new Map();
  for (const grp of fbGroups(J)) {
    const ks = grp.map(j => kOf(j, P));
    const ts = grp.map(j => (1 - (P.postTx ?? 0)) * Math.min(tEnd, j.end / Math.max(v, 1e-6)) + (P.postTx ?? 0) * tEnd);
    const r = fbRace(ks, ts, horizon);
    grp.forEach((j, n) => { win.set(j.id, r.p[n]); none.set(j.id, r.none); });
  }
  const fractions = {}, junctionUse = {}, exonUse = {};
  for (const j of J) junctionUse[j.id] = win.get(j.id) ?? 0;
  for (const iso of g.isoforms ?? []) {
    const set = new Set(iso.exon_ids ?? []);
    let w = 1;
    for (const grp of fbGroups(J)) {
      const used = grp.find(j => set.size && (g.exons ?? []).some(e => set.has(e.id) && e.end === j.start) && (g.exons ?? []).some(e => set.has(e.id) && e.start === j.end)
        && !(g.exons ?? []).some(e => set.has(e.id) && e.start >= j.start && e.end <= j.end));
      w *= used ? (win.get(used.id) ?? 0) : (none.get(grp[0].id) ?? 0);
    }
    fractions[iso.id] = w;
  }
  const tot = Object.values(fractions).reduce((a, b) => a + b, 0) || 1;
  for (const k of Object.keys(fractions)) fractions[k] /= tot;
  for (const e of g.exons ?? []) {
    exonUse[e.id] = (g.isoforms ?? []).filter(i => (i.exon_ids ?? []).includes(e.id)).reduce((a, i) => a + (fractions[i.id] ?? 0), 0);
  }
  const psi = (g.exons ?? []).filter(e => e.type === 'cassette').map(e => ({ event: `${e.id} inclusion`, exon: e.id, value: exonUse[e.id] }));
  return { fractions, isoforms: (g.isoforms ?? []).map(i => ({ id: i.id, name: i.name, value: fractions[i.id] ?? 0 })), junctionUse, exonUse, psi, fallback: true };
}
// fallback stochastic sim, same rules sampled
class FbSim {
  constructor(g, J) { this.g = g; this.J = J; this.span = Math.max(...J.map(j => j.end), ...(g.exons ?? []).map(e => e.end), 1); this.reset(); }
  reset() { this.t = 0; this.pols = []; this.counts = {}; this.events = []; this.nextId = 1; }
  step(dt, P) {
    this.t += dt;
    if (Math.random() < 1 - Math.exp(-(P.initRate ?? 0.02) * dt) && this.pols.length < 10)
      this.pols.push({ id: this.nextId++, pos: 0, tssPos: 0, released: !(Math.random() < 0.9), pauseAt: P.pausePos ?? 50, committed: [], removed: [], skipped: [], pending: [], txDone: false, postTx: Math.random() < (P.postTx ?? 0.2) });
    for (const r of this.pols) {
      if (!r.txDone) {
        if (!r.released) { if (r.pos < r.pauseAt) r.pos = Math.min(r.pauseAt, r.pos + (P.elongRate ?? 25) * dt); else if (Math.random() < 1 - Math.exp(-(P.pauseEscape ?? 0.002) * dt)) r.released = true; }
        else r.pos += (P.elongRate ?? 25) * dt;
        if (r.pos >= this.span) { r.pos = this.span; r.txDone = true; }
      }
      for (const j of this.J) {
        if (r.committed.includes(j.idx) || r.skipped.includes(j.idx)) continue;
        if (!(r.postTx ? r.txDone : r.pos >= j.end)) continue;
        if (Math.random() < 1 - Math.exp(-kOf(j, P) * dt)) {
          r.committed.push(j.idx); r.pending.push({ idx: j.idx, at: this.t + (P.removalTime ?? 480) });
          for (const o of this.J) if (conflict(o, j) && !r.skipped.includes(o.idx)) r.skipped.push(o.idx);
          this.events.push({ kind: 'commit', intron: j.id, pol: r.id });
          break;
        }
      }
      for (let i = r.pending.length - 1; i >= 0; i--) if (r.pending[i].at <= this.t) {
        r.removed.push(r.pending[i].idx);
        this.events.push({ kind: 'spliced', intron: this.J[r.pending[i].idx].id, pol: r.id });
        r.pending.splice(i, 1);
      }
      if (r.txDone && !r.pending.length && this.J.every(j => r.committed.includes(j.idx) || r.skipped.includes(j.idx) || (r.postTx ? false : false))) {
        const ids = new Set(r.committed.map(i => this.J[i].key));
        const iso = (this.g.isoforms ?? []).find(o => {
          const set = new Set(o.exon_ids ?? []);
          const chain = (this.g.exons ?? []).filter(e => set.has(e.id)).sort(ntSort);
          const sig = new Set(chain.slice(0, -1).map((e, k) => `${e.end}_${chain[k + 1].start}`));
          return sig.size === ids.size && [...sig].every(x => ids.has(x));
        });
        const key = iso?.id ?? 'other';
        this.counts[key] = (this.counts[key] ?? 0) + 1;
        this.events.push({ kind: 'finished', pol: r.id, isoform: key });
        r.done = true;
      }
    }
    this.pols = this.pols.filter(r => !r.done);
    if (this.events.length > 300) this.events.splice(0, this.events.length - 300);
  }
  view() { return this.pols.map(r => ({ id: r.id, pos: r.pos, tssPos: r.tssPos, paused: !r.released, txDone: r.txDone, committed: r.committed.map(i => this.J[i].id), removed: r.removed.map(i => this.J[i].id), skipped: r.skipped.map(i => this.J[i].id) })); }
}

// ---------------------------------------------------------------- x mapping (introns compressed)
function makeMap(g, J, W, padL, padR) {
  const ex = (g.exons ?? []).slice().sort(ntSort);
  const s0 = Math.min(...ex.map(e => e.start), ...(g.tss ?? []).map(t => t.pos ?? Infinity));
  const s1 = Math.max(...ex.map(e => e.end), ...J.map(j => j.end));
  const lead = Math.max(60, Math.round((s1 - s0) * 0.02));
  const segs = [{ s: s0 - lead, e: s0, lead: true }];
  for (let i = 0; i < ex.length; i++) {
    segs.push({ s: ex[i].start, e: ex[i].end, exon: ex[i] });
    if (i + 1 < ex.length && ex[i + 1].start > ex[i].end) segs.push({ s: ex[i].end, e: ex[i + 1].start, gap: true });
  }
  if (segs[segs.length - 1].e < s1) segs.push({ s: segs[segs.length - 1].e, e: s1, gap: true });
  const exNt = segs.filter(s => s.exon).reduce((a, s) => a + (s.e - s.s), 0) || 1;
  const gaps = segs.filter(s => s.gap);
  const gw = gaps.reduce((a, s) => a + Math.log10(1 + (s.e - s.s)), 0) || 1;
  const avail = W - padL - padR, fExon = gaps.length ? 0.5 : 0.95;
  for (const s of segs) {
    s.px = s.lead ? avail * 0.04
      : s.exon ? (s.e - s.s) / exNt * avail * fExon
        : Math.log10(1 + (s.e - s.s)) / gw * avail * (0.96 - fExon);
  }
  // every exon gets a readable minimum width; the gaps give up the difference
  const MINEX = 15;
  const short = segs.filter(s => s.exon && s.px < MINEX);
  if (short.length) {
    const need = short.reduce((a, s) => a + (MINEX - s.px), 0);
    const pool = segs.filter(s => !s.exon || s.px >= MINEX).reduce((a, s) => a + s.px, 0) || 1;
    const take = Math.min(need, pool * 0.55);
    for (const s of segs) {
      if (s.exon && s.px < MINEX) s.px = MINEX;
      else s.px -= take * s.px / pool;
    }
  }
  const tot = segs.reduce((a, s) => a + s.px, 0) || 1;
  let x = padL;
  for (const s of segs) { s.px *= avail / tot; s.x0 = x; x += s.px; s.x1 = x; }
  const map = (pos) => {
    if (pos <= segs[0].s) return segs[0].x0;
    for (const s of segs) if (pos <= s.e) return s.x0 + (pos - s.s) / Math.max(1, s.e - s.s) * s.px;
    return segs[segs.length - 1].x1;
  };
  map.segs = segs; map.s0 = s0; map.s1 = s1;
  map.compressed = gaps.some(s => (s.e - s.s) > 400);
  return map;
}

// How much faster the bead moves through the most expanded segment of the compressed x axis than
// it would on a linear one (1 = linear). Clamped: ESR1 would otherwise read ~50x.
function burstOf(m) {
  if (!m?.segs?.length) return 1;
  if (m._burst != null) return m._burst;
  const totPx = m.segs.reduce((a, s) => a + s.px, 0) || 1;
  const totNt = m.segs.reduce((a, s) => a + Math.max(1, s.e - s.s), 0) || 1;
  const mean = totPx / totNt;
  let worst = mean;
  for (const s of m.segs) if (!s.lead) worst = Math.max(worst, s.px / Math.max(1, s.e - s.s));
  m._burst = Math.min(Math.max(worst / mean, 1), 60);
  return m._burst;
}

// ---------------------------------------------------------------- svg helpers
const el = (tag, attrs, parent) => {
  const n = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs ?? {})) n.setAttribute(k, v);
  if (parent) parent.appendChild(n);
  return n;
};
const COL = { exon: '#e9e4d8', alt: '#ffb547', intron: '#6a7183', pol: '#5ee0a0', rna: '#e85d9c', don: '#4fc3f7', acc: '#c39bff' };
const VIEW_W = 1000, VIEW_H = 520;
const SUBSTEP = 0.5;       // model seconds per sub-step (the race needs a fine dt)

// Genes whose predicted splicing contradicts the measured data shipped in this repo. They are off
// the default demo path and only appear behind the opt-in below, each with its reason shown.
// Default gene for the demo path: one where the model is internally coherent (the leftover
// "unannotated combinations" bucket is 0% and the canonical transcript is not at zero).
const DEFAULT_GENES = ['PPP1R15A', 'CDKN2A'];
// Parameters that src/splice/model_exact.js never reads: they only drive the bead animation.
// (verified by grep: pauseEscape / pausePos / pauseProb / initRate appear 0 times in the model)
const ANIM_ONLY = new Set(['initRate', 'pauseProb', 'pausePos', 'pauseEscape']);
const CONTRADICTED = {
  MAPT: 'the model predicts 2N4R tau at 82.3% of MAPT output, but 2N tau is the minority form (~10%) in human brain; predicted exon-3 inclusion 0.952 vs 0.062 in VastDB.',
  MDM2: 'the model predicts exon-10 inclusion PSI 0.856, while the measured value is 0.023 across 56 GTEx tissues (ling2020) — the prediction is sign-inverted.',
};

export class SpliceTab {
  constructor() {
    this.PARAMS = FALLBACK_PARAMS; this.SOURCES = FALLBACK_SOURCES;
    this.P = {}; for (const [k, d] of Object.entries(FALLBACK_PARAMS)) this.P[k] = d.value;
    // RELATIVE animation speed: 1 = the watchable pace pace.js solves for this gene and these
    // rates. The absolute compression it works out to is this.pace.ts, printed on the page.
    this.mul = 1; this.animRange = { min: 0.25, max: 4, step: 0.25 };
    this.pace = { ts: 150, crossSim: 1, crossWall: 1 };
    this.view = 'diagram'; this.simple = true; this.active = false;
    this.data = null; this.gene = null; this.sim = null; this.exp = null;
    this.v3 = null; this.raf = null; this.last = 0; this.t = 0;
    this.anim = []; this.drops = []; this.evAt = 0;
    // set by main.js: open one of these isoforms' 5′ leader on the ② Translation page.
    // Returns null on success, or the reason it is not there (which is printed on the rows).
    this.onTranslate = null; this.xnote = null;
    // --- isoform diversity: combinatorial space, live per-transcript identity, divergence watch
    this.sp = null;                      // DIV.space() for the current gene
    this.sigs = [];                      // one signature per live polymerase (this tick)
    this.dw = new DIV.DivergenceWatch();
    this.divPulse = [];                  // one-shot "same site, different choice" pulses
    try { if (localStorage.getItem('splice-view') === 'struct3d') this.view = 'struct3d'; } catch {}
  }

  // ---------- lifecycle ----------
  async init() {
    this.svg = $('ssvg');
    document.querySelectorAll('#spliceViewToggle button').forEach(b => b.onclick = () => this.setView(b.dataset.sv));
    $('slowBtn').onclick = () => this.setSlow(!this.isSlow());
    $('sgene').onchange = () => { this.noMatch = null; this.setGene(this.genes[+$('sgene').value]); };
    // "→ translate this isoform": delegated, because the headline re-renders several times a second
    $('sheadline').addEventListener('click', (e) => {
      const b = e.target.closest?.('button.sxl');
      if (!b) return;
      const why = this.onTranslate?.(this.gene?.gene, b.dataset.iso);
      this.xnote = why ?? null;
      this.renderPanel();
    });
    const all = $('sAllGenes');
    if (all) all.onchange = () => this.fillGenes();
    this.mountDiversity();
    await this.loadMods();
    this.renderControls(); this.renderSources();
    await this.loadData();
    this.setView(this.view);
  }
  async loadMods() {
    const get = async (k) => { try { return MODS['./' + k] ? await MODS['./' + k]() : null; } catch (e) { console.warn('splice: ' + k + ' not usable', e); return null; } };
    const [m, x, p] = await Promise.all([get('model.js'), get('model_exact.js'), get('params.js')]);
    this.Gene = m?.Gene ?? null;
    this.junctionsFn = typeof m?.junctions === 'function' ? m.junctions : null;
    this.tssListFn = typeof m?.tssList === 'function' ? m.tssList : null;
    this.expectedFn = typeof x?.expected === 'function' ? x.expected : null;
    this.calibrateFn = typeof x?.calibrate === 'function' ? x.calibrate : null;
    if (p?.PARAMS) {
      this.PARAMS = p.PARAMS;
      this.SOURCES = { ...FALLBACK_SOURCES, ...(p.SOURCES ?? {}) };
      this.P = typeof p.values === 'function' ? p.values() : {};
      for (const [k, d] of Object.entries(p.PARAMS)) if (this.P[k] == null) this.P[k] = d.value;
    }
    this.role = {};
    for (const [r, keys] of Object.entries(ROLES)) this.role[r] = keys.find(k => this.PARAMS[k]) ?? null;
    this.P0 = { ...this.P };
  }
  async loadData() {
    $('sheadline').innerHTML = '<div class="shl-load">loading…</div>';
    let d = null;
    try { const r = await fetch('./splice_data.json'); if (r.ok) d = await r.json(); } catch {}
    this.placeholder = !d || !Array.isArray(d.genes) || !d.genes.length;
    this.data = this.placeholder ? PLACEHOLDER : d;
    this.allGenes = this.data.genes.filter(g => (g.exons ?? []).length >= 2 && (g.isoforms ?? []).length);
    if (!this.allGenes.length) this.allGenes = this.data.genes;
    if (this.calibrateFn && !this.placeholder) { try { this.calibrateFn(this.data); } catch (e) { console.warn('splice: calibrate failed', e); } }
    this.fillGenes(true);
    const meta = this.data.meta ?? {};
    $('sprov').innerHTML = this.placeholder
      ? '<b class="warn">Placeholder gene geometry — <code>splice_data.json</code> was not found, so the layout below is an illustration and the splice-site scores are blank.</b>'
      : `Transcripts: ${esc(meta.release ?? 'see sources')}${meta.methods?.ss_score ? ` · splice-site scores: ${esc(meta.methods.ss_score)}` : ''}${meta.methods?.psi ? ` · PSI: ${esc(meta.methods.psi)}` : ''}.
         Percentages are expected values from the ${this.expectedFn ? 'exact splicing model' : 'fallback model'}; the beads are a stochastic simulation of the same rules, run at the compressed wall-clock tempo stated above the gene diagram. Model, not measurement.`;
    this.setGene(this.genes.find(g => DEFAULT_GENES.includes(g.gene)) ?? this.genes[0]);
  }

  // Genes where the prediction contradicts the measured data are only offered behind the opt-in.
  fillGenes(initial = false) {
    const all = !!$('sAllGenes')?.checked;
    const cur = this.gene;
    this.genes = (this.allGenes ?? []).filter(g => all || !CONTRADICTED[g.gene]);
    if (!this.genes.length) this.genes = this.allGenes ?? [];
    $('sgene').innerHTML = this.genes.map((g, i) => {
      const bad = CONTRADICTED[g.gene];
      return `<option value="${i}" title="${esc(bad ?? '')}">${esc(g.gene)} · ${(g.isoforms ?? []).length} isoforms${bad ? ' · contradicts measurement' : ''}</option>`;
    }).join('');
    this.renderGeneNotes();
    if (initial) return;
    const keep = this.genes.indexOf(cur);
    if (keep >= 0) $('sgene').value = keep; else this.setGene(this.genes[0]);
  }

  // one place for the two things the gene picker has to admit: what it hides, and what it could not show
  renderGeneNotes() {
    const warn = $('sgeneWarn'); if (!warn) return;
    const all = !!$('sAllGenes')?.checked;
    const hidden = (this.allGenes ?? []).filter(g => CONTRADICTED[g.gene]).map(g => g.gene);
    const miss = this.noMatch
      ? `<div class="sbad"><b>No splicing data for ${esc(this.noMatch)}.</b> <code>splice_data.json</code> has no exon/intron model for it, so this tab cannot show where its isoforms come from
         — it is still showing <b>${esc(this.gene?.gene ?? '—')}</b>, which is <b>not</b> the gene selected in the translation tab.</div>` : '';
    const hid = !hidden.length ? ''
      : all ? `<div class="sbad"><b>Shown with a known contradiction:</b> ${hidden.map(h => `${esc(h)} — ${esc(CONTRADICTED[h])}`).join(' · ')}</div>`
        : `<div class="sgap">${hidden.length} gene${hidden.length > 1 ? 's' : ''} (${hidden.map(esc).join(', ')}) ${hidden.length > 1 ? 'are' : 'is'} hidden here because the model's predicted splicing contradicts the measured data in this repo. Tick the box above to inspect ${hidden.length > 1 ? 'them' : 'it'} with the reason shown.</div>`;
    warn.innerHTML = miss + hid;
  }
  setActive(on) {
    this.active = on;
    if (on) { this.resize(); this.last = performance.now(); if (!this.raf) this.loop(); this.v3?.setActive(this.view === 'struct3d'); }
    else { if (this.raf) cancelAnimationFrame(this.raf); this.raf = null; this.v3?.setActive(false); }
  }
  setSimple(s) { this.simple = s; this.renderControls(); if (this.gene) { this.renderPanel(); this.renderDiversity(); } }

  // ---------- the slow-Pol II moment ----------
  slowVal() {
    const k = this.role.elong; if (!k) return null;
    const d = this.PARAMS[k];
    return String(d.unit).includes('nt/s') ? Math.max(d.min, 8.3) : Math.max(d.min, this.P0[k] / 3);
  }
  isSlow() { const k = this.role.elong; return !!k && this.P[k] <= (this.slowVal() ?? 0) * 1.05; }
  setSlow(on) {
    const k = this.role.elong; if (!k) return;
    this.setP(k, on ? this.slowVal() : this.P0[k]);
    this.syncSlowUI();
    this.renderControls(); this.recalc(); this.renderPace();
  }
  syncSlowUI() {
    const on = this.isSlow(), k = this.role.elong;
    document.body.classList.toggle('slowpol', on);
    $('slowBtn').textContent = on ? '⏩ Normal Pol II' : '⏳ Slow Pol II';
    if (k) $('slowPill').textContent = `⏳ SLOW Pol II: ${fmt(kbmin(this.P[k], this.PARAMS[k].unit), 2)} kb/min`;
  }

  setView(v) {
    this.view = v;
    try { localStorage.setItem('splice-view', v); } catch {}
    document.querySelectorAll('#spliceViewToggle button').forEach(b => b.classList.toggle('on', b.dataset.sv === v));
    const is3 = v === 'struct3d';
    document.body.classList.toggle('splice3d', is3);   // CSS picks the stage (see style.css)
    $('shint').textContent = is3
      ? 'true molecular scale · drag = orbit · wheel = zoom · ⏳ slow Pol II changes what the spliceosome has time to do · 1/2 = view'
      : 'gene diagram · move a slider, or hit ⏳ Slow Pol II · 1/2 = view · D = detail';
    if (is3) this.ensure3d(); else this.note();
    this.resize();
  }
  async ensure3d() {
    if (!this.v3) {
      try {
        const mod = await import('./view3d.js');
        this.v3 = new mod.Splice3DView($('splice3d'));
        this.v3.onReady = () => this.note();   // the PDB assets land a moment later
      } catch (e) { console.warn('splice: 3D view failed', e); $('sviewNote').innerHTML = '<b class="warn">3D view failed to start — see the console.</b>'; return; }
      if (this.gene) { this.v3.build(this.gene, this.J, this.P); this.v3.setSpace?.(this.sp); }
    }
    this.v3.setActive(this.active && this.view === 'struct3d');
    this.resize(); this.note();
  }
  note() {
    const n = $('sviewNote');
    n.innerHTML = this.view === 'struct3d' && this.v3 ? this.v3.noteHTML()
      : `<b>Real: the exon/intron structure and the splice-site scores${this.placeholder ? ' (none in this placeholder)' : ''}. Illustrative: the moving Pol II, the loops and the bins are an animation of the kinetic model.</b>
         Arc thickness above the gene = predicted usage of that junction. Introns are ${this.map?.compressed ? 'log-compressed' : 'drawn to scale'}; percentages are expected values from the ${this.expectedFn ? 'exact model' : 'fallback model'}, not the bead counts.`;
    n.style.display = '';
  }

  // ---------- isoform diversity: the combinatorial space, and the live identities ----------
  // The panel block is created from here rather than index.html so this feature is additive.
  mountDiversity() {
    const host = $('splicePanel'); if (!host || $('sdiv')) return;
    // headline counts + what is alive on the gene right now, directly under the fractions…
    const top = document.createElement('div');
    top.id = 'sdivtop';
    top.innerHTML = '<div class="dvtopbody"></div>';
    (($('sheadline') ?? null) ? host.insertBefore(top, $('sheadline')) : host.appendChild(top));
    // …and the fan-out tree itself lower down, next to the simulated bins it explains
    const d = document.createElement('div');
    d.id = 'sdiv';
    d.innerHTML = '<div class="ftitle">Where those isoforms come from — promoter × splice-site choices</div><div class="dvbody"><div class="shl-load">…</div></div>';
    (($('sbins') ?? null) ? host.insertBefore(d, $('sbins')) : host.appendChild(d));
    // "→ translate this isoform" on the tree leaves: same action as the rows above, delegated
    top.addEventListener('click', (e) => {
      const h = e.target.closest?.('.dvnode'); if (h) h.classList.toggle('open');
    });
    d.addEventListener('click', (e) => {
      const b = e.target.closest?.('button.sxl');
      if (b) {
        const why = this.onTranslate?.(this.gene?.gene, b.dataset.iso);
        this.xnote = why ?? null;
        this.renderPanel(); this.renderDiversity();
        return;
      }
      const h = e.target.closest?.('.dvnode');
      if (h) h.classList.toggle('open');
    });
    // the ② Translation join (protein per mRNA per isoform) and the measured Isopedia usage both
    // load in the background; each re-renders when it lands
    DIV.translation().then(() => { if (this.gene) { this.renderDiversity(); this.renderPanel(); } });
    DIV.usage().then(() => { if (this.gene) { this.renderDiversity(); this.renderPanel(); } });
  }
  // Full re-render of the fan-out tree: only on gene / parameter change (it is a big DOM).
  renderDiversity() {
    const host = $('sdiv'); if (!host || !this.gene) return;
    const body = host.querySelector('.dvbody'); if (!body) return;
    const top = $('sdivtop')?.querySelector('.dvtopbody');
    if (!this.sp || !this.exp) { body.innerHTML = '<div class="shl-load">…</div>'; return; }
    try {
      if (top) {
        // The measured validation comes FIRST: it is the only thing on this page that was tested
        // against data, and it does not flatter the model.
        top.innerHTML = DIV.validationHTML(this.gene.gene)
          + DIV.denominatorHTML(this.gene.gene)
          + `<div class="ftitle">How many real choices does this gene make?</div>`
          + DIV.countHTML(this.gene, this.J, this.sp, this.exp)
          + `<div id="sdivlive"></div>`;
      }
      body.innerHTML = `<div class="dvlegend">Every % here is an expected value from the ${this.expectedFn ? 'exact splicing model' : 'fallback model'}, never the bead counts:
             promoter branches use the model's own promoter weights, junction branches its <code>junctionUse</code>, leaves its isoform fractions.
             <b>No transcript is an average:</b> each polymerase on the gene picks one promoter and resolves each competing junction one way, so the chip above each nascent RNA is that transcript's isoform as far as it has been decided.</div>`
        + DIV.treeHTML(this.gene, this.J, this.sp, this.exp, { simple: this.simple, onTranslate: !!this.onTranslate });
    } catch (e) {
      console.warn('splice: diversity render failed', e);
      body.innerHTML = '<div class="sbad"><b>The isoform-space diagram failed to render for this gene</b> — see the console. The numbers above are unaffected.</div>';
    }
    this.renderDivLive();
  }
  // Cheap, per-frame part: what is live on the gene right now and the divergences seen.
  renderDivLive() {
    const el = $('sdivlive'); if (!el) return;
    const sigs = this.sigs ?? [];
    const rows = sigs.slice(0, 6).map(s => `<div class="dvlive" style="--idc:${s.colour}"><span class="dvchip">${esc(s.chip)}</span>
      <span class="mute">${s.resolved ? 'decided' : `${s.open} decision${s.open > 1 ? 's' : ''} still open`}</span></div>`).join('');
    const log = (this.dw.log ?? []).map(d => `<div class="dvdiv"><b>same site, different choice</b> at ${esc(this.siteName(d.at))}: <b>${esc(d.a)}</b> vs <b>${esc(d.b)}</b></div>`).join('');
    el.innerHTML = `<div class="dvliveh">${sigs.length} transcript${sigs.length === 1 ? '' : 's'} being written right now${this.dw.count ? ` · ${this.dw.count} divergence${this.dw.count === 1 ? '' : 's'} seen` : ''}</div>
      ${rows || '<div class="mute dvnone">waiting for the first polymerase to initiate…</div>'}${log}`;
  }
  // A competing group starts at one 5′ss when every junction in it shares that donor (a plain
  // 5′-site choice); otherwise it is a race between junctions with different donors, and saying
  // "the 5′ss" would be wrong.
  siteName(at) {
    const g = (this.sp?.groups ?? []).find(x => x.at === at);
    const same = g && g.members.every(i => (this.J[i].donorPos ?? this.J[i].start) === at);
    const nt = Number(at).toLocaleString();
    return same ? `the 5′ss at ${nt} nt` : `the junction race that opens at ${nt} nt`;
  }
  // one signature per live polymerase, plus the divergence watch (runs in both views)
  tickIdentity(pols, dt = 1 / 60) {
    if (!this.sp) { this.sigs = []; return; }
    try {
      this.sigs = pols.map(p => {
        const s = DIV.polSignature(this.gene, this.J, this.sp, p);
        s.polId = p.id; s.pos = p.pos; s.tssPos = p.tssPos; s.paused = p.paused;
        return s;
      });
      for (const d of this.dw.tick(this.sigs, this.t)) this.divPulse.push({ ...d, age: 0 });
    } catch (e) { console.warn('splice: identity tick failed', e); this.sigs = []; }
    for (const p of this.divPulse) p.age += dt;
    this.divPulse = this.divPulse.filter(p => p.age < 2.6).slice(-3);
  }

  // ---------- gene ----------
  setGene(g) {
    if (!g) return;
    this.gene = g; this.xnote = null;
    const i = this.genes.indexOf(g); if (i >= 0) $('sgene').value = i;
    // junctions: the model's list (it adds the skipping junctions implied by the isoforms)
    let J = null;
    if (this.junctionsFn) { try { J = this.junctionsFn(g); } catch (e) { console.warn('splice: junctions() failed', e); } }
    if (!J?.length) J = (g.introns ?? []).map((i2, idx) => ({ id: i2.id, key: `${i2.start}_${i2.end}`, idx, donorPos: i2.start, acceptorPos: i2.end, donor: i2.donor, acceptor: i2.acceptor, bp: i2.bp }));
    this.J = J.map(j => ({ ...j, start: j.donorPos ?? j.start, end: j.acceptorPos ?? j.end }));
    this.J.forEach((j, k) => { if (!Number.isFinite(j.idx)) j.idx = k; });
    this.bpOf = new Map((g.introns ?? []).filter(i2 => Number.isFinite(i2.bp?.pos)).map(i2 => [`${i2.start}_${i2.end}`, i2.bp.pos]));
    this.map = makeMap(g, this.J, VIEW_W, 58, 26);
    this.sim = this.makeSim();
    this.anim = []; this.drops = []; this.evAt = 0;
    this.sigs = []; this.divPulse = []; this.dw.reset();   // identity/divergence state is per gene
    this.buildSVG();
    this.v3?.build(g, this.J, this.P);
    this.recalc(); this.note(); this.syncSlowUI(); this.renderGeneNotes(); this.renderPace();
  }
  makeSim() {
    if (this.Gene) {
      try {
        const s = new this.Gene(this.gene, this.placeholder ? null : this.data);
        s.step(0.5, this.P);
        if (typeof s.view === 'function' && Array.isArray(s.pols)) return s;
      } catch (e) { console.warn('splice: Gene() failed, using the fallback sim', e); }
    }
    return new FbSim(this.gene, this.J);
  }
  // Never silently substitute a different gene: if this tab has no model for `sym`, keep what is
  // shown and say so (the caller can also act on the false return).
  matchGene(sym) {
    const g = this.genes?.find(x => x.gene === sym);
    this.noMatch = g ? null : sym;
    if (g && g !== this.gene) this.setGene(g); else this.renderGeneNotes();
    return !!g;
  }

  // ---------- numbers ----------
  recalc() {
    if (!this.gene) return;
    this.exp = this.expectedOf(this.P);
    const k = this.role.elong;
    this.cmpFast = k ? this.expectedOf({ ...this.P, [k]: this.P0[k] }) : this.exp;
    this.cmpSlow = k ? this.expectedOf({ ...this.P, [k]: this.slowVal() }) : this.exp;
    // the combinatorial space (cheap: ~1 ms/gene) is rebuilt whenever the numbers move, because the
    // branch labels and the near-deterministic flags are expected()'s own values
    try { this.sp = DIV.space(this.gene, this.J, this.tssListFn ? this.tssListFn(this.gene) : (this.gene.tss ?? []), this.exp); }
    catch (e) { console.warn('splice: diversity space failed', e); this.sp = null; }
    this.renderPanel();
    this.renderDiversity();
    this.v3?.setSpace?.(this.sp);
    this.v3?.setParams(this.P);
  }
  expectedOf(P) {
    if (this.expectedFn) {
      try {
        const r = this.expectedFn(this.gene, P);
        if (r && (r.isoforms || r.fractions)) return r;
      } catch (e) { console.warn('splice: expected() threw', e); }
    }
    try { return fbExpected(this.gene, this.J, P); } catch (e) { console.warn('splice: fallback expected failed', e); return { isoforms: [], psi: [], junctionUse: {}, fallback: true }; }
  }
  isoRows(r) {
    const rows = (r?.isoforms ?? []).map(i => ({ id: i.id, name: i.name ?? i.id, value: i.value }));
    const other = r?.fractions?.other ?? r?.other ?? 0;
    // NOT intron retention: the model has no retained-intron state. Every available non-conflicting
    // junction commits with probability 1, so this bucket is junction combinations that match no
    // annotated transcript.
    if (other > 0.0005) rows.push({ id: 'other', name: 'unannotated junction combinations (model artifact)', value: other, artifact: true });
    return rows.sort((a, b) => b.value - a.value);
  }

  // The model has no intron-retention output and no TSS weights. Both limits show up directly in the
  // numbers above, so they are stated next to them rather than buried in the sources list.
  modelLimitHTML(rows) {
    const g = this.gene;
    const zero = (rows ?? []).filter(r => r.id !== 'other' && !(r.value > 0.0005));
    const tss = (g?.tss ?? []).filter(t => Number.isFinite(t.pos ?? NaN));
    const weighted = tss.some(t => Number.isFinite(t.weight ?? t.frac ?? t.usage ?? NaN));
    const bits = [];
    bits.push(`<b>This model has no intron-retention output.</b> Every available non-conflicting junction commits with probability 1,
      so the leftover bucket is <b>unannotated junction combinations</b>, not retained introns.`);
    if (zero.length) bits.push(`${zero.length} annotated isoform${zero.length > 1 ? 's' : ''} sit${zero.length > 1 ? '' : 's'} at <b>0%</b>
      (${zero.slice(0, 6).map(z => esc(z.name)).join(', ')}${zero.length > 6 ? ', …' : ''}) because every junction they use is a strict subset of another isoform's set —
      <b>they may be real but are unreachable in this model</b>, not predicted to be absent.`);
    if (tss.length > 1 && !weighted) bits.push(`<b>TSS weights assumed equal:</b> ${esc(g.gene)} has ${tss.length} annotated TSSs and
      <code>splice_data.json</code> carries no usage weights, so they are mixed ${tss.length} × ${(100 / tss.length).toFixed(0)}% —
      which caps any single-TSS isoform at ${(100 / tss.length).toFixed(0)}% whatever the splicing does.`);
    const bad = CONTRADICTED[g?.gene];
    return `${bad ? `<div class="sbad"><b>${esc(g.gene)} is off the default demo path:</b> ${esc(bad)} It is shown only because the opt-in above is ticked.</div>` : ''}
      <div class="smodlim">${bits.join(' ')}</div>`;
  }

  // abundance × protein-per-mRNA: this isoform's share of the gene's protein output, over the
  // isoforms that HAVE a 5′ leader in utr_data.json (the others say so instead of being counted).
  trShare(rows, f) {
    if (!DIV.translationReady() || !this.gene) return null;
    let den = 0;
    for (const r of rows) {
      if (r.artifact) continue;
      const e = DIV.trEntry(this.gene.gene, r.name);
      if (e) den += r.value * e.perMin;
    }
    const e = DIV.trEntry(this.gene.gene, f.name);
    return e && den > 0 ? f.value * e.perMin / den : null;
  }

  // ---------- panel ----------
  renderPanel() {
    const g = this.gene;
    if (!g || !this.exp) { $('sheadline').innerHTML = '<div class="shl-load">no data yet</div>'; return; }
    const rows = this.isoRows(this.exp);
    $('sheadline').innerHTML = `<div class="ftitle">Predicted isoform fractions ${this.exp.fallback ? '<span class="warn">(fallback model)</span>' : ''}</div>
      ${rows.length ? rows.map(f => `<div class="sfr${f.artifact ? ' art' : ''}${f.value > 0.0005 ? '' : ' zero'}" title="${f.artifact ? 'junction combinations that match no annotated transcript — the model has no intron-retention state' : f.value > 0.0005 ? esc(f.name) : esc(f.name) + ' — 0%: every junction it uses is a strict subset of another isoform’s set, so this model can never reach it. It may still be a real transcript.'}"><span class="n">${esc(f.name)}</span><span class="b"><i style="width:${(100 * clamp(f.value, 0, 1)).toFixed(1)}%"></i></span><span class="v">${pct(f.value)}</span>${this.onTranslate && !f.artifact
        ? `<button type="button" class="sxl" data-iso="${esc(f.name)}" title="run the ribosome-scanning model on this isoform's 5′ leader (② Translation page)">→ translate this isoform</button>` : ''}${f.artifact ? '' : DIV.measChip(g.gene, f.id)}${f.artifact ? '' : DIV.trChip(g.gene, f.name, f.value, this.trShare(rows, f))}</div>`).join('')
        : '<div class="shl-load">this gene has no competing junctions in the data</div>'}
      ${this.xnote ? `<div class="sbad">${this.xnote}</div>` : ''}
      ${this.modelLimitHTML(rows)}
      ${this.keyMomentHTML()}`;
    const n = this.finished();
    // Below ~30 finished transcripts the sampling error is larger than anything it could show, so
    // the bars are withheld rather than captioned as if they had converged.
    const ENOUGH = 30;
    $('sbins').innerHTML = `<div class="ftitle">Simulated transcripts finished (n = ${n.total})</div>
      ${n.total >= ENOUGH
        ? rows.map(f => `<div class="sfr sim" title="${esc(f.name)}"><span class="n">${esc(f.name)}</span><span class="b"><i style="width:${(100 * n.frac(f.id)).toFixed(1)}%"></i></span><span class="v">${pct(n.frac(f.id))}</span></div>`).join('')
          + '<div class="fine">stochastic simulation of the same rules — it converges on the numbers above once a few hundred transcripts have finished</div>'
        : `<div class="fine">warming up — ${n.total} of ${ENOUGH} transcripts needed before the sampled fractions are worth showing (${pace.label(this.pace.ts)}; splicing a gene takes minutes).</div>`}`;
    const ev = (this.sim?.events ?? []).slice(-7).reverse();
    $('sevents').innerHTML = ev.length
      ? `<div class="ftitle">Last events</div><div class="sev">${ev.map(x => `<span>${esc(x.kind)}${x.intron ? ' ' + esc(x.intron) : ''}${x.name ? ' → ' + esc(x.name) : x.isoform ? ' → ' + esc(x.isoform) : ''}</span>`).join('')}</div>` : '';
  }
  keyMomentHTML() {
    const psi = this.exp.psi ?? [];
    const pick = psi.find(p => /inclusion/.test(p.event)) ?? psi[0];
    if (!pick) return '';
    const f = (this.cmpFast?.psi ?? []).find(p => p.event === pick.event);
    const s = (this.cmpSlow?.psi ?? []).find(p => p.event === pick.event);
    const k = this.role.elong, d = k ? this.PARAMS[k] : null;
    const meas = (this.gene.measured?.psi ?? []).find(p => p.event === pick.event)
      ?? (this.gene.measured?.psi ?? []).find(p => pick.exon && String(p.event).includes(pick.exon));
    const dd = f && s ? s.value - f.value : null;
    // The pair above is a FIXED pair of elongation rates, so it must not look like a live readout:
    // the live value (which the slider and the ⏳ button do move) is printed under it.
    const live = pct(pick.value);
    const liveKb = d ? fmt(kbmin(this.P[k], d.unit), 2) : null;
    return `<div class="kmoment"><div class="ftitle">Kinetic coupling — ${esc(pick.event)}</div>
      <div class="kbig"><b>${pct(f?.value)}</b><span>normal Pol II<br>${d ? fmt(kbmin(this.P0[k], d.unit), 1) + ' kb/min' : ''}</span><em>→</em><b class="slow">${pct(s?.value)}</b><span>slow Pol II<br>${d ? fmt(kbmin(this.slowVal(), d.unit), 2) + ' kb/min' : ''}</span></div>
      <div class="klive">At the <b>current</b> setting${liveKb ? ` (<b>${liveKb} kb/min</b>)` : ''}: <b>${live}</b>.
        <span class="mute">The two numbers above are a fixed ${d ? `${fmt(kbmin(this.P0[k], d.unit), 1)} vs ${fmt(kbmin(this.slowVal(), d.unit), 2)} kb/min` : 'normal vs slow'} comparison and do not follow the slider; this line does.</span></div>
      <div class="fine">${dd == null ? '' : Math.abs(dd) < 0.01
        ? 'this event barely moves — slowing Pol II does not change every exon (Fong 2014)'
        : `slowing Pol II <b>${dd > 0 ? 'raises' : 'lowers'}</b> it by ${(100 * Math.abs(dd)).toFixed(0)} points (de la Mata 2003)`}</div>
      ${meas ? `<div class="meas">Measured ${fmt(meas.value, 2)}${meas.tissue ? ` · ${esc(meas.tissue)}` : ''}${meas.n ? ` · n=${meas.n}` : ''}${meas.cite ? ` · ${esc(meas.cite)}` : ''}</div>`
        : '<div class="honest">no measured PSI for this event in the data</div>'}</div>`;
  }
  finished() {
    const c = this.sim?.counts ?? {};
    const total = Object.values(c).reduce((a, b) => a + b, 0);
    return { total, get: (k) => (c[k] ?? 0), frac: (k) => (total ? (c[k] ?? 0) / total : 0) };
  }
  renderControls() {
    // ANIM_ONLY sliders are kept out of the headline: the exact model never reads them, so moving
    // one there would look like it changes the prediction when it only changes the animation.
    const head = [this.role.elong, this.role.pause, this.role.strength, this.role.cotxn]
      .filter(k => k && !ANIM_ONLY.has(k));
    const row = (k) => {
      const d = this.PARAMS[k]; if (!d) return '';
      const extra = k === this.role.elong ? ` <small>= ${fmt(kbmin(this.P[k], d.unit), 2)} kb/min</small>` : '';
      const anim = ANIM_ONLY.has(k) ? ' <small class="warn">(animation only — does not affect the predicted fractions)</small>' : '';
      return `<label class="sctl${anim ? ' animonly' : ''}"><span class="l">${esc(d.label)}${d.cite ? '' : ' <small>(assumed)</small>'}${extra}${anim}</span>
        <input type="range" data-pk="${k}" min="${d.min}" max="${d.max}" step="${d.step}" value="${this.P[k]}">
        <span class="v" data-pv="${k}">${this.P[k]}${d.unit ? ' ' + d.unit : ''}</span></label>`;
    };
    $('sctrls').innerHTML = head.map(row).join('');
    const ar = this.animRange;
    $('sctrlsAll').innerHTML = `<label class="sctl"><span class="l">animation speed × the watchable pace <small>(the absolute compression is worked out per gene and shown above the gene diagram)</small></span><input type="range" data-pk="_speed" min="${ar.min}" max="${ar.max}" step="${ar.step}" value="${this.mul}"><span class="v" data-pv="_speed">${this.mul}×</span></label>`
      + Object.keys(this.PARAMS).filter(k => !head.includes(k)).map(row).join('');
    for (const inp of document.querySelectorAll('#splicePanel [data-pk]')) {
      const k = inp.dataset.pk;
      inp.oninput = () => {
        const v = +inp.value;
        if (k === '_speed') this.mul = v; else this.setP(k, v);
        const out = document.querySelector(`[data-pv="${k}"]`);
        if (out) out.textContent = k === '_speed' ? v + '×' : v + (this.PARAMS[k]?.unit ? ' ' + this.PARAMS[k].unit : '');
        if (k !== '_speed') { this.syncSlowUI(); this.recalc(); }
        this.renderPace();
      };
    }
  }
  setP(k, v) { this.P[k] = v; }
  renderSources() {
    $('ssources').innerHTML = `<div class="ftitle">Sources</div><ul>${Object.entries(this.SOURCES).map(([k, s]) =>
      `<li><b>${esc(k)}</b>: ${s.url ? `<a href="${esc(s.url)}" target="_blank" rel="noopener">${esc(s.text)}</a>` : esc(s.text)}</li>`).join('')}</ul>
      <div class="fine">Assumed (no citation yet): ${Object.values(this.PARAMS).filter(d => !d.cite).map(d => esc(d.label)).join(' · ')}</div>`;
  }

  // ---------- SVG scaffold ----------
  buildSVG() {
    const g = this.gene, m = this.map;
    this.svg.setAttribute('viewBox', `0 0 ${VIEW_W} ${VIEW_H}`);
    this.svg.innerHTML = '';
    const defs = el('defs', {}, this.svg);
    const mk = el('marker', { id: 'arr', viewBox: '0 0 10 10', refX: 8, refY: 5, markerWidth: 6, markerHeight: 6, orient: 'auto' }, defs);
    el('path', { d: 'M0 0 L10 5 L0 10 z', fill: COL.exon }, mk);
    // the ② Translation page's bright-dot look, done cheaply in SVG: one blur+merge filter, applied
    // ONLY to moving actors (Pol II, committed spliceosome, lariat, finished transcripts) — never to
    // text, exon boxes or splice-site scores, which the layout audit found crowded already.
    const fl = el('filter', { id: 'sglow', x: '-120%', y: '-120%', width: '340%', height: '340%' }, defs);
    el('feGaussianBlur', { in: 'SourceGraphic', stdDeviation: 3.2, result: 'b1' }, fl);
    el('feGaussianBlur', { in: 'SourceGraphic', stdDeviation: 7, result: 'b2' }, fl);
    const fm = el('feMerge', {}, fl);
    el('feMergeNode', { in: 'b2' }, fm);
    el('feMergeNode', { in: 'b1' }, fm);
    el('feMergeNode', { in: 'SourceGraphic' }, fm);
    const Y = this.Y = { dna: 150, rna: 290, bin: 414 };
    const root = el('g', {}, this.svg);
    const alt = (t) => t && !['constitutive', 'first', 'last'].includes(t);
    // a label is drawn only when it will not collide with the last one on its row
    const rows = {};
    const room = (row, x, w) => {
      const last = rows[row] ?? -1e9;
      if (x - w / 2 < last) return false;
      rows[row] = x + w / 2 + 4; return true;
    };
    // promoter + TSS (several TSSs may sit within a few px of each other)
    const tss = (g.tss ?? [{ pos: m.s0 }]).slice().sort((a, b) => (a.pos ?? 0) - (b.pos ?? 0));
    for (const t of tss) {
      const xT = m(t.pos ?? m.s0);
      el('path', { d: `M${xT} ${Y.dna - 6} L${xT} ${Y.dna - 30} L${xT + 24} ${Y.dna - 30}`, fill: 'none', stroke: COL.exon, 'stroke-width': 2.2, 'marker-end': 'url(#arr)' }, root);
      if (!room('tss', xT + 18, 58)) continue;
      const lab = el('text', { x: xT + 2, y: Y.dna - 36, class: 'slab on' }, root);
      lab.textContent = (tss.length > 1 && t.id ? 'TSS ' + t.id : 'TSS');
      if (t.support) el('title', {}, lab).textContent = `${t.support}${t.cite ? ' · ' + t.cite : ''}`;
    }
    const xP = m(m.s0);
    el('line', { x1: Math.max(4, xP - 48), y1: Y.dna, x2: xP, y2: Y.dna, stroke: COL.intron, 'stroke-width': 5 }, root);
    el('text', { x: Math.max(4, xP - 48), y: Y.dna - 12, class: 'slab' }, root).textContent = 'promoter';
    // intron lines
    for (const j of (g.introns ?? [])) {
      const a = m(j.start), b = m(j.end);
      el('line', { x1: a, y1: Y.dna, x2: b, y2: Y.dna, stroke: COL.intron, 'stroke-width': 1.6 }, root);
      const n = Math.max(1, Math.floor((b - a) / 26));
      for (let k = 1; k <= n; k++) {
        const x = a + (b - a) * k / (n + 1);
        el('path', { d: `M${x - 3} ${Y.dna - 4} L${x + 2} ${Y.dna} L${x - 3} ${Y.dna + 4}`, fill: 'none', stroke: COL.intron, 'stroke-width': 1.2 }, root);
      }
    }
    // junction arcs above, thickness = predicted usage (set in frame())
    this.arcs = new Map();
    const multi = this.J.filter(j => this.J.some(o => o !== j && (o.start === j.start || o.end === j.end)));
    for (const j of this.J) {
      const a = m(j.start), b = m(j.end), h = clamp((b - a) * 0.22, 16, 54) + (multi.includes(j) ? 10 : 0);
      const p = el('path', { d: `M${a} ${Y.dna - 15} Q${(a + b) / 2} ${Y.dna - 15 - 2 * h} ${b} ${Y.dna - 15}`, fill: 'none', stroke: COL.rna, 'stroke-width': 1.4, opacity: 0.85 }, root);
      const t = el('text', { x: (a + b) / 2, y: Y.dna - 18 - h, class: 'slab', 'text-anchor': 'middle', fill: COL.rna }, root);
      el('title', {}, p).textContent = `${j.id}: donor ${fmt(j.donor?.score, 1)} / acceptor ${fmt(j.acceptor?.score, 1)}${j.donor?.scale ? ' (' + j.donor.scale + ')' : ''}`;
      this.arcs.set(j.id, { path: p, lab: t });
    }
    // exons
    for (const e of g.exons ?? []) {
      const a = m(e.start), b = m(e.end), w = Math.max(6, b - a);
      const r = el('rect', { x: a, y: Y.dna - 13, width: w, height: 26, rx: 3, fill: alt(e.type) ? COL.alt : COL.exon, stroke: '#0b0d12', 'stroke-width': 1 }, root);
      el('title', {}, r).textContent = `${e.id} · ${e.type ?? 'exon'} · ${e.end - e.start} nt`;
      if (w >= 20) el('text', { x: a + w / 2, y: Y.dna + 5, class: 'exlab', 'text-anchor': 'middle' }, root).textContent = e.id;
      if (alt(e.type) && room('exalt', a + w / 2, 8 * String(e.type).length))
        el('text', { x: a + w / 2, y: Y.dna + 68, class: 'slab alt', 'text-anchor': 'middle' }, root).textContent = e.type;
    }
    // splice sites with their strength scores
    const seen = new Set();
    for (const j of this.J) {
      for (const [side, s, pos] of [['donor', j.donor, j.start], ['acceptor', j.acceptor, j.end]]) {
        const key = side + pos;
        if (seen.has(key)) continue; seen.add(key);
        const x = m(pos), up = side === 'donor';
        const ln = el('line', { x1: x, y1: Y.dna + 13, x2: x, y2: Y.dna + 26, stroke: up ? COL.don : COL.acc, 'stroke-width': 1.6 }, root);
        const txt = `${up ? '5′ss' : '3′ss'} ${s && Number.isFinite(s.score) ? s.score.toFixed(1) : '—'}`;
        const tip = s && Number.isFinite(s.score)
          ? `${side} ${s.seq ?? ''} · ${s.score.toFixed(2)} ${s.scale ?? ''}`.trim()
          : `${side}: no score in the data`;
        el('title', {}, ln).textContent = tip;
        const row = up ? 'ssD' : 'ssA';
        if (!room(row, x, 7 * txt.length)) continue;
        const t = el('text', { x, y: Y.dna + (up ? 38 : 50), class: 'sslab', 'text-anchor': 'middle', fill: up ? COL.don : COL.acc }, root);
        t.textContent = txt;
        el('title', {}, t).textContent = tip;
      }
      const bp = this.bpOf.get(`${j.start}_${j.end}`);
      if (Number.isFinite(bp)) {
        const x = m(bp);
        el('circle', { cx: x, cy: Y.dna + 20, r: 3, fill: COL.acc }, root);
        el('text', { x, y: Y.dna + 32, class: 'sslab', 'text-anchor': 'middle', fill: COL.acc }, root).textContent = 'BP';
      }
    }
    // live layers
    this.gRNA = el('g', {}, root); this.gPol = el('g', {}, root); this.gDrop = el('g', {}, root);
    // bins
    this.gBin = el('g', {}, root);
    this.bins = new Map();
    const list = (g.isoforms ?? []).map(i => ({ id: i.id, name: i.name ?? i.id, iso: i }));
    // not "retained": the model has no intron-retention state (see modelLimitHTML)
    list.push({ id: 'other', name: 'unannotated combos', iso: null });
    const bw = Math.min(200, (VIEW_W - 100) / list.length - 14);
    list.forEach((o, k) => {
      const x = 58 + k * (bw + 14), gb = el('g', {}, this.gBin);
      el('rect', { x, y: Y.bin, width: bw, height: 66, rx: 7, fill: 'rgba(233,228,216,.05)', stroke: 'rgba(233,228,216,.18)' }, gb);
      const fill = el('rect', { x: x + 1, y: Y.bin + 65, width: bw - 2, height: 0, rx: 6, fill: o.iso ? COL.rna : COL.intron, opacity: 0.32 }, gb);
      el('text', { x: x + 8, y: Y.bin + 17, class: 'binlab' }, gb).textContent = String(o.name).slice(0, 24);
      const v = el('text', { x: x + 8, y: Y.bin + 48, class: 'binpct' }, gb);
      const c = el('text', { x: x + bw - 8, y: Y.bin + 48, class: 'bincnt', 'text-anchor': 'end' }, gb);
      if (o.iso) {
        const set = new Set(o.iso.exon_ids ?? []);
        const kept = (g.exons ?? []).filter(e => set.has(e.id));
        const tot = kept.reduce((a, e) => a + (e.end - e.start), 0) || 1;
        let gx = x + 8;
        el('line', { x1: x + 8, y1: Y.bin + 29, x2: x + bw - 8, y2: Y.bin + 29, stroke: COL.intron, 'stroke-width': 1 }, gb);
        for (const e of kept) {
          const w = Math.max(3, (e.end - e.start) / tot * (bw - 16));
          el('rect', { x: gx, y: Y.bin + 24, width: w, height: 10, rx: 2, fill: alt(e.type) ? COL.alt : COL.exon }, gb);
          gx += w + 1;
        }
      }
      this.bins.set(o.id, { x, bw, fill, v, c, cy: Y.bin });
    });
    el('text', { x: 58, y: Y.bin - 12, class: 'slab' }, root).textContent = 'finished transcripts drop into bins →';
    el('text', { x: VIEW_W - 26, y: VIEW_H - 6, class: 'slab', 'text-anchor': 'end' }, root).textContent =
      `${esc(g.gene)} · ${(m.s1 - m.s0).toLocaleString()} nt${m.compressed ? ' · introns log-compressed' : ''}`;
  }

  // ---------- per-frame ----------
  pump(dt) {
    this.pace = this.solvePace();
    const want = dt * this.pace.ts;
    // more sub-steps than before: a 1000x+ compression on a 478 kb gene must not coarsen the race
    const n = clamp(Math.ceil(want / SUBSTEP), 1, 160);
    for (let i = 0; i < n; i++) { try { this.sim.step(want / n, this.P); } catch (e) { console.warn('splice: step failed', e); this.sim = new FbSim(this.gene, this.J); break; } }
    // consume events for the one-shot animations
    const ev = this.sim.events ?? [];
    if (this.evAt > ev.length) this.evAt = 0;
    for (; this.evAt < ev.length; this.evAt++) {
      const e = ev[this.evAt];
      if (e.kind === 'spliced') this.anim.push({ kind: 'lariat', pol: e.pol, intron: e.intron, age: 0 });
      // the finished transcript falls in the colour of the isoform it became (same key the tree
      // leaves use), so the drop, the leaf and the bin all agree
      else if (e.kind === 'finished') this.drops.push({ key: e.isoform ?? 'other', age: 0, colour: e.isoform && e.isoform !== 'other' ? DIV.idColour(e.isoform) : COL.intron });
    }
    for (const a of this.anim) a.age += dt;
    this.anim = this.anim.filter(a => a.age < 1.5).slice(-12);
    for (const d of this.drops) d.age += dt;
    this.drops = this.drops.filter(d => d.age < 1.1).slice(-20);
    // per-transcript identity + the "same site, different choice" watch: done here, not in frame(),
    // so it also runs while the 3D view is on screen
    this.tickIdentity((typeof this.sim.view === 'function' ? this.sim.view() : this.sim.pols) ?? [], dt);
  }
  frame(dt) {
    const g = this.gene, m = this.map, Y = this.Y;
    if (!g || !m) return;
    const pols = (typeof this.sim.view === 'function' ? this.sim.view() : this.sim.pols) ?? [];
    this.gPol.innerHTML = ''; this.gRNA.innerHTML = '';
    const byId = new Map(this.J.map(j => [j.id, j]));
    // each nascent transcript's own emerging isoform (promoter + the choices it has made so far)
    const sigBy = new Map((this.sigs ?? []).map(s => [s.polId, s]));
    this.chipRow = {};
    for (const p of pols.slice(0, 8)) {
      const pos = p.pos ?? 0, x = m(pos);
      const rem = new Set(p.removed ?? []), com = new Set(p.committed ?? []), dead = new Set(p.skipped ?? []);
      // nascent RNA, 5' end anchored at the TSS
      const wig = (x0, x1, amp) => {
        let s = '';
        const n = clamp(Math.round((x1 - x0) / 16), 1, 80);
        for (let k = 1; k <= n; k++) {
          const xx = x0 + (x1 - x0) * k / n;
          s += ` Q${xx - (x1 - x0) / n / 2} ${Y.rna + Math.sin(this.t * 2.6 + xx * 0.2 + p.id) * amp} ${xx} ${Y.rna}`;
        }
        return s;
      };
      let d = `M${m(p.tssPos ?? m.s0)} ${Y.rna}`, cur = p.tssPos ?? m.s0;
      const live = this.J.filter(j => j.start < pos && !dead.has(j.id) && j.start >= (p.tssPos ?? 0)).sort(ntSort);
      for (const j of live) {
        if (j.start < cur) continue;
        const a = m(j.start), b = m(Math.min(j.end, pos));
        d += wig(m(cur), a, 2.4);
        if (rem.has(j.id)) { d += ` L${b} ${Y.rna}`; }
        else {
          const h = clamp(26 + (b - a) * 0.45, 24, 92);
          d += ` C${a + (b - a) * 0.08} ${Y.rna + h}, ${b - (b - a) * 0.08} ${Y.rna + h}, ${b} ${Y.rna}`;
          if (com.has(j.id)) { // spliceosome committed to this intron: amber bead on the loop
            const c = el('circle', { cx: (a + b) / 2, cy: Y.rna + h * 0.78, r: 6, fill: COL.alt, opacity: .9, class: 'glow' }, this.gRNA);
            el('title', {}, c).textContent = `${j.id}: spliceosome committed, waiting to cut`;
            if (b - a > 90) el('text', { x: (a + b) / 2, y: Y.rna + h * 0.78 + 17, class: 'slab alt', 'text-anchor': 'middle' }, this.gRNA).textContent = 'committed';
          }
        }
        cur = Math.min(j.end, pos);
      }
      d += wig(m(cur), x, 3.2);
      // the RNA is drawn in the colour of the isoform signature this transcript has committed to so
      // far, so two polymerases that have made different choices are different colours on screen
      const sig = sigBy.get(p.id) ?? null;
      el('path', { d, fill: 'none', stroke: sig?.colour ?? COL.rna, 'stroke-width': 2.2, 'stroke-linecap': 'round', opacity: sig && !sig.resolved ? .8 : .95 }, this.gRNA);
      // lariats released by this polymerase
      for (const a of this.anim.slice(-2)) {
        if (a.pol !== p.id) continue;
        const j = byId.get(a.intron); if (!j) continue;
        const xa = m(j.start), xb = m(j.end), cx = (xa + xb) / 2, k = a.age / 1.5;
        const r = clamp((xb - xa) / 5, 7, 24) * (1 - .45 * k);
        el('circle', { cx, cy: Y.rna + 64 + k * 74, r, fill: 'none', stroke: COL.intron, 'stroke-width': 2, opacity: (1 - k).toFixed(2), class: 'glow' }, this.gRNA);
        el('text', { x: cx, y: Y.rna + 64 + k * 74 - r - 5, class: 'slab', 'text-anchor': 'middle', opacity: (1 - k).toFixed(2) }, this.gRNA).textContent = 'intron lariat';
      }
      // Pol II
      const gp = el('g', {}, this.gPol);
      el('line', { x1: x, y1: Y.dna + 3, x2: x, y2: Y.rna, stroke: COL.pol, 'stroke-width': 1, opacity: .3, 'stroke-dasharray': '3 3' }, gp);
      el('circle', { cx: x, cy: Y.dna, r: 11, fill: COL.pol, opacity: p.paused ? .5 : 1, class: p.paused ? '' : 'glow' }, gp);
      el('text', { x, y: Y.dna + 4, class: 'pollab', 'text-anchor': 'middle' }, gp).textContent = 'Pol';
      if (p.paused) el('text', { x, y: Y.dna - 50, class: 'slab warn2', 'text-anchor': 'middle' }, gp).textContent = 'paused';
      if (p.postTx) el('text', { x, y: Y.rna + 4, class: 'slab warn2', 'text-anchor': 'start' }, gp).textContent = ' post-txn';
      // the identity chip: this transcript's isoform as far as it has been decided ("?" = still racing)
      if (sig) this.drawChip(x, sig);
    }
    this.drawDivergence(m, Y);
    // falling transcripts
    this.gDrop.innerHTML = '';
    for (const dr of this.drops) {
      const b = this.bins.get(dr.key) ?? this.bins.get('other'); if (!b) continue;
      const k = clamp(dr.age / 1.0, 0, 1), xe = m(m.s1);
      const x = xe + (b.x + b.bw / 2 - xe) * k, y = Y.rna + (b.cy + 30 - Y.rna) * k * k;
      el('rect', { x: x - 9, y: y - 3, width: 18, height: 6, rx: 3, fill: dr.colour ?? COL.rna, opacity: (1 - .3 * k).toFixed(2), class: 'glow' }, this.gDrop);
    }
    // bins + arcs from the expected model
    const n = this.finished();
    const rows = new Map(this.isoRows(this.exp).map(r => [r.id, r.value]));
    for (const [id, b] of this.bins) {
      const f = n.total ? n.frac(id) : 0, h = f * 64;
      b.fill.setAttribute('height', h.toFixed(1)); b.fill.setAttribute('y', (b.cy + 65 - h).toFixed(1));
      b.v.textContent = rows.has(id) ? pct(rows.get(id)) : (n.total ? pct(f) : '—');
      b.c.textContent = `${n.get(id)} made`;
    }
    const use = this.exp?.junctionUse ?? {};
    for (const [id, a] of this.arcs) {
      const u = use[id];
      a.path.setAttribute('stroke-width', (0.6 + 6 * clamp(u ?? 0, 0, 1)).toFixed(2));
      a.path.setAttribute('opacity', (0.25 + 0.7 * clamp(u ?? 0, 0, 1)).toFixed(2));
      a.lab.textContent = Number.isFinite(u) ? pct(u) : '';
    }
  }
  // ---------- live isoform identity on the gene diagram ----------
  // One chip per nascent transcript: "T2 · −E3 · ?" — the promoter it started from, then one entry
  // per competing junction group, filled in as that race is decided. Two chips that differ are two
  // transcripts of the same gene becoming different isoforms.
  drawChip(x, sig) {
    const Y = this.Y;
    const txt = sig.chip;
    const w = Math.max(34, 6.1 * txt.length + 12);
    // stagger over three rows so chips do not overlap when several polymerases are close together
    let row = 0;
    for (; row < 3; row++) { const last = this.chipRow[row] ?? -1e9; if (x - w / 2 > last) break; }
    if (row > 2) return;
    this.chipRow[row] = x + w / 2 + 6;
    const y = Y.rna - 18 - row * 17;
    const g = el('g', {}, this.gRNA);
    el('rect', {
      x: (x - w / 2).toFixed(1), y: y - 11, width: w.toFixed(1), height: 15, rx: 7,
      fill: 'rgba(11,13,18,.82)', stroke: sig.colour, 'stroke-width': sig.resolved ? 1.6 : 1,
      'stroke-dasharray': sig.resolved ? '' : '3 2',
    }, g);
    const t = el('text', { x, y: y, class: 'dvchiplab', 'text-anchor': 'middle', fill: sig.colour }, g);
    t.textContent = txt;
    el('title', {}, t).textContent = sig.resolved
      ? 'this transcript has made every choice it can: this is the isoform it will be'
      : '“?” = this transcript has not yet resolved that competing junction — the race is still running on this molecule';
    el('line', { x1: x, y1: y + 4, x2: x, y2: Y.rna - 3, stroke: sig.colour, 'stroke-width': 1, opacity: .45 }, g);
  }
  // A divergence moment: two live transcripts resolved the SAME junction group differently.
  drawDivergence(m, Y) {
    for (const d of this.divPulse) {
      const k = clamp(d.age / 2.6, 0, 1), x = m(d.at);
      const g = el('g', { opacity: (1 - k * k).toFixed(2) }, this.gRNA);
      el('circle', { cx: x, cy: Y.dna + 13, r: (8 + 26 * k).toFixed(1), fill: 'none', stroke: COL.alt, 'stroke-width': (2.4 * (1 - k)).toFixed(2) }, g);
      el('circle', { cx: x, cy: Y.dna + 13, r: 5, fill: COL.alt, class: 'glow' }, g);
      const y = Y.dna + 92;
      el('line', { x1: x, y1: Y.dna + 18, x2: x, y2: y - 13, stroke: COL.alt, 'stroke-width': 1.2, 'stroke-dasharray': '3 2' }, g);
      const txt = `same site, different choice: ${d.a} vs ${d.b}`;
      const w = 6.4 * txt.length + 16;
      el('rect', { x: (x - w / 2).toFixed(1), y: y - 13, width: w.toFixed(1), height: 18, rx: 9, fill: 'rgba(255,181,71,.14)', stroke: COL.alt, 'stroke-width': 1.2 }, g);
      el('text', { x, y: y, class: 'dvpulse', 'text-anchor': 'middle', fill: COL.alt }, g).textContent = txt;
    }
  }
  loop() {
    this.raf = requestAnimationFrame((now) => {
      const dt = Math.min((now - this.last) / 1000, 0.05); this.last = now;
      this.t += dt;
      if (this.active && this.sim) {
        this.pump(dt);
        if (this.view === 'struct3d') this.v3?.update(dt, this.sim, this.J);
        else this.frame(dt);
        this.uiAcc = (this.uiAcc ?? 0) + dt;
        if (this.uiAcc > 0.4) { this.uiAcc = 0; this.renderPanel(); this.renderPace(); this.renderDivLive(); }
      }
      this.loop();
    });
  }
  resize() { this.v3?.resize(); }

  // ---------- animation pace (view only: no printed number depends on it) ----------
  // Pol II really takes 4.6 min to cross PPP1R15A and 5.3 h to cross ESR1; at one fixed
  // "× real time" one of those is a blur and the other is frozen. pace.js maps the real rate
  // onto a watchable tempo per gene, keeping (a compressed version of) the slow<->fast ratio.
  solvePace(mul = this.mul, rate = null) {
    const m = this.map, ek = this.role.elong;
    const nt = m ? m.s1 - m.s0 : 1;
    const r = rate ?? (ek ? this.P[ek] : 25);
    // the two things that can leave the gene looking dead: no new polymerase arrives, and a
    // promoter-proximal pause nobody escapes
    const init = this.P.initRate ?? 0.02, pe = this.P.pauseEscape ?? 0.002;
    const gapSim = Math.max(init > 0 ? 1 / init : 0, (this.P.pauseProb ?? 0) > 0.02 && pe > 0 ? 1 / pe : 0);
    // the x axis is log-compressed, so the bead dashes through the exons and crawls through a long
    // intron: the blur ceiling has to be set by the FASTEST stretch, not by the average
    return pace.solve({ nt, rate: r, mul, refNt: 20000, refRate: (ek ? this.P0?.[ek] : null) ?? r, gapSim, fastRate: r * burstOf(m), blurSoft: true });
  }
  renderPace() {
    let n = $('spaceNote');
    if (!n) {
      n = document.createElement('div');
      n.id = 'spaceNote'; n.className = 'pacenote';
      const host = $('sctrls') ?? $('sbins');
      if (!host) return;                       // panel not mounted yet; try again next frame
      host.before(n);
    }
    const p = this.pace = this.solvePace();
    const m = this.map;
    if (!m || !(p.ts > 0)) { n.textContent = ''; return; }
    const ek = this.role.elong, d = ek ? this.PARAMS[ek] : null;
    const ratio = d ? pace.ratioNote('Pol II elongation', this.solvePace(this.mul, d.min), this.solvePace(this.mul, d.max)) : '';
    n.innerHTML = `<b>⏱ ${esc(pace.label(p.ts))}</b>
      <span>Pol II crosses ${esc(this.gene?.gene ?? 'this gene')} (${(m.s1 - m.s0).toLocaleString()} nt) in
      <b>${fmt(p.crossWall)} s</b> of your time — ${esc(dur(p.crossSim))} in the cell at ${fmt(ek ? this.P[ek] : 25, 0)} nt/s.</span>
      <span>${esc(ratio)}${p.limit ? ` · here the picture is held at the <b>${p.limit}</b> watchable pace, so the on-screen difference is smaller than the real one` : ''}</span>
      <span class="fine">Only the wall-clock tempo is compressed. Every percentage on this page is the exact model at the real rates.</span>`;
  }

  // ---------- pacing diagnostics (read-only) ----------
  // Used by scratchpad/pace/measure.mjs to MEASURE the Pol II bead's px/s and traversal time.
  // Reports positions in CSS pixels; changes nothing.
  paceInfo() {
    const m = this.map, g = this.gene;
    const o = { gene: g?.gene ?? null, nt: m ? m.s1 - m.s0 : 0, timeScale: this.pace?.ts, mul: this.mul, anim: this.animRange, view: this.view };
    const ctm = this.svg?.getScreenCTM?.();
    if (m && ctm) {
      const px = (x) => ctm.a * x + ctm.e;
      o.lane = { x0: px(m(m.s0)), x1: px(m(m.s1)) };
      const pols = (typeof this.sim?.view === 'function' ? this.sim.view() : this.sim?.pols) ?? [];
      o.marks = pols.map(p => ({ id: p.id, px: px(m(p.pos ?? 0)), nt: p.pos ?? 0, paused: !!p.paused }));
    }
    o.pace = this.pace ?? null;
    return o;
  }
}

// ---------------------------------------------------------------- placeholder (clearly labelled)
// Geometry only, so the tab is demoable before pipeline/splice writes splice_data.json.
// No scores and no measured values: nothing here is a number we would show as data.
const PLACEHOLDER = {
  meta: { release: 'placeholder (not real data)', sources: [], methods: {} },
  genes: [{
    gene: 'PLACEHOLDER (cassette exon)', strand: '+',
    tss: [{ id: 'T1', pos: 0, support: 'placeholder', cite: null }],
    exons: [
      { id: 'E1', start: 0, end: 210, type: 'first' },
      { id: 'E2', start: 2100, end: 2370, type: 'cassette' },
      { id: 'E3', start: 5200, end: 5900, type: 'last' },
    ],
    introns: [
      { id: 'I1', start: 210, end: 2100, donor: { pos: 210, score: null }, acceptor: { pos: 2100, score: null }, bp: { pos: null } },
      { id: 'I2', start: 2370, end: 5200, donor: { pos: 2370, score: null }, acceptor: { pos: 5200, score: null }, bp: { pos: null } },
    ],
    isoforms: [
      { id: 'inc', name: 'E2 included', tss_id: 'T1', exon_ids: ['E1', 'E2', 'E3'], utr5_id: null },
      { id: 'skp', name: 'E2 skipped', tss_id: 'T1', exon_ids: ['E1', 'E3'], utr5_id: null },
    ],
    measured: { psi: [], isoform_fraction: [] },
  }],
};
