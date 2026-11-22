// "Target: therapeutic levers" — NOT a page of its own. The same module is mounted twice, once per
// page, and each instance shows only the levers that act on that page's biology:
//   ② Translation  (eIF4A / ISR / mTOR / start-codon levers — parameters from src/params.js)
//   ① Transcription (splice-site, SF3B1, CLK-SRPK and elongation-rate levers — src/splice/params.js)
// The split is computed from each lever's `param` (see leverPage); a lever may also state
// `"page": "transcription" | "translation"` in app/public/targets.json to override it.
//
// Self-contained module (own DOM, own CSS, no edits to style.css / index.html):
//   mountTarget(hooks) -> { mount(rootEl, ctx), update(ctx), open(), close(), el }
//
// What it does: a lever maps a *model parameter* onto a real drug class. The dose slider moves
// that parameter from its baseline (src/params.js defaults) to its extreme, and the predicted
// output ratio (treated / untreated) is recomputed LIVE for every gene with the same exact
// scanning model the rest of the app uses (src/analytic.js `expected`). Nothing is precomputed.
//
// The readout is a THERAPEUTIC WINDOW: how much more the oncogene set loses than the reference set.
// It is a model prediction, not a measurement — see the honesty block, which is always visible.
import { PARAMS, values } from '../params.js';
import { expected } from '../analytic.js';
import { motifStyle, chipText, motifsOf } from '../motifs.js';
import './panel.css';

// optional splicing model: present only if SPLICE-MODEL has landed it. import.meta.glob resolves to
// {} when the file does not exist, so this build never breaks on its absence.
const SPLICE_GLOB = import.meta.glob('../splice/model_exact.js');
let spliceMod = null, spliceTried = false;
async function loadSplice() {
  if (spliceTried) return spliceMod;
  spliceTried = true;
  const k = Object.keys(SPLICE_GLOB)[0];
  if (!k) return null;
  try { spliceMod = await SPLICE_GLOB[k](); } catch { spliceMod = null; }
  return spliceMod;
}

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const clamp = (v, a, b) => Math.max(a, Math.min(b, v));
const f2 = (x, d = 2) => (Number.isFinite(x) ? Number(x).toFixed(d) : '—');
// "8.3" style fold change
const fmtX = (x) => (!Number.isFinite(x) ? '—' : x >= 100 ? x.toFixed(0) : x >= 10 ? x.toFixed(1) : x.toFixed(2));
const RMIN = 1e-4, RMAX = 1e4;

// ---------- honesty: the model's dynamic range is not a fold change ----------
// Predicted output spans 4.3e6-fold across these leaders, so a printed response ratio is clamped
// to >100x / <0.01x and means rank order only.
const CAP_HI = 100, CAP_LO = 1 / CAP_HI;
const RANK_TIP = "the model's dynamic range is unphysically wide; treat these as rank order, not fold-change";
// "↑ 8.3×" / "↓ 4.1×", or "↑ >100×" / "↓ >100×" past the cap
const fmtArrow = (r) => (!Number.isFinite(r) ? '—'
  : Math.abs(Math.log10(r)) < 0.01 ? '= no change'
    : r < 1 ? (r < CAP_LO ? '↓ &gt;100×' : `↓ ${fmtX(1 / r)}×`)
      : (r > CAP_HI ? '↑ &gt;100×' : `↑ ${fmtX(r)}×`));
// a bare window number: ">100×" past the cap, with the uncapped value kept for the tooltip
const fmtWin = (w) => (!Number.isFinite(w) ? '—' : w > CAP_HI ? '&gt;100×' : w < CAP_LO ? '&lt;0.01×' : `${fmtX(w)}×`);
const capNote = (w) => (Number.isFinite(w) && (w > CAP_HI || w < CAP_LO) ? ` <span class="mute">(uncapped ${fmtX(w > CAP_HI ? w : 1 / w)}× — rank order only)</span>` : '');

const ROLE_COL = { oncogene: '#e85d9c', 'tumour suppressor': '#ffb547', 'tumor suppressor': '#ffb547', reference: '#5ee0d0' };
const OTHER_COL = '#8a8f9c';
const roleKey = (r) => String(r || '').toLowerCase().trim();
const roleCol = (r) => ROLE_COL[roleKey(r)] ?? OTHER_COL;
const isOnco = (r) => roleKey(r) === 'oncogene' || roleKey(r).includes('oncogene');
const isRef = (r) => roleKey(r) === 'reference' || roleKey(r).includes('reference') || roleKey(r).includes('housekeep');

const isMane = (iso) => iso?.tags?.includes('MANE_Select');
const repIso = (g) => g.isoforms.find(isMane) ?? g.isoforms[0];

// ---------- lever -> parameter mapping ----------
// targets.json gives {param, direction}. `param` may be a name, a list of names, or {name: extreme}.
// The primary knob. `also_params` are listed by targets.json as secondary and are only sweept when
// the user asks: moving every one of them to its extreme is OUR guess, not the curator's statement.
function leverParams(lv, withSecondary = false) {
  const p = lv?.param;
  let out = [];
  if (!p) out = [];
  else if (typeof p === 'string') out = p.split(/[\s,+]+/).filter(Boolean);
  else if (Array.isArray(p)) out = p.filter(x => typeof x === 'string');
  else if (typeof p === 'object') out = Object.keys(p);
  if (withSecondary && Array.isArray(lv.also_params)) out = out.concat(lv.also_params.filter(x => typeof x === 'string'));
  return [...new Set(out)];
}
const secondaryParams = (lv) => (Array.isArray(lv?.also_params) ? lv.also_params.filter(k => PARAMS[k]) : []);
// `direction` is free text in targets.json ("down", "both", "down at the masked site (up for …)").
// Rule: the extreme is the parameter's max only when the direction says up and never says down.
function dirUp(lv) {
  const s = String(lv?.direction || '').toLowerCase();
  return /\b(up|rais|increas|higher)/.test(s) && !/\b(down|lower|decreas)/.test(s);
}
function extremeOf(lv, k) {
  if (lv.param && typeof lv.param === 'object' && !Array.isArray(lv.param) && Number.isFinite(lv.param[k])) return lv.param[k];
  const d = PARAMS[k];
  if (!d) return null;
  return dirUp(lv) ? d.max : d.min;
}
const isSpliceLever = (lv) => /splice/i.test(lv.model || '') || /splice/i.test(String(lv.param || '')) ||
  leverParams(lv).some(k => !PARAMS[k]);

// ---------- which page does a lever act on? ----------
// Translation parameters live in src/params.js, splicing parameters in src/splice/params.js. A lever
// is placed by its own `page` field if targets.json states one, else by where its parameters live.
// The splice parameter file is read through import.meta.glob so this build never breaks on its
// absence (same rule as the splice model above); SPLICE_PARAM_NAMES is the fallback.
const SPLICE_PARAMS_GLOB = import.meta.glob('../splice/params.js', { eager: true });
const SPLICE_PARAM_NAMES = new Set([
  ...Object.keys(Object.values(SPLICE_PARAMS_GLOB)[0]?.PARAMS ?? {}),
  // fallback list, used only if src/splice/params.js is not in the build
  'initRate', 'pauseProb', 'pausePos', 'pauseEscape', 'elongRate', 'elongCV', 'commitScale',
  'scoreSlope', 'removalTime', 'ctLag', 'exonWindow', 'exonDefGain', 'postTx', 'site_score',
]);
export function leverPage(lv) {
  const p = String(lv?.page || '').toLowerCase();
  if (p === 'transcription' || p === 'translation') return p;
  // the PRIMARY knob decides (a secondary one may well live in the other model)
  for (const ks of [leverParams(lv), leverParams(lv, true)]) {
    if (ks.some(k => PARAMS[k])) return 'translation';            // a translation parameter is moved
    if (ks.some(k => SPLICE_PARAM_NAMES.has(k))) return 'transcription';
  }
  return isSpliceLever(lv) ? 'transcription' : 'translation';
}

// treated parameter set at dose d (0..1)
function treatedParams(lv, base, d, withSecondary) {
  const out = { ...base };
  for (const k of leverParams(lv, withSecondary)) {
    const x = extremeOf(lv, k);
    if (!PARAMS[k] || x == null) continue;
    out[k] = base[k] + d * (x - base[k]);
  }
  return out;
}

// ---------- live model calls ----------
const perMin = (iso, P, opts) => {
  try { return expected(iso, P, opts).perMin; } catch { return NaN; }
};
function ratioOf(iso, Pt, Pb, opts) {
  const a = perMin(iso, Pb, opts), b = perMin(iso, Pt, opts);
  if (!Number.isFinite(a) || !Number.isFinite(b) || a <= 1e-12) return NaN;
  return clamp(b / a, RMIN, RMAX);
}

const allIdx = (arr) => new Set((arr || []).map((_, i) => i));

// Which feature drives this isoform's sensitivity: remove the feature, recompute the ratio,
// and report how much of the log-fold response disappears. All from the same exact model.
function drivers(iso, Pt, Pb) {
  const r0 = ratioOf(iso, Pt, Pb);
  if (!Number.isFinite(r0)) return { r0, list: [] };
  const l0 = Math.log10(r0);
  const tests = [];
  const hp = iso.hairpins || [], uo = iso.uorfs || [];
  if (hp.length) tests.push({ name: 'structured 5′ leader', detail: `${hp.length} hairpin${hp.length > 1 ? 's' : ''} melted`, opts: { melted: allIdx(hp) } });
  const capHP = hp.map((h, i) => [h, i]).filter(([h]) => h.start < Pb.capZone).map(([, i]) => i);
  if (capHP.length && capHP.length < hp.length) tests.push({ name: 'cap-proximal structure', detail: `${capHP.length} hairpin(s) within ${Pb.capZone} nt of the cap melted`, opts: { melted: new Set(capHP) } });
  if (uo.length) tests.push({ name: 'uORFs / upstream AUGs', detail: `${uo.length} uAUG${uo.length > 1 ? 's' : ''} mutated`, opts: { disabledU: allIdx(uo) } });
  const over = uo.map((u, i) => [u, i]).filter(([u]) => u.type === 'overlapping' || u.end > iso.len).map(([, i]) => i);
  if (over.length && over.length < uo.length) tests.push({ name: 'overlapping uORFs', detail: `${over.length} CDS-overlapping uORF(s) mutated`, opts: { disabledU: new Set(over) } });
  const list = [];
  for (const t of tests) {
    const r = ratioOf(iso, Pt, Pb, t.opts);
    if (!Number.isFinite(r)) continue;
    list.push({ ...t, r, d: Math.log10(r) - l0 });
  }
  // start-codon context: force the main start to a strong Kozak
  if (iso.main_kozak?.class && iso.main_kozak.class !== 'strong') {
    const clone = { ...iso, main_kozak: { ...iso.main_kozak, class: 'strong' } };
    const r = ratioOf(clone, Pt, Pb);
    if (Number.isFinite(r)) list.push({ name: 'weak main-start Kozak', detail: `forced to strong (was ${iso.main_kozak.class})`, r, d: Math.log10(r) - l0 });
  }
  list.sort((a, b) => Math.abs(b.d) - Math.abs(a.d));
  return { r0, list };
}

const geoMean = (xs) => (xs.length ? Math.pow(10, xs.reduce((s, x) => s + Math.log10(x), 0) / xs.length) : NaN);

// ---------- chart ----------
function chartSVG(rows) {
  if (!rows.length) return '<div class="tgt-wait">no genes in this set pass the evidence bar.</div>';
  // the SVG scales to the column width, so a phone gets a narrower viewBox (bigger effective text)
  const narrow = innerWidth <= 760;
  const W = narrow ? 376 : 580, lab = narrow ? 60 : 92, padR = narrow ? 56 : 74,
    rowH = narrow ? 20 : 21, top = 16, H = top + rows.length * rowH + 6;
  const vals = rows.map(r => r.ratio).concat([1]);
  let lo = Math.min(...vals), hi = Math.max(...vals);
  lo = Math.min(lo / 1.35, 1 / 1.05); hi = Math.max(hi * 1.35, 1.05);
  // the axis stops at the cap as well: a bar running out to 1e-4 would assert a fold change the
  // model cannot support, so anything past the cap is drawn at the edge
  lo = Math.max(lo, CAP_LO); hi = Math.min(hi, CAP_HI);
  const plotW = W - lab - padR;
  const X = (v) => lab + (Math.log10(clamp(v, CAP_LO, CAP_HI)) - Math.log10(lo)) / (Math.log10(hi) - Math.log10(lo)) * plotW;
  const x1 = X(1);
  // ticks: powers of ten inside the range, plus 1
  const ticks = [];
  for (let e = Math.ceil(Math.log10(lo)); e <= Math.floor(Math.log10(hi)); e++) ticks.push(Math.pow(10, e));
  if (!ticks.includes(1)) ticks.push(1);
  const grid = ticks.map(t => `<line class="${t === 1 ? 'ax' : 'gr'}" x1="${X(t).toFixed(1)}" x2="${X(t).toFixed(1)}" y1="${top - 4}" y2="${H - 4}"/>
    <text class="tk" x="${X(t).toFixed(1)}" y="${top - 7}" text-anchor="middle">${t === 1 ? (narrow ? '1×' : '1× (no change)') : t < 1 ? '÷' + fmtX(1 / t) : '×' + fmtX(t)}</text>`).join('');
  const body = rows.map((r, i) => {
    const y = top + i * rowH, col = roleCol(r.role);
    const xr = X(r.ratio), x = Math.min(x1, xr), w = Math.max(Math.abs(xr - x1), 1.5);
    const txt = fmtArrow(r.ratio);
    const tip = `${r.gene} (${r.role || 'unclassified'}) · predicted output ratio ${f2(r.ratio, 3)} = treated / untreated · ${r.isoName} (${r.iso?.len ?? '?'} nt leader) · ${RANK_TIP}`;
    return `<g data-gene="${esc(r.gene)}"><title>${esc(tip)}</title>
      <rect class="hit" x="0" y="${y}" width="${W}" height="${rowH}"/>
      <text class="gl${isRef(r.role) ? ' ref' : ''}" x="${lab - 7}" y="${y + rowH / 2 + 4}" text-anchor="end">${esc(r.gene)}</text>
      <rect x="${x.toFixed(1)}" y="${y + 4}" width="${w.toFixed(1)}" height="${rowH - 9}" rx="2" fill="${col}"/>
      <text class="gv" x="${W - padR + 8}" y="${y + rowH / 2 + 4}" fill="${col}">${txt}</text></g>`;
  }).join('');
  return `<svg viewBox="0 0 ${W} ${H}" class="tgt-chart">${grid}${body}</svg>
    <div class="tgt-capnote">Axis and labels are capped at <b>×100 / ×0.01</b>: beyond that the model's range (4.3×10⁶-fold across these leaders) is unphysical, so bars at the edge mean <b>rank order</b>, not fold change.</div>`;
}

function driverSVG(list) {
  const cs = list.filter(c => Number.isFinite(c.d) && Math.abs(c.d) > 1e-4).slice(0, 5);
  if (!cs.length) return '<div class="tgt-why">no single leader feature changes this gene’s response to this lever.</div>';
  const narrow = innerWidth <= 760;
  const W = 420, lab = narrow ? 104 : 150, padR = 74, rowH = 18, H = cs.length * rowH + 4;
  const max = Math.max(...cs.map(c => Math.abs(c.d)), 1e-9);
  const plotW = W - lab - padR, x0 = lab + plotW / 2, sc = (plotW / 2 - 4) / max;
  const rows = cs.map((c, i) => {
    const y = 2 + i * rowH, v = c.d, w = Math.max(Math.abs(v) * sc, 1.5), x = v >= 0 ? x0 : x0 - w;
    const col = v >= 0 ? '#e85d9c' : '#4fc3f7';
    const lbl = (v >= 0 ? '+' : '−') + f2(Math.abs(v)) + ' log₁₀';
    return `<g><title>${esc(c.name)}: ${esc(c.detail || '')} → output ratio becomes ${f2(c.r, 3)} instead of the full-model value; ${v >= 0 ? 'this feature drives the response' : 'this feature buffers the response'}</title>
      <rect x="0" y="${y}" width="${W}" height="${rowH}" fill="transparent"/>
      <text x="${lab - 4}" y="${y + 12}" text-anchor="end" class="wl">${esc(c.name)}</text>
      <rect x="${x.toFixed(1)}" y="${y + 3}" width="${w.toFixed(1)}" height="${rowH - 6}" rx="2" fill="${col}"/>
      <text x="${(v >= 0 ? x + w + 3 : x - 3).toFixed(1)}" y="${y + 12}" text-anchor="${v >= 0 ? 'start' : 'end'}" class="wv">${lbl}</text></g>`;
  }).join('');
  return `<svg viewBox="0 0 ${W} ${H}" class="wsvg"><line x1="${x0}" x2="${x0}" y1="0" y2="${H}" class="wax"/>${rows}</svg>`;
}

const HONEST = `<div class="tgt-honest"><div class="hh">What this is, and what it is not</div><ul>
  <li><b>Predicted differential sensitivity from a kinetic model</b> — the same ribosome-scanning model as the ② Translation page, not a measured drug response.</li>
  <li><b>The “dose” is a model parameter change, not a drug concentration.</b> 100% = that parameter's extreme in <code>src/params.js</code>; nothing here is calibrated to any compound's IC<sub>50</sub>.</li>
  <li><b>Parameters are assumed</b> unless cited in <code>src/params.js</code> (see the Sources list on the ② Translation page); the baseline is the app's default set, calm cell.</li>
  <li><b>The model has no IRES, no RBP binding and no ribosome stalling / collision mechanism</b>, so any transcript whose regulation runs through those is modelled wrongly.</li>
  <li><b>It showed no ability to rank isoforms against measured TE</b> — within genes, performance was at chance. Between-gene structure is all this ranking can claim.</li>
  <li>Per-gene output here uses one leader per gene (MANE Select where annotated); <b>isoform abundance is not modelled in this panel</b>.</li>
  <li><b>Hypothesis generation for experiments, not evidence of efficacy.</b> Nothing on this page says a drug works.</li>
</ul></div>`;

// ---------- the panel ----------
export function mountTarget(hooks = {}) {
  // one instance per page: ① transcription (splicing levers) / ② translation (scanning levers)
  const PAGE = hooks.page === 'transcription' ? 'transcription' : 'translation';
  const PAGE_LABEL = PAGE === 'transcription' ? '① Transcription' : '② Translation';
  const OTHER = PAGE === 'transcription' ? 'translation' : 'transcription';
  const OTHER_LABEL = PAGE === 'transcription' ? '② Translation' : '① Transcription';
  const WHERE = PAGE === 'transcription'
    ? 'Levers that act while the gene is being transcribed and spliced — splice sites, SF3b, CLK/SRPK, Pol II elongation rate.'
    : 'Levers that act on the scanning ribosome — eIF4A, the integrated stress response, mTOR / cap dependence, start-codon context.';
  const el = document.createElement('div');
  el.id = 'tgtPane-' + PAGE; el.className = 'tgt-pane';
  const btn = document.createElement('button');
  btn.id = 'tgtToggle-' + PAGE; btn.type = 'button'; btn.className = 'tgt-toggle';
  const btnHTML = (on) => `<i>◎</i><b class="blong">${on ? 'close the levers panel'
    : `Therapeutic levers${PAGE === 'transcription' ? ' (splicing)' : ' (translation)'}`}</b><b class="bshort">${on ? 'close' : 'levers'}</b>`;
  btn.innerHTML = btnHTML(false);

  const st = { targets: null, spliceData: null, leverId: null, dose: 0.6, gene: null, allGenes: false, secondary: false, open: false, err: null };

  const getData = () => (hooks.getData ? hooks.getData() : null);
  const getEv = () => (hooks.getEvidence ? hooks.getEvidence() : null);

  // Same evidence bar as tab 1: data/evidence.json decides for any gene it audits; for a gene only
  // targets.json knows about, its own tier is used (targets.json states it applies the same bar).
  function passesBar(t) {
    if (st.allGenes || isRef(t.role)) return true;   // comparators carry their own citation
    const ev = getEv();
    const e = ev?.genes?.[t.gene];
    if (e) return e.show_as_example === true;
    if (!ev?.genes) return true;                     // no audit file: same fallback as tab 1
    return t.tier === 'A' && t.isoform_resolved !== false;
  }
  function evBadge(t) {
    const name = typeof t === 'string' ? t : t.gene;
    const e = getEv()?.genes?.[name];
    const tier = e?.tier ?? (typeof t === 'object' ? t.tier : null);
    const cav = e?.caveat ?? (typeof t === 'object' ? t.caveat : '');
    if (!tier) return `<span class="evb tB" title="${esc(cav || '')}">no evidence tier</span>`;
    const weak = e ? e.show_as_example !== true : (typeof t === 'object' && t.isoform_resolved === false);
    return `<span class="evb t${esc(tier)}" title="${esc(cav || '')}">${tier === 'A' ? '● experimentally characterised leader' : '◐ indirect / transcriptome-wide evidence'}${weak ? ' · no isoform-level evidence' : ''}</span>`;
  }

  // rows for the current lever + dose, computed live
  function compute() {
    const data = getData(), T = st.targets;
    if (!data || !T) return null;
    const lv = (T.levers || []).find(l => l.id === st.leverId);
    if (!lv) return null;
    const base = values();
    const Pt = treatedParams(lv, base, st.dose, st.secondary);
    const byName = new Map(data.genes.map(g => [g.gene, g]));
    // the comparator set lives in its own block in targets.json; score it alongside the targets
    const refs = (T.reference_set || []).map(r => (typeof r === 'string' ? { gene: r } : r))
      .filter(r => r?.gene).map(r => ({ ...r, role: r.role || 'reference' }));
    const seen = new Set();
    const cand = [...(T.genes || []), ...refs].filter(t => t?.gene && !seen.has(t.gene) && seen.add(t.gene));
    const rows = [], missing = [], barred = [];
    for (const t of cand) {
      const g = byName.get(t.gene);
      if (!g) { missing.push(t.gene); continue; }     // no leader in utr_data.json yet
      if (!passesBar(t)) { barred.push(t.gene); continue; }
      const iso = repIso(g);
      const ratio = ratioOf(iso, Pt, base);
      if (!Number.isFinite(ratio)) { missing.push(t.gene); continue; }
      rows.push({ ...t, role: t.role || '', gene: t.gene, g, iso, isoName: iso.name, ratio });
    }
    rows.sort((a, b) => a.ratio - b.ratio);
    const refR = rows.filter(r => isRef(r.role)).map(r => r.ratio);
    const onco = rows.filter(r => isOnco(r.role));
    // how many oncogenes targets.json lists at all, vs how many could be scored
    const nOncoListed = (T.genes || []).filter(t => isOnco(t.role)).length;
    const oncoNames = new Set(onco.map(r => r.gene));
    // ---- length-matched controls -------------------------------------------------------------
    // The single biggest driver of the response in this model is leader length / structure, so an
    // oncogene-vs-housekeeping comparison is confounded: the references are 50-84 nt and the
    // oncogenes 175-424 nt. Each scored oncogene is therefore also matched to the nearest-length
    // NON-oncogene leader in the same data set, and the window recomputed against those.
    const oncoGenes = new Set((T.genes || []).filter(t => isOnco(t.role)).map(t => t.gene));
    const pool = [];
    for (const g of data.genes) {
      if (oncoGenes.has(g.gene)) continue;
      for (const iso of g.isoforms) {
        const r = ratioOf(iso, Pt, base);
        if (Number.isFinite(r) && Number.isFinite(iso.len)) pool.push({ gene: g.gene, name: iso.name, len: iso.len, ratio: r });
      }
    }
    const used = new Set(), matched = [];
    for (const o of onco) {
      let best = null;
      for (const c of pool) {
        if (used.has(c.gene + '|' + c.name)) continue;
        if (!best || Math.abs(c.len - o.iso.len) < Math.abs(best.len - o.iso.len)) best = c;
      }
      if (best) { used.add(best.gene + '|' + best.name); matched.push({ ...best, forGene: o.gene, forLen: o.iso.len }); }
    }
    // the gene whose single-gene extremum used to be printed as "the window", and its closest
    // non-oncogene twin by leader length: if that twin scores the same, the extremum is about length
    const ext = (onco.length ? onco : rows).slice().sort((a, b) => a.ratio - b.ratio)[0];
    let twin = null;
    if (ext) for (const c of pool) if (!twin || Math.abs(c.len - ext.iso.len) < Math.abs(twin.len - ext.iso.len)) twin = c;
    // ---- is the response just structure? r(log10 response, most stable hairpin dG) -------------
    const xs = [], ys = [];
    for (const g of data.genes) for (const iso of g.isoforms) {
      const hp = iso.hairpins || []; if (!hp.length) continue;
      const r = ratioOf(iso, Pt, base);
      if (!Number.isFinite(r) || r <= 0) continue;
      xs.push(Math.min(...hp.map(h => h.dG))); ys.push(Math.log10(r));
    }
    let corr = null;
    if (xs.length > 3) {
      const n = xs.length, mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
      let sxy = 0, sxx = 0, syy = 0;
      for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2; }
      const r = sxy / Math.sqrt(sxx * syy || 1);
      if (Number.isFinite(r)) corr = { r, n };
    }
    return {
      lv, base, Pt, rows, missing, barred, onco, nOncoListed,
      refGeo: geoMean(refR), nRef: refR.length,
      oncoGeo: geoMean(onco.map(r => r.ratio)),
      matched, matchGeo: geoMean(matched.map(m => m.ratio)), twin, ext, corr,
      oncoNames,
    };
  }

  // `disabled` + `disabled_reason` in targets.json take a lever off the demo path entirely: it stays
  // visible and selectable-looking, but it cannot be scored and the reason is printed on it.
  const leverOK = (lv) => { const ps = leverParams(lv); return !lv?.disabled && !isSpliceLever(lv) && ps.length > 0 && ps.every(k => PARAMS[k]); };

  // why a splicing lever cannot be scored here — stated precisely, not hand-waved
  function spliceWhy(lv) {
    if (!spliceMod) return 'splicing lever — the splicing model (src/splice/model_exact.js) is not in this build';
    const have = (st.spliceData?.genes || []).map(g => g.gene);
    if (!have.length) return 'splicing lever — the splicing model is loaded but splice_data.json has no genes yet';
    const want = (st.targets?.genes || []).filter(g => (g.levers || []).includes(lv.id)).map(g => g.gene);
    const ov = want.filter(w => have.includes(w));
    if (!ov.length) return `splicing lever — the splicing model is loaded, but splice_data.json covers only ${have.join(', ')}; the genes this lever targets (${want.join(', ') || 'see targets.json'}) have no transcription model yet, so their isoform abundance cannot be moved`;
    return `splicing lever — ${ov.join(', ')} could be scored, but output = abundance × protein-per-mRNA is not wired up in this build`;
  }
  const citeLink = (c, u, fb = 'source') => (u ? `<a href="${esc(u)}" target="_blank" rel="noopener">${esc(c || fb)}</a>` : esc(c || ''));

  function leverLI(lv) {
    const on = lv.id === st.leverId, ok = leverOK(lv), ps = leverParams(lv, st.secondary);
    const pdesc = ok ? ps.map(k => `${k} ${f2(values()[k])} → ${f2(extremeOf(lv, k))}`).join(' · ') : String(lv.param ?? '—');
    const stat = String(lv.status || lv.clinical_status || '').split(/[;(]/)[0].trim();
    return `<li class="${on ? 'on' : ''}${ok ? '' : ' off'}" data-lever="${esc(lv.id)}" ${ok ? '' : 'aria-disabled="true"'}>
      <div class="lv-t">${esc(lv.label || lv.id)}${stat ? `<span class="lv-s" title="${esc(lv.status || lv.clinical_status)}">${esc(stat.length > 46 ? stat.slice(0, 45) + '…' : stat)}</span>` : ''}</div>
      <div class="lv-p">model parameter · ${esc(pdesc)}</div>
      ${lv.agent ? `<div class="lv-a">${esc(lv.agent)}</div>` : ''}
      ${lv.mechanism ? `<div class="lv-m clamp">${esc(lv.mechanism)}</div>` : ''}
      ${lv.cite ? `<div class="lv-m">${citeLink(lv.cite, lv.url)}</div>` : ''}
      ${ok ? '' : `<div class="lv-x">${lv.disabled ? `<b>disabled:</b> ${esc(lv.disabled_reason || 'this lever does not reproduce its own cited mechanism in this model')}`
        : isSpliceLever(lv) ? esc(spliceWhy(lv)) : 'no model parameter mapped'}</div>`}</li>`;
  }

  // everything targets.json knows about the selected lever: agent, mechanism, status, quote, citations
  function leverDossier(lv) {
    if (!lv) return '';
    const sup = (lv.supporting || []).map(s => `<li>${citeLink(s.cite, s.url)}${s.quote ? ` — “${esc(s.quote)}”` : ''}</li>`).join('');
    return `<div class="tgt-dos"><div class="ftitle">The lever · ${esc(lv.label || lv.id)}</div>
      ${lv.disabled ? `<div class="tgt-dis"><b>Disabled — not scored.</b> ${esc(lv.disabled_reason || '')}</div>` : ''}
      ${lv.direction_warning ? `<div class="tgt-dir"><b>Which way the dose goes:</b> ${esc(lv.direction_warning)}</div>` : ''}
      ${lv.agent ? `<div class="dr"><span>agent</span><div>${esc(lv.agent)}</div></div>` : ''}
      ${lv.mechanism ? `<div class="dr"><span>mechanism</span><div>${esc(lv.mechanism)}</div></div>` : ''}
      ${lv.status || lv.clinical_status ? `<div class="dr"><span>clinical status</span><div>${esc(lv.status || lv.clinical_status)}</div></div>` : ''}
      <div class="dr"><span>model parameter</span><div><code>${esc(String(lv.param ?? '—'))}</code> · direction: ${esc(lv.direction || '—')}${lv.also_params?.length ? `. Secondary knobs listed by targets.json: <code>${lv.also_params.map(esc).join('</code>, <code>')}</code> — <b>${st.secondary ? 'moved with the dose' : 'not moved'}</b>, since taking each to its extreme would be the app's guess, not the curator's claim.` : '.'}
        <b>The drug → parameter mapping is the app's, not the papers'</b> — the cited work does not mention these parameters.</div></div>
      ${lv.quote ? `<div class="dr"><span>what the paper says</span><div>“${esc(lv.quote)}” <span class="mute">(${esc(lv.quote_source || 'source')}, ${citeLink(lv.cite, lv.url)})</span></div></div>` : `<div class="dr"><span>source</span><div>${citeLink(lv.cite, lv.url)}</div></div>`}
      ${sup ? `<div class="dr"><span>supporting</span><div><ul>${sup}</ul></div></div>` : ''}
      ${lv.caveat ? `<div class="dr"><span>caveat</span><div class="cav">${esc(lv.caveat)}</div></div>` : ''}</div>`;
  }

  // what could not be scored, and why — stated, not silently dropped
  function gapNote(c) {
    const req = st.targets?.pipeline_requests;
    const bits = [];
    if (c.missing.length) bits.push(`<b>${c.missing.length} gene${c.missing.length > 1 ? 's' : ''} in targets.json ${c.missing.length > 1 ? 'have' : 'has'} no 5′ leader in <code>utr_data.json</code> yet</b> (${c.missing.map(esc).join(', ')}) — the pipeline has to be rebuilt with them before the model can score them.`);
    if (c.barred.length) bits.push(`${c.barred.length} did not pass the evidence bar: ${c.barred.map(esc).join(', ')}.`);
    if (!c.nRef) bits.push('No comparator gene could be scored, so <b>no therapeutic window can be computed</b> — the window needs the reference set in the leader data.');
    if (req?.note) bits.push(`<span class="mute">${esc(req.note)}</span>`);
    return bits.length ? `<div class="tgt-gap">${bits.join(' ')}</div>` : '';
  }

  function headline(c) {
    const { lv, rows } = c;
    const pct = Math.round(st.dose * 100);
    const ps = leverParams(lv, st.secondary).filter(k => PARAMS[k]);
    const pline = ps.map(k => `<code>${k}</code> ${f2(c.base[k])} → <b>${f2(c.Pt[k])}</b>`).join(' · ');
    const cand = (c.onco.length ? c.onco : rows).slice().sort((a, b) => a.ratio - b.ratio)[0];
    if (!cand || !Number.isFinite(c.refGeo) || !Number.isFinite(c.oncoGeo)) {
      const refNames = (st.targets?.reference_set || []).map(r => (typeof r === 'string' ? r : r?.gene)).filter(Boolean);
      return `<div class="tgt-win"><div class="wk">Therapeutic window · ${esc(lv.label || lv.id)} at ${pct}% of the parameter range</div>
        <div class="wt">${cand ? `${esc(cand.gene)} output moves <b title="${esc(RANK_TIP)}">${fmtArrow(cand.ratio)}</b> at this dose <span class="mute">(one gene, not a window)</span>.` : 'No gene in this set could be scored.'}
        ${Number.isFinite(c.refGeo) ? '' : ` <b>No window yet:</b> the comparator set${refNames.length ? ` (${refNames.map(esc).join(', ')})` : ''} has no 5′ leader in <code>utr_data.json</code>, so there is nothing to measure the drop against.`}</div>
        <div class="wk" style="margin-top:4px">${pline} · predicted, model only</div></div>`;
    }
    // SET vs SET is the window. The old headline divided the geometric mean of the references by
    // the single most-affected oncogene — an extremum over one set against an average of the other,
    // which inflates it (eIF4A at full dose: 360x that way, 189x set-vs-set, ~5x length-matched).
    const win = Number.isFinite(c.oncoGeo) ? c.refGeo / c.oncoGeo : NaN;
    const lenWin = Number.isFinite(c.matchGeo) && Number.isFinite(c.oncoGeo) ? c.matchGeo / c.oncoGeo : NaN;
    const dirO = c.oncoGeo < 1 ? 'falls' : 'rises', dirR = c.refGeo < 1 ? 'falls' : 'rises';
    const wn = win >= 1.05 ? fmtWin(win)
      : `<span class="none">no window${win < 0.95 ? ` · ${fmtWin(1 / win)} the wrong way` : ''}</span>`;
    const lens = c.onco.map(r => r.iso.len), refLens = c.rows.filter(r => isRef(r.role)).map(r => r.iso.len);
    return `<div class="tgt-win"><div class="wk">Therapeutic window · ${esc(lv.label || lv.id)} at ${pct}% of the parameter range</div>
      <div class="wn" title="${esc(RANK_TIP)}">${wn}</div>
      <div class="wt"><b>Set vs set</b> (the only comparison this headline makes): the oncogene set (n=${c.onco.length}, geometric mean) ${dirO}
        <b>${fmtWin(c.oncoGeo < 1 ? 1 / c.oncoGeo : c.oncoGeo)}</b>${capNote(c.oncoGeo < 1 ? 1 / c.oncoGeo : c.oncoGeo)}
        while the reference set (n=${c.nRef}) ${dirR} <b>${fmtWin(c.refGeo < 1 ? 1 / c.refGeo : c.refGeo)}</b>.
        <span class="mute">n = ${c.onco.length} of ${c.nOncoListed} oncogenes listed in <code>targets.json</code> have 5′ leader data${c.onco.length < c.nOncoListed ? ` — the other ${c.nOncoListed - c.onco.length} are not scored at all` : ''}.</span></div>
      ${Number.isFinite(lenWin) ? `<div class="wlen"><b>${lenWin >= 0.95 ? 'Most of this separation is leader length / structure, not oncogene biology.' : 'There is no favourable window here.'}</b>
        Against <b>length-matched non-oncogene controls</b> from the same data (nearest leader length, n=${c.matched.length}) the window is only <b>${fmtWin(lenWin)}</b>${capNote(lenWin)}
        ${st.leverId === 'eif4a' && st.dose >= 0.99 ? ' — an independent audit measured ~4.4× for eIF4A at full dose.' : ''}
        ${lenWin < 0.95 ? '<b>This is below 1, i.e. the length-matched controls are hit <i>harder</i> than the oncogenes — the window is inverted, not favourable.</b>' : ''} The oncogenes here are ${Math.min(...lens)}–${Math.max(...lens)} nt
        and the reference set ${refLens.length ? `${Math.min(...refLens)}–${Math.max(...refLens)} nt` : 'much shorter'}, so the two sets differ in length as much as in biology.
        ${c.corr ? `Across all ${c.corr.n} leaders, r(log₁₀ response, most-stable-hairpin ΔG) = <b>${f2(c.corr.r)}</b>.` : ''}
        ${c.twin && c.ext ? `The single most-affected gene (${esc(c.ext.gene)}, ${c.ext.iso.len} nt) is matched by a non-oncogene leader of nearly the same length —
          ${esc(c.twin.gene)}/${esc(c.twin.name)} (${c.twin.len} nt), which scores ${fmtArrow(c.twin.ratio)}.` : ''}
        <span class="matchlist">matched pairs: ${c.matched.map(m => `${esc(m.forGene)} ${m.forLen} nt → ${esc(m.gene)}/${esc(m.name)} ${m.len} nt ${fmtArrow(m.ratio)}`).join(' · ')}</span></div>` : ''}
      <div class="wk" style="margin-top:5px">${pline} · predicted, model only</div></div>`;
  }

  function detailHTML(c) {
    const row = c.rows.find(r => r.gene === st.gene);
    if (!row) return '';
    const g = row.g;
    const isoRows = g.isoforms.map(iso => ({ iso, ratio: ratioOf(iso, c.Pt, c.base) })).filter(x => Number.isFinite(x.ratio));
    const lomax = Math.max(...isoRows.map(x => Math.abs(Math.log10(x.ratio))), 0.05);
    const bars = isoRows.map(({ iso, ratio }) => {
      const l = Math.log10(ratio), half = 50, w = Math.min(50, Math.abs(l) / lomax * 50);
      return `<div class="tgt-iso"><span class="nm" title="${esc(iso.name)} · ${iso.len} nt${isMane(iso) ? ' · MANE Select' : ''}">${esc(iso.name)}${isMane(iso) ? ' ·MANE' : ''}</span>
        <span class="tr"><i style="left:${l < 0 ? half - w : half}%;width:${w}%;background:${l < 0 ? '#e85d9c' : '#5ee0a0'}"></i></span>
        <span title="${esc(RANK_TIP)}">${fmtArrow(ratio)}</span></div>`;
    }).join('');
    const { list } = drivers(row.iso, c.Pt, c.base);
    const chips = motifsOf(row.iso).filter(m => m.type !== 'uORF_stop').slice(0, 10).map(m =>
      `<span class="mchip${m.strength === 'weak' ? ' weak' : ''}" style="--mc:${motifStyle(m.type).color}" title="${esc(motifStyle(m.type).name)} at nt ${m.start}–${m.end}${m.strength ? ' · ' + m.strength : ''}">${esc(chipText(m))}</span>`).join('');
    const nU = (row.iso.uorfs || []).length, nH = (row.iso.hairpins || []).length;
    return `<div class="tgt-det"><div class="dh"><b>${esc(row.gene)}</b> ${evBadge(row)}
        <span class="mute">${esc(row.role || 'unclassified')}${row.disease ? ' · ' + esc(row.disease) : ''}</span>
        <button type="button" class="dgo">open on the ${esc(PAGE_LABEL)} page →</button></div>
      ${row.why_targetable || row.why ? `<div class="tgt-why">${esc(row.why_targetable || row.why)}</div>` : ''}
      ${row.mechanism ? `<div class="tgt-why"><b>Leader / event:</b> ${esc(row.mechanism)}${row.event ? ` <span class="mute">(${esc(row.event)})</span>` : ''}</div>` : ''}
      ${row.quote ? `<div class="tgt-why">“${esc(row.quote)}” <span class="mute">(${esc(row.quote_source || 'source')}${row.measured ? '; measured: ' + esc(row.measured) : ''}, ${citeLink(row.cite, row.url)})</span></div>`
        : row.cite ? `<div class="tgt-why">${citeLink(row.cite, row.url)}</div>` : ''}
      ${row.caveat ? `<div class="tgt-why cav">${esc(row.caveat)}</div>` : ''}
      <div class="ftitle" style="margin-top:7px">What drives ${esc(row.gene)}'s sensitivity to this lever</div>
      ${driverSVG(list)}
      <div class="tgt-key"><span><i style="background:#e85d9c"></i>removing this feature removes the response (it drives it)</span><span><i style="background:#4fc3f7"></i>removing it makes the response bigger (it buffers)</span></div>
      <div class="ftitle" style="margin-top:9px">Per-isoform predicted ratio · ${esc(row.iso.name)} leader: ${row.iso.len} nt · ${nU} uORF${nU === 1 ? '' : 's'} · ${nH} hairpin${nH === 1 ? '' : 's'} · Kozak ${esc(row.iso.main_kozak?.class ?? '—')}</div>
      ${bars}
      <div class="tgt-chips">${chips || '<span class="mute" style="font-size:11px">no annotated motifs</span>'}</div></div>`;
  }

  function render() {
    const data = getData(), T = st.targets;
    const nOther = (T?.levers || []).filter(lv => leverPage(lv) === OTHER).length;
    const head = `<div class="tgt-h"><h2>Therapeutic levers · ${esc(PAGE_LABEL)}</h2>
        <span class="mute" style="font-size:11.5px">the model's parameters are the targets of real drugs</span>
        <button type="button" class="tgt-close">← back to ${esc(PAGE_LABEL)}</button></div>
      <p class="tgt-sub">${esc(WHERE)} Pick a lever, move the dose, and the predicted output ratio (treated / untreated) is recomputed for every gene
        with the same exact scanning model the translation page uses. Sorted, log scale. Click a gene for the per-isoform breakdown.</p>
      ${nOther ? `<p class="tgt-sub tgt-elsewhere">${nOther} lever${nOther > 1 ? 's act' : ' acts'} at the other stage and ${nOther > 1 ? 'are' : 'is'} on the
        <b>${esc(OTHER_LABEL)}</b> page instead. <button type="button" class="tgt-other">go to ${esc(OTHER_LABEL)} →</button></p>` : ''}
      <div class="tgt-flag">Predicted differential sensitivity from a kinetic model — <b>not</b> measured drug response.
        The “dose” is a change to a model parameter, <b>not</b> a drug concentration. Hypothesis generation, not evidence of efficacy.
        <span class="mute">Full caveats at the bottom of this tab.</span></div>`;
    if (!data) { el.innerHTML = head + '<div class="tgt-wait">loading transcripts…</div>' + HONEST; return; }
    if (!T) {
      el.innerHTML = head + `<div class="tgt-wait">loading <code>targets.json</code>…${st.err ? ` <span class="mute">(${esc(st.err)} — the levers + oncogene set are written by the TARGET-DATA agent into <code>app/public/targets.json</code>; this tab polls for it)</span>` : ''}</div>` + HONEST;
      bind(); return;
    }
    // only this page's levers (split by where the lever acts — see leverPage)
    const levers = (T.levers || []).filter(lv => leverPage(lv) === PAGE);
    if (!st.leverId || !levers.some(l => l.id === st.leverId)) st.leverId = (levers.find(leverOK) ?? levers[0])?.id ?? null;
    const c = compute();
    // A lever that cannot be scored must NOT get a window: with no parameter moved every ratio is
    // exactly 1.00x, which would read as "this drug does nothing" instead of "this build cannot say".
    const lvSel = levers.find(l => l.id === st.leverId) ?? null;
    const scorable = !!lvSel && leverOK(lvSel);
    const notScorable = (lv) => (!lv ? 'No lever is selected.'
      : `<b>${esc(lv.label || lv.id)} cannot be scored in this build</b> — ${lv.disabled ? esc(lv.disabled_reason || 'this lever does not reproduce its own cited mechanism in this model')
        : isSpliceLever(lv) ? esc(spliceWhy(lv)) : 'no model parameter is mapped to it'}.
        <b>No therapeutic window and no per-gene ratios are shown for it</b>; the dose slider moves nothing. What <code>targets.json</code> knows about the lever is below.`);
    el.innerHTML = head + `<div class="tgt-grid">
      <div><div class="ftitle">Lever${PAGE === 'transcription' ? 's acting during transcription / splicing' : 's acting on translation'}</div>
        <ul class="tgt-levers">${levers.map(leverLI).join('') || `<li class="off">no lever in <code>targets.json</code> acts at this stage</li>`}</ul>
        <div class="tgt-dose"><div class="ftitle">Dose</div>
          <input type="range" min="0" max="100" step="1" value="${Math.round(st.dose * 100)}" id="tgtDose" aria-label="dose"${scorable ? '' : ' disabled'}>
          <div class="dval"><b>${Math.round(st.dose * 100)}%</b> of the way from baseline to the parameter's extreme</div>
          <div class="dnote">This is a model parameter change, not a drug concentration.</div>
          ${(() => { const lv = levers.find(l => l.id === st.leverId);
            return lv?.direction_warning ? `<div class="tgt-dir">${esc(lv.direction_warning)}</div>` : ''; })()}
          ${(() => { const lv = levers.find(l => l.id === st.leverId), sec = secondaryParams(lv);
            return sec.length ? `<label class="tgt-allg"><input type="checkbox" id="tgtSec" ${st.secondary ? 'checked' : ''}> also move the secondary knobs targets.json lists (<code>${sec.map(esc).join('</code>, <code>')}</code>)</label>` : ''; })()}</div>
        <label class="tgt-allg"><input type="checkbox" id="tgtAll" ${st.allGenes ? 'checked' : ''}> show genes without isoform-level evidence</label>
      </div>
      <div>${c && scorable ? headline(c) : `<div class="tgt-gap">${notScorable(lvSel)}</div>`}
        <div class="tgt-chartbox"><div class="ftitle">Predicted output ratio, treated / untreated (model)</div>
          ${c && scorable ? chartSVG(c.rows) : '<div class="tgt-wait">not computed — this lever moves no model parameter in this build.</div>'}
          <div class="tgt-key"><span><i style="background:${ROLE_COL.oncogene}"></i>oncogene</span><span><i style="background:${ROLE_COL['tumour suppressor']}"></i>tumour suppressor</span><span><i style="background:${ROLE_COL.reference}"></i>reference / comparator set</span><span><i style="background:${OTHER_COL}"></i>unclassified</span></div>
          ${c && scorable ? gapNote(c) : ''}</div>
        ${c && scorable ? detailHTML(c) : ''}
        ${leverDossier(c?.lv ?? levers.find(l => l.id === st.leverId))}</div></div>` + HONEST;
    bind(c);
  }

  function bind(c) {
    el.querySelector('.tgt-close')?.addEventListener('click', leave);
    el.querySelector('.tgt-other')?.addEventListener('click', () => { setOpen(false); hooks.setPage?.(OTHER, { levers: true }); });
    el.querySelector('#tgtDose')?.addEventListener('input', (e) => { st.dose = clamp(+e.target.value / 100, 0, 1); render(); });
    el.querySelector('#tgtAll')?.addEventListener('change', (e) => { st.allGenes = e.target.checked; render(); });
    el.querySelector('#tgtSec')?.addEventListener('change', (e) => { st.secondary = e.target.checked; render(); });
    el.querySelectorAll('.tgt-levers li[data-lever]').forEach(li => li.addEventListener('click', () => {
      st.leverId = li.dataset.lever; render();
      // one column: the chart is above the lever list, so bring the new answer back into view
      if (innerWidth <= 1100) el.querySelector('.tgt-win, .tgt-chartbox')?.scrollIntoView({ block: 'start', behavior: 'smooth' });
    }));
    el.querySelectorAll('.tgt-chart g[data-gene]').forEach(g => g.addEventListener('click', () => {
      st.gene = g.dataset.gene === st.gene ? null : g.dataset.gene; render();
    }));
    const go = el.querySelector('.dgo');
    if (go && c) go.addEventListener('click', () => { if (hooks.openGene?.(st.gene)) leave(); });
  }

  // --- page wiring. The panel is no longer a page of its own: each page carries its own button,
  // and opening it overlays that page (the page switcher stays visible and keeps working).
  // hidden until the app says which page is showing (setPage -> setPageActive)
  function mountButton() { btn.style.display = 'none'; document.body.appendChild(btn); }
  const enter = () => setOpen(true);
  const leave = () => setOpen(false);          // back to this page's own stage

  function setOpen(on) {
    st.open = on;
    el.classList.toggle('open', on);
    // body.tgt-open hides the stage chrome; it must stay on only while SOME pane is open
    document.body.classList.toggle('tgt-open', !!document.querySelector('.tgt-pane.open'));
    btn.classList.toggle('on', on);
    btn.innerHTML = btnHTML(on);
    if (on) {
      const r = document.getElementById('tabs')?.getBoundingClientRect();
      el.style.top = r && r.height ? `${Math.round(r.bottom + 6)}px` : '';
      render(); el.scrollTop = 0;
    }
  }
  // this page became (in)active: only the active page shows its levers button, and leaving a page
  // closes its pane so the next page opens on its own stage
  function setPageActive(on) {
    btn.style.display = on ? '' : 'none';
    if (!on && st.open) setOpen(false);
  }

  async function fetchTargets() {
    for (let i = 0; i < 160 && !st.targets; i++) {
      try {
        const r = await fetch('./targets.json', { cache: 'no-store' });
        if (r.ok) {
          const j = await r.json();
          if (j && (j.levers || j.genes)) {
            st.targets = j; st.err = null;
            await loadSplice();
            try { const s = await fetch('./splice_data.json'); if (s.ok) { const sj = await s.json(); if (sj?.genes) st.spliceData = sj; } } catch {}
            if (st.open) render();
            return;
          }
          st.err = 'targets.json has no levers or genes yet';
        } else st.err = `not published yet (HTTP ${r.status})`;
      } catch (e) { st.err = 'not published yet'; }
      if (st.open) render();
      await new Promise(res => setTimeout(res, 4000));
    }
  }

  btn.addEventListener('click', () => (st.open ? leave() : enter()));
  document.addEventListener('keydown', (e) => {
    if (['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
    if (e.code === 'Escape' && st.open) leave();
  });

  return {
    el, btn, page: PAGE,
    setActive: setOpen,
    setPageActive,
    isOpen: () => st.open,
    mount(rootEl) {
      (rootEl ?? document.body).appendChild(el);
      mountButton();
      fetchTargets();
      return this;
    },
    update() { if (st.open) render(); },
    open: () => setOpen(true),
    close: () => setOpen(false),
  };
}

// Simple contract for main.js: mount(rootEl, ctx) / update(ctx) / setPage(page).
// ctx = { data, evidence, openGene(name, page) -> bool, setPage(page) }
// One panel per page, same module, levers filtered by leverPage().
const panels = [];
let CTX = {};
export function mount(rootEl, ctx = {}) {
  CTX = { ...CTX, ...ctx };
  for (const page of ['translation', 'transcription']) {
    panels.push(mountTarget({
      page,
      getData: () => CTX.data,
      getEvidence: () => CTX.evidence,
      openGene: (name) => !!CTX.openGene?.(name, page),
      setPage: CTX.setPage ? ((p, o) => CTX.setPage(p, o)) : null,
    }).mount(rootEl));
  }
  return panels;
}
export function update(ctx = {}) { CTX = { ...CTX, ...ctx }; panels.forEach(p => p.update()); }
// which page is showing: only that page's levers button is live, and the other page's pane closes
export function setPage(page, opts = {}) {
  for (const p of panels) {
    p.setPageActive(p.page === page);
    if (p.page === page && opts.levers) p.open();
  }
}
export function anyOpen() { return panels.some(p => p.isOpen()); }
