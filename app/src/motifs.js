// Motif presentation: colours, chip text, per-lane summary, sequence strip, tooltip HTML.
// Motifs come precomputed from the pipeline (iso.motifs, see data/SPEC_v2.md); every field may be missing.
import { SOURCES } from './params.js';

export const MOTIF = {
  kozak_main: { color: '#e85d9c', name: 'main start Kozak', prio: 0 },
  uAUG:       { color: '#ffb547', name: 'upstream AUG', prio: 1 },
  TISU:       { color: '#b8f36b', name: 'TISU', prio: 2 },
  TOP:        { color: '#5ee0d0', name: '5′ TOP', prio: 3 },
  G4:         { color: '#a990ff', name: 'G-quadruplex', prio: 4 },
  uORF_stop:  { color: '#ff6b4a', name: 'uORF stop', prio: 5 },
  nearcog:    { color: '#ffe08a', name: 'near-cognate start', prio: 6 },
  PRTE:       { color: '#6fb7ff', name: 'PRTE', prio: 7 },
  polyU:      { color: '#f7a1c4', name: 'poly-U', prio: 8 },
};
const DEF = { color: '#c8c8c8', name: 'motif', prio: 9 };
export const motifStyle = (t) => MOTIF[t] ?? DEF;
const S1 = { strong: 'S', adequate: 'A', weak: 'W' };
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

export const motifsOf = (iso) => (Array.isArray(iso?.motifs) ? iso.motifs.filter(m => m && Number.isFinite(m.start)) : []);

export function chipText(m) {
  const s = m.strength ? S1[m.strength] ?? '' : '';
  switch (m.type) {
    case 'uAUG': return `uAUG${s ? '·' + s : ''}`;
    case 'kozak_main': return `Kozak${s ? '·' + s : ''}`;
    case 'uORF_stop': return 'stop';
    case 'nearcog': return (m.seq || '').toUpperCase().replace(/T/g, 'U').slice(0, 3) || 'near-cog';
    case 'polyU': return 'U≥6';
    default: return m.type;
  }
}
// visual weight by strength: strong > adequate > weak (unknown = adequate)
export const strengthScale = (m) => (m.strength === 'strong' ? 1 : m.strength === 'weak' ? 0.78 : 0.88);

// Which motifs get an in-scene chip (keep the scene readable): all but weak ones first, at most `max`.
export function chipMotifs(iso, max = 12, keep = null) {
  const ms = motifsOf(iso).filter(m => m.type !== 'uORF_stop' && (!keep || keep(m))); // stops are already drawn as uORF bars
  ms.sort((a, b) => (a.strength === 'weak') - (b.strength === 'weak') || motifStyle(a.type).prio - motifStyle(b.type).prio || a.start - b.start);
  return ms.slice(0, max).sort((a, b) => a.start - b.start);
}

export function summaryLine(iso) {
  const ms = motifsOf(iso);
  if (!Array.isArray(iso?.motifs)) return '';
  if (!ms.length) return '<div class="msum">no annotated motifs</div>';
  const by = {};
  for (const m of ms) (by[m.type] ??= []).push(m);
  const parts = [];
  const dot = (t) => `<i style="background:${motifStyle(t).color}"></i>`;
  if (by.uAUG) {
    const st = by.uAUG.filter(m => m.strength === 'strong').length;
    parts.push(`${dot('uAUG')}${by.uAUG.length} uAUG${st ? ` (${st} strong)` : ''}`);
  }
  for (const t of Object.keys(by)) {
    if (t === 'uAUG' || t === 'kozak_main') continue;
    const n = by[t].length;
    parts.push(`${dot(t)}${n > 1 ? n + ' ' : ''}${t === 'uORF_stop' ? 'uORF stop' + (n > 1 ? 's' : '') : t === 'nearcog' ? 'near-cognate' : t === 'polyU' ? 'poly-U' : t}`);
  }
  if (by.kozak_main) parts.push(`${dot('kozak_main')}Kozak ${by.kozak_main[0].strength ?? ''}`);
  return `<div class="msum">${parts.join(' · ')}</div>`;
}

export function sourceFor(key, data) {
  if (!key) return null;
  const m = (data?.meta?.sources || []).find(s => s.key === key);
  if (m) return { text: m.citation, url: m.url };
  const p = SOURCES[key];
  return p ? { text: p.text, url: p.url } : { text: key, url: null };
}

export function tipHTML(m, iso, data) {
  const st = motifStyle(m.type);
  const seq = esc(String(m.seq || '').toUpperCase().replace(/T/g, 'U')).replace(/AUG/g, '<b class="aug">AUG</b>');
  const up = iso?.len != null ? iso.len - m.start : null;
  const src = sourceFor(m.cite, data);
  return `<div class="tt-h"><i style="background:${st.color}"></i>${esc(m.label || st.name)}${m.strength ? ` <span class="mute">· ${esc(m.strength)}</span>` : ''}</div>
    ${seq ? `<div class="tt-seq">${seq}</div>` : ''}
    <div class="tt-pos">nt ${m.start}–${m.end} from cap${up != null && m.type !== 'kozak_main' ? ` · ${up} nt before the start codon` : ''}</div>
    ${m.note ? `<div class="tt-note">${esc(m.note)}</div>` : ''}
    ${src ? `<div class="tt-src">Source: ${src.url ? `<a href="${esc(src.url)}" target="_blank" rel="noopener">${esc(src.text)}</a>` : esc(src.text)}</div>` : ''}
    <div class="tt-hint">click chip to pin</div>`;
}

export function legendHTML() {
  return `<div class="mlegend"><div class="ftitle">Motifs</div>${Object.entries(MOTIF).map(([k, v]) =>
    `<span><i style="background:${v.color}"></i>${k === 'nearcog' ? 'near-cognate (CUG…)' : v.name}</span>`).join('')}
    <span class="mute">S/A/W = strong / adequate / weak Kozak context; bigger chip = stronger</span></div>`;
}

// Monospace strip: first 120 nt + last 30 nt, motif spans tinted in their colour.
export function stripHTML(iso) {
  const seq = String(iso?.utr5 || '');
  if (!seq) return '';
  const L = seq.length, ms = motifsOf(iso).filter(m => m.type !== 'kozak_main');
  const col = new Array(L).fill(null), pr = new Array(L).fill(99);
  for (const m of ms) {
    const s = motifStyle(m.type);
    for (let i = Math.max(0, m.start); i < Math.min(L, m.end); i++) if (s.prio < pr[i]) { pr[i] = s.prio; col[i] = s.color; }
  }
  const seg = (a, b) => {
    let out = '', cur = undefined, buf = '';
    const flush = () => { if (!buf) return; out += cur ? `<span style="background:${cur}40;box-shadow:inset 0 -2px ${cur}">${buf}</span>` : buf; buf = ''; };
    for (let i = a; i < b; i++) { if (col[i] !== cur) { flush(); cur = col[i]; } buf += seq[i]; }
    flush(); return out;
  };
  const body = L <= 150 ? seg(0, L) : `${seg(0, 120)}<span class="gap"> …${L - 150} nt… </span>${seg(L - 30, L)}`;
  const aug = String(iso.cds_context || '').slice(6, 9) || 'AUG';
  return `<details class="seqd"><summary>sequence (5′ cap → start)</summary><div class="strip">${body}<span class="startc">${esc(aug)}</span></div></details>`;
}
