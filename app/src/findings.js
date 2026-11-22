// Published findings as constraints on ONE shared parameter set.
// Each finding is a direction check computed with the exact model (analytic.js), never tuned per finding.
// status: 'mechanism' = should be reproducible by scanning physics; 'missing' = needs a mechanism the model lacks
import { expected } from './analytic.js';

const STRESS_TC = 0.15;
const find = (data, gene, pred) => {
  const g = data.genes.find(x => x.gene === gene); if (!g) return null;
  return pred ? g.isoforms.find(pred) : g.isoforms.find(i => i.tags?.includes('MANE_Select'));
};
const out = (iso, p, opts) => expected(iso, p, opts).perMin;
const calm = p => ({ ...p, tc: 1 }), stress = p => ({ ...p, tc: STRESS_TC });
const uIdx = (iso, f) => new Set(iso.uorfs.map((u, i) => (f(u) ? i : -1)).filter(i => i >= 0));

export const FINDINGS = [
  {
    id: 'bulk', title: 'Stress lowers translation of most mRNAs', cite: 'eIF2α-P general repression',
    status: 'mechanism',
    test(data, p) {
      const r = data.genes.flatMap(g => g.isoforms).map(i => out(i, stress(p)) / Math.max(out(i, calm(p)), 1e-12)).filter(Number.isFinite).sort((a, b) => a - b);
      const med = r[Math.floor(r.length / 2)];
      return { pass: med < 0.6, detail: `median stress/calm = ${med.toFixed(2)}` };
    },
  },
  {
    id: 'atf4_up', title: 'ATF4 output rises under stress', cite: 'vattem2004', status: 'mechanism',
    test(data, p) {
      const iso = find(data, 'ATF4'); if (!iso) return null;
      const a = out(iso, calm(p)), b = out(iso, stress(p));
      return { pass: b > 1.5 * a, detail: `${a.toFixed(2)} → ${b.toFixed(2)} (×${(b / a).toFixed(1)})` };
    },
  },
  {
    id: 'atf4_u1', title: 'ATF4: removing uORF1 abolishes stress induction', cite: 'vattem2004', status: 'mechanism',
    test(data, p) {
      const iso = find(data, 'ATF4'); if (!iso) return null;
      const o = { disabledU: uIdx(iso, u => u.type === 'uORF' && u.aa > 1 && u.aa <= 5) };
      if (!o.disabledU.size) return null;
      const a = out(iso, calm(p), o), b = out(iso, stress(p), o);
      return { pass: b <= 1.2 * a, detail: `${a.toFixed(2)} → ${b.toFixed(2)}` };
    },
  },
  {
    id: 'atf4_u2', title: 'ATF4: removing uORF2 raises basal output', cite: 'vattem2004', status: 'mechanism',
    test(data, p) {
      const iso = find(data, 'ATF4'); if (!iso) return null;
      const o = { disabledU: uIdx(iso, u => u.type === 'overlapping') };
      const wt = out(iso, calm(p)), mut = out(iso, calm(p), o);
      return { pass: mut > 3 * wt, detail: `WT ${wt.toFixed(2)} vs mutant ${mut.toFixed(2)}` };
    },
  },
  {
    id: 'chop_up', title: 'CHOP (DDIT3) output rises under stress', cite: 'palam2011', status: 'missing',
    note: 'Needs ribosome stalling inside the CHOP uORF (not in the model). A failure here points at that mechanism.',
    test(data, p) {
      const iso = find(data, 'DDIT3'); if (!iso) return null;
      const a = out(iso, calm(p)), b = out(iso, stress(p));
      return { pass: b > 1.5 * a, detail: `${a.toFixed(2)} → ${b.toFixed(2)}` };
    },
  },
];

// Pairwise isoform findings are appended at runtime from data.findings (written by the pipeline), shape:
// { id, title, cite, gene, better: transcript id, worse: transcript id }
export function allFindings(data) {
  const extra = (data.findings || []).map(f => ({
    id: f.id, title: f.title, cite: f.cite, status: 'mechanism',
    test(d, p) {
      const g = d.genes.find(x => x.gene === f.gene); if (!g) return null;
      const by = id => g.isoforms.find(i => i.id.split('.')[0] === id.split('.')[0]);
      const hi = by(f.better), lo = by(f.worse); if (!hi || !lo) return null;
      const a = out(hi, calm(p)), b = out(lo, calm(p));
      return { pass: a > b, detail: `${hi.name} ${a.toFixed(2)} vs ${lo.name} ${b.toFixed(2)}` };
    },
  }));
  return [...FINDINGS, ...extra];
}

export function score(data, p) {
  const res = allFindings(data).map(f => ({ ...f, result: f.test(data, p) })).filter(f => f.result);
  const mech = res.filter(f => f.status === 'mechanism');
  return { res, passed: mech.filter(f => f.result.pass).length, total: mech.length };
}
