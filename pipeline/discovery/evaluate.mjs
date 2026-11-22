// Run the exact expected-value model (app/src/analytic.js) on every discovery isoform,
// calm (tc=1) and stress (tc=STRESS_TC=0.15, same as findings.js), at the default params plus every
// feasible grid point (feasible.json), for both fold methods (global / local hairpins).
// Writes pipeline/discovery/work/model_out.json: per isoform {calm[fold][k], stress[fold][k]} with k = point index.
import fs from 'node:fs';
import { expected } from '../../app/src/analytic.js';
import { values } from '../../app/src/params.js';

const STRESS_TC = 0.15;
const root = new URL('../../', import.meta.url).pathname;
const W = root + 'pipeline/discovery/work/';
const { groups } = JSON.parse(fs.readFileSync(W + 'isoforms.json', 'utf8'));
const feas = JSON.parse(fs.readFileSync(W + 'feasible.json', 'utf8'));
const base = values();
const points = [{ label: 'default', ...base }, ...feas.feasible.map(q => ({
  label: `v${q.scanSpeed}_tc${q.tcRate}_b${q.loadTCexp}`, ...base, scanSpeed: q.scanSpeed, tcRate: q.tcRate, loadTCexp: q.loadTCexp }))];
// drop exact duplicates of the default
const seen = new Set(); const P = points.filter(p => { const k = `${p.scanSpeed}|${p.tcRate}|${p.loadTCexp}`; if (seen.has(k)) return false; seen.add(k); return true; });

const out = {};
for (const g of groups) for (const iso of g.isoforms) {
  const r = { calm: {}, stress: {}, lost_calm_default: {} };
  for (const fm of ['global', 'local']) {
    const m = { ...iso, hairpins: iso.fold[fm].hairpins };
    r.calm[fm] = P.map(p => expected(m, { ...p, tc: 1 }).perMin);
    r.stress[fm] = P.map(p => expected(m, { ...p, tc: STRESS_TC }).perMin);
    const e = expected(m, { ...P[0], tc: 1 });
    r.lost_calm_default[fm] = { pMain: e.pMain, flux: e.flux, Tblock: e.Tblock, ...e.lost };
    // ablations at every point: all upstream AUGs disabled / all hairpins melted
    const noU = { disabledU: new Set(iso.uorfs.map((_, i) => i)) };
    const noH = { melted: new Set(m.hairpins.map((_, i) => i)) };
    r.abl = r.abl || {};
    const noUH = { ...noU, ...noH };
    r.abl[fm] = {
      noUH_calm: P.map(p => expected(m, { ...p, tc: 1 }, noUH).perMin), noUH_stress: P.map(p => expected(m, { ...p, tc: STRESS_TC }, noUH).perMin),
      noU_calm: P.map(p => expected(m, { ...p, tc: 1 }, noU).perMin), noU_stress: P.map(p => expected(m, { ...p, tc: STRESS_TC }, noU).perMin),
      noH_calm: P.map(p => expected(m, { ...p, tc: 1 }, noH).perMin), noH_stress: P.map(p => expected(m, { ...p, tc: STRESS_TC }, noH).perMin) };
  }
  out[iso.id] = r;
}
fs.writeFileSync(W + 'model_out.json', JSON.stringify({ points: P.map(p => ({ label: p.label, scanSpeed: p.scanSpeed, tcRate: p.tcRate, loadTCexp: p.loadTCexp })), stress_tc: STRESS_TC, iso: out }));
console.log('points', P.length, 'isoforms', Object.keys(out).length);
