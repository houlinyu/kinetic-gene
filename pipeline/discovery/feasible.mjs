// Feasible parameter region: grid over lambda = tcRate/scanSpeed (via tcRate at several scanSpeeds) and loadTCexp,
// everything else at params.js defaults. A point is feasible when every 'mechanism' finding passes
// (findings.score on data/utr_data.json). Writes pipeline/discovery/work/feasible.json.
import fs from 'node:fs';
import { values } from '../../app/src/params.js';
import { score } from '../../app/src/findings.js';

const root = new URL('../../', import.meta.url).pathname;
const data = JSON.parse(fs.readFileSync(root + 'data/utr_data.json', 'utf8'));
const base = values();
const pts = [];
for (const scanSpeed of [10, 20, 40])
  for (const tcRate of [0.1, 0.2, 0.3, 0.5, 0.75, 1, 1.5, 2, 3])
    for (const loadTCexp of [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5]) {
      const p = { ...base, scanSpeed, tcRate, loadTCexp };
      const s = score(data, p);
      pts.push({ scanSpeed, tcRate, loadTCexp, lambda: tcRate / scanSpeed, passed: s.passed, total: s.total,
                 failed: s.res.filter(f => f.status === 'mechanism' && !f.result.pass).map(f => f.id) });
    }
const d = score(data, base);
const out = { default: { ...base, passed: d.passed, total: d.total,
  failed: d.res.filter(f => f.status === 'mechanism' && !f.result.pass).map(f => f.id) },
  findings: d.res.map(f => ({ id: f.id, status: f.status, pass: f.result.pass, detail: f.result.detail })),
  grid: pts, feasible: pts.filter(q => q.passed === q.total) };
fs.mkdirSync(root + 'pipeline/discovery/work', { recursive: true });
fs.writeFileSync(root + 'pipeline/discovery/work/feasible.json', JSON.stringify(out, null, 1));
console.log('default', d.passed, '/', d.total, out.default.failed);
console.log('feasible', out.feasible.length, '/', pts.length);
for (const q of out.feasible) console.log(q.scanSpeed, q.tcRate, q.loadTCexp, q.lambda.toFixed(4));
