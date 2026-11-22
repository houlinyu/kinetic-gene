// Reads isoform-schema JSON lines on stdin, writes "pMain\tflux\tperMin" per line using app/src/analytic.js
// expected(iso, values()) at the app's default parameters (or overrides from argv[2] JSON).
import { expected } from '../../app/src/analytic.js';
import { values } from '../../app/src/params.js';
import { readFileSync } from 'fs';
const p = { ...values(), ...(process.argv[2] ? JSON.parse(process.argv[2]) : {}) };
const lines = readFileSync(0, 'utf8').split('\n').filter(Boolean);
const out = [];
for (const l of lines) { const r = expected(JSON.parse(l), p); out.push(`${r.pMain}\t${r.flux}\t${r.perMin}`); }
process.stdout.write(out.join('\n') + '\n');
