// Mechanistic prediction of isoform fractions for every gene in data/splice_data.json,
// using the app's exact co-transcriptional splicing model at its default parameters.
// Writes TSV: gene<TAB>isoform_id<TAB>pred_fraction  (+ an "other" row per gene).
import { expected, calibrate } from '../../app/src/splice/model_exact.js';
import { values } from '../../app/src/splice/params.js';
import { readFileSync } from 'fs';
const data = JSON.parse(readFileSync(new URL('../../data/splice_data.json', import.meta.url), 'utf8'));
calibrate(data);
const p = { ...values(), ...(process.argv[2] ? JSON.parse(process.argv[2]) : {}) };
const out = ['gene\tisoform\tpred'];
for (const g of data.genes) {
  const r = expected(g, p);
  for (const i of r.isoforms) out.push(`${g.gene}\t${i.id}\t${i.value}`);
  out.push(`${g.gene}\tother\t${r.other}`);
}
process.stdout.write(out.join('\n') + '\n');
