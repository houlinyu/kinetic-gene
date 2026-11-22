// Predicted isoform fractions for every gene in splice_data.json, from the exact kinetic model,
// plus two null baselines. Writes JSON to stdout.
//   node pipeline/isopedia/predict.mjs > data/figures/isopedia_pred.json
import { readFileSync } from 'node:fs';
import { expected, calibrate } from '../../app/src/splice/model_exact.js';
import { values } from '../../app/src/splice/params.js';

const data = JSON.parse(readFileSync(new URL('../../app/public/splice_data.json', import.meta.url)));
calibrate(data);
const p = values();
const out = {};
for (const gene of data.genes) {
  const e = expected(gene, p);
  // model
  const model = {};
  for (const i of e.isoforms) model[i.id] = i.value;
  // baseline 1: uniform over the gene's annotated isoforms
  const uni = {};
  for (const i of gene.isoforms) uni[i.id] = 1 / gene.isoforms.length;
  // baseline 2: promoter weights only (CAGE share of each TSS, split equally among the isoforms
  // that start at that TSS; no splicing kinetics at all). Falls back to uniform TSS weights.
  const tssW = new Map();
  let msum = 0;
  for (const t of gene.tss) if (Number.isFinite(t.weight) && t.weight >= 0) msum += t.weight;
  const useMeasured = msum > 0;
  for (const t of gene.tss) tssW.set(t.id, useMeasured ? (Number.isFinite(t.weight) ? t.weight : 0) / msum : 1 / gene.tss.length);
  const perTss = new Map();
  for (const i of gene.isoforms) perTss.set(i.tss_id, (perTss.get(i.tss_id) || 0) + 1);
  const prom = {};
  for (const i of gene.isoforms) prom[i.id] = (tssW.get(i.tss_id) || 0) / (perTss.get(i.tss_id) || 1);
  out[gene.gene] = {
    model, uniform: uni, promoter_only: prom,
    other: e.other,
    tss_mode: e.tssWeighting.mode, tss_matched: e.tssWeighting.matched, tss_total: e.tssWeighting.total,
  };
}
process.stdout.write(JSON.stringify(out, null, 1));
