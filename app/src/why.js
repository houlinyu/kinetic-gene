// "Why" panels: per-isoform ridge contributions + global fit-to-real-data card. SVG, no libraries.
// Everything is optional in the data (iso.ridge, meta.ridge); missing -> placeholder text.
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const f2 = (x, d = 2) => (Number.isFinite(x) ? Number(x).toFixed(d) : '—');
const UP = '#5ee0a0', DOWN = '#ff4d5e';
const trunc = (s, n) => (s.length > n ? s.slice(0, n - 1) + '…' : s);

// Horizontal diverging bars of the top-k |contribution| for one isoform.
export function contribSVG(ridge, k = 6) {
  const cs = (ridge?.contrib || []).filter(c => Number.isFinite(c?.contribution))
    .sort((a, b) => Math.abs(b.contribution) - Math.abs(a.contribution)).slice(0, k);
  if (!cs.length) return '<div class="mute small">no contributions</div>';
  const W = 300, rowH = 17, lab = 128, pad = 34, H = cs.length * rowH + 6;
  const max = Math.max(...cs.map(c => Math.abs(c.contribution)), 1e-9);
  const plotW = W - lab - pad * 2, x0 = lab + pad + plotW / 2, sc = (plotW / 2) / max;
  const rows = cs.map((c, i) => {
    const y = 3 + i * rowH, v = c.contribution, w = Math.max(Math.abs(v) * sc, 1.5), x = v >= 0 ? x0 : x0 - w;
    const name = c.label || c.name;
    const val = (v >= 0 ? '+' : '−') + f2(Math.abs(v));
    return `<g><title>${esc(name)}: feature value ${f2(c.value)} → ${val} to predicted MRL</title>
      <rect x="0" y="${y}" width="${W}" height="${rowH}" fill="transparent"/>
      <text x="${lab - 4}" y="${y + 12}" text-anchor="end" class="wl">${esc(trunc(name, 24))}</text>
      <rect x="${x}" y="${y + 3}" width="${w}" height="${rowH - 6}" rx="2" fill="${v >= 0 ? UP : DOWN}"/>
      <text x="${v >= 0 ? x + w + 3 : x - 3}" y="${y + 12}" text-anchor="${v >= 0 ? 'start' : 'end'}" class="wv">${val}</text></g>`;
  }).join('');
  return `<svg viewBox="0 0 ${W} ${H}" class="wsvg"><line x1="${x0}" x2="${x0}" y1="0" y2="${H}" class="wax"/>${rows}</svg>`;
}

// The MRL models (ridge here, Optimus 5′ for comparison) were fitted on 25–100 nt MPRA reporter
// UTRs. A mean ribosome load is a count of ribosomes, so a negative value is not a small number —
// it is out of the model's range. 18 of 82 isoforms in this data set predict below 0.
const MRL_TRAIN = [25, 100];
const MRL_OOR = 'the ridge / Optimus 5′ MRL models were fitted on 25–100 nt reporter UTRs and extrapolate badly to long endogenous leaders; a mean ribosome load cannot be negative, so no number is shown';
function mrlCell(iso, r) {
  const v = r?.mrl_pred;
  if (!Number.isFinite(v)) return '';
  const len = iso?.len;
  const far = Number.isFinite(len) && (len < MRL_TRAIN[0] || len > 3 * MRL_TRAIN[1]);
  if (v < 0 || far) {
    return `MRL <b class="oor" title="${esc(`raw prediction ${f2(v)}${Number.isFinite(len) ? `, leader ${len} nt vs ${MRL_TRAIN[0]}–${MRL_TRAIN[1]} nt training range` : ''} — ${MRL_OOR}`)}">out of model range</b>`;
  }
  return `MRL <b title="predicted mean ribosome load (clamped at 0; models trained on ${MRL_TRAIN[0]}–${MRL_TRAIN[1]} nt reporter UTRs)">${f2(Math.max(0, v))}</b>`;
}

export function whyHTML(isoA, isoB, nameA, nameB) {
  const any = isoA?.ridge || isoB?.ridge;
  const one = (tag, iso, name) => {
    const r = iso?.ridge;
    if (!r) return `<div class="wlane"><div class="whead"><span class="tag sm">${tag}</span> ${esc(name)} <span class="mute">· ridge model not fitted for this gene</span></div></div>`;
    return `<div class="wlane"><div class="whead"><span class="tag sm">${tag}</span> ${esc(name)}
      <span class="wpred">${mrlCell(iso, r)}${Number.isFinite(r.te_pred) ? ` · TE <b>${f2(r.te_pred)}</b>` : ''}</span></div>${contribSVG(r)}</div>`;
  };
  return `<div class="ftitle">What the data-fitted model says drives initiation here</div>
    ${any ? '' : '<div class="mute small">ridge model not fitted for this gene — it was fitted on the earlier gene set, so these leaders have no contributions.</div>'}
    ${one('A', isoA, nameA)}${one('B', isoB, nameB)}
    ${any ? `<div class="fkey"><span><i style="background:${UP}"></i>raises predicted ribosome load</span><span><i style="background:${DOWN}"></i>lowers it</span></div>
      <div class="oornote">Trained on ${MRL_TRAIN[0]}–${MRL_TRAIN[1]} nt MPRA reporter UTRs. Endogenous leaders are much longer, and the fit extrapolates badly —
        18 of 82 isoforms here predict a <b>negative</b> ribosome load, which is impossible; those are shown as “out of model range”.</div>` : ''}`;
}

// ---- global fit card ----
function modelsOf(meta) {
  const R = meta?.ridge; if (!R) return null;
  const M = R.models ?? R;
  return { mrl: M.mrl, te: M.te, physics: R.physics_fit ?? M.physics_fit ?? null };
}
const MCOL = { ridge: '#4fc3f7', optimus: '#8a8f9c', physics: '#ffb547' };

function r2Bars(ms) {
  const sets = [];
  if (ms.mrl) sets.push({ name: 'MRL reporter (MPRA)', sub: ms.mrl.dataset, n: ms.mrl.n, bars: [['ridge', ms.mrl.r2_cv], ['optimus', ms.mrl.r2_optimus_same_split]] });
  if (ms.te) sets.push({ name: 'endogenous TE', sub: ms.te.dataset, n: ms.te.n, bars: [['ridge', ms.te.r2_cv], ['optimus', ms.te.r2_optimus_same_split]] });
  const ph = ms.physics;
  const mp = ms.mrl?.r2_physics_fitted_same_split ?? ms.mrl?.r2_physics_same_split;
  if (Number.isFinite(mp)) sets[0]?.bars.push(['physics', mp]);
  else if (ph && Number.isFinite(ph.r2_cv)) {
    const tgt = /te|cao|endog/i.test(ph.dataset || '') && ms.te ? sets.find(s => s.name.includes('TE')) : sets[0];
    if (tgt) tgt.bars.push(['physics', ph.r2_cv]);
  }
  for (const s of sets) s.bars = s.bars.filter(b => Number.isFinite(b[1]));
  const vals = sets.flatMap(s => s.bars.map(b => b[1]));
  if (!vals.length) return '<div class="mute small">fitting…</div>';
  const vmax = Math.max(0.1, ...vals), vmin = Math.min(0, ...vals);
  const W = 300, lab = 64, padR = 34, rowH = 14;
  const X = (v) => lab + (v - vmin) / (vmax - vmin) * (W - lab - padR);
  let y = 0, body = '';
  for (const s of sets) {
    body += `<text x="0" y="${y + 11}" class="wl2">${esc(s.name)}${s.n ? ` · n=${s.n}` : ''}</text>`; y += 15;
    for (const [k, v] of s.bars) {
      const xa = X(Math.min(0, v)), xb = X(Math.max(0, v));
      body += `<g><title>${k}: cross-validated r² = ${f2(v, 3)}${s.sub ? ` (${esc(s.sub)})` : ''}</title>
        <text x="${lab - 4}" y="${y + 10}" text-anchor="end" class="wl">${k}</text>
        <rect x="${xa}" y="${y + 2}" width="${Math.max(xb - xa, 1.5)}" height="${rowH - 4}" rx="2" fill="${MCOL[k]}"/>
        <text x="${xb + 3}" y="${y + 10}" class="wv">${f2(v, 2)}</text></g>`;
      y += rowH;
    }
    y += 6;
  }
  const zx = X(0);
  return `<svg viewBox="0 0 ${W} ${y}" class="wsvg"><line x1="${zx}" x2="${zx}" y1="14" y2="${y - 4}" class="wax"/>${body}</svg>`;
}

function forest(model, title, k = 7) {
  const fs = (model?.features || []).filter(f => Number.isFinite(f?.coef))
    .sort((a, b) => Math.abs(b.coef) - Math.abs(a.coef)).slice(0, k);
  if (!fs.length) return '';
  const lo = Math.min(0, ...fs.map(f => Array.isArray(f.ci) && Number.isFinite(f.ci[0]) ? f.ci[0] : f.coef));
  const hi = Math.max(0, ...fs.map(f => Array.isArray(f.ci) && Number.isFinite(f.ci[1]) ? f.ci[1] : f.coef));
  const W = 300, lab = 140, rowH = 15, pad = 8, H = fs.length * rowH + 4;
  const X = (v) => lab + pad + (v - lo) / Math.max(hi - lo, 1e-9) * (W - lab - 2 * pad);
  const rows = fs.map((f, i) => {
    const y = 2 + i * rowH + rowH / 2, col = f.coef >= 0 ? UP : DOWN;
    const ci = Array.isArray(f.ci) && f.ci.every(Number.isFinite) ? f.ci : null;
    return `<g><title>${esc(f.label || f.name)}: coef ${f2(f.coef, 3)}${ci ? ` [${f2(ci[0], 3)}, ${f2(ci[1], 3)}]` : ''} (standardized)</title>
      <rect x="0" y="${y - rowH / 2}" width="${W}" height="${rowH}" fill="transparent"/>
      <text x="${lab - 2}" y="${y + 3.5}" text-anchor="end" class="wl">${esc(trunc(f.label || f.name, 25))}</text>
      ${ci ? `<line x1="${X(ci[0])}" x2="${X(ci[1])}" y1="${y}" y2="${y}" stroke="${col}" stroke-width="2" stroke-linecap="round" opacity=".6"/>` : ''}
      <circle cx="${X(f.coef)}" cy="${y}" r="4" fill="${col}" stroke="#14171f" stroke-width="2"/></g>`;
  }).join('');
  return `<div class="wsub">${esc(title)} <span class="mute">· top coefficients ± 95% CI</span></div>
    <svg viewBox="0 0 ${W} ${H}" class="wsvg"><line x1="${X(0)}" x2="${X(0)}" y1="0" y2="${H}" class="wax"/>${rows}</svg>`;
}

export function fitCardHTML(meta) {
  const ms = modelsOf(meta);
  const head = `<h2 class="fith">Fit to real data</h2>
    <p class="sub2">Cross-validated r² of a ridge model on interpretable UTR features, against the deep-learning baseline on the same split.</p>`;
  if (!ms || (!ms.mrl && !ms.te)) return head + '<div class="mute small">fitting…</div>';
  return `${head}${r2Bars(ms)}
    <div class="fkey"><span><i style="background:${MCOL.ridge}"></i>ridge (features)</span><span><i style="background:${MCOL.optimus}"></i>Optimus 5′ (CNN)</span>${ms.physics || Number.isFinite(ms.mrl?.r2_physics_same_split) ? `<span><i style="background:${MCOL.physics}"></i>scanning physics (this app's model)</span>` : ''}</div>
    ${forest(ms.mrl, 'MRL model')}${forest(ms.te, 'TE model')}`;
}
