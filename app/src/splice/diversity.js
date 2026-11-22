// Isoform DIVERSITY: how one gene makes many transcripts.
//
// The point this module exists to make visible: nothing here is a per-gene "isoform mix" being
// dialled. Each individual polymerase picks ONE promoter and then resolves each competing junction
// ONE way, independently of every other polymerase, and the isoform it becomes is that list of
// choices. Promoter choice × splice-site choices is therefore combinatorial, and most introns are
// not choices at all.
//
// Everything displayed as a PROBABILITY here comes from model_exact.expected():
//   - promoter branch weights: the same normalisation expected() applies (measured CAGE share of
//     the annotated TSSs, else uniform) — recomputed here from the same gene.tss[].weight fields,
//     and labelled with expected()'s own tssWeighting.note;
//   - junction branch labels: expected().junctionUse (marginal usage of that junction);
//   - leaf values: expected().fractions.
// The bead simulation (model.js) is never used for a displayed probability — only for what the
// live transcripts on screen are doing right now.
//
// Owner: ISOFORM-DIVERSITY. Pure functions + HTML/SVG builders; no DOM state of its own except the
// lazily-loaded translation table.

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const pct = (x) => (Number.isFinite(x) ? (100 * x).toFixed(x > 0 && x < 0.005 ? 2 : 0) + '%' : '—');
const pct1 = (x) => (Number.isFinite(x) ? (100 * x).toFixed(1) + '%' : '—');
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

// A junction set is identified exactly as model.js does it: the sorted list of donor_acceptor keys.
const keyOf = (j) => j.key ?? `${j.donorPos ?? j.start}_${j.acceptorPos ?? j.end}`;
const sigOf = (J, idxs) => idxs.map(i => keyOf(J[i])).sort().join('|');

// Five-colour identity palette (the app's existing palette; see style.css :root). Cream is left
// out on purpose: it is the exon fill, so an RNA drawn in it would read as part of the gene.
export const IDCOL = ['#e85d9c', '#4fc3f7', '#ffb547', '#5ee0a0', '#c39bff'];
export function idColour(sig) {
  if (!sig) return IDCOL[0];
  let h = 2166136261;
  for (let i = 0; i < sig.length; i++) { h ^= sig.charCodeAt(i); h = Math.imul(h, 16777619); }
  return IDCOL[Math.abs(h) % IDCOL.length];
}

// ---------------------------------------------------------------- geometry of the choices
// conflicts: overlapping introns cannot both survive, so they compete (model.js's rule). Recomputed
// here only when the junction list we were handed has no conflicts array (fallback path).
function conflictsOf(J) {
  if (J.every(j => Array.isArray(j.conflicts))) return J.map(j => j.conflicts);
  const s = (j) => j.donorPos ?? j.start, e = (j) => j.acceptorPos ?? j.end;
  return J.map((j, i) => J.map((o, k) => k).filter(k => k !== i && s(j) < e(J[k]) && s(J[k]) < e(j)));
}

// Maximal independent sets of one competing component = the outcomes that component can have
// (every junction that is still available and unconflicted commits with probability 1 in this
// model, so a final state is always maximal). Capped so a pathological component cannot hang the UI.
const MIS_CAP = 17;
function maximalSets(grp, conf) {
  if (grp.length > MIS_CAP) return { alts: [[grp[0]]], capped: true };
  const set = new Set(grp);
  const adj = new Map(grp.map(i => [i, conf[i].filter(c => set.has(c))]));
  const all = [];
  const go = (i, cur) => {
    if (i >= grp.length) { all.push(cur.slice()); return; }
    const v = grp[i];
    go(i + 1, cur);
    if (!cur.some(u => adj.get(v).includes(u))) { cur.push(v); go(i + 1, cur); cur.pop(); }
  };
  go(0, []);
  const alts = all.filter(s => { const t = new Set(s); return grp.every(v => t.has(v) || adj.get(v).some(u => t.has(u))); });
  return { alts: alts.length ? alts : [[]], capped: false };
}

// connected components of the conflict graph among `avail`
function components(avail, conf) {
  const set = new Set(avail), seen = new Set(), out = [];
  for (const i of avail) {
    if (seen.has(i)) continue;
    const st = [i], grp = [];
    while (st.length) {
      const v = st.pop();
      if (seen.has(v)) continue;
      seen.add(v); grp.push(v);
      for (const c of conf[v]) if (set.has(c) && !seen.has(c)) st.push(c);
    }
    out.push(grp.sort((a, b) => a - b));
  }
  return out;
}

// A junction is readable as what it does to the exons: "−E3" removes E3, "+E3" splices straight
// into E3 (nothing skipped). Falls back to the junction id.
export function junctionLabel(gene, J, i) {
  const j = J[i], a = j.donorPos ?? j.start, b = j.acceptorPos ?? j.end;
  const ex = (gene?.exons ?? []);
  const inside = ex.filter(e => e.start >= a && e.end <= b);
  if (inside.length) return '−' + inside.map(e => e.id).join(',');
  const to = ex.find(e => e.start === b);
  return to ? '+' + to.id : (j.id ?? `${a}→${b}`);
}
const altLabel = (gene, J, alt) => (alt.length ? alt.map(i => junctionLabel(gene, J, i)).join(' ') : 'none of these');

/**
 * The combinatorial space of one gene, as the model can actually produce it.
 * tssL = model.tssList(gene) (or gene.tss); exp = model_exact.expected(gene, P).
 */
export function space(gene, J, tssL, exp) {
  const conf = conflictsOf(J);
  const tss = (tssL ?? gene?.tss ?? []).filter(t => Number.isFinite(t.pos));
  // promoter weights: exactly expected()'s normalisation (measured CAGE share, else uniform)
  const measured = tss.filter(t => Number.isFinite(t.weight) && t.weight >= 0);
  const msum = measured.reduce((a, t) => a + t.weight, 0);
  const useMeasured = measured.length > 0 && msum > 0;
  const promoters = tss.map(t => ({
    id: t.id ?? `tss@${t.pos}`, pos: t.pos, support: t.support ?? null,
    w: useMeasured ? (Number.isFinite(t.weight) ? t.weight : 0) / msum : 1 / (tss.length || 1),
  }));
  const used = promoters.filter(p => p.w > 0);

  const use = exp?.junctionUse ?? {};
  const sigs = new Set();
  let capped = false;
  const byPos = new Map();
  for (const p of (used.length ? used : promoters)) {
    if (byPos.has(p.pos)) continue;
    const avail = J.map((j, i) => i).filter(i => (J[i].donorPos ?? J[i].start) >= p.pos);
    const comps = components(avail, conf);
    const solo = comps.filter(c => c.length === 1).map(c => c[0]);
    const groups = [];
    for (const c of comps.filter(c => c.length > 1)) {
      const m = maximalSets(c, conf);
      if (m.capped) capped = true;
      const uses = c.map(i => use[J[i].id]).filter(Number.isFinite);
      groups.push({
        members: c, alts: m.alts, capped: m.capped,
        at: Math.min(...c.map(i => J[i].donorPos ?? J[i].start)),
        // "near-deterministic": one alternative takes essentially all the usage
        maxUse: uses.length ? Math.max(...uses) : null,
      });
    }
    groups.sort((a, b) => a.at - b.at);
    // distinct junction sets this promoter can produce (the cross product of the groups' outcomes)
    let states = [solo];
    for (const g of groups) {
      const next = [];
      for (const s of states) for (const alt of g.alts) next.push(s.concat(alt));
      states = next;
      if (states.length > 4000) { capped = true; states = states.slice(0, 4000); }
    }
    for (const s of states) sigs.add(sigOf(J, s));
    byPos.set(p.pos, { groups, solo, combos: states.length });
  }

  // annotated transcripts, and which of them this model can reach at all
  const annSig = new Map();
  const exMap = new Map((gene?.exons ?? []).map(e => [e.id, e]));
  for (const iso of gene?.isoforms ?? []) {
    const chain = (iso.exon_ids ?? []).map(id => exMap.get(id)).filter(Boolean).sort((a, b) => a.start - b.start);
    const s = [];
    for (let k = 0; k + 1 < chain.length; k++) s.push(`${chain[k].end}_${chain[k + 1].start}`);
    annSig.set(iso.id, s.sort().join('|'));
  }
  const reachable = [...annSig.values()].filter(s => sigs.has(s)).length;

  // which junctions are a choice, which are uncontested, and which are not on any transcript the
  // model makes at all (their 5′ss sits upstream of every promoter it uses)
  const competing = new Set(), available = new Set();
  const distinct = new Map();       // one entry per distinct competing group, however many promoters see it
  for (const v of byPos.values()) {
    for (const i of v.solo) available.add(i);
    for (const g of v.groups) {
      for (const i of g.members) { competing.add(i); available.add(i); }
      const k = g.members.join(',');
      if (!distinct.has(k)) distinct.set(k, g);
    }
  }
  const groups = [...distinct.values()].sort((a, b) => a.at - b.at);
  // alternatives of each group that some annotated transcript actually uses (checked gene-wide, not
  // per promoter: this model names a transcript by its junction set, which several promoters share)
  const altUsed = new Map();
  for (const g of groups) {
    const s = new Set();
    for (const sig of annSig.values()) {
      const have = new Set(sig.split('|').filter(Boolean));
      s.add(g.members.filter(i => have.has(keyOf(J[i]))).map(i => keyOf(J[i])).sort().join('|'));
    }
    altUsed.set(g.members.join(','), s);
  }

  return {
    promoters, used, byPos, useMeasured, groups, altUsed,
    tssNote: exp?.tssWeighting?.note ?? null,
    tssMode: exp?.tssWeighting?.mode ?? null,
    combos: sigs.size, capped,
    competingJ: competing.size,
    constitutiveJ: available.size - competing.size,
    offPromoterJ: J.length - available.size,
    totalJ: J.length,
    groupCount: groups.length,
    realGroups: groups.filter(g => !(g.maxUse >= 0.95)).length,
    annotated: (gene?.isoforms ?? []).length, reachable,
    other: exp?.fractions?.other ?? exp?.other ?? 0,
    annSig,
  };
}

// ---------------------------------------------------------------- the live per-transcript identity
/**
 * What ONE polymerase has decided so far: its promoter, plus one entry per competing group it can
 * still see — resolved to an alternative, or still open ('?'). This is the transcript's identity as
 * it is being written; it is read from the stochastic simulation, not from expected().
 */
export function polSignature(gene, J, sp, pol) {
  // byPos holds one entry per promoter the model actually uses; a transcript from a zero-weight
  // promoter (possible only in the bead sim) is read against the nearest promoter at or upstream.
  let part = sp.byPos.get(pol.tssPos);
  if (!part) {
    let best = null;
    for (const [pos, v] of sp.byPos) if (pos <= pol.tssPos && (best == null || pos > best)) { best = pos; part = v; }
    part = part ?? [...sp.byPos.values()][0];
  }
  const prom = sp.promoters.find(p => p.pos === pol.tssPos) ?? null;
  const com = new Set(pol.committed ?? []), dead = new Set(pol.skipped ?? []);
  const idOf = (i) => J[i].id;
  const steps = (part?.groups ?? []).map(g => {
    const won = g.members.filter(i => com.has(idOf(i)));
    const lost = g.members.filter(i => dead.has(idOf(i)));
    if (won.length) {
      // resolved as far as this group goes once no member is still undecided
      const open = g.members.filter(i => !com.has(idOf(i)) && !dead.has(idOf(i)));
      return { group: g, state: open.length ? 'partial' : 'done', label: altLabel(gene, J, won), won };
    }
    if (lost.length === g.members.length) return { group: g, state: 'none', label: 'none of these', won: [] };
    return { group: g, state: 'open', label: '?', won: [] };
  });
  const sig = (prom?.id ?? '?') + '|' + steps.map(s => (s.state === 'open' ? '?' : s.label)).join('·');
  const open = steps.filter(s => s.state === 'open' || s.state === 'partial').length;
  return {
    promoter: prom, steps, sig, open,
    chip: [prom ? prom.id : 'TSS?', ...steps.map(s => s.label)].join(' · '),
    colour: idColour(sig),
    resolved: open === 0,
  };
}

/**
 * Divergence watch: two transcripts that are alive at the same time and resolved the SAME competing
 * junction group DIFFERENTLY. Call once per animation tick with the live signatures; it returns the
 * new divergences only (so the caller can pulse them once).
 */
export class DivergenceWatch {
  constructor() { this.prev = new Map(); this.log = []; this.count = 0; }
  reset() { this.prev = new Map(); this.log = []; this.count = 0; }
  tick(sigs, t = 0) {
    const fresh = [];
    const now = new Map();
    for (const s of sigs) {
      const m = new Map();
      s.sig0 = s.sig;
      for (const st of s.steps) if (st.state === 'done' || st.state === 'none') m.set(st.group.at, st.label);
      now.set(s.polId, m);
    }
    for (const [id, m] of now) {
      const was = this.prev.get(id) ?? new Map();
      for (const [at, label] of m) {
        if (was.get(at) === label) continue;         // not a new decision
        for (const [other, om] of now) {
          if (other === id) continue;
          const ol = om.get(at);
          if (ol && ol !== label) {
            fresh.push({ at, a: label, b: ol, t, pols: [id, other] });
            break;
          }
        }
      }
    }
    this.prev = now;
    for (const f of fresh) {
      this.count++;
      this.log.unshift(f);
    }
    if (this.log.length > 4) this.log.length = 4;
    return fresh;
  }
}

// ---------------------------------------------------------------- translation join (② page model)
// protein per mRNA for one isoform, from the SAME scanning model the ② Translation page runs
// (src/analytic.js on the real ViennaRNA-folded 5′ leader in utr_data.json, default parameters).
let TR = null, TRfail = null;
export async function translation() {
  if (TR || TRfail) return TR;
  try {
    const [ana, par, res] = await Promise.all([
      import('../analytic.js'), import('../params.js'), fetch('./utr_data.json'),
    ]);
    if (!res.ok) throw new Error('utr_data.json ' + res.status);
    const d = await res.json();
    const p = typeof par.values === 'function' ? par.values() : {};
    const byGene = new Map();
    for (const g of d.genes ?? []) {
      const m = new Map();
      for (const iso of g.isoforms ?? []) {
        try {
          const r = ana.expected(iso, p, {});
          m.set(iso.name, { perMin: r.perMin, pMain: r.pMain, flux: r.flux, id: iso.id });
        } catch { /* one bad isoform must not take the table down */ }
      }
      byGene.set(g.gene, m);
    }
    TR = { byGene, release: d.meta?.release ?? null };
  } catch (e) {
    TRfail = String(e);
    console.warn('diversity: translation join unavailable', e);
  }
  return TR;
}
export const translationReady = () => !!TR;
export const translationFail = () => TRfail;

/** per-isoform protein-per-mRNA chip; explicit when there is no entry (never silently blank). */
export function trEntry(geneSym, isoName) {
  if (!TR) return null;
  return TR.byGene.get(geneSym)?.get(isoName) ?? null;
}
export function trChip(geneSym, isoName, frac, share) {
  if (!TR) return '';
  const e = trEntry(geneSym, isoName);
  if (!e) return `<span class="dvtr none" title="utr_data.json carries no 5′ leader for this isoform, so the ② Translation page's scanning model cannot be run on it">no 5′ leader in <code>utr_data.json</code> — protein per mRNA not computable</span>`;
  const out = Number.isFinite(frac) && Number.isFinite(e.perMin) ? frac * e.perMin : null;
  return `<span class="dvtr" title="② Translation page model (src/analytic.js) on this isoform's real folded 5′ leader, default parameters, no stress: ${e.perMin.toFixed(2)} protein/min per mRNA">
    ${e.perMin.toFixed(2)} protein/min per mRNA${out == null ? '' : ` · abundance × that = <b>${out.toFixed(2)}</b>${Number.isFinite(share) ? ` (${pct(share)} of this gene's protein output)` : ''}`}</span>`;
}

// ---------------------------------------------------------------- the fan-out tree (HTML + CSS elbows)
/**
 * treeHTML(gene, J, sp, exp, opts) — promoter choice branching into splice-site choices, leaves =
 * the annotated isoforms plus the unannotated-combination bucket. Branch labels are expected()'s
 * numbers; leaf values are expected().fractions.
 */
export function treeHTML(gene, J, sp, exp, opts = {}) {
  const simple = !!opts.simple, onTranslate = !!opts.onTranslate, geneSym = gene?.gene;
  const frac = (id) => exp?.fractions?.[id] ?? (exp?.isoforms ?? []).find(i => i.id === id)?.value ?? 0;
  const use = exp?.junctionUse ?? {};
  // protein-output shares across the isoforms that HAVE a translation entry (said so explicitly)
  let shareDen = 0;
  for (const iso of gene?.isoforms ?? []) {
    const e = trEntry(geneSym, iso.name);
    if (e) shareDen += frac(iso.id) * e.perMin;
  }
  const shareOf = (iso) => {
    const e = trEntry(geneSym, iso.name);
    return e && shareDen > 0 ? frac(iso.id) * e.perMin / shareDen : null;
  };

  // group every annotated isoform by its promoter, then by its choice at each competing group
  const isoByProm = new Map();
  for (const iso of gene?.isoforms ?? []) {
    const p = sp.promoters.find(x => x.id === iso.tss_id) ?? null;
    const k = p ? p.pos : '?';
    if (!isoByProm.has(k)) isoByProm.set(k, []);
    isoByProm.get(k).push(iso);
  }
  const sigSet = (iso) => new Set((sp.annSig.get(iso.id) ?? '').split('|').filter(Boolean));

  const leafHTML = (iso) => {
    const v = frac(iso.id), zero = !(v > 0.0005);
    const sh = shareOf(iso);
    return `<li class="dvleaf${zero ? ' zero' : ''}"><div class="dvnode leaf" style="--idc:${idColour(iso.id)}">
      <span class="dvname" title="${esc(iso.name)}">${esc(iso.name)}</span>
      <span class="dvbar"><i style="width:${(100 * clamp(v, 0, 1)).toFixed(1)}%"></i></span>
      <span class="dvpct">${pct1(v)}</span>
      ${onTranslate ? `<button type="button" class="sxl" data-iso="${esc(iso.name)}" title="run the ribosome-scanning model on this isoform's 5′ leader (② Translation page)">→ translate this isoform</button>` : ''}
      ${measChip(geneSym, iso.id)}
      ${trChip(geneSym, iso.name, v, sh)}
      ${zero ? `<span class="dvzero">0%: every junction it uses is a strict subset of another isoform’s set, so this model can never reach it. It may still be a real transcript.</span>` : ''}
    </div></li>`;
  };

  // one promoter branch
  const promHTML = (p, isos) => {
    const part = p ? sp.byPos.get(p.pos) : null;
    const groups = part?.groups ?? [];
    const w = p?.w ?? 0;
    // recursive split of `isos` by their choice at groups[gi]
    const split = (list, gi) => {
      if (gi >= groups.length || !list.length) return list.map(leafHTML).join('');
      const g = groups[gi];
      const buckets = new Map();
      for (const iso of list) {
        const s = sigSet(iso);
        const won = g.members.filter(i => s.has(keyOf(J[i])));
        const exact = g.alts.find(a => a.length === won.length && a.every(i => won.includes(i)));
        const key = won.map(i => keyOf(J[i])).sort().join('|') || '∅';
        if (!buckets.has(key)) buckets.set(key, { won, exact: !!exact, isos: [] });
        buckets.get(key).isos.push(iso);
      }
      // Alternatives of this group that NO annotated transcript of this gene uses — checked
      // gene-wide (sp.altUsed), not against the transcripts in this branch, because several
      // promoters reach the same junction set and the annotation files it under one of them.
      const used = sp.altUsed.get(g.members.join(',')) ?? new Set();
      const unused = g.alts.filter(a => !used.has(a.map(i => keyOf(J[i])).sort().join('|')));
      const branch = (b) => {
        const u = b.won.length === 1 ? use[J[b.won[0]].id] : (b.won.length ? Math.max(...b.won.map(i => use[J[i].id] ?? 0)) : null);
        const det = Number.isFinite(u) && (u >= 0.95);
        return `<li class="dvbranch"><div class="dvnode ${b.won.length ? 'choice' : 'nochoice'}${det ? ' det' : ''}">
          <span class="dvsite">${esc(altLabel(gene, J, b.won))}</span>
          ${b.won.length ? `<span class="dvuse" title="expected().junctionUse — the fraction of transcripts of this gene that use this junction">${pct(u)} of transcripts use this junction</span>`
            : `<span class="dvuse warnc">no junction here — unreachable: every available junction with no surviving competitor commits with probability 1</span>`}
          ${!b.exact && b.won.length ? '<span class="dvuse warnc">a strict subset of a competitor’s set</span>' : ''}
        </div><ul>${split(b.isos, gi + 1)}</ul></li>`;
      };
      const extra = unused.length && !simple
        ? `<li class="dvbranch other"><div class="dvnode choice art"><span class="dvsite">${unused.length} further outcome${unused.length > 1 ? 's' : ''} of this junction group</span>
             <span class="dvuse">${unused.map(a => esc(altLabel(gene, J, a))).join(' / ')} — no annotated transcript of this gene uses ${unused.length > 1 ? 'these' : 'this'}, so ${unused.length > 1 ? 'they' : 'it'} land${unused.length > 1 ? '' : 's'} in the unannotated bucket below</span></div></li>` : '';
      return [...buckets.values()].sort((a, b2) => b2.isos.length - a.isos.length).map(branch).join('') + extra;
    };
    const real = groups.filter(g => !(g.maxUse >= 0.95)).length;
    const detr = groups.length - real;
    const body = part
      ? `${groups.length ? `<span class="dvuse">${groups.length} competing junction group${groups.length > 1 ? 's' : ''}${detr ? ` (${detr} near-deterministic)` : ''} · ${part.combos} junction-set combination${part.combos > 1 ? 's' : ''}</span>`
        : '<span class="dvuse">no competing junctions downstream of this promoter — one outcome</span>'}
         ${part.solo.length ? `<span class="dvconst">${part.solo.length} junction${part.solo.length > 1 ? 's' : ''} with no competitor: used by every transcript from this promoter (p = 1) — not a choice</span>` : ''}`
      : `<span class="dvuse warnc">no transcript starts here in this model, so its splice-site choices are not drawn. The transcript${isos.length > 1 ? 's' : ''} below ${isos.length > 1 ? 'are' : 'is'} filed under this promoter by the annotation; whatever percentage ${isos.length > 1 ? 'they' : 'it'} carr${isos.length > 1 ? 'y' : 'ies'} comes from the promoters above, which reach the same junction set.</span>`;
    return `<li class="dvprom"><div class="dvnode prom${w > 0 ? '' : ' zero'}">
      <span class="dvsite">promoter ${esc(p?.id ?? '—')}</span>
      <span class="dvuse">${w > 0 ? `${pct(w)} of transcripts start here` : 'weight 0 in this model — no transcript starts here'}</span>
      ${body}
    </div><ul>${split(isos, 0)}</ul></li>`;
  };

  const order = [...sp.promoters].sort((a, b) => b.w - a.w || a.pos - b.pos);
  const proms = order.filter(p => p.w > 0);
  // every promoter the model never fires is collapsed into one branch: repeating the same caveat
  // per promoter buries the tree, but the caveat itself is kept in full, once.
  const dead = order.filter(p => !(p.w > 0) && isoByProm.has(p.pos));
  const deadIsos = dead.flatMap(p => isoByProm.get(p.pos) ?? []);
  const orphan = isoByProm.get('?') ?? [];
  const other = sp.other;
  return `<ul class="dvtree">
    <li class="dvroot"><div class="dvnode root"><span class="dvsite">every transcript of ${esc(gene?.gene ?? '—')}</span>
      <span class="dvuse">each polymerase makes these choices on its own — promoter first, then one per competing junction</span>
      <span class="dvuse warnc">This model names a finished transcript by its <b>junction set alone</b>, so two promoters that lead to the same junction set land on the same leaf: a leaf’s % is <b>summed over every promoter</b> that reaches it, not just the branch it is drawn under.</span></div>
      <ul>${proms.map(p => promHTML(p, isoByProm.get(p.pos) ?? [])).join('')}
      ${dead.length ? `<li class="dvprom"><div class="dvnode prom zero"><span class="dvsite">${dead.length} promoter${dead.length > 1 ? 's' : ''} with weight 0 (${dead.map(p => esc(p.id)).join(', ')})</span>
        <span class="dvuse warnc">no transcript starts at ${dead.length > 1 ? 'these' : 'this one'} in this model, so ${dead.length > 1 ? 'their' : 'its'} splice-site choices are not drawn.
          The transcript${deadIsos.length > 1 ? 's' : ''} below ${deadIsos.length > 1 ? 'are' : 'is'} filed under ${dead.length > 1 ? 'them' : 'it'} by the annotation; whatever percentage ${deadIsos.length > 1 ? 'they' : 'it'} carr${deadIsos.length > 1 ? 'y' : 'ies'} comes from the promoters above, which reach the same junction set.</span></div>
        <ul>${deadIsos.map(leafHTML).join('')}</ul></li>` : ''}
      ${orphan.length ? `<li class="dvprom"><div class="dvnode prom zero"><span class="dvsite">promoter not in the TSS list</span>
        <span class="dvuse">these annotated transcripts name a TSS that <code>splice_data.json</code> does not list, so no promoter branch can be drawn for them</span></div>
        <ul>${orphan.map(leafHTML).join('')}</ul></li>` : ''}
      ${other > 0.0005 ? `<li class="dvbranch other"><div class="dvnode leaf art">
        <span class="dvname">unannotated junction combinations (model artifact)</span>
        <span class="dvbar"><i style="width:${(100 * clamp(other, 0, 1)).toFixed(1)}%"></i></span>
        <span class="dvpct">${pct1(other)}</span>
        <span class="dvzero">junction combinations that match no annotated transcript — the model has no intron-retention state</span></div></li>` : ''}
      </ul></li></ul>`;
}

// ---------------------------------------------------------------- "count the space"
export function countHTML(gene, J, sp, exp) {
  const nz = Object.values(exp?.fractions ?? {}).filter(v => v > 0.0005).length;
  const prom = sp.used.length;
  const perProm = [...sp.byPos.entries()].map(([pos, v]) => {
    const p = sp.promoters.find(x => x.pos === pos);
    return `${esc(p?.id ?? pos)}: ${v.groups.length ? v.groups.map(g => g.alts.length).join(' × ') : '1'} = ${v.combos}`;
  });
  return `<div class="dvcount">
    <div class="dvrow"><b>${sp.competingJ}</b> of ${sp.totalJ} junctions are genuinely competing (in ${sp.groupCount} group${sp.groupCount === 1 ? '' : 's'}${sp.realGroups < sp.groupCount ? `, ${sp.groupCount - sp.realGroups} of them near-deterministic` : ''});
      <b>${sp.constitutiveJ}</b> have no competitor at all, so every transcript that reaches them uses them
      (probability 1 in this model) — <b>most introns are not choices</b>.${sp.offPromoterJ ? ` ${sp.offPromoterJ} more sit upstream of every promoter this model uses, so no transcript here contains ${sp.offPromoterJ > 1 ? 'them' : 'it'}.` : ''}</div>
    <div class="dvrow"><b>${prom}</b> promoter${prom === 1 ? '' : 's'} with non-zero weight × the junction groups
      (${perProm.join(' · ')}) ⇒ <b>${sp.combos}</b> distinct junction-set combination${sp.combos === 1 ? '' : 's'} this model can produce${sp.capped ? ' (enumeration capped — at least this many)' : ''},
      against <b>${sp.annotated}</b> annotated transcript${sp.annotated === 1 ? '' : 's'} in the data,
      of which <b>${sp.reachable}</b> ${sp.reachable === 1 ? 'is' : 'are'} reachable here and <b>${nz}</b> outcome${nz === 1 ? '' : 's'} carry more than 0.05% of the output.</div>
    ${sp.annotated > sp.reachable ? `<div class="dvrow mute">${sp.annotated - sp.reachable} annotated transcript${sp.annotated - sp.reachable === 1 ? '' : 's'} sit at <b>0%</b> because every junction they use is a strict subset of another isoform’s set — <b>they may be real but are unreachable in this model</b>, not predicted to be absent.</div>` : ''}
    ${sp.tssNote ? `<div class="dvrow mute">Promoter weights: ${esc(sp.tssNote)}.</div>` : ''}
  </div>`;
}

export { esc as dvEsc, pct as dvPct };

// ---------------------------------------------------------------- measured isoform usage (Isopedia)
// app/public/isoform_usage.json: per-isoform share of ALL reads Isopedia assigned to the gene, in
// 821 public long-read samples. Key = UNVERSIONED Ensembl transcript id, so the app's versioned
// isoforms[].id must be split on '.'. Two app transcripts that share an intron chain share one
// measured value (the join is structural, not by id).
let US = null, USfail = null;
export async function usage() {
  if (US || USfail) return US;
  try {
    const r = await fetch('./isoform_usage.json');
    if (!r.ok) throw new Error('isoform_usage.json ' + r.status);
    US = await r.json();
  } catch (e) { USfail = String(e); console.warn('diversity: measured usage unavailable', e); }
  return US;
}
export const usageReady = () => !!US;
export const usageMeta = () => US?.meta ?? null;
export const usageGene = (geneSym) => US?.genes?.[geneSym] ?? null;

/** measured usage of one annotated isoform, or why there is none. Never returns 0 for "missing". */
export function measuredOf(geneSym, isoId) {
  if (!US) return { kind: 'loading' };
  const g = usageGene(geneSym);
  if (!g) return { kind: 'nogene' };
  if (g.measured === false) return { kind: 'nopanel', note: g.note ?? null };
  const base = String(isoId ?? '').split('.')[0];
  const v = g.overall?.[base];
  if (!Number.isFinite(v)) return { kind: 'below', n: g.n_samples };
  // a value shared by more than one of this app's isoforms means they share an intron chain
  const shared = Object.values(g.overall).filter(x => x === v).length > 1;
  return { kind: 'value', value: v, median: g.overall_median?.[base], n: g.n_samples, shared };
}

/** the measured bar that sits beside the predicted one, on a row or a leaf. */
export function measChip(geneSym, isoId) {
  const m = measuredOf(geneSym, isoId);
  if (m.kind === 'loading' || m.kind === 'nogene') return '';
  if (m.kind === 'nopanel') return `<span class="dvmeas none" title="${esc(m.note ?? '')}">measured: <b>not in the Isopedia panel</b></span>`;
  if (m.kind === 'below') return `<span class="dvmeas below" title="Isopedia pools every splice pattern under 5% usage into 'other', so this transcript has no value of its own — that is a detection limit, not a measurement of zero">measured: <b>&lt;5% (below detection)</b> · n=${m.n} samples</span>`;
  return `<span class="dvmeas" title="Isopedia long-read usage: this isoform's share of all reads assigned to the gene, mean over ${m.n} samples (median ${pct1(m.median)})${m.shared ? '. Shared with another app transcript that has the same intron chain.' : ''}">
    <span class="dvmbar"><i style="width:${(100 * clamp(m.value, 0, 1)).toFixed(1)}%"></i></span>
    measured <b>${pct1(m.value)}</b> · n=${m.n}${m.shared ? ' · shared intron chain' : ''}</span>`;
}

// Gene-specific joins worth naming, read out of data/figures/isopedia_join.csv (the structural
// join behind isoform_usage.json). Only facts that file states; nothing inferred.
const GENE_NOTES = {
  ATF4: 'The dominant splice pattern here (68.8% of ATF4 reads) is catalogued in GENCODE v49 as ENST00000960426.1 and is matched to this app\'s ATF4-206 (ENST00000674920) <b>by intron chain, not by transcript id</b> (±10 bp flank merge); the app\'s other four ATF4 transcripts have no catalogued pattern of their own and fall into “other”.',
};

/** what the measured numbers are a fraction OF — stated, not assumed. */
export function denominatorHTML(geneSym) {
  const g = usageGene(geneSym);
  const meta = usageMeta();
  if (!meta) return '';
  if (!g) return `<div class="dvmeasbox"><b>No measured usage for ${esc(geneSym)}</b> in <code>isoform_usage.json</code>.</div>`;
  if (g.measured === false) return `<div class="dvmeasbox"><b>${esc(geneSym)} is not in the Isopedia panel.</b> ${esc(g.note ?? '')}</div>`;
  const share = g.app_isoform_share_of_gene, other = g.pooled_other_fraction, unmatched = g.unmatched_catalog_fraction;
  return `<div class="dvmeasbox">
    <b>Measured usage (Isopedia long-read, n=${g.n_samples} samples).</b>
    Each measured % is that isoform's share of <b>all ${esc(geneSym)} reads</b>, so the app's isoforms do not sum to 1:
    our ${g.n_app_isoforms} isoform${g.n_app_isoforms === 1 ? '' : 's'} (${g.n_app_isoforms_matched} matched to a measured splice pattern) cover <b>${pct1(share)}</b> of this gene's reads;
    ${pct1(unmatched)} goes to catalogued splice patterns this app does not model and ${pct1(other)} to the pooled “other” bin below Isopedia's 5% threshold.
    ${g.n_app_isoforms_matched < g.n_app_isoforms ? `The ${g.n_app_isoforms - g.n_app_isoforms_matched} unmatched app transcript${g.n_app_isoforms - g.n_app_isoforms_matched === 1 ? ' is' : 's are'} below that threshold in ≥95% of samples, so ${g.n_app_isoforms - g.n_app_isoforms_matched === 1 ? 'it has' : 'they have'} no value of ${g.n_app_isoforms - g.n_app_isoforms_matched === 1 ? 'its' : 'their'} own — <b>not</b> a measured zero.` : ''}
    ${GENE_NOTES[geneSym] ? `<div class="dvgnote">${GENE_NOTES[geneSym]}</div>` : ''}
    <details class="dvcav"><summary>how this was measured, and what it is not</summary>
      <div>${esc(meta.how_quantified ?? '')}</div>
      <div><b>Denominator:</b> ${esc(meta.denominator ?? '')}</div>
      <ul>${(meta.caveats ?? []).map(c => `<li>${esc(c)}</li>`).join('')}</ul>
      <div class="mute">${esc(meta.source ?? '')} · ${esc(meta.citation_or_internal ?? '')}</div>
    </details></div>`;
}

// ---------------------------------------------------------------- the validation headline
// Straight out of data/figures/isopedia_validation_summary.csv (one row per predictor, the same
// 32 isoforms in 11 genes and the same 36 within-gene pairs). Printed, not rounded away.
export const VALIDATION = {
  n_isoforms: 32, n_genes: 11, n_pairs: 36,
  rows: [
    { key: 'model', label: 'kinetic splicing model (this page)', pearson: 0.552, spearman: 0.442, mae: 0.263, pair: 0.486 },
    { key: 'promoter_only', label: 'CAGE promoter weights only, no splicing', pearson: 0.776, spearman: 0.709, mae: 0.193, pair: 0.681 },
    { key: 'uniform', label: 'uniform over annotated isoforms', pearson: 0.555, spearman: 0.469, mae: 0.197, pair: 0.528 },
  ],
  inverted: ['MAPT', 'MDM2', 'CDKN2A', 'TP63'],
  source: 'data/figures/isopedia_validation_summary.csv — predictions from this page\'s model vs Isopedia long-read usage (821 samples)',
};

/** three predictors, two metrics, as bars — never a table. */
export function validationHTML(geneSym) {
  const V = VALIDATION;
  const bar = (v, cls, lo = 0, hi = 1) => `<span class="dvvbar ${cls}"><i style="width:${(100 * clamp((v - lo) / (hi - lo), 0, 1)).toFixed(1)}%"></i></span>`;
  const row = (r) => `<div class="dvvrow${r.key === 'model' ? ' me' : ''}">
    <span class="dvvlab">${esc(r.label)}</span>
    <span class="dvvm">Pearson r ${bar(r.pearson, 'r')} <b>${r.pearson.toFixed(2)}</b></span>
    <span class="dvvm">within-gene pair accuracy ${bar(r.pair, 'p')} <b>${r.pair.toFixed(2)}</b></span>
    <span class="dvvm mute">Spearman ${r.spearman.toFixed(2)} · MAE ${r.mae.toFixed(2)}</span></div>`;
  const hit = V.inverted.includes(geneSym);
  return `<div class="dvvbox">
    <div class="ftitle">Does this model predict real isoform usage? — measured answer</div>
    <div class="dvvlead">Tested against Isopedia long-read usage over <b>${V.n_isoforms} isoforms in ${V.n_genes} genes</b> (${V.n_pairs} within-gene pairs).
      Chance on the pair test is 0.50.</div>
    ${V.rows.map(row).join('')}
    <div class="dvvsum"><b>Promoter choice predicts measured isoform usage better than this page's splicing kinetics do, and the kinetics sit at the coin-flip baseline</b>
      (${V.rows[0].pair.toFixed(2)} vs ${V.rows[1].pair.toFixed(2)} on the pair test).
      With n=${V.n_isoforms} isoforms and ${V.n_pairs} pairs that gap is <b>suggestive, not significant</b>; what is solid is that the kinetic model is <b>not above chance</b> here.</div>
    <div class="dvvinv">It gets the dominant isoform <b>inverted</b> in ${V.inverted.map(esc).join(', ')}${hit ? ` — including <b>${esc(geneSym)}</b>, the gene on screen` : ''}.</div>
    <div class="dvvsrc mute">${esc(V.source)}</div>
  </div>`;
}
