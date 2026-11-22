// Stochastic kinetic model of transcription + co-transcriptional splicing for one gene.
// Positions are nt along the pre-mRNA in genomic order (same frame as splice_data.json).
//
// Rules (one sentence each, cf. data/SPEC_SPLICE.md):
//  - Pol II initiates at a TSS as a Poisson process (rate initRate), pauses at pausePos with
//    probability pauseProb, and escapes the pause at rate pauseEscape.
//  - Each polymerase elongates at its own speed, drawn once from a lognormal with mean elongRate
//    and coefficient of variation elongCV (polymerases differ; one polymerase is steady).
//  - A junction (one 5'ss paired with one 3'ss) can only be committed once Pol II has transcribed
//    its 3'ss plus ctLag nt, so upstream junctions become available earlier.
//  - Commitment is a Poisson event whose rate is commitScale x strength(5'ss) x strength(3'ss),
//    with strength a logistic function of the site's score (slope scoreSlope) so any score scale works.
//  - Exon definition: a junction's 3'ss gets a boost up to exonDefGain for a downstream exon
//    shorter than exonWindow, fading as the exon grows.
//  - Junctions whose introns overlap are mutually exclusive, so they race: committing one marks
//    every overlapping competitor 'skipped'.
//  - The intron is physically removed removalTime seconds after commitment (commitment, not removal,
//    decides the isoform), and the transcript is finished when transcription and all races are done.
//  - A fraction postTx of transcripts splice post-transcriptionally: every site is available at once,
//    so there is no kinetic coupling for them.
import { values } from './params.js';

// ---------------------------------------------------------------- score handling
// Splice-site scores arrive on an arbitrary (MaxEntScan-like) scale, so centre and spread are
// estimated from the data itself and `scoreSlope` is dimensionless.
const median = (a) => { const s = [...a].sort((x, y) => x - y); return s.length ? s[Math.floor(s.length / 2)] : 0; };

export function scoreStats(data) {
  const xs = [];
  for (const g of data?.genes || []) for (const i of g.introns || []) {
    for (const s of [i.donor, i.acceptor]) if (s && Number.isFinite(s.score)) xs.push(s.score);
  }
  if (xs.length < 4) return { centre: 0, spread: 1, n: xs.length };
  const c = median(xs);
  const spread = median(xs.map((x) => Math.abs(x - c))) * 1.4826 || 1; // MAD -> sigma
  return { centre: c, spread: Math.max(spread, 1e-6), n: xs.length };
}

// strength in (0,2), 1 at a median-scoring site; a missing score is treated as median (=1).
export function strength(site, p, stats) {
  if (!site || !Number.isFinite(site.score)) return 1;
  const z = (site.score - stats.centre) / stats.spread;
  return 2 / (1 + Math.exp(-p.scoreSlope * z));
}

// ---------------------------------------------------------------- gene geometry
function exonMap(gene) { return new Map((gene.exons || []).map((e) => [e.id, e])); }

// Every junction any isoform (or the intron list) uses, with its two sites and its competitors.
export function junctions(gene) {
  const ex = exonMap(gene);
  const sites = { d: new Map(), a: new Map() };
  for (const i of gene.introns || []) {
    if (i.donor?.pos != null) sites.d.set(i.donor.pos, i.donor);
    if (i.acceptor?.pos != null) sites.a.set(i.acceptor.pos, i.acceptor);
  }
  const map = new Map();
  const add = (d, a, id) => {
    const key = `${d}_${a}`;
    if (!map.has(key)) map.set(key, {
      key, id: id || key, donorPos: d, acceptorPos: a,
      donor: sites.d.get(d) || null, acceptor: sites.a.get(a) || null, isoforms: [],
    });
    const j = map.get(key);
    if (id && j.id === key) j.id = id;
    return j;
  };
  for (const i of gene.introns || []) {
    if (i.donor?.pos != null && i.acceptor?.pos != null) add(i.donor.pos, i.acceptor.pos, i.id);
  }
  // Junctions implied by the isoform exon chains (this is what gives us the skipping junctions).
  for (const iso of gene.isoforms || []) {
    const chain = (iso.exon_ids || []).map((id) => ex.get(id)).filter(Boolean).sort((x, y) => x.start - y.start);
    for (let k = 0; k + 1 < chain.length; k++) add(chain[k].end, chain[k + 1].start, null).isoforms.push(iso.id);
  }
  const J = [...map.values()].sort((x, y) => x.donorPos - y.donorPos || x.acceptorPos - y.acceptorPos);
  J.forEach((j, idx) => { j.idx = idx; });
  // overlapping introns cannot both survive -> they compete
  for (const j of J) j.conflicts = J.filter((o) => o !== j && j.donorPos < o.acceptorPos && o.donorPos < j.acceptorPos).map((o) => o.idx);
  // first exon starting at or after the 3'ss: the exon this acceptor defines
  for (const j of J) {
    const e = (gene.exons || []).filter((x) => x.start === j.acceptorPos).sort((a, b) => (a.end - a.start) - (b.end - b.start))[0];
    j.exonLen = e ? e.end - e.start : null;
  }
  return J;
}

export function geneSpan(gene) {
  let end = 0;
  for (const e of gene.exons || []) end = Math.max(end, e.end);
  for (const i of gene.introns || []) end = Math.max(end, i.end ?? 0, i.acceptor?.pos ?? 0);
  return Math.max(end, 1);
}

export function tssList(gene) {
  const t = (gene.tss || []).filter((x) => Number.isFinite(x.pos));
  return t.length ? t : [{ id: 'tss0', pos: 0 }];
}

// Commitment rate of one junction, /s.
export function rateOf(j, p, stats) {
  const sd = strength(j.donor, p, stats);
  const sa = strength(j.acceptor, p, stats);
  const def = j.exonLen == null ? 1 : 1 + p.exonDefGain / (1 + j.exonLen / Math.max(p.exonWindow, 1));
  return Math.max(p.commitScale * sd * sa * def, 1e-9);
}

// Isoform lookup: a mature transcript is identified by the set of junctions it used.
export function isoformIndex(gene) {
  const ex = exonMap(gene);
  const m = new Map();
  for (const iso of gene.isoforms || []) {
    const chain = (iso.exon_ids || []).map((id) => ex.get(id)).filter(Boolean).sort((x, y) => x.start - y.start);
    const sig = [];
    for (let k = 0; k + 1 < chain.length; k++) sig.push(`${chain[k].end}_${chain[k + 1].start}`);
    const key = sig.sort().join('|');
    if (!m.has(key)) m.set(key, iso);
  }
  return m;
}

export const sigOf = (J, committed) => committed.map((i) => J[i].key).sort().join('|');

// ---------------------------------------------------------------- stochastic simulation
const expDraw = (rate) => -Math.log(1 - Math.random()) / Math.max(rate, 1e-12);
function lognormalDraw(mean, cv) {
  if (!(cv > 0)) return mean;
  const s2 = Math.log(1 + cv * cv);
  const z = Math.sqrt(-2 * Math.log(1 - Math.random())) * Math.cos(2 * Math.PI * Math.random());
  return mean * Math.exp(Math.sqrt(s2) * z - s2 / 2);
}

export class Gene {
  constructor(gene, data = null) {
    this.gene = gene;
    this.J = junctions(gene);
    this.stats = scoreStats(data || { genes: [gene] });
    this.isoIndex = isoformIndex(gene);
    this.span = geneSpan(gene);
    this.tss = tssList(gene);
    this.reset();
  }

  reset() {
    this.t = 0;
    this.pols = [];        // nascent transcripts, each with its own polymerase
    this.counts = {};      // finished transcripts tallied by isoform id (or 'other')
    this.finished = 0;
    this.events = [];      // consumed by the view: {kind, intron, isoform, pos, pol}
    this.nextId = 1;
  }

  // rates are shared across polymerases, so cache them per parameter set
  rates(p) {
    if (!this._rk || this._rp !== p.commitScale + 1e3 * p.scoreSlope + 1e6 * p.exonDefGain + 1e9 * p.exonWindow) {
      this._rk = this.J.map((j) => rateOf(j, p, this.stats));
      this._rp = p.commitScale + 1e3 * p.scoreSlope + 1e6 * p.exonDefGain + 1e9 * p.exonWindow;
    }
    return this._rk;
  }

  _drawTss() {
    if (this._cum === undefined) {
      const w = this.tss.map(t => (typeof t.weight === 'number' && t.weight > 0 ? t.weight : 0));
      const tot = w.reduce((s, x) => s + x, 0);
      if (tot > 0) { let a = 0; this._cum = w.map(x => (a += x / tot)); }
      else this._cum = null; // no measured weights -> uniform
    }
    if (!this._cum) return this.tss[Math.floor(Math.random() * this.tss.length)];
    const r = Math.random();
    for (let i = 0; i < this._cum.length; i++) if (r <= this._cum[i]) return this.tss[i];
    return this.tss[this.tss.length - 1];
  }

  spawn(p) {
    // Draw the promoter by measured CAGE weight when the data has them; uniform otherwise.
    // Keeps the animation consistent with expected(), which uses the same weights.
    const tss = this._drawTss();
    const paused = Math.random() < p.pauseProb;
    this.pols.push({
      id: this.nextId++, tssId: tss.id, tssPos: tss.pos,
      pos: tss.pos, len: 0,
      speed: lognormalDraw(p.elongRate, p.elongCV),
      paused, pauseAt: tss.pos + p.pausePos, released: !paused,
      postTx: Math.random() < p.postTx,
      committed: [], removed: [], pending: [], skipped: [],
      txDone: false, done: false, isoform: null,
    });
  }

  // junctions this transcript could still commit: compatible, in range, and transcribed
  candidates(r, p) {
    const out = [];
    for (const j of this.J) {
      if (r.committed.includes(j.idx)) continue;
      if (j.donorPos < r.tssPos) continue;                                  // not on this transcript
      if (r.committed.some((c) => this.J[c].conflicts.includes(j.idx))) continue; // lost the race
      const ready = r.postTx ? r.txDone : r.pos >= j.acceptorPos + p.ctLag;
      if (ready) out.push(j);
    }
    return out;
  }

  // junctions that could still become available later (used to decide "the races are over")
  pendingAny(r) {
    return this.J.some((j) => !r.committed.includes(j.idx) && j.donorPos >= r.tssPos
      && !r.committed.some((c) => this.J[c].conflicts.includes(j.idx)));
  }

  step(dt, p = values()) {
    this.t += dt;
    const k = this.rates(p);
    if (Math.random() < 1 - Math.exp(-p.initRate * dt)) this.spawn(p);

    for (const r of this.pols) {
      // --- move
      if (!r.txDone) {
        if (!r.released) {
          if (r.pos < r.pauseAt) r.pos = Math.min(r.pauseAt, r.pos + r.speed * dt);
          else if (Math.random() < 1 - Math.exp(-p.pauseEscape * dt)) r.released = true;
        } else {
          r.pos += r.speed * dt;
        }
        if (r.pos >= this.span) { r.pos = this.span; r.txDone = true; this.events.push({ kind: 'txEnd', pol: r.id, pos: r.pos }); }
        r.len = r.pos - r.tssPos;
      }
      // --- commit (race)
      for (const j of this.candidates(r, p)) {
        if (Math.random() < 1 - Math.exp(-k[j.idx] * dt)) {
          r.committed.push(j.idx);
          r.pending.push({ idx: j.idx, at: this.t + p.removalTime });
          this.events.push({ kind: 'commit', intron: j.id, pol: r.id, pos: j.donorPos, to: j.acceptorPos });
          for (const c of j.conflicts) {
            if (!r.committed.includes(c) && !r.skipped.includes(c)) {
              r.skipped.push(c);
              this.events.push({ kind: 'skipped', intron: this.J[c].id, pol: r.id, pos: this.J[c].donorPos });
            }
          }
          break; // one commitment per step keeps the race fair at coarse dt
        }
      }
      // --- removal
      for (let i = r.pending.length - 1; i >= 0; i--) {
        if (r.pending[i].at <= this.t) {
          const idx = r.pending[i].idx;
          r.removed.push(idx); r.pending.splice(i, 1);
          this.events.push({ kind: 'spliced', intron: this.J[idx].id, pol: r.id, pos: this.J[idx].donorPos, to: this.J[idx].acceptorPos });
        }
      }
      // --- finish
      if (r.txDone && !r.pending.length && !this.candidates(r, p).length && !this.pendingAny(r)) {
        const iso = this.isoIndex.get(sigOf(this.J, r.committed));
        r.isoform = iso ? iso.id : 'other';
        this.counts[r.isoform] = (this.counts[r.isoform] || 0) + 1;
        this.finished++;
        this.events.push({ kind: 'finished', pol: r.id, isoform: r.isoform, name: iso?.name || 'other' });
        r.done = true;
      }
    }
    this.pols = this.pols.filter((r) => !r.done);
    if (this.events.length > 400) this.events.splice(0, this.events.length - 400);
  }

  fractions() {
    const tot = Object.values(this.counts).reduce((a, b) => a + b, 0);
    const f = {};
    for (const [k2, v] of Object.entries(this.counts)) f[k2] = tot ? v / tot : 0;
    return f;
  }

  // nascent-RNA view state: one entry per polymerase
  view() {
    return this.pols.map((r) => ({
      id: r.id, pos: r.pos, len: r.len, tssPos: r.tssPos, paused: !r.released,
      txDone: r.txDone, postTx: r.postTx,
      committed: r.committed.map((i) => this.J[i].id),
      removed: r.removed.map((i) => this.J[i].id),
      skipped: r.skipped.map((i) => this.J[i].id),
    }));
  }
}
