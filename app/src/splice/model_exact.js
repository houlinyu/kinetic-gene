// Deterministic expected-value version of the co-transcriptional splicing model in model.js
// (same rules, same parameters). This is what the UI displays; model.js is the animation.
//
// The commitment race is a continuous-time Markov chain over "which junctions are committed".
// Rates are constant; only availability changes with time, and only at the moments Pol II passes a
// 3' splice site. So:
//   - Phase 1, up to the last availability breakpoint: integrate the state distribution forward.
//   - Phase 2, after the last breakpoint: every remaining junction is available forever, so the
//     winner of each race is exactly k_j / sum(k) -- solved by recursion, no time integration.
// Transcripts that splice post-transcriptionally (fraction postTx) skip phase 1 entirely.
//
// PROMOTER MIXING. Multiple TSSs are mixed with the measured `tss[i].weight` from the data file
// (FANTOM5 CAGE share of the gene's 5'-cap tags over a 41-tissue adult panel; see
// data/splice_data.json meta.methods.tss_weight and pipeline/tss/). The rule, reported back in
// `result.tssWeighting`, is: if any TSS of the gene carries a weight, use the weights and give the
// unmatched promoters 0 (no CAGE peak in the panel = not detectably used; renormalised over the
// matched ones); if none does, fall back to UNIFORM mixing and say so with
// tssWeighting.mode = 'uniform (assumed)'. A weight is never invented here.
//
// PROMOTER-PROXIMAL PAUSING. With probability pauseProb a polymerase stalls at pausePos nt past the
// TSS for an Exp(pauseEscape) time D before resuming (model.js draws exactly this). Everything the
// polymerase has NOT yet reached is delayed by the same D, so the availability schedule of those
// junctions is translated in time -- and the commitment race only sees the *gaps* between
// availability times, so a translation of all of them cannot change the answer. Pausing therefore
// changes the isoform mixture only through splice sites whose 3'ss (+ctLag) is already transcribed
// when the polymerase stalls, i.e. within pausePos nt of the TSS: those sites get the whole pause
// duration to themselves. That is implemented exactly (`pauseResolvent` below integrates the chain
// over the Exp(pauseEscape) stall in closed form), and `result.pauseEffective` says whether any
// gene/TSS actually has such a site -- if not, the pause branch is provably identical and skipped.
//
// ---------------------------------------------------------------------------------------------
// Why this is fast enough to run on every slider move (it used to be exponential in the number
// of introns -- 11 s for MDM2, no answer at all for ERBB2):
//
//  1. FACTORISATION. A junction is only ever interfered with by a junction it *conflicts* with
//     (overlapping intron). Commitment clocks are independent exponentials, so the chain factorises
//     exactly over the connected components of the conflict graph. A junction with no conflicts
//     commits with probability 1 and carries no state at all -- so the ~24 constitutive introns of
//     ERBB2 contribute a factor of 1 each instead of doubling the state space. Real genes decompose
//     into a handful of components of size <= 12, i.e. <= ~1000 states instead of 2^36.
//  2. EXACT PER-COMPONENT PROPAGATOR. Inside a component the old fixed-grid integrator is replaced
//     by uniformisation (exp(Q t) by a Poisson-weighted series), which is exact to machine
//     precision and needs ~1 matrix-vector product per unit of k*t instead of 40.
//  3. EARLY TERMINATION. Availability intervals are usually many commitment half-lives long, so the
//     propagation stops as soon as the probability still able to move drops below `EXACT_OPTS.eps`;
//     the leftover is carried into phase 2 (which resolves everything anyway) and reported.
//  4. PRUNING. Joint branches below `EXACT_OPTS.eps` are dropped. `result.prunedMass` adds up
//     everything either dropped or left for phase 2 by (3), so the loss can be audited; it is
//     <= 1e-8 on every gene in the shipped data, and the isoform fractions still sum to 1 to 1e-12.
//
// Phase 1 is now the exact CTMC rather than the old explicit one-event-per-dt grid, so on genes
// with real competition the answers move by that grid's discretisation error (~1e-3 at the old
// step size); refining the old grid converges onto these numbers. Genes without competing
// junctions are bit-identical.
import { junctions, scoreStats, rateOf, isoformIndex, geneSpan, tssList, sigOf, strength } from './model.js';
import { values } from './params.js';

const PRUNE = 1e-10;

// Solver tolerances. `eps` is the only one worth touching: it is both the smallest joint branch
// kept and the smallest amount of probability the time integration bothers to chase.
export const EXACT_OPTS = {
  eps: 1e-9,        // prune/stop threshold on probability mass (reported as prunedMass)
  stateEps: 1e-13,  // drop single states below this during propagation
  chunkMu: 24,      // uniformisation chunk length, in units of 1/Lambda (doubles if more is needed)
  maxMu: 200,       // keep exp(-mu) representable
  tail: 1e-15,      // Poisson truncation tolerance
  maxCompStates: 1 << 18,
};

// 5-point Gauss-Hermite quadrature over the between-polymerase elongation-rate spread
// (model.js draws one lognormal speed per polymerase; averaging over it here is what makes the
// deterministic and stochastic versions agree rather than differ by Jensen's inequality).
const GH_Z = [0, 1.3556262, -1.3556262, 2.8569700, -2.8569700];
const GH_W = [0.5333333, 0.2220759, 0.2220759, 0.0112574, 0.0112574];
export function speedQuadrature(mean, cv) {
  if (!(cv > 0)) return [{ v: mean, w: 1 }];
  const s = Math.sqrt(Math.log(1 + cv * cv));
  return GH_Z.map((z, i) => ({ v: mean * Math.exp(s * z - s * s / 2), w: GH_W[i] }));
}

// ---------------------------------------------------------------- one competing component
// A component is a set of mutually reachable (via conflicts) junctions. State = bitmask of the
// members already committed; `free` = members that may still commit (not committed, not knocked
// out by a committed conflicting neighbour). `free` depends on the mask alone, so the state graph
// is a DAG and can be enumerated once per gene geometry.
function compGraph(J, grp) {
  const m = grp.length;
  // the shipped data's largest component is 12 members / 961 states, so this is ~250x of headroom;
  // it is a loud error rather than a silent approximation if some future annotation exceeds it.
  if (m > 26) throw new Error(`splice: ${m} mutually competing junctions is beyond the exact solver`);
  const at = new Map(); grp.forEach((idx, i) => at.set(idx, i));
  const conf = new Int32Array(m);
  for (let i = 0; i < m; i++) {
    for (const c of J[grp[i]].conflicts) { const q = at.get(c); if (q !== undefined) conf[i] |= 1 << q; }
  }
  const FULL = (1 << m) - 1;
  const masks = [0], free = [FULL], idOf = new Map([[0, 0]]);
  for (let s = 0; s < masks.length; s++) {
    const f = free[s];
    for (let i = 0; i < m; i++) {
      if (!(f & (1 << i))) continue;
      const nm = masks[s] | (1 << i);
      if (idOf.has(nm)) continue;
      idOf.set(nm, masks.length);
      masks.push(nm);
      free.push(f & ~(1 << i) & ~conf[i]);
      if (masks.length > EXACT_OPTS.maxCompStates) throw new Error('splice: competing component too large');
    }
  }
  const n = masks.length;
  const succ = new Int32Array(n * m).fill(-1);
  for (let s = 0; s < n; s++) {
    for (let i = 0; i < m; i++) if (free[s] & (1 << i)) succ[s * m + i] = idOf.get(masks[s] | (1 << i));
  }
  return {
    m, n, grp, masks, free: Int32Array.from(free), succ,
    buf: {
      st: { v: new Float64Array(n), act: [] },
      cv: new Float64Array(n), nv: new Float64Array(n), acc: new Float64Array(n),
      cAct: [], nAct: [], aAct: [],
      inA: new Int32Array(n), inN: new Int32Array(n), stamp: 0,
      er: new Float64Array(n), erS: new Int32Array(n), es: 0,
      rk: new Float64Array(m), memo: new Map(), token: -1,
    },
  };
}

// total rate out of a state: sum of k over the set bits of `f` (the still-choosable, available members)
function exitRate(f, rk) {
  let r = 0, t = f;
  while (t) { const b = t & -t; r += rk[31 - Math.clz32(b)]; t ^= b; }
  return r;
}

// exp(Q*theta) applied to the current distribution, by uniformisation:
//   pi(theta) = sum_n Poisson(n; Lambda*theta) * pi * P^n,   P = I + Q/Lambda
// Exact up to the Poisson truncation; Lambda is the largest exit rate still reachable (exit rates
// only decrease as junctions commit, so the current maximum is a valid bound for the whole chunk).
function uniformChunk(C, openMask, st, Lam, mu, diag) {
  const m = C.m, free = C.free, succ = C.succ, B = C.buf, rk = B.rk;
  let cv = B.cv, nv = B.nv, cAct = B.cAct, nAct = B.nAct;
  const acc = B.acc, aAct = B.aAct, inA = B.inA, inN = B.inN;
  const er = B.er, erS = B.erS, es = B.es;
  for (const s of cAct) cv[s] = 0; cAct.length = 0;
  aAct.length = 0;
  const sa = ++B.stamp;
  for (const s of st.act) { cv[s] = st.v[s]; cAct.push(s); }
  let w = Math.exp(-mu), rem = 1 - w;
  const addA = (s, val) => { if (inA[s] !== sa) { inA[s] = sa; acc[s] = 0; aAct.push(s); } acc[s] += val; };
  for (const s of cAct) addA(s, w * cv[s]);
  const N = Math.ceil(mu + 6 * Math.sqrt(mu) + 14);
  for (let it = 1; it <= N; it++) {
    const sn = ++B.stamp;
    nAct.length = 0;
    for (const s of cAct) {
      const pv = cv[s];
      if (!(pv > 0)) continue;
      const f = free[s] & openMask;
      let r;
      if (erS[s] === es) r = er[s];
      else { r = exitRate(f, rk); er[s] = r; erS[s] = es; }
      if (inN[s] !== sn) { inN[s] = sn; nv[s] = 0; nAct.push(s); }
      nv[s] += pv * (1 - r / Lam);
      let t = f;
      while (t) {
        const b = t & -t, i = 31 - Math.clz32(b);
        t ^= b;
        const dst = succ[s * m + i];
        if (inN[dst] !== sn) { inN[dst] = sn; nv[dst] = 0; nAct.push(dst); }
        nv[dst] += pv * rk[i] / Lam;
      }
    }
    const tv = cv; cv = nv; nv = tv;
    const ta = cAct; cAct = nAct; nAct = ta;
    w *= mu / it; rem -= w;
    for (const s of cAct) if (cv[s] > 0) addA(s, w * cv[s]);
    if (rem < EXACT_OPTS.tail && it > mu) break;
  }
  if (rem > 0) for (const s of cAct) if (cv[s] > 0) addA(s, rem * cv[s]);
  B.cv = cv; B.nv = nv; B.cAct = cAct; B.nAct = nAct;
  for (const s of st.act) st.v[s] = 0;
  st.act.length = 0;
  for (const s of aAct) {
    const val = acc[s];
    if (val > EXACT_OPTS.stateEps) { st.v[s] = val; st.act.push(s); }
    else if (val > 0) diag.pruned += val;
  }
}

// Advance the chain over a stall whose duration is Exp(lam) instead of a known length: the
// promoter-proximal pause. Averaging the propagator over the escape time is a resolvent,
//   y = integral_0^inf lam e^{-lam d} pi e^{Q d} dd = lam pi (lam I - Q)^{-1},
// and because a junction can only ever be ADDED to the committed set, the state ids are a
// topological order of Q (succ(s,i) > s always), so lam I - Q is triangular and y comes out of one
// forward substitution -- exact, no series, O(states x members). Row-wise:
//   y_u = (lam pi_u + sum_{s->u} y_s k_{s->u}) / (lam + r_u),   r_u = exit rate of u.
// `openMask` is the set available during the stall (by construction only the pre-pause members).
function pauseResolvent(C, openMask, st, lam, diag) {
  const m = C.m, n = C.n, free = C.free, succ = C.succ, B = C.buf, rk = B.rk;
  if (!st.act.length) return;
  const inflow = B.cv, out = B.nv;              // reuse the propagation scratch
  inflow.fill(0); out.fill(0);
  let lo = n, hi = 0;
  for (const s of st.act) { inflow[s] = lam * st.v[s]; if (s < lo) lo = s; if (s > hi) hi = s; }
  for (const s of st.act) st.v[s] = 0;
  st.act.length = 0;
  for (let u = lo; u < n; u++) {
    const inf = inflow[u];
    if (!(inf > 0)) continue;
    const f = free[u] & openMask;
    const y = inf / (lam + exitRate(f, rk));
    out[u] = y;
    let t = f;
    while (t) {
      const b = t & -t, i = 31 - Math.clz32(b);
      t ^= b;
      inflow[succ[u * m + i]] += y * rk[i];
    }
  }
  for (let u = lo; u < n; u++) {
    const y = out[u];
    if (y > EXACT_OPTS.stateEps) { st.v[u] = y; st.act.push(u); }
    else if (y > 0) diag.pruned += y;
  }
  inflow.fill(0); out.fill(0);
  B.cAct.length = 0; B.nAct.length = 0;   // cv/nv are now all-zero, so their active lists are empty
}

// integrate one availability interval of length L, stopping once the mass that can still move is
// negligible (phase 2 resolves whatever is left over, so stopping early costs < eps)
function propagate(C, openMask, st, L, diag) {
  const free = C.free, B = C.buf, rk = B.rk, er = B.er, erS = B.erS;
  const es = ++B.es;   // the cached exit rates are only valid for this interval's `openMask`
  let remain = L, muT = EXACT_OPTS.chunkMu;
  while (remain > 1e-9) {
    let Lam = 0, mov = 0;
    for (const s of st.act) {
      const pv = st.v[s];
      if (!(pv > 0)) continue;
      const f = free[s] & openMask;
      if (!f) continue;
      mov += pv;
      let r;
      if (erS[s] === es) r = er[s];
      else { r = exitRate(f, rk); er[s] = r; erS[s] = es; }
      if (r > Lam) Lam = r;
    }
    if (!(Lam > 0)) return;
    if (mov < EXACT_OPTS.eps) { diag.deferred += mov; return; }
    const mu = Math.min(muT, Lam * remain, EXACT_OPTS.maxMu);
    if (!(mu > 1e-12)) return;
    uniformChunk(C, openMask, st, Lam, mu, diag);
    remain -= mu / Lam;
    muT = Math.min(muT * 2, EXACT_OPTS.maxMu);
  }
}

// phase 2 inside one component: every member is available, so each race is won with probability
// k_j / sum(k). Memoised on the state id -> polynomial in the number of states, not orderings.
function raceTail(C, s, memo) {
  const hit = memo.get(s);
  if (hit) return hit;
  const m = C.m, f = C.free[s], rk = C.buf.rk;
  let out;
  if (!f) out = new Map([[C.masks[s], 1]]);
  else {
    const K = exitRate(f, rk);
    out = new Map();
    let t = f;
    while (t) {
      const b = t & -t, i = 31 - Math.clz32(b);
      t ^= b;
      const w = rk[i] / K;
      for (const [tm, tv] of raceTail(C, C.succ[s * m + i], memo)) out.set(tm, (out.get(tm) || 0) + w * tv);
    }
  }
  memo.set(s, out);
  return out;
}

// terminal distribution of one component: Map<committed bitmask, probability>
// `pause` = null, or { tP, rate }: the polymerase stalls at time tP for an Exp(rate) interval, so
// the schedule is the unpaused one with one extra, random-length stall inserted at tP (everything
// the polymerase has not yet reached is pushed back by the same amount, which is why nothing else
// in the schedule changes).
function solveComp(C, k, availT, coTx, token, diag, pause = null) {
  const B = C.buf, m = C.m;
  for (let i = 0; i < m; i++) B.rk[i] = k[C.grp[i]];
  if (B.token !== token) { B.memo.clear(); B.token = token; }
  const st = B.st;
  for (const s of st.act) st.v[s] = 0;
  st.act.length = 0;
  st.v[0] = 1; st.act.push(0);
  if (coTx) {
    const marks = [...availT, 0];
    if (pause) marks.push(pause.tP);
    const bps = [...new Set(marks)].sort((a, b) => a - b);
    for (let b = 0; b < bps.length; b++) {
      const t0 = bps[b], L = b + 1 < bps.length ? bps[b + 1] - t0 : 0;
      let openMask = 0;
      for (let i = 0; i < m; i++) if (availT[i] <= t0 + 1e-12) openMask |= 1 << i;
      if (pause && Math.abs(t0 - pause.tP) < 1e-12 && openMask) pauseResolvent(C, openMask, st, pause.rate, diag);
      if (!(L > 0) || !openMask) continue;
      propagate(C, openMask, st, L, diag);
    }
  }
  const out = new Map();
  for (const s of st.act) {
    const pv = st.v[s];
    if (!(pv > 0)) continue;
    for (const [tm, tv] of raceTail(C, s, B.memo)) out.set(tm, (out.get(tm) || 0) + pv * tv);
  }
  return out;
}

// ---------------------------------------------------------------- gene geometry (cached)
// Everything that depends only on the gene's exon/intron geometry, not on the parameters: the
// junction list, the isoform lookup, which junctions cut which exon, and -- per TSS -- the
// decomposition of the still-reachable junctions into competing components plus the conflict-free
// junctions, which always commit. Keyed weakly on the gene object, so the findings' perturbed
// copies get their own entry and nothing leaks.
const geoCache = new WeakMap();
let TOKEN = 0;

function geometry(gene) {
  let G = geoCache.get(gene);
  if (G) return G;
  const J = junctions(gene);
  const exons = gene.exons || [];
  const cuts = exons.map((e) => J.filter((j) => j.donorPos < e.end && e.start < j.acceptorPos).map((j) => j.idx));
  const tss = tssList(gene);
  const parts = new Map();
  // different TSSs usually see the same competing group; share one graph (and so one phase-2
  // memo, since the rates do not depend on the TSS) instead of rebuilding it per TSS.
  const graphs = new Map();
  for (const t of tss) {
    if (parts.has(t.pos)) continue;
    const on = J.filter((j) => j.donorPos >= t.pos);
    const live = new Set(on.map((j) => j.idx));
    const seen = new Set(), comps = [], solo = [];
    for (const j of on) {
      if (seen.has(j.idx)) continue;
      const stack = [j.idx], grp = [];
      while (stack.length) {
        const i = stack.pop();
        if (seen.has(i)) continue;
        seen.add(i); grp.push(i);
        for (const c of J[i].conflicts) if (live.has(c) && !seen.has(c)) stack.push(c);
      }
      if (grp.length > 1) {
        const gk = grp.sort((a, b) => a - b).join(',');
        let C = graphs.get(gk);
        if (!C) { C = compGraph(J, grp); graphs.set(gk, C); }
        comps.push(C);
      } else solo.push(grp[0]);
    }
    solo.sort((a, b) => a - b);
    parts.set(t.pos, { comps, solo });
  }
  G = { J, index: isoformIndex(gene), span: geneSpan(gene), tss, exons, cuts, parts };
  geoCache.set(gene, G);
  return G;
}

/**
 * expected(gene, params) -> deterministic expected outcome for one gene.
 * { fractions, isoforms, junctionUse, exonUse, psi, meanRemovalDistance, prunedMass }
 */
export function expected(gene, params = values()) {
  const p = params;
  const G = geometry(gene);
  const J = G.J;
  const stats = expected.stats || scoreStats({ genes: [gene] });
  const k = J.map((j) => rateOf(j, p, stats));
  const index = G.index;
  const token = ++TOKEN;
  const diag = { pruned: 0, deferred: 0 };

  // --- promoter mixing weights (measured CAGE share, else uniform; see the header note)
  const measured = G.tss.filter((t) => Number.isFinite(t.weight) && t.weight >= 0);
  let msum = measured.reduce((a, t) => a + t.weight, 0);
  const useMeasured = measured.length > 0 && msum > 0;
  if (!useMeasured) msum = G.tss.length;
  const tssWeighting = {
    mode: useMeasured ? 'measured' : 'uniform (assumed)',
    matched: useMeasured ? measured.length : 0,
    total: G.tss.length,
    source: useMeasured ? (measured[0].weight_source?.dataset || null) : null,
    note: useMeasured
      ? `promoter weights from ${measured.length}/${G.tss.length} annotated TSSs with a measured CAGE peak; unmatched promoters get 0 and the rest are renormalised`
      : `no measured promoter weight for this gene: its ${G.tss.length} TSSs are mixed uniformly (assumed)`,
  };

  // --- promoter-proximal pause (see the header note: only pre-pause 3'ss can feel it)
  const pauseOn = p.pauseProb > 0 && p.pauseEscape > 0 && Number.isFinite(p.pauseEscape);
  let pauseEffective = false;

  const final = new Map(); // committedKey -> {c, p, tssPos}
  for (const t of G.tss) {
    const w = (useMeasured ? (Number.isFinite(t.weight) ? t.weight : 0) : 1) / msum;
    if (!(w > 0)) continue;   // promoter not used (measured 0 / unmatched): contributes nothing
    const { comps, solo } = G.parts.get(t.pos);
    for (const coTx of [true, false]) {
      const mix = coTx ? 1 - p.postTx : p.postTx;
      if (mix <= 0) continue;
      // average the co-transcriptional branch over the spread of polymerase speeds. The mixture is
      // over whole transcripts, so it stays OUTSIDE the product over components.
      const quad = coTx ? speedQuadrature(p.elongRate, p.elongCV) : [{ v: p.elongRate, w: 1 }];
      for (const q of quad) {
        const v = Math.max(q.v, 1e-6);
        const tP = p.pausePos / v;                       // when the polymerase reaches the pause site
        const availT = [];
        let anyPre = false;
        for (const C of comps) {
          let a = null;
          if (coTx) {
            a = new Float64Array(C.m);
            for (let i = 0; i < C.m; i++) {
              a[i] = Math.max(0, (J[C.grp[i]].acceptorPos + p.ctLag - t.pos) / v);
              if (a[i] <= tP + 1e-12) anyPre = true;     // a 3'ss inside the paused window
            }
          }
          availT.push(a);
        }
        // Only split into paused / unpaused branches when the pause can actually be felt; when no
        // 3'ss is available during the stall the two branches are provably identical.
        const usePause = coTx && pauseOn && anyPre;
        if (usePause) pauseEffective = true;
        const branches = usePause
          ? [{ pause: { tP, rate: p.pauseEscape }, w: p.pauseProb }, { pause: null, w: 1 - p.pauseProb }]
          : [{ pause: null, w: 1 }];
        for (const br of branches) {
          if (!(br.w > 0)) continue;
          let states = [{ c: solo, p: 1 }];
          for (let ci = 0; ci < comps.length; ci++) {
            const dist = solveComp(comps[ci], k, availT[ci], coTx, token, diag, br.pause);
            const next = [];
            for (const st of states) {
              for (const [tm, tv] of dist) {
                const pr = st.p * tv;
                if (pr < EXACT_OPTS.eps) { diag.pruned += pr; continue; }
                const c = st.c.slice();
                for (let i = 0; i < comps[ci].m; i++) if (tm & (1 << i)) c.push(comps[ci].grp[i]);
                next.push({ c, p: pr });
              }
            }
            states = next;
          }
          for (const st of states) {
            const c = st.c === solo ? solo.slice() : st.c.sort((a, b) => a - b);
            const kk = `${t.pos}|${c.join(',')}`;
            const pr = w * mix * q.w * br.w * st.p;
            const e = final.get(kk);
            if (e) e.p += pr; else final.set(kk, { c, p: pr, tssPos: t.pos });
          }
        }
      }
    }
  }

  // --- isoform fractions
  const fractions = {};
  const junctionUse = {};
  const exonUse = {};
  for (const j of J) junctionUse[j.id] = 0;
  for (const e of G.exons) exonUse[e.id] = 0;
  const on = new Uint8Array(J.length);
  for (const st of final.values()) {
    const iso = index.get(sigOf(J, st.c));
    const id = iso ? iso.id : 'other';
    fractions[id] = (fractions[id] || 0) + st.p;
    for (const i of st.c) { junctionUse[J[i].id] += st.p; on[i] = 1; }
    for (let ei = 0; ei < G.exons.length; ei++) {
      const e = G.exons[ei];
      if (e.start < st.tssPos) continue;
      let cut = false;
      for (const i of G.cuts[ei]) if (on[i]) { cut = true; break; }
      if (!cut) exonUse[e.id] += st.p;
    }
    for (const i of st.c) on[i] = 0;
  }
  const tot = Object.values(fractions).reduce((a, b) => a + b, 0) || 1;
  for (const key of Object.keys(fractions)) fractions[key] /= tot;
  for (const key of Object.keys(junctionUse)) junctionUse[key] /= tot;
  for (const key of Object.keys(exonUse)) exonUse[key] /= tot;

  // --- per-event PSI: alternative exons, plus junction usage for alt 5'/3' choices
  const psi = [];
  for (const e of G.exons) {
    // only whole-exon skipping has a meaningful "inclusion" PSI; alt 5'/3' ends are reported
    // as junction usage below (both alt exons share a start, so neither is ever fully removed).
    if (e.type === 'cassette') psi.push({ event: `${e.id} inclusion`, exon: e.id, value: exonUse[e.id] });
  }
  for (const j of J) {
    if (J.some((o) => o !== j && (o.donorPos === j.donorPos || o.acceptorPos === j.acceptorPos))) {
      psi.push({ event: `${j.id} usage`, junction: j.id, value: junctionUse[j.id] });
    }
  }

  // mean nt Pol II travels past a 3'ss before that intron is physically gone (Drexler 2020 check)
  const dists = J.map((j) => p.elongRate * (1 / k[j.idx] + p.removalTime));
  const meanRemovalDistance = dists.length ? dists.reduce((a, b) => a + b, 0) / dists.length : 0;

  return {
    fractions,
    isoforms: (gene.isoforms || []).map((i) => ({ id: i.id, name: i.name, value: fractions[i.id] || 0 })),
    other: fractions.other || 0,
    junctionUse, exonUse, psi, meanRemovalDistance,
    tssWeighting,
    // true only if some 3'ss of this gene is already transcribed when the polymerase stalls, i.e.
    // if pauseProb/pausePos/pauseEscape can change this gene's numbers at all (see header).
    pauseEffective,
    prunedMass: diag.pruned + diag.deferred,
    junctions: J.map((j, i) => ({ id: j.id, key: j.key, donorPos: j.donorPos, acceptorPos: j.acceptorPos, rate: k[i], use: junctionUse[j.id] })),
  };
}

// Calibrate the score scale on the whole dataset once (so `scoreSlope` means the same thing
// for every gene). Safe to skip: expected() then calibrates per gene.
export function calibrate(data) { expected.stats = scoreStats(data); return expected.stats; }

// ---------------------------------------------------------------- findings
// Direction checks on ONE shared parameter set, computed with the exact model (never tuned per finding).
// status: 'mechanism' = should fall out of the kinetics; 'check' = internal consistency, not a measurement.
const SLOW = 8, FAST = 67;  // nt/s = 0.5 and 4 kb/min, the ends of the Jonkers 2014 measured range

// Groups of mutually exclusive junctions (connected components of the conflict graph).
export function competingGroups(J) {
  const seen = new Set(), out = [];
  for (const j of J) {
    if (seen.has(j.idx) || !j.conflicts.length) continue;
    const stack = [j.idx], grp = [];
    while (stack.length) {
      const i = stack.pop();
      if (seen.has(i)) continue;
      seen.add(i); grp.push(i);
      for (const c of J[i].conflicts) if (!seen.has(c)) stack.push(c);
    }
    if (grp.length > 1) out.push(grp);
  }
  return out;
}

// The de la Mata / Dujardin read-out, uniform across event types: in each competing group, how
// often does the junction whose 3'ss appears FIRST win? (For a cassette exon that is inclusion.)
export const upstreamUse = (g, p) => {
  const J = junctions(g);
  const groups = competingGroups(J);
  if (!groups.length) return null;
  const r = expected(g, p);
  const xs = groups.map((grp) => {
    const first = grp.reduce((a, b) => (J[b].acceptorPos < J[a].acceptorPos ? b : a));
    return r.junctionUse[J[first].id];
  }).filter(Number.isFinite);
  return xs.length ? xs.reduce((a, b) => a + b, 0) / xs.length : null;
};
const withAlt = (data) => (data.genes || []).filter((g) => competingGroups(junctions(g)).length);

export const FINDINGS_SPLICE = [
  {
    id: 'kinetic_coupling',
    title: 'Slow Pol II raises inclusion of alternative exons',
    cite: 'dujardin2014',
    status: 'mechanism',
    note: 'Direction and wording from Dujardin 2014 ("In most cases, slow Pol II elongation allows weak splice sites to be recognized, leading to higher inclusion of alternative exons"); the effect itself was shown by de la Mata 2003. Slow/fast = 0.5 and 4 kb/min, the ends of the Jonkers 2014 measured range. CFTR exon 9 (Dujardin 2014) is a known exception the model does not contain.',
    test(data, p) {
      const gs = withAlt(data);
      if (!gs.length) return null;
      const d = gs.map((g) => {
        const s = upstreamUse(g, { ...p, elongRate: SLOW }), f = upstreamUse(g, { ...p, elongRate: FAST });
        return s == null || f == null ? null : s - f;
      }).filter((x) => x != null);
      if (!d.length) return null;
      const mean = d.reduce((a, b) => a + b, 0) / d.length;
      const up = d.filter((x) => x > 0.005).length;
      return { pass: mean > 0.005, detail: `mean ΔPSI slow−fast = ${(mean * 100).toFixed(1)} pp (up in ${up}/${d.length} genes)` };
    },
  },
  {
    id: 'site_strength',
    title: 'The stronger of two competing splice sites is used more',
    cite: 'dujardin2014',
    status: 'mechanism',
    note: 'Dujardin 2014 state the perturbation directly: "Weakening of E9 5′ ss increases ... E9 skipping, whereas strengthening of the 5′ ss usage has the opposite effect." So this is run as a perturbation (strengthen one competing 5′ss by one MAD of the score scale and re-solve), not as a correlation — raw score order does not have to predict usage, because availability timing and exon definition also matter.',
    test(data, p) {
      const stats = expected.stats || scoreStats(data);
      let ok = 0, n = 0;
      const bump = (g, pos, d) => ({
        ...g,
        introns: (g.introns || []).map((i) => (i.donor?.pos === pos
          ? { ...i, donor: { ...i.donor, score: (i.donor.score ?? stats.centre) + d } } : i)),
      });
      for (const g of data.genes || []) {
        const J = junctions(g);
        for (const grp of competingGroups(J)) {
          for (const idx of grp) {
            const j = J[idx];
            if (!j.donor || !Number.isFinite(j.donor.score)) continue;
            // only sites whose 5'ss belongs to one junction: strengthening a shared donor would
            // raise every junction that uses it, which tests nothing.
            if (J.filter((o) => o.donorPos === j.donorPos).length !== 1) continue;
            const before = expected(g, p).junctionUse[j.id];
            const after = expected(bump(g, j.donorPos, stats.spread), p).junctionUse[j.id];
            if (before == null || after == null) continue;
            n++;
            if (after > before - 1e-9) ok++;
          }
        }
      }
      if (!n) return null;
      return { pass: n > 0 && ok === n, detail: `strengthening a 5′ss raised its own use in ${ok}/${n} competing sites` };
    },
  },
  {
    id: 'cotx_timing',
    title: 'Introns are removed minutes after synthesis, not instantly',
    cite: 'singh2009',
    status: 'mechanism',
    note: 'Singh & Padgett 2009: U2-dependent introns spliced "within 5-10 min of synthesis"; Drexler 2020 reports mammalian splicing half-lives of 7-14 min. Pass window 5-14 min is the union of the two.',
    test(data, p) {
      const ts = [];
      const stats = expected.stats || scoreStats(data);
      for (const g of data.genes || []) for (const j of junctions(g)) ts.push((1 / rateOf(j, p, stats) + p.removalTime) / 60);
      if (!ts.length) return null;
      ts.sort((a, b) => a - b);
      const med = ts[Math.floor(ts.length / 2)];
      return { pass: med >= 5 && med <= 14, detail: `median ${med.toFixed(1)} min from 3′ss synthesis to removal` };
    },
  },
  {
    id: 'coupling_needs_cotx',
    title: 'The elongation-rate effect disappears if splicing is post-transcriptional',
    cite: 'delamata2003',
    status: 'check',
    note: 'Internal consistency, not a measurement: de la Mata 2003 showed elongation controls alternative splicing in vivo, which requires splicing to see sites in the order they are made. With postTx = 1 the model must lose the effect.',
    test(data, p) {
      const gs = withAlt(data);
      if (!gs.length) return null;
      const q = { ...p, postTx: 1 };
      const d = gs.map((g) => Math.abs(upstreamUse(g, { ...q, elongRate: SLOW }) - upstreamUse(g, { ...q, elongRate: FAST }))).filter(Number.isFinite);
      if (!d.length) return null;
      const worst = Math.max(...d);
      return { pass: worst < 1e-6, detail: `max |ΔPSI| with postTx=1 = ${worst.toExponential(1)}` };
    },
  },
  {
    id: 'distal_splicing',
    title: 'Human introns are spliced kilobases behind Pol II',
    cite: 'drexler2020',
    status: 'mechanism',
    note: 'Drexler 2020: "the vast majority of human introns were not spliced until after Pol II transcribed 4 kb downstream". Checked as a model readout (>= 2 kb) because the value is emergent, not a parameter.',
    test(data, p) {
      const ds = (data.genes || []).map((g) => expected(g, p).meanRemovalDistance).filter(Number.isFinite);
      if (!ds.length) return null;
      const mean = ds.reduce((a, b) => a + b, 0) / ds.length;
      return { pass: mean >= 2000, detail: `mean ${(mean / 1000).toFixed(1)} kb past the 3′ss at removal` };
    },
  },
];

// Measured isoform/PSI comparisons appended at runtime from the data file, shape:
// data.splice_findings = [{ id, title, cite, gene, higher: isoform id, lower: isoform id }]
export function allSpliceFindings(data) {
  const extra = (data.splice_findings || []).map((f) => ({
    id: f.id, title: f.title, cite: f.cite, status: 'mechanism', note: f.note || null,
    test(d, p) {
      const g = (d.genes || []).find((x) => x.gene === f.gene);
      if (!g) return null;
      const r = expected(g, p);
      const pick = (id) => r.isoforms.find((i) => i.id.split('.')[0] === String(id).split('.')[0]);
      const hi = pick(f.higher), lo = pick(f.lower);
      if (!hi || !lo) return null;
      return { pass: hi.value > lo.value, detail: `${hi.name || hi.id} ${(hi.value * 100).toFixed(0)}% vs ${lo.name || lo.id} ${(lo.value * 100).toFixed(0)}%` };
    },
  }));
  return [...FINDINGS_SPLICE, ...extra];
}

export function scoreSplice(data, p = values()) {
  const res = allSpliceFindings(data).map((f) => ({ ...f, result: f.test(data, p) })).filter((f) => f.result);
  const mech = res.filter((f) => f.status === 'mechanism');
  return { res, passed: mech.filter((f) => f.result.pass).length, total: mech.length };
}
