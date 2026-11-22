// Wall-clock pacing for the two animations. VIEW LAYER ONLY.
//
// The models keep their real units (nt/s, /s, s) and still produce every number the pages
// print; nothing here is read by sim.js / analytic.js / splice/model*.js. This module answers
// one question for the renderers: how many MODEL seconds should pass per second of WALL time,
// so that the moving marker (ribosome bead, Pol II bead, nascent RNA) is neither frozen nor a
// blur — at both ends of every speed slider, and for the shortest and the longest gene / 5'
// leader in the data.
//
// It is a compression, never a misrepresentation:
//   - label() returns the factor, which the page prints live ("1 s on screen = 12 s in the cell");
//   - the slow<->fast ratio is only COMPRESSED (screen ratio = real ratio ^ GAMMA), never
//     flattened, so slowing Pol II still LOOKS slower;
//   - ratioNote() states that compression, with the real and the on-screen ratio side by side,
//     and the `limit` field says when an extreme hit the edge of the watchable window,
//     which both pages print as "held at the slowest/fastest watchable pace".

export const PACE = {
  MID: 6,                       // wall seconds one traversal should take at the DEFAULT settings
  MIN: 2.5, MAX: 20,            // preferred window for one traversal (the user's 2-20 s target)
  HARD_MIN: 1.5, HARD_MAX: 26,  // never outside this, whatever else asks
  MAX_GAP: 12,                  // longest wall gap between two visible events before it looks dead
  GAMMA: 0.45,                  // share of the real slow<->fast ratio kept on screen
  LGAMMA: 0.18,                 // share of the real short<->long gene ratio kept on screen
  MAX_WPS: 1.2,                 // blur ceiling: lane widths per second for the FASTEST marker
  MIN_STEP: 0.03,               // freeze floor: the marker must move this much of its lane ...
  MIN_STEP_S: 30,               // ... within this many wall seconds
};

const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

// nt       : nucleotides the marker has to cross (whole lane / whole gene)
// rate     : the model's real speed for that marker, nt/s
// refNt    : reference length (so a long gene still reads as slower, just not 100x slower)
// refRate  : the parameter's own DEFAULT value (that setting lands on PACE.MID)
// gapSim   : model seconds between two visible events (optional; keeps the lane from going empty)
// fastRate : the fastest marker in the same picture (ribosomes elongating through a uORF, say),
//            so a second control cannot blur it even when `rate` is slow
// mul      : the user's relative animation-speed control, 1 = the watchable pace
// blurSoft : the fast stretch is brief (a bead dashing through the exons of a log-compressed
//            axis), so meet the blur ceiling half way instead of head on — a hard cap there is
//            exactly proportional to 1/rate and would flatten the slow<->fast ratio to nothing
export function solve({ nt, rate, refNt, refRate, gapSim = 0, mul = 1, fastRate = 0, blurSoft = false }) {
  const n = Math.max(1, nt || 1), r = Math.max(1e-6, rate || 1);
  const crossSim = n / r;
  const want = PACE.MID
    * Math.pow(n / Math.max(1, refNt || n), PACE.LGAMMA)
    * Math.pow(Math.max(1e-6, refRate || r) / r, PACE.GAMMA);
  const ideal = crossSim / clamp(want, PACE.MIN, PACE.MAX);
  let ts = ideal * (mul > 0 ? mul : 1);
  const empty = gapSim > 0 && gapSim / ts > PACE.MAX_GAP;
  if (empty) ts = Math.max(ts, gapSim / PACE.MAX_GAP);         // lane would be empty -> speed up
  const lo = crossSim / PACE.HARD_MAX, hi = crossSim / PACE.HARD_MIN;
  let out = clamp(ts, lo, hi);
  // hard blur ceiling, highest priority: nothing on screen may cross more than MAX_WPS lane
  // widths per second, whichever control made it fast
  const blur = PACE.MAX_WPS * n / Math.max(r, fastRate || 0);
  if (blur < out) out = blurSoft ? Math.sqrt(out * blur) : blur;
  // hard freeze floor, last word: the marker must still visibly move
  out = Math.max(out, PACE.MIN_STEP * crossSim / PACE.MIN_STEP_S);
  return {
    ts: out,
    crossSim,
    crossWall: crossSim / out,
    gapWall: gapSim > 0 ? gapSim / out : 0,
    // which bound is biting, for the honest on-screen note.
    // 'slowest' = the real kinetics would have been slower still and we sped the picture up;
    // 'fastest' = they would have been faster and we held the picture back so it does not blur.
    limit: out < ts * 0.999 ? 'fastest' : out > ts * 1.001 ? 'slowest'
      : want > PACE.MAX ? 'slowest' : want < PACE.MIN ? 'fastest' : null,
    emptyBump: empty,
  };
}

// "1 s on screen = 12 s in the cell", or the other way round when the picture is slowed down.
export function label(ts) {
  const f = (x) => (x >= 100 ? Math.round(x).toLocaleString() : x >= 10 ? x.toFixed(0) : x.toFixed(x >= 1 ? 1 : 2));
  if (!(ts > 0)) return '';
  if (ts >= 1.02) return `1 s on screen = ${f(ts)} s in the cell`;
  if (ts <= 0.98) return `1 s in the cell is stretched over ${f(1 / ts)} s on screen`;
  return 'on screen in real time (1 s = 1 s)';
}

// The honesty line: how much of the real slow<->fast difference survives on screen.
// lo/hi are two solve() results for the same gene at the two ends of one control.
export function ratioNote(name, lo, hi) {
  const real = lo.crossSim / hi.crossSim;                 // slowest / fastest, in the cell
  const screen = lo.crossWall / hi.crossWall;             // slowest / fastest, on screen
  if (!(real > 1.05)) return '';
  const f = (x) => (x >= 10 ? x.toFixed(0) : x.toFixed(1));
  return `${name} slow↔fast is ×${f(real)} in the cell, shown as ×${f(screen)} on screen`
    + (screen < real / 1.15 ? ' (ratio compressed so neither end freezes or blurs)' : '');
}
