// Deterministic expected-value version of the scanning model in sim.js (same rules, same params).
// Used for the findings scorecard and parameter sweeps: no noise, ~microseconds per UTR.
//
// Scanning mass is split into a = with eIF2·TC, b = without. We sweep positions 5'->3':
//  - between points, b -> a at rate tcRate*tc over the travel time (d / scanSpeed);
//  - a hairpin adds an exponential stall W with rate k; TC can be acquired while stalled:
//      P(no acquisition during W) = k / (k + r);
//  - at an AUG, fraction pInit of a initiates; overlapping uORFs (and in-frame extensions) are lost,
//    others reinitiate with P = reinitMax*exp(-aa/reinitAA) and re-enter as b at the uORF end;
//  - at the main start, a*pInit makes protein; the rest is "scanned past".
// Output flux J = 1 / (1/loadRateEff + T_block), T_block = longest mean hairpin stall
// (with footprint exclusion one stalled scanner gates everyone behind it).

function kozakP(cls, p) {
  return cls === 'strong' ? p.kozakStrong : cls === 'adequate' ? p.kozakAdeq : p.kozakWeak;
}
import { hpRate, loadRateEff } from './sim.js';

export function expected(iso, p, opts = {}) {
  const disabledU = opts.disabledU || new Set(), melted = opts.melted || new Set();
  const len = iso.len;
  const nonAUG = (iso.tags || []).some(t => t.startsWith('non_AUG_start'));
  const pts = [];
  (iso.uorfs || []).forEach((u, i) => { if (!disabledU.has(i)) pts.push({ pos: u.start, kind: 'aug', u }); });
  (iso.hairpins || []).forEach((h, i) => { if (!melted.has(i)) pts.push({ pos: h.start, kind: 'hp', h }); });
  pts.push({ pos: len, kind: 'main' });
  const pending = []; // {pos, mass} reinitiating ribosomes re-entering as b
  for (const q of pts) if (q.kind === 'aug') q.order = 1; else if (q.kind === 'hp') q.order = 0; else q.order = 2;
  pts.sort((x, y) => x.pos - y.pos || x.order - y.order);

  const r = p.tcRate * p.tc;
  let a = 1, b = 0, pos = 0, lost = { uorfDrop: 0, overlapLoss: 0 }, Tblock = 0;
  const travel = (to) => {
    // absorb pending re-entries that occur before `to`, in order
    pending.sort((x, y) => x.pos - y.pos);
    while (pending.length && pending[0].pos <= to) {
      const e = pending.shift();
      const keep = Math.exp(-r * Math.max(0, e.pos - pos) / p.scanSpeed);
      a += b * (1 - keep); b *= keep; pos = Math.max(pos, e.pos); b += e.mass;
    }
    const keep = Math.exp(-r * Math.max(0, to - pos) / p.scanSpeed);
    a += b * (1 - keep); b *= keep; pos = Math.max(pos, to);
  };
  let main = 0, ext = 0;
  for (const q of pts) {
    travel(q.pos);
    if (q.kind === 'hp') {
      const k = hpRate(q.h, p);
      Tblock = Math.max(Tblock, 1 / Math.max(k, 1e-9));
      const keep = k / (k + r);
      a += b * (1 - keep); b *= keep;
    } else if (q.kind === 'aug') {
      const m = a * kozakP(q.u.kozak_class, p) * (q.pos < p.capAUGnt ? p.capAUG : 1);
      a -= m;
      if (q.u.type === 'in-frame-extension') { ext += m; }
      else if (q.u.type === 'overlapping' || q.u.end > len) { lost.overlapLoss += m; }
      else {
        const pr = p.reinitMax * Math.exp(-Math.max(0, q.u.aa || 0) / p.reinitAA);
        pending.push({ pos: q.u.end, mass: m * pr }); lost.uorfDrop += m * (1 - pr);
      }
    } else {
      main = ext + a * kozakP(nonAUG ? 'weak' : iso.main_kozak?.class, p) * (len < p.capAUGnt ? p.capAUG : 1);
    }
  }
  const hps = (iso.hairpins || []).map((h, i) => ({ ...h, melted: melted.has(i) }));
  const L = loadRateEff(iso, p, hps);
  const J = L > 0 ? 1 / (1 / L + Tblock) : 0;
  return { pMain: main, flux: J, perMin: 60 * J * main, Tblock, lost };
}
