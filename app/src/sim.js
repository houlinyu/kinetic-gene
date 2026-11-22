// Stochastic scanning model of one 5' UTR ("lane").
// Positions in nt from the cap (0). Main AUG at iso.len. CDS stub drawn after it.
//
// Rules (one sentence each):
//  - 43S complexes load at the cap and scan 5'->3'; they can't overlap (footprint exclusion).
//  - A hairpin stalls the scanner until it is unwound: rate = min(5, exp((dG + helicase*hpFree) / hpKT)).
//  - Loading is TC-limited (rate ∝ tc) and blocked by stable hairpins within capZone nt of the cap.
//  - At an AUG, a scanner carrying eIF2·TC initiates with a Kozak-dependent probability; otherwise it scans past (leaky).
//  - After a short uORF the ribosome may resume scanning without TC (reinitiation) and must re-acquire TC before it can start again.
//  - An overlapping uORF runs past the main AUG, so that ribosome never makes the protein.

export const hpRate = (h, p) => Math.min(5, Math.exp((h.dG + p.helicase * p.hpFree) / p.hpKT));
export function loadRateEff(iso, p, hairpins) {
  let m = 1;
  for (const h of hairpins || iso.hairpins || []) if (!h.melted && h.start < p.capZone) m = Math.min(m, Math.min(1, Math.exp((h.dG + p.capFree) / p.hpKT)));
  return p.loadRate * Math.pow(p.tc, p.loadTCexp ?? 1) * m;
}

export const CDS_STUB = 90; // nt of CDS drawn after the main AUG

export const FATES = ['main', 'uorfDrop', 'overlapLoss', 'leakyEnd'];

export class Lane {
  constructor(iso) {
    this.iso = iso;
    this.uorfs = (iso.uorfs || []).map((u, i) => ({ ...u, idx: i, enabled: true }));
    this.hairpins = (iso.hairpins || []).map((h, i) => ({ ...h, idx: i, melted: false }));
    this.reset();
  }

  reset() {
    this.ribos = [];
    this.t = 0;
    this.loadClock = 0;
    this.proteins = 0;
    this.fates = { main: 0, uorfDrop: 0, overlapLoss: 0, leakyEnd: 0 };
    this.recent = []; // timestamps of protein events (sim seconds)
    this.events = [];  // visual events for the renderer: {kind, pos}
  }

  // Sorted list of things a scanner can hit.
  buildEvents() {
    const ev = [];
    for (const u of this.uorfs) if (u.enabled) ev.push({ pos: u.start, kind: 'aug', u });
    ev.push({ pos: this.iso.len, kind: 'aug', u: null }); // main AUG
    for (const h of this.hairpins) if (!h.melted) ev.push({ pos: h.start, kind: 'hp', h });
    ev.sort((a, b) => a.pos - b.pos || (a.kind === 'hp' ? -1 : 1));
    this.ev = ev;
  }

  kozakP(cls, p) {
    return cls === 'strong' ? p.kozakStrong : cls === 'adequate' ? p.kozakAdeq : p.kozakWeak;
  }

  spawn(p) {
    // cap must be free
    const blocking = this.ribos.some(r => r.pos < p.footprint);
    if (blocking) return;
    // a newly loaded 43S always carries eIF2·TC (loading itself is TC-limited, see step())
    this.ribos.push({
      pos: 0, mode: 'scan', hasTC: true, evIdx: 0,
      waitHp: null, target: 0, done: false, id: Math.random(),
    });
  }

  step(dt, p) {
    if (!this.ev) this.buildEvents();
    this.t += dt;
    // Poisson loading; rate scales with TC availability and drops for cap-proximal hairpins
    if (Math.random() < 1 - Math.exp(-loadRateEff(this.iso, p, this.hairpins) * dt)) this.spawn(p);

    // process front-most first for exclusion
    this.ribos.sort((a, b) => b.pos - a.pos);
    let aheadPos = Infinity;
    const end = this.iso.len + CDS_STUB;
    for (const r of this.ribos) {
      const limit = aheadPos - p.footprint;
      this.advance(r, dt, p, limit);
      if (!r.done) aheadPos = r.pos;
      if (r.pos > end) { r.done = true; }
    }
    this.ribos = this.ribos.filter(r => !r.done);
    const win = 60;
    while (this.recent.length && this.recent[0] < this.t - win) this.recent.shift();
  }

  advance(r, dt, p, limit) {
    if (r.mode === 'scan') {
      if (!r.hasTC && Math.random() < 1 - Math.exp(-p.tcRate * p.tc * dt)) r.hasTC = true;
      if (r.waitHp) {
        const k = hpRate(r.waitHp, p);
        if (Math.random() < 1 - Math.exp(-k * dt)) { r.waitHp = null; r.evIdx++; }
        else return;
      }
      let move = p.scanSpeed * dt;
      while (move > 0) {
        const e = this.ev[r.evIdx];
        const next = e ? e.pos : Infinity;
        const room = Math.min(next, limit) - r.pos;
        if (room <= 0 && limit <= next) return; // blocked by ribosome ahead
        if (move < room) { r.pos += move; return; }
        r.pos += Math.max(room, 0); move -= Math.max(room, 0);
        if (r.pos >= limit - 1e-9 && limit <= next) return;
        // hit event e
        if (e.kind === 'hp') { r.waitHp = e.h; this.events.push({ kind: 'stall', pos: e.pos, h: e.h }); return; }
        r.evIdx++;
        // non-AUG main starts (e.g. CUG) are treated as weak starts
        const nonAUG = !e.u && (this.iso.tags || []).some(t => t.startsWith('non_AUG_start'));
        const cls = e.u ? e.u.kozak_class : nonAUG ? 'weak' : this.iso.main_kozak?.class;
        const capPen = e.pos < p.capAUGnt ? p.capAUG : 1; // AUGs very close to the cap are mostly skipped
        if (r.hasTC && Math.random() < this.kozakP(cls, p) * capPen) {
          if (!e.u) { // main ORF
            r.mode = 'cds'; this.proteins++; this.fates.main++; this.recent.push(this.t);
            this.events.push({ kind: 'protein', pos: e.pos });
          } else {
            r.mode = 'uorf'; r.target = e.u.end; r.u = e.u;
            this.events.push({ kind: 'uorfInit', pos: e.pos });
          }
          return;
        }
        if (!e.u) { // scanned past the main AUG without starting
          this.fates.leakyEnd++; r.done = true; return;
        }
      }
    } else if (r.mode === 'uorf') {
      r.pos = Math.min(r.pos + p.elongSpeed * dt, limit, r.target);
      if (r.pos >= r.target) {
        const u = r.u;
        if (u.type === 'in-frame-extension') { // makes an N-terminally extended protein
          this.proteins++; this.fates.main++; this.recent.push(this.t); r.mode = 'cds'; return;
        }
        if (u.type === 'overlapping' || u.end > this.iso.len) {
          this.fates.overlapLoss++; r.done = true; return;
        }
        const pr = p.reinitMax * Math.exp(-Math.max(0, u.aa) / p.reinitAA);
        if (Math.random() < pr) {
          r.mode = 'scan'; r.hasTC = false;
          // skip events behind us
          while (this.ev[r.evIdx] && this.ev[r.evIdx].pos <= r.pos) r.evIdx++;
        } else { this.fates.uorfDrop++; r.done = true; }
      }
    } else if (r.mode === 'cds') {
      r.pos = Math.min(r.pos + p.elongSpeed * dt, limit === Infinity ? Infinity : Math.max(limit, r.pos));
    }
  }

  ratePerMin() {
    const win = Math.min(60, this.t);
    return win > 1 ? (this.recent.length / win) * 60 : 0;
  }

  setUorf(i, enabled) { this.uorfs[i].enabled = enabled; this.ev = null; this.reset(); }
  setHairpin(i, melted) { this.hairpins[i].melted = melted; this.ev = null; this.reset(); }
}
