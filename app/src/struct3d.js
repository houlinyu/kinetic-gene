// "3D structure" view: a molecular movie at true physical scale (world unit = u per Ångström, set from the 2D scaffold).
// - Each stem of the predicted secondary structure becomes an ideal A-form helix (rise / bp per turn / diameter from
//   struct/constants.json when present, else textbook A-form values, labelled as such).
// - Single-stranded segments follow fold.layout2d, bent smoothly onto the helix ends.
// - Driven by the live sim (lane.ribos): duplex melts ahead of each 43S, re-folds behind it; a stalled 43S makes the
//   blocking stem strain and its junction pairs flicker; breathing amplitude ∝ fold.punp.
// - The 43S is drawn from struct/43s.json (PDB-derived coarse-grained points) when present, else an ellipsoid blob.
// NOT a predicted tertiary structure.
import * as THREE from 'three';
import { CDS_STUB } from './sim.js';
import { COLORS, parsePairs } from './world.js';
import { FoldView } from './folded.js';

const MAX_RIBO = 60;
const BASE = { A: 0x9ad18b, C: 0x7fb8ff, G: 0xf2d06b, U: 0xff8f8f };
export const PAIRCOL = { GC: 0x4fc3f7, AU: 0x9ad18b, GU: 0xc39bff };
const Y = new THREE.Vector3(0, 1, 0), Zv = new THREE.Vector3(0, 0, 1);
const pairType = (a, b) => { const s = [a, b].sort().join(''); return s === 'CG' ? 'GC' : s === 'AU' ? 'AU' : s === 'GU' ? 'GU' : 'GC'; };
const clamp01 = (x) => Math.max(0, Math.min(1, x));
const ease = (x) => x * x * (3 - 2 * x);

// ---------------- assets (set by main.js after fetching; all optional) ----------------
export const ASSETS = { constants: null, ribo: null };

// textbook fallbacks — used only when constants.json is missing; the UI says so
const FALLBACK = {
  rise: { value: 2.8, unit: 'Å', label: 'A-form rise per bp', source: 'textbook A-form value (constants file missing)' },
  bpturn: { value: 11, unit: 'bp', label: 'A-form bp per turn', source: 'textbook A-form value (constants file missing)' },
  diam: { value: 19.5, unit: 'Å', label: 'A-form helix diameter (P–P)', source: 'Tian 2023 (phosphate-to-phosphate); fallback — constants file missing' },
  ribo: { value: 250, unit: 'Å', label: '40S max dimension', source: 'approximate (constants file missing)' },
  chan: { value: 12, unit: 'nt', label: 'nt in the 40S mRNA channel', source: 'approximate (constants file missing)' },
  foot: { value: 28, unit: 'nt', label: 'ribosome footprint', source: 'approximate (constants file missing)' },
};
const PATTERNS = {
  rise: /rise/i, bpturn: /(bp|base.?pairs?).*(per|_)turn|per_turn/i, diam: /helix.*diam|diam.*helix|a.?form.*diam/i,
  ribo: /40s.*(max|dim|size|diam)|(max|dim).*40s/i, chan: /channel/i, foot: /footprint/i,
};
function toAngstrom(v, unit, key) {
  if (!Number.isFinite(v)) return null;
  const u = String(unit || '').toLowerCase();
  if (key === 'bpturn' || key === 'chan' || key === 'foot') return v;
  if (u.startsWith('nm')) return v * 10;
  if (u.startsWith('pm')) return v / 100;
  return v; // Å assumed
}
// walk any JSON shape and pick the first numeric entry whose key/label matches
function findConst(obj, key) {
  const re = PATTERNS[key];
  let hit = null;
  const visit = (o, path) => {
    if (hit || o == null) return;
    if (Array.isArray(o)) { o.forEach((x, i) => visit(x, path + '.' + (x?.key ?? x?.name ?? i))); return; }
    if (typeof o !== 'object') return;
    const name = `${path} ${o.key ?? ''} ${o.name ?? ''} ${o.label ?? ''}`;
    const val = Number.isFinite(o.value) ? o.value : Number.isFinite(o.val) ? o.val : null;
    if (val != null && re.test(name)) {
      hit = { value: toAngstrom(val, o.unit ?? o.units, key), unit: key === 'bpturn' ? 'bp' : (key === 'chan' || key === 'foot') ? 'nt' : 'Å',
        label: o.label ?? o.name ?? path.split('.').pop(), source: o.source ?? o.cite ?? o.citation ?? o.ref ?? '', raw: `${val} ${o.unit ?? o.units ?? ''}` };
      return;
    }
    for (const [k, v] of Object.entries(o)) {
      if (Number.isFinite(v) && re.test(path + '.' + k)) {
        const unit = o[k + '_unit'] ?? o.unit ?? (/_nm$/.test(k) ? 'nm' : /_a$|_ang/i.test(k) ? 'Å' : '');
        hit = { value: toAngstrom(v, unit, key), unit: key === 'bpturn' ? 'bp' : 'Å', label: k, source: o.source ?? o[k + '_source'] ?? '', raw: `${v} ${unit}` };
        return;
      }
      visit(v, path + '.' + k);
    }
  };
  visit(ASSETS.constants, '');
  return hit && Number.isFinite(hit.value) ? { ...hit, found: true } : { ...FALLBACK[key], found: false };
}
export function constants() {
  const out = {};
  for (const k of Object.keys(FALLBACK)) out[k] = findConst(ASSETS.constants, k);
  return out;
}

// parse the coarse-grained 43S: returns {pos: Float32Array (Å, centred at the channel/P-site), col, n, chan: [[x,y,z]...], meta}
const PARTCOL = [[/18s|rrna/i, 0x9fb4d8], [/protein|rps|40s/i, 0xb7a6e8], [/trna|met/i, 0x5ee0a0], [/eif|factor/i, 0xffb547], [/mrna/i, 0xe85d9c]];
let RIBO_CACHE = null;
function parseRibo() {
  const d = ASSETS.ribo;
  if (!d) return null;
  if (RIBO_CACHE?.src === d) return RIBO_CACHE.val;
  const pts = [], cols = [];
  const asPts = (a) => (Array.isArray(a) ? a.filter(p => Array.isArray(p) && p.length >= 3 && p.slice(0, 3).every(Number.isFinite)) : []);
  const pickPts = (p) => asPts(p?.points ?? p?.coords ?? p?.spheres ?? p?.xyz ?? p?.centers ?? p?.beads ?? (Array.isArray(p) ? p : null));
  let parts = d.parts ?? d.components ?? null;
  if (parts && !Array.isArray(parts)) parts = Object.entries(parts).map(([name, v]) => ({ name, ...(Array.isArray(v) ? { points: v } : v) }));
  if (!parts && (d.points || d.coords)) parts = [{ name: 'all', points: d.points ?? d.coords }];
  if (!parts) return null;
  let chan = asPts(d.mrna_channel ?? d.channel_path ?? d.channel ?? d.mrna_path ?? d.mRNA_channel?.points ?? d.channel?.points);
  for (const p of parts) {
    const name = String(p.name ?? p.id ?? p.part ?? '');
    const ps = pickPts(p);
    if (/mrna/i.test(name) && !chan.length) chan = ps.map(q => q.slice(0, 3));
    const c = new THREE.Color((PARTCOL.find(([re]) => re.test(name)) ?? [null, 0xc8c8d8])[1]);
    for (const q of ps) { pts.push(q); cols.push(c); }
  }
  if (pts.length < 10) return null;
  // decimate to <= 5000 points
  const stride = Math.max(1, Math.ceil(pts.length / 5000));
  const sel = pts.filter((_, i) => i % stride === 0), sc = cols.filter((_, i) => i % stride === 0);
  // centre: channel midpoint if known, else centroid
  const cen = new THREE.Vector3();
  if (chan.length) chan.forEach(q => cen.add(new THREE.Vector3(q[0], q[1], q[2]))), cen.multiplyScalar(1 / chan.length);
  else sel.forEach(q => cen.add(new THREE.Vector3(q[0], q[1], q[2]))), cen.multiplyScalar(1 / sel.length);
  const pos = new Float32Array(sel.length * 3), col = new Float32Array(sel.length * 3);
  let rad = 0, maxR = 0;
  sel.forEach((q, i) => {
    pos[3 * i] = q[0] - cen.x; pos[3 * i + 1] = q[1] - cen.y; pos[3 * i + 2] = q[2] - cen.z;
    col[3 * i] = sc[i].r; col[3 * i + 1] = sc[i].g; col[3 * i + 2] = sc[i].b;
    rad += Number.isFinite(q[3]) ? q[3] : 0;
    maxR = Math.max(maxR, Math.hypot(pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]));
  });
  const chanDir = chan.length >= 2
    ? new THREE.Vector3(chan[chan.length - 1][0] - chan[0][0], chan[chan.length - 1][1] - chan[0][1], chan[chan.length - 1][2] - chan[0][2]).normalize()
    : new THREE.Vector3(1, 0, 0);
  const val = { pos, col, n: sel.length, chanDir, beadR: rad / sel.length || 0, maxR,
    meta: { pdb: d.pdb_id ?? d.pdb ?? d.meta?.pdb_id ?? d.source?.pdb_id ?? null, citation: d.citation ?? d.meta?.citation ?? d.source?.citation ?? null,
      dims: d.dimensions_nm ?? d.dims_nm ?? d.dimensions ?? null } };
  RIBO_CACHE = { src: d, val };
  return val;
}
export function riboMeta() { return parseRibo()?.meta ?? null; }

export class Struct3DView extends FoldView {
  build(lane) {
    this.lane = lane;
    this.content.clear();
    this.pickables = [];
    const iso = lane.iso, L = iso.len, seq = iso.utr5 || '';
    this.len = L;
    const K = this.K = constants();
    const RISE = K.rise.value, TWIST = 2 * Math.PI / K.bpturn.value, RAD = K.diam.value / 2, PHI = 150 * Math.PI / 180;
    this.computeCoords(iso);
    const p2 = this.pts.map(p => new THREE.Vector3(p.x, p.y, 0));
    const st = this.step, u = this.u = st / 6.3;   // world units per Å: one layout step ≈ one single-strand P–P distance (~6.3 Å)
    const rise = RISE * u, R = RAD * u;
    this.R = R;
    const half = Math.ceil(K.foot.value / 2);
    this.LEAD = half; this.TRAIL = half; this.MELT = 6; this.REFOLD = 10;
    this.CHAN = K.chan.value;

    // ---- stems ----
    const pairs = parsePairs(iso.fold?.dotbracket || '').filter(([i, j]) => i != null && j < L).sort((a, b) => a[0] - b[0]);
    const stems = [];
    for (const [i, j] of pairs) {
      const s = stems[stems.length - 1];
      if (s) { const [li, lj] = s.bp[s.bp.length - 1]; if (i === li + 1 && j === lj - 1) { s.bp.push([i, j]); continue; } }
      stems.push({ bp: [[i, j]] });
    }
    const p3 = p2.map(p => p.clone()), isPaired = new Uint8Array(L);
    const pairInfo = [], pairOf = new Int32Array(L).fill(-1);
    const hpOf = (i, j) => lane.hairpins.findIndex(h => i >= h.start && j < h.end + 1);
    for (const s of stems) {
      const n = s.bp.length, [i0, j0] = s.bp[0], [i1, j1] = s.bp[n - 1];
      const m0 = p2[i0].clone().add(p2[j0]).multiplyScalar(0.5), m1 = p2[i1].clone().add(p2[j1]).multiplyScalar(0.5);
      let ax = m1.clone().sub(m0);
      if (n < 2 || ax.length() < st * 0.3) {
        const inner = new THREE.Vector3(); let c = 0;
        for (let k = i1 + 1; k < j1; k++) { inner.add(p2[k]); c++; }
        ax = c ? inner.multiplyScalar(1 / c).sub(m0) : new THREE.Vector3(-(p2[j0].y - p2[i0].y), p2[j0].x - p2[i0].x, 0);
      }
      if (ax.lengthSq() < 1e-12) ax.set(1, 0, 0);
      ax.normalize();
      const v = new THREE.Vector3().crossVectors(Zv, ax).normalize();
      const side = Math.sign(p2[i0].clone().sub(p2[j0]).dot(v)) || 1;
      const a0 = (Math.PI - PHI) / 2 + (side > 0 ? 0 : Math.PI);
      const mid = m0.clone().add(m1).multiplyScalar(0.5);
      s.bp.forEach(([i, j], b) => {
        const c = mid.clone().addScaledVector(ax, (b - (n - 1) / 2) * rise);
        const ti = a0 + b * TWIST, tj = ti + PHI;
        const ri = v.clone().multiplyScalar(Math.cos(ti)).addScaledVector(Zv, Math.sin(ti));
        const rj = v.clone().multiplyScalar(Math.cos(tj)).addScaledVector(Zv, Math.sin(tj));
        p3[i].copy(c).addScaledVector(ri, R); p3[j].copy(c).addScaledVector(rj, R);
        isPaired[i] = isPaired[j] = 1;
        pairOf[i] = pairOf[j] = pairInfo.length;
        // relaxed single-stranded position: strands swing apart and extend sideways (~1 nt step outward)
        const away = ri.clone().sub(rj).normalize();
        const oi = ri.clone().multiplyScalar(R * 0.9).addScaledVector(away, st * 1.1).addScaledVector(ax, -rise * b * 0.15);
        const oj = rj.clone().multiplyScalar(R * 0.9).addScaledVector(away, -st * 1.1).addScaledVector(ax, -rise * b * 0.15);
        pairInfo.push({ i, j, c, oi, oj, hp: hpOf(i, j), type: pairType(seq[i], seq[j]), open: 0, b, n });
      });
    }
    // single-stranded nts: 2D position + blend of neighbouring anchors' 3D offsets
    let prev = -1;
    const off = (k) => p3[k].clone().sub(p2[k]);
    for (let k = 0; k <= L; k++) {
      if (k < L && !isPaired[k]) continue;
      if (k - prev > 1) {
        const da = prev >= 0 ? off(prev) : null, db = k < L ? off(k) : null;
        for (let t = prev + 1; t < k; t++) {
          const f = (t - prev) / (k - prev);
          p3[t].copy(p2[t]).add(da && db ? da.clone().lerp(db, f) : (da ?? db ?? new THREE.Vector3()));
        }
      }
      prev = k;
    }
    this.base = p3; this.cur = p3.map(p => p.clone()); this.pts = p3;
    this.pairInfo = pairInfo; this.pairOf = pairOf;
    // breathing amplitude per nt ∝ predicted unpaired probability (fold.punp), up to ~2 Å
    const punp = Array.isArray(iso.fold?.punp) && iso.fold.punp.length >= L ? iso.fold.punp : null;
    this.breathe = L <= 800 && punp ? Float32Array.from({ length: L }, (_, k) => (Number(punp[k]) || 0) * 2.2 * u) : null;
    this.phase = Float32Array.from({ length: L }, (_, k) => (k * 2.399) % (2 * Math.PI));

    const cen = new THREE.Vector3(); p3.forEach(p => cen.add(p)); cen.multiplyScalar(1 / Math.max(L, 1));
    const last = p3[L - 1] ?? new THREE.Vector3();
    const out = last.clone().sub(cen); out.z = 0; if (out.lengthSq() < 1e-9) out.set(1, 0, 0); out.normalize();
    this.cen = cen; this.dir = out;
    this.aug = last.clone().addScaledVector(out, st * 2);
    this.tailLen = CDS_STUB * st * 0.15;                // CDS stub drawn compressed (not to scale)
    this.tailEnd = this.aug.clone().addScaledVector(out, this.tailLen);
    const nr = st * 0.28;
    this.nr = nr;

    // ---- instanced meshes ----
    const M = Math.max(L, 1);
    this.ntMesh = new THREE.InstancedMesh(new THREE.SphereGeometry(nr, 8, 6), new THREE.MeshStandardMaterial({ roughness: 0.45 }), M);
    this.segMesh = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 6, 1, true), new THREE.MeshStandardMaterial({ roughness: 0.5, metalness: 0.1 }), M);
    this.rungMesh = new THREE.InstancedMesh(new THREE.CylinderGeometry(1, 1, 1, 6), new THREE.MeshStandardMaterial({ roughness: 0.4 }), Math.max(pairInfo.length, 1));
    for (const im of [this.ntMesh, this.segMesh, this.rungMesh]) { im.instanceMatrix.setUsage(THREE.DynamicDrawUsage); im.frustumCulled = false; this.content.add(im); }
    this.rungMesh.count = pairInfo.length;
    const c = new THREE.Color();
    for (let k = 0; k < L; k++) this.ntMesh.setColorAt(k, c.set(BASE[seq[k]] ?? 0x888888).multiplyScalar(0.85));
    pairInfo.forEach((p, n) => this.rungMesh.setColorAt(n, c.set(PAIRCOL[p.type])));
    this.paintBackbone();

    // uORF pick sleeves (closed model) + hairpin pick spheres + halos
    this.uorfBars = [];
    lane.uorfs.forEach((uo, idx) => {
      const a = Math.max(0, uo.start), b = Math.min(uo.end, L + CDS_STUB), n = Math.max(2, Math.min(120, Math.ceil(b - a)));
      const sp = []; for (let s = 0; s <= n; s++) sp.push(this.P(a + (b - a) * s / n, new THREE.Vector3(), true));
      const mesh = new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3(sp), n, Math.max(nr * 2.5, R * 0.7), 6), new THREE.MeshBasicMaterial({ visible: false }));
      mesh.userData = { kind: 'uorf', idx, lane: this };
      this.content.add(mesh); this.pickables.push(mesh); this.uorfBars[idx] = mesh;
    });
    this.hp = lane.hairpins.map((h, k) => {
      const hc = new THREE.Vector3(); let n = 0, r = 0;
      for (let i = h.start; i <= Math.min(h.end, L - 1); i++) { hc.add(p3[i]); n++; }
      hc.multiplyScalar(1 / Math.max(n, 1));
      for (let i = h.start; i <= Math.min(h.end, L - 1); i++) r = Math.max(r, p3[i].distanceTo(hc));
      const hit = new THREE.Mesh(new THREE.SphereGeometry(Math.max(r * 0.8, R * 1.5), 12, 8), new THREE.MeshBasicMaterial({ visible: false }));
      hit.position.copy(hc); hit.userData = { kind: 'hairpin', idx: k, lane: this };
      this.content.add(hit); this.pickables.push(hit);
      const halo = new THREE.Mesh(new THREE.SphereGeometry(Math.max(r * 0.8, R * 1.8), 20, 14),
        new THREE.MeshBasicMaterial({ color: 0xff2a3d, transparent: true, opacity: 0, depthWrite: false, blending: THREE.AdditiveBlending }));
      halo.position.copy(hc); halo.visible = false; this.content.add(halo);
      return { cen: hc, r, halo, stalled: 0, tag: null };
    });
    this.arcGroups = [];
    this.heli = new THREE.InstancedMesh(new THREE.OctahedronGeometry(Math.max(R * 0.6, nr * 2)), new THREE.MeshBasicMaterial({ color: new THREE.Color(0x5ee0a0).multiplyScalar(2) }), 16);
    this.heli.count = 0; this.content.add(this.heli);

    // CDS tail + start cone + cap
    const tail = new THREE.Mesh(new THREE.TubeGeometry(new THREE.LineCurve3(this.aug, this.tailEnd), 8, nr * 0.8, 10),
      new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.9 }));
    this.cds = tail; this.cdsFlash = 0; this.content.add(tail);
    const cone = new THREE.Mesh(new THREE.ConeGeometry(nr * 2, nr * 5, 14), new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.7 }));
    cone.position.copy(this.aug); cone.quaternion.setFromUnitVectors(Y, this.dir); this.content.add(cone);
    const cap = new THREE.Mesh(new THREE.SphereGeometry(nr * 2.4, 16, 12), new THREE.MeshStandardMaterial({ color: 0xe9e4d8, emissive: 0xe9e4d8, emissiveIntensity: 0.4 }));
    cap.position.copy(p3[0] ?? this.aug); this.content.add(cap);

    // ---- 43S: PDB point cloud (true size) or ellipsoid fallback ----
    const rb = parseRibo();
    this.riboReal = !!rb;
    this.riboObjs = [];
    this.riboQ = new Map();
    if (rb) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(rb.pos, 3));
      g.setAttribute('color', new THREE.Float32BufferAttribute(rb.col, 3));
      const size = Math.max(rb.beadR * 2 || 9, 6) * u;
      for (let k = 0; k < MAX_RIBO; k++) {
        const o = new THREE.Points(g, new THREE.PointsMaterial({ size, vertexColors: true, transparent: true, opacity: 0.38, depthWrite: false, sizeAttenuation: true }));
        o.scale.setScalar(u); o.visible = false; o.renderOrder = 2; o.frustumCulled = false;
        this.content.add(o); this.riboObjs.push(o);
      }
      this.chanDir = rb.chanDir.clone();
    } else {
      const sph = new THREE.SphereGeometry(1, 24, 16);
      const Rr = K.ribo.value / 2 * u;
      for (let k = 0; k < MAX_RIBO; k++) {
        const o = new THREE.Mesh(sph, new THREE.MeshStandardMaterial({ color: 0xdfe6ff, transparent: true, opacity: 0.16, roughness: 0.3, depthWrite: false }));
        o.scale.set(Rr, Rr * 0.8, Rr * 0.8); o.visible = false; o.renderOrder = 2;
        this.content.add(o); this.riboObjs.push(o);
      }
      this.chanDir = new THREE.Vector3(1, 0, 0);
    }
    this.riboMesh.geometry = new THREE.SphereGeometry(1, 12, 10);
    this.coreR = Math.max(R * 1.2, nr * 3);

    // scale bar: 10 nm
    let ymin = Infinity, xmin = Infinity; for (const p of p3) { ymin = Math.min(ymin, p.y); xmin = Math.min(xmin, p.x); }
    const bar = new THREE.Mesh(new THREE.BoxGeometry(100 * u, nr * 0.8, nr * 0.8), new THREE.MeshBasicMaterial({ color: 0xe9e4d8 }));
    bar.position.set(xmin + 50 * u, ymin - st * 4, 0); this.content.add(bar);
    this.scaleBarPos = bar.position.clone();

    this.writeAll(true);
    this.refreshToggles();
  }

  // main.js hooks for CSS2D labels
  attachStallTags(make) {
    this.lane.hairpins.forEach((h, k) => { const o = make(`stalled · hairpin ΔG ${h.dG.toFixed(1)} kcal/mol`, this.hp[k].cen); this.hp[k].tag = o; o.visible = false; });
  }
  get scaleLabel() { return { text: '10 nm', pos: this.scaleBarPos }; }

  paintBackbone() {
    const L = this.len, c = new THREE.Color(), inU = new Array(L).fill(null);
    this.lane.uorfs.forEach(uo => { for (let k = Math.max(0, uo.start); k < Math.min(uo.end, L); k++) inU[k] = uo; });
    for (let k = 0; k < L; k++) {
      const uo = inU[k];
      if (uo && uo.enabled) c.set(uo.type === 'overlapping' || uo.end > L ? COLORS.overlap : COLORS.uorf).multiplyScalar(1.1);
      else c.set(COLORS.strand).multiplyScalar(0.75);
      this.segMesh.setColorAt(k, c);
    }
    if (this.segMesh.instanceColor) this.segMesh.instanceColor.needsUpdate = true;
  }

  P(pos, target = new THREE.Vector3(), closed = false) {
    const pts = closed || !this.cur ? (this.base ?? this.pts) : this.cur, L = pts.length;
    if (L === 0) return target.set(0, 0, 0);
    if (pos <= 0) return target.copy(pts[0]);
    if (pos < L - 1) { const i = Math.floor(pos), f = pos - i; return target.copy(pts[i]).lerp(pts[i + 1], f); }
    if (pos < this.len) { const f = (pos - (L - 1)) / Math.max(this.len - (L - 1), 1e-6); return target.copy(pts[L - 1]).lerp(this.aug, f); }
    const f = Math.min((pos - this.len) / CDS_STUB, 1.2);
    return target.copy(this.aug).addScaledVector(this.dir, f * this.tailLen);
  }
  motifAnchor(m) {
    const p = this.P(Math.min(m.start, this.len - 0.01), new THREE.Vector3(), true);
    const o = p.clone().sub(this.cen); o.z = 0; if (o.lengthSq() < 1e-9) o.set(0, 1, 0); o.normalize();
    return { pos: p.clone().addScaledVector(o, this.R * 1.6), center: [o.x < 0 ? 1 : 0, o.y > 0 ? 1 : 0] };
  }

  writeAll() {
    const L = this.len, m = new THREE.Matrix4(), q = new THREE.Quaternion(), s = new THREE.Vector3(), d = new THREE.Vector3(), mid = new THREE.Vector3();
    const cur = this.cur, br = this.nr * 0.42;
    for (let k = 0; k < L; k++) { m.makeTranslation(cur[k].x, cur[k].y, cur[k].z); this.ntMesh.setMatrixAt(k, m); }
    const seg = (n, a, b) => {
      d.subVectors(b, a); const len = d.length(); mid.addVectors(a, b).multiplyScalar(0.5);
      q.setFromUnitVectors(Y, len > 1e-9 ? d.multiplyScalar(1 / len) : Y);
      m.compose(mid, q, s.set(br, len, br)); this.segMesh.setMatrixAt(n, m);
    };
    for (let k = 0; k < L - 1; k++) seg(k, cur[k], cur[k + 1]);
    if (L) seg(L - 1, cur[L - 1], this.aug);
    this.pairInfo.forEach((p, n) => {
      const a = cur[p.i], b = cur[p.j];
      d.subVectors(b, a); const len = d.length(); mid.addVectors(a, b).multiplyScalar(0.5);
      q.setFromUnitVectors(Y, len > 1e-9 ? d.multiplyScalar(1 / len) : Y);
      const keep = p.open < 0.25 ? 1 - p.open * 4 : 0;      // the H-bond rung breaks first, then strands separate
      m.compose(mid, q, s.set(this.nr * 0.5 * keep + 1e-6, len * 0.85 * keep + 1e-6, this.nr * 0.5 * keep + 1e-6));
      this.rungMesh.setMatrixAt(n, m);
    });
    this.ntMesh.instanceMatrix.needsUpdate = this.segMesh.instanceMatrix.needsUpdate = this.rungMesh.instanceMatrix.needsUpdate = true;
  }

  refreshToggles() { if (this.pairInfo) this.paintBackbone(); }

  update(dt) {
    const lane = this.lane; if (!lane || !this.pairInfo) return;
    const t = performance.now() / 1000, ribos = lane.ribos;
    const stalled = new Set();
    for (const r of ribos) if (r.mode === 'scan' && r.waitHp) stalled.add(r.waitHp.idx);
    const { LEAD, TRAIL, MELT, REFOLD } = this;
    const k = 1 - Math.exp(-dt * 10);
    // 1) base-pair opening targets from the ribosome positions
    for (const p of this.pairInfo) {
      const h = p.hp >= 0 ? lane.hairpins[p.hp] : null;
      let tgt = 0;
      if (h && h.melted) tgt = 1;
      else if (p.hp >= 0 && stalled.has(p.hp)) {
        // blocked: junction pairs flicker open/closed, the rest of the stem holds
        if (p.b < 3) tgt = Math.sin(t * 17 + p.i * 2.1) > 0.35 ? 0.55 : 0;
      } else {
        for (const r of ribos) {
          const ahead = p.i - (r.pos + LEAD);                 // > 0: pair still ahead of the leading edge
          let a = ahead > 0 ? clamp01(1 - ahead / MELT) : 1;   // melts progressively as the edge approaches
          const behind = r.pos - TRAIL - p.j;                  // > 0: trailing edge has passed the pair
          if (behind > 0) a *= clamp01(1 - behind / REFOLD);   // re-folds behind the ribosome
          if (a > tgt) tgt = a;
          if (tgt >= 1) break;
        }
      }
      p.open += (tgt - p.open) * k;
    }
    // 2) positions: closed model + opening + breathing
    const cur = this.cur, base = this.base, L = this.len, B = this.breathe, ph = this.phase;
    for (let n = 0; n < L; n++) {
      cur[n].copy(base[n]);
      if (B) { const a = B[n] * Math.sin(t * 2.3 + ph[n]); cur[n].x += a; cur[n].y += a * 0.6; cur[n].z += B[n] * Math.cos(t * 1.9 + ph[n]); }
    }
    for (const p of this.pairInfo) {
      if (p.open < 1e-3) continue;
      const e = ease(clamp01(p.open));
      cur[p.i].addScaledVector(p.oi, e); cur[p.j].addScaledVector(p.oj, e);
    }
    // strain on a blocked stem: whole stem vibrates
    this.hp.forEach((H, idx) => {
      if (!stalled.has(idx)) return;
      const h = lane.hairpins[idx], amp = this.R * 0.12;
      for (let n = h.start; n <= Math.min(h.end, L - 1); n++) { const w = Math.sin(t * 40 + n); cur[n].x += w * amp; cur[n].z += Math.cos(t * 33 + n) * amp; }
    });
    // 3) thread the mRNA through each 43S channel: nts within ±CHAN/2 of the P-site lie on a straight, extended path
    const tv = new THREE.Vector3(), c0 = new THREE.Vector3(), a1 = new THREE.Vector3(), a2 = new THREE.Vector3();
    const tangents = [];
    for (const r of ribos) {
      if (r.pos >= L + 2) { tangents.push(this.dir); continue; }
      this.P(r.pos, c0, true);
      this.P(Math.max(0, r.pos - 4), a1, true); this.P(r.pos + 4, a2, true);
      tv.subVectors(a2, a1); if (tv.lengthSq() < 1e-12) tv.copy(this.dir); tv.normalize();
      tangents.push(tv.clone());
      const half = this.CHAN / 2, spacing = this.u * 6.3;
      for (let n = Math.max(0, Math.floor(r.pos - half - 3)); n <= Math.min(L - 1, Math.ceil(r.pos + half + 3)); n++) {
        const dn = n - r.pos, w = clamp01(1 - (Math.abs(dn) - half) / 3);
        if (w <= 0) continue;
        cur[n].lerp(a1.copy(c0).addScaledVector(tv, dn * spacing), ease(w));
      }
    }
    this.writeAll();

    // 4) halos, helicase, tags
    const m = new THREE.Matrix4(), v = new THREE.Vector3(), q = new THREE.Quaternion(), s = new THREE.Vector3(), c = new THREE.Color();
    let hn = 0;
    this.hp.forEach((H, idx) => {
      const on = stalled.has(idx);
      H.stalled += ((on ? 1 : 0) - H.stalled) * (1 - Math.exp(-dt * 8));
      H.halo.visible = H.stalled > 0.02;
      H.halo.material.opacity = 0.25 * H.stalled * (0.7 + 0.3 * Math.sin(t * 9));
      if (H.tag) H.tag.visible = on;
      if (on && hn < 16) {
        this.P(lane.hairpins[idx].start + 1, v);
        const a = t * 7 + idx;
        v.x += Math.cos(a) * this.R * 1.2; v.y += Math.sin(a) * this.R * 1.2; v.z += this.R;
        m.makeRotationY(t * 6).setPosition(v); this.heli.setMatrixAt(hn++, m);
      }
    });
    this.heli.count = hn; this.heli.instanceMatrix.needsUpdate = true;

    // 5) ribosomes: 43S body oriented with its mRNA channel along the local backbone + bright P-site core (bloom)
    let n = 0;
    for (let ri = 0; ri < ribos.length && n < MAX_RIBO; ri++) {
      const r = ribos[ri];
      this.P(r.pos, v, true);
      const stall = r.mode === 'scan' && r.waitHp;
      const o = this.riboObjs[n];
      o.visible = true; o.position.copy(v);
      const want = q.setFromUnitVectors(this.chanDir, tangents[ri] ?? this.dir);
      let cq = this.riboQ.get(r.id); if (!cq) { cq = want.clone(); this.riboQ.set(r.id, cq); }
      cq.slerp(want, 1 - Math.exp(-dt * 4));
      o.quaternion.copy(cq);
      const col = stall ? 0xff8090 : r.mode === 'cds' ? COLORS.riboCds : r.mode === 'uorf' ? COLORS.riboUorf : r.hasTC ? 0xffffff : COLORS.riboNoTC;
      o.material.color.set(col);
      m.compose(v, q.identity(), s.setScalar(this.coreR)); this.riboMesh.setMatrixAt(n, m);
      const cc = r.mode === 'cds' ? COLORS.riboCds : r.mode === 'uorf' ? COLORS.riboUorf : r.hasTC ? COLORS.riboTC : COLORS.riboNoTC;
      this.riboMesh.setColorAt(n, c.set(stall ? 0xff4d5e : cc).multiplyScalar(r.mode === 'scan' ? (r.hasTC ? 1.4 : 0.9) : 1.9));
      n++;
    }
    for (let k2 = n; k2 < this.riboObjs.length; k2++) this.riboObjs[k2].visible = false;
    if (this.riboQ.size > 200) { const live = new Set(ribos.map(r => r.id)); for (const id of this.riboQ.keys()) if (!live.has(id)) this.riboQ.delete(id); }
    this.riboMesh.count = n;
    this.riboMesh.instanceMatrix.needsUpdate = true;
    if (this.riboMesh.instanceColor) this.riboMesh.instanceColor.needsUpdate = true;

    // 6) protein bursts
    for (const e of lane.events) if (e.kind === 'protein') {
      this.cdsFlash = 1;
      for (let j = 0; j < 7; j++) {
        const a = (Math.random() - 0.5) * 2.2, sp = 0.5 + Math.random() * 0.9;
        const d = this.dir.clone().applyAxisAngle(Zv, a);
        this.bursts.push({ x: this.tailEnd.x, y: this.tailEnd.y, vx: d.x * sp, vy: d.y * sp, life: 2.4 });
      }
    }
    lane.events.length = 0;
    this.cdsFlash = Math.max(0, this.cdsFlash - dt * 2.5);
    this.cds.material.emissiveIntensity = 0.9 + 4 * this.cdsFlash;
    let bn = 0;
    this.bursts = this.bursts.filter(b => (b.life -= dt) > 0);
    const bs = Math.max(this.nr * 1.5, 0.05);
    for (const b of this.bursts) {
      if (bn >= 300) break;
      b.x += b.vx * dt; b.y += b.vy * dt;
      m.makeRotationZ(b.life * 3).setPosition(b.x, b.y, this.tailEnd.z);
      m.scale(s.setScalar(bs * Math.min(1, b.life / 1.2)));
      this.burstMesh.setMatrixAt(bn++, m);
    }
    this.burstMesh.count = bn; this.burstMesh.instanceMatrix.needsUpdate = true;
  }
}

// panel note (honest label + the physical numbers with sources)
export function structNoteHTML() {
  const K = constants(), rb = riboMeta();
  const esc = (x) => String(x ?? '').replace(/[&<>"]/g, ch => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[ch]));
  const riboTxt = rb ? `48S scanning complex, PDB ${esc(rb.pdb ?? '?')}${rb.citation ? ` (${esc(rb.citation)})` : ''}` : '43S drawn as a simple blob of approximate 40S size (structure file not loaded)';
  const row = (k) => `<span title="${esc(K[k].source || 'no source given')}">${esc(K[k].label)}: <b>${K[k].unit === 'Å' ? (K[k].value / 10).toFixed(K[k].value < 100 ? 2 : 0) + ' nm' : K[k].value + ' ' + K[k].unit}</b>${K[k].found ? '' : '*'}</span>`;
  const anyFallback = Object.values(K).some(x => !x.found);
  return `<b>3D model built from the predicted secondary structure with ideal A-form helices; ${riboTxt}. Not a predicted tertiary structure; dynamics are illustrative, driven by the scanning model.</b><br>
    True scale (bar = 10 nm). Helices melt ahead of each 43S and re-fold behind it; <span class="hot">red</span> = blocked by a stable hairpin (green = eIF4A helicase); breathing ∝ predicted unpaired probability.
    <span class="pk"><i style="background:#4fc3f7"></i>G·C <i style="background:#9ad18b"></i>A·U <i style="background:#c39bff"></i>G·U</span>
    <span class="kn">${['diam', 'rise', 'bpturn', 'ribo', 'chan'].map(row).join(' · ')}${anyFallback ? '<br>* textbook/approximate value; constants file not loaded' : ''}</span>`;
}
