// Tab 2, "3D real scale": the same transcription/splicing model drawn at true molecular scale.
// World unit = 1 nm. Everything geometric here is measured, not remembered:
//  - the DNA duplex is built from app/public/struct/txn_constants.json (rise, bp per turn,
//    diameter, bubble size — each measured from the cited structure),
//  - Pol II is the deposited elongation complex (PDB 6GMH) drawn from its surface spheres,
//  - U1 snRNP is PDB 4PJO and the spliceosome is the human C* complex PDB 5XJC,
//  - the nascent-RNA rise per nucleotide is measured from 6GMH's own transcript trace.
// The dynamics (bubble travel, strand folding, snRNP binding order, lariat release) are an
// illustration of the kinetic model, and only the lead polymerase is drawn.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

// Fallbacks, used only when struct/txn_constants.json is missing. Labelled as approximate in the UI.
const FALLBACK = {
  rise:   { value: 0.34, unit: 'nm', label: 'DNA rise per bp', src: 'approximate B-form value (constants file not loaded)', found: false },
  bpturn: { value: 10.5, unit: 'bp', label: 'DNA bp per turn', src: 'approximate B-form value (constants file not loaded)', found: false },
  diam:   { value: 2.0,  unit: 'nm', label: 'DNA duplex diameter', src: 'approximate B-form value (constants file not loaded)', found: false },
  bubble: { value: 11,   unit: 'nt', label: 'Transcription bubble', src: 'approximate (constants file not loaded)', found: false },
  pol:    { value: 24.6, unit: 'nm', label: 'Pol II elongation complex', src: 'approximate (constants file not loaded)', found: false },
  spl:    { value: 36,   unit: 'nm', label: 'Spliceosome', src: 'approximate (constants file not loaded)', found: false },
  u1:     { value: 12.6, unit: 'nm', label: 'U1 snRNP', src: 'approximate (constants file not loaded)', found: false },
  rna:    { value: 0.6,  unit: 'nm', label: 'ssRNA rise per nt', src: 'approximate; no sourced value available', found: false },
};
// exact constant names in struct/txn_constants.json
const WANT = {
  rise: 'dna_duplex_rise_per_bp', bpturn: 'dna_duplex_bp_per_turn',
  diam: 'dna_duplex_diameter_phosphate_convention', bubble: 'transcription_bubble_size',
  pol: 'polii_elongation_complex_max_dimension', spl: 'spliceosome_max_dimension',
  u1: 'u1_snrnp_max_dimension',
};
const toNm = (v, unit) => {
  const u = String(unit ?? '').toLowerCase();
  if (u === 'a' || u === 'å' || u.startsWith('ang')) return v / 10;
  if (u.startsWith('nm')) return v;
  return v;
};

export const ASSETS = { consts: null, polii: null, spliceosome: null };
let LOADED = null;
export function loadTxnAssets() {
  if (LOADED) return LOADED;
  const get = async (u) => { try { const r = await fetch(u); if (!r.ok) return null; const j = await r.json(); return typeof j === 'object' ? j : null; } catch { return null; } };
  LOADED = Promise.all([get('./struct/txn_constants.json'), get('./struct/polii.json'), get('./struct/spliceosome.json')])
    .then(([c, p, s]) => { ASSETS.consts = c; ASSETS.polii = p; ASSETS.spliceosome = s; return ASSETS; });
  return LOADED;
}
function constants() {
  const list = ASSETS.consts?.constants ?? [];
  const by = new Map(list.map(c => [c.name, c]));
  const out = {};
  for (const [k, name] of Object.entries(WANT)) {
    const c = by.get(name);
    out[k] = c && Number.isFinite(c.value)
      ? { value: +toNm(c.value, c.unit).toFixed(3), unit: String(c.unit).includes('/') || String(c.unit) === 'nt' ? c.unit : 'nm',
          label: c.label ?? name, src: c.source?.citation ?? '', url: c.source?.url ?? null, note: c.note ?? '', found: true }
      : { ...FALLBACK[k] };
  }
  // ssRNA rise per nt: measured from the nascent transcript of the Pol II structure itself
  const rna = rnaRise();
  out.rna = rna ?? { ...FALLBACK.rna };
  return out;
}
// mean P–P distance along the modelled nascent RNA in the Pol II asset
function rnaRise() {
  const tr = (ASSETS.polii?.geometry?.traces ?? []).filter(t => /nascent_RNA/i.test(t.part ?? ''));
  let sum = 0, n = 0;
  for (const t of tr) {
    const a = t.xyz ?? [];
    for (let i = 3; i < a.length; i += 3) {
      const d = Math.hypot(a[i] - a[i - 3], a[i + 1] - a[i - 2], a[i + 2] - a[i - 1]);
      if (d > 2 && d < 12) { sum += d; n++; }
    }
  }
  if (n < 3) return null;
  const pdb = ASSETS.polii?.meta?.pdb_id ?? '?';
  return { value: +(sum / n / 10).toFixed(3), unit: 'nm', label: 'nascent RNA rise per nt', found: true,
    src: `measured here from the ${n + 1} modelled transcript phosphates of PDB ${pdb}`, url: null, note: '' };
}
// surface spheres -> {xyzr (nm, centred), n, maxdim, meta}
function spheres(doc, geom) {
  const g = geom ?? doc?.geometry;
  const raw = g?.surface_spheres?.xyzr;
  if (!Array.isArray(raw) || raw.length < 16) return null;
  const n = Math.floor(raw.length / 4);
  const stride = Math.max(1, Math.ceil(n / 1400));
  const pts = [];
  for (let i = 0; i < n; i += stride) pts.push([raw[4 * i] / 10, raw[4 * i + 1] / 10, raw[4 * i + 2] / 10, raw[4 * i + 3] / 10]);
  const cen = [0, 1, 2].map(a => pts.reduce((s, q) => s + q[a], 0) / pts.length);
  let maxdim = 0;
  for (const a of [0, 1, 2]) {
    const v = pts.map(q => q[a]);
    maxdim = Math.max(maxdim, Math.max(...v) - Math.min(...v));
  }
  return { pts: pts.map(q => [q[0] - cen[0], q[1] - cen[1], q[2] - cen[2], q[3]]), n: pts.length, cen, maxdim,
    meta: { pdb: doc?.meta?.pdb_id ?? null, title: doc?.meta?.title ?? null, cite: doc?.meta?.citation ?? doc?.meta?.source?.citation ?? null } };
}
// the asset's own DNA axis, so the deposited complex can be aligned to our helix
function dnaAxis(doc, cen) {
  const tr = (doc?.geometry?.traces ?? []).filter(t => /template_DNA/i.test(t.part ?? ''));
  let best = null;
  for (const t of tr) { const a = t.xyz ?? []; if (a.length >= 6 && (!best || a.length > best.length)) best = a; }
  if (!best) return null;
  const v = new THREE.Vector3(best[best.length - 3] - best[0], best[best.length - 2] - best[1], best[best.length - 1] - best[2]);
  return v.lengthSq() > 1 ? v.normalize() : null;
}
function exitPoint(doc, cen) {
  const tr = (doc?.geometry?.traces ?? []).filter(t => /nascent_RNA/i.test(t.part ?? ''));
  const a = tr[0]?.xyz;
  if (!a || a.length < 3) return null;
  return new THREE.Vector3(a[0] / 10 - cen[0], a[1] / 10 - cen[1], a[2] / 10 - cen[2]);
}

const COL = { back: 0x7c88a3, exon: 0xe9e4d8, intron: 0x8a8f9c, alt: 0xffb547, pol: 0x5ee0a0,
  rna: 0xe85d9c, u1: 0x4fc3f7, u2: 0xc39bff, spl: 0xff8f8f, bg: 0x0b0d12 };
const WIN = 200;           // bp of duplex kept on screen
const RNA_PTS = 240;

export class Splice3DView {
  constructor(host) {
    this.host = host;
    this.renderer = new THREE.WebGLRenderer({ antialias: true });
    this.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
    this.renderer.outputColorSpace = THREE.SRGBColorSpace;
    host.appendChild(this.renderer.domElement);
    this.scene = new THREE.Scene();
    this.scene.background = new THREE.Color(COL.bg);
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.5, 4000);
    this.camera.position.set(6, 10, 165);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true;
    this.controls.minDistance = 20; this.controls.maxDistance = 900;
    this.scene.add(new THREE.HemisphereLight(0xdde6ff, 0x20242e, 1.15));
    const dl = new THREE.DirectionalLight(0xffffff, 1.0); dl.position.set(40, 70, 60); this.scene.add(dl);
    this.world = new THREE.Group(); this.scene.add(this.world);
    this.controls.target.set(0, -16, 0);
    this.ov = document.createElement('div'); this.ov.className = 's3-ov'; host.appendChild(this.ov);
    this.bar = document.createElement('div'); this.bar.className = 's3-bar'; this.bar.innerHTML = '<i></i><span>10 nm</span>';
    host.appendChild(this.bar);
    this.labels = new Map();
    this.tmp = new THREE.Object3D();
    this.K = constants();
    this.buildDNA(); this.buildRNA(); this.buildParticles();
    this.active = false;
    this.resize();
    loadTxnAssets().then(() => {           // real geometry arrives -> rebuild with it
      this.K = constants();
      this.rebuildAll();
      this.onReady?.();
    });
  }
  rebuildAll() {
    this.world.clear();
    this.buildDNA(); this.buildRNA(); this.buildParticles();
    this.n0 = null;
  }
  label(key, text, cls) {
    let d = this.labels.get(key);
    if (!d) { d = document.createElement('div'); d.className = 's3-lab ' + (cls ?? ''); this.ov.appendChild(d); this.labels.set(key, d); }
    if (d.textContent !== text) d.textContent = text;
    return d;
  }

  // ---------- duplex ----------
  buildDNA() {
    const R = this.K.diam.value / 2;
    this.back = new THREE.InstancedMesh(new THREE.SphereGeometry(0.4, 8, 6),
      new THREE.MeshStandardMaterial({ color: COL.back, roughness: .5 }), 2 * WIN);
    const rg = new THREE.CylinderGeometry(0.15, 0.15, Math.max(0.4, 2 * R - 0.55), 5);
    rg.rotateZ(Math.PI / 2);
    this.rung = new THREE.InstancedMesh(rg, new THREE.MeshStandardMaterial({ roughness: .6 }), WIN);
    this.rung.instanceColor = new THREE.InstancedBufferAttribute(new Float32Array(WIN * 3), 3);
    for (const m of [this.back, this.rung]) { m.instanceMatrix.setUsage(THREE.DynamicDrawUsage); m.frustumCulled = false; this.world.add(m); }
    this.n0 = null;
  }
  // ---------- nascent RNA ----------
  buildRNA() {
    const g = new THREE.BufferGeometry();
    this.rnaPos = new Float32Array(RNA_PTS * 3);
    this.rnaCol = new Float32Array(RNA_PTS * 3);
    g.setAttribute('position', new THREE.BufferAttribute(this.rnaPos, 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute('color', new THREE.BufferAttribute(this.rnaCol, 3).setUsage(THREE.DynamicDrawUsage));
    g.setDrawRange(0, 0);
    this.rna = new THREE.Line(g, new THREE.LineBasicMaterial({ vertexColors: true }));
    this.rna.frustumCulled = false; this.world.add(this.rna);
    this.rnaBeads = new THREE.InstancedMesh(new THREE.SphereGeometry(0.35, 7, 5),
      new THREE.MeshStandardMaterial({ color: COL.rna, roughness: .5 }), RNA_PTS);
    this.rnaBeads.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.rnaBeads.count = 0; this.rnaBeads.frustumCulled = false;
    this.world.add(this.rnaBeads);
  }
  // ---------- Pol II, U1, spliceosome ----------
  particle(sph, col, fallbackNm, opacity = 0.9) {
    const g = new THREE.Group();
    if (sph) {
      const m = new THREE.InstancedMesh(new THREE.SphereGeometry(1, 8, 6),
        new THREE.MeshStandardMaterial({ color: col, roughness: .5, transparent: opacity < 1, opacity, flatShading: true }), sph.n);
      sph.pts.forEach((q, i) => {
        this.tmp.position.set(q[0], q[1], q[2]);
        this.tmp.scale.setScalar(Math.max(q[3], 0.2));
        this.tmp.rotation.set(0, 0, 0); this.tmp.updateMatrix();
        m.setMatrixAt(i, this.tmp.matrix);
      });
      m.instanceMatrix.needsUpdate = true; m.frustumCulled = false;
      g.add(m);
      g.userData.real = true; g.userData.meta = sph.meta; g.userData.maxdim = sph.maxdim;
    } else {
      const m = new THREE.Mesh(new THREE.IcosahedronGeometry(fallbackNm / 2, 2),
        new THREE.MeshStandardMaterial({ color: col, roughness: .5, transparent: true, opacity: .7, flatShading: true }));
      g.add(m); g.userData.real = false; g.userData.maxdim = fallbackNm;
    }
    this.world.add(g);
    return g;
  }
  buildParticles() {
    const pol = spheres(ASSETS.polii);
    this.polSph = pol;
    this.polG = this.particle(pol, COL.pol, this.K.pol.value, 0.55);
    // align the deposited complex so its DNA runs along our helix axis (+X)
    if (pol) {
      const ax = dnaAxis(ASSETS.polii, pol.cen);
      if (ax) this.polG.quaternion.setFromUnitVectors(ax, new THREE.Vector3(1, 0, 0));
      const ep = exitPoint(ASSETS.polii, pol.cen);
      this.exitLocal = ep ? ep.clone().applyQuaternion(this.polG.quaternion) : new THREE.Vector3(2, -4, 7);
    } else this.exitLocal = new THREE.Vector3(2, -4, 7);
    const spl = spheres(ASSETS.spliceosome);
    this.splSph = spl;
    this.splG = this.particle(spl, COL.spl, this.K.spl.value, 0.85);
    const comp = (ASSETS.spliceosome?.companion_structures ?? [])[0];
    const u1 = comp ? spheres(comp, comp.geometry) : null;
    this.u1Sph = u1; this.u1Meta = comp?.meta ?? null;
    this.u1G = this.particle(u1, COL.u1, this.K.u1.value, 0.9);
    if (u1) { this.u1G.userData.meta = { pdb: comp.meta?.pdb_id ?? null, title: comp.meta?.title ?? null }; }
    // U2 at the branch point: no structure of U2 snRNP alone here -> an approximate blob
    this.u2G = this.particle(null, COL.u2, 14);
    this.u2G.userData.approx = true;
    this.lariat = new THREE.Mesh(new THREE.TorusGeometry(6, 0.45, 6, 28),
      new THREE.MeshStandardMaterial({ color: COL.intron, roughness: .6, transparent: true }));
    this.world.add(this.lariat);
    for (const g of [this.splG, this.u1G, this.u2G, this.lariat]) g.visible = false;
  }

  build(gene, J, P) {
    this.gene = gene; this.P = P;
    this.J = (J ?? []).slice().sort((a, b) => a.start - b.start);
    const ex = gene.exons ?? [];
    this.s0 = ex.length ? Math.min(...ex.map(e => e.start)) : 0;
    this.s1 = Math.max(...ex.map(e => e.end), ...this.J.map(j => j.end), 1);
    this.n0 = null; this.cutSeen = new Map(); this.lastCut = null;
    this.rna.geometry.setDrawRange(0, 0); this.rnaBeads.count = 0;
  }
  setParams(P) { this.P = P; }
  setActive(on) { this.active = on; if (on) this.resize(); }
  resize() {
    const w = this.host.clientWidth || 800, h = this.host.clientHeight || 600;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h; this.camera.updateProjectionMatrix();
  }
  kindAt(pos) {
    for (const e of this.gene?.exons ?? []) if (pos >= e.start && pos < e.end)
      return e.type && !['constitutive', 'first', 'last'].includes(e.type) ? 'alt' : 'exon';
    return 'intron';
  }

  // duplex window centred on the polymerase; the bubble melts the bp inside Pol II
  dna(pos) {
    const rise = this.K.rise.value, R = this.K.diam.value / 2, tw = 2 * Math.PI / this.K.bpturn.value;
    const bub = this.K.bubble.value;
    const bp = Math.round(pos - this.s0);
    const n0 = Math.max(0, Math.round(bp - WIN / 2));
    if (this.n0 === null || Math.abs(n0 - this.n0) >= 8) {
      this.n0 = n0;
      const c = new THREE.Color();
      for (let i = 0; i < WIN; i++) {
        const k = this.kindAt(this.s0 + n0 + i);
        c.set(k === 'intron' ? COL.intron : k === 'alt' ? COL.alt : COL.exon);
        this.rung.setColorAt(i, c);
      }
      if (this.rung.instanceColor) this.rung.instanceColor.needsUpdate = true;
    }
    const bi = bp - this.n0;
    for (let i = 0; i < WIN; i++) {
      const n = this.n0 + i, th = n * tw, x = i * rise, d = i - bi;
      const open = Math.abs(d) <= bub / 2 ? 1 - Math.abs(d) / (bub / 2 + 0.6) : 0;
      this.tmp.position.set(x, 0, 0); this.tmp.rotation.set(th, 0, 0);
      this.tmp.scale.set(1, 1, Math.max(0.001, 1 - open)); this.tmp.updateMatrix();
      this.rung.setMatrixAt(i, this.tmp.matrix);
      for (const [s, ph] of [[0, 0], [1, 2.45]]) {
        const rr = R * (1 + 0.5 * open);
        this.tmp.position.set(x, rr * Math.cos(th + ph), rr * Math.sin(th + ph));
        this.tmp.rotation.set(0, 0, 0); this.tmp.scale.setScalar(1); this.tmp.updateMatrix();
        this.back.setMatrixAt(2 * i + s, this.tmp.matrix);
      }
    }
    this.rung.instanceMatrix.needsUpdate = true;
    this.back.instanceMatrix.needsUpdate = true;
    const px = bi * rise;
    this.world.position.x = -px;         // slide the duplex so Pol II stays at the origin
    this.polG.position.set(px, 0, 0);
    return px;
  }
  // folded, wiggling single strand out of the exit channel
  rnaPath(nt, px, t) {
    const step = this.K.rna.value;
    const nPts = clamp(Math.round(nt / 8), 2, RNA_PTS), per = Math.max(1, nt / nPts);
    const seg = step * per * 0.4;
    const e = this.exitLocal;
    const ax = px + e.x, ay = e.y - 6, az = e.z + 6;   // the strand coils around this anchor
    const R0 = 16 + Math.min(20, nt / 400);            // nm, how far the coil spreads
    let x = px + e.x, y = e.y, z = e.z, dx = 0.1, dy = -0.5, dz = 0.85;
    const col = new THREE.Color();
    for (let i = 0; i < nPts; i++) {
      const h = Math.sin(i * 12.9898) * 43758.5453;
      const r1 = h - Math.floor(h), r2 = (h * 7.13) % 1, r3 = (h * 3.71) % 1;
      dx = dx * .55 + (r1 - .5) * 1.4; dy = dy * .55 + (r2 - .5) * 1.4; dz = dz * .55 + (r3 - .5) * 1.4;
      // keep it folded near the exit channel instead of running off screen
      const rx = x - ax, ry = y - ay, rz = z - az, rr = Math.hypot(rx, ry, rz);
      if (rr > R0) { const k = 0.9 * (rr - R0) / R0; dx -= rx / rr * k; dy -= ry / rr * k; dz -= rz / rr * k; }
      const L = Math.hypot(dx, dy, dz) || 1; dx /= L; dy /= L; dz /= L;
      const w = Math.sin(t * 2.2 + i * .7) * .5;
      this.rnaPos[3 * i] = x + w * .3; this.rnaPos[3 * i + 1] = y + w; this.rnaPos[3 * i + 2] = z + w * .6;
      const k = this.kindAt(this.s0 + Math.max(0, nt - i * per));
      col.set(k === 'intron' ? COL.intron : k === 'alt' ? COL.alt : COL.rna);
      this.rnaCol[3 * i] = col.r; this.rnaCol[3 * i + 1] = col.g; this.rnaCol[3 * i + 2] = col.b;
      x += dx * seg; y += dy * seg; z += dz * seg;
    }
    this.rna.geometry.setDrawRange(0, nPts);
    this.rna.geometry.attributes.position.needsUpdate = true;
    this.rna.geometry.attributes.color.needsUpdate = true;
    let b = 0;
    for (let i = 0; i < nPts; i += 3) {
      this.tmp.position.set(this.rnaPos[3 * i], this.rnaPos[3 * i + 1], this.rnaPos[3 * i + 2]);
      this.tmp.rotation.set(0, 0, 0); this.tmp.scale.setScalar(1); this.tmp.updateMatrix();
      this.rnaBeads.setMatrixAt(b++, this.tmp.matrix);
    }
    this.rnaBeads.count = b; this.rnaBeads.instanceMatrix.needsUpdate = true;
    this.nPts = nPts; this.per = per;
  }

  update(dt, sim, J) {
    if (!this.active || !this.gene) return;
    if (J && J !== this.J) this.J = J.slice().sort((a, b) => a.start - b.start);
    this.t = (this.t ?? 0) + dt;
    const pols = (typeof sim?.view === 'function' ? sim.view() : sim?.pols) ?? [];
    const p = pols.filter(x => !x.paused).sort((a, b) => (b.pos ?? 0) - (a.pos ?? 0))[0] ?? pols[0] ?? null;
    const t0 = p?.tssPos ?? this.s0;
    const pos = p ? (p.pos ?? this.s0) : this.s0;
    const px = this.dna(pos);
    const nt = Math.max(1, pos - t0);
    this.rnaPath(nt, px, this.t);
    const rem = new Set(p?.removed ?? []), com = new Set(p?.committed ?? []), dead = new Set(p?.skipped ?? []);
    for (const j of this.J ?? []) {
      const key = (p?.id ?? 0) + ':' + j.id;
      if (rem.has(j.id) && !this.cutSeen.get(key)) { this.cutSeen.set(key, true); this.lastCut = { j, t: this.t }; }
    }
    const live = (this.J ?? []).find(j => j.start < pos && !dead.has(j.id) && !rem.has(j.id));
    const onRNA = (genPos) => {
      const i = clamp(Math.round((nt - (genPos - t0)) / Math.max(1, this.per ?? 1)), 0, (this.nPts ?? 1) - 1);
      return new THREE.Vector3(this.rnaPos[3 * i], this.rnaPos[3 * i + 1], this.rnaPos[3 * i + 2]);
    };
    const place = (g, v, k, lift) => {
      g.position.copy(v).add(new THREE.Vector3(Math.sin(this.t * 1.1 + k) * .8, Math.cos(this.t * .9 + k) * .8 - lift, Math.sin(this.t * .7 + k) * .8));
      g.rotation.y = this.t * 0.05 + k;
    };
    const r1 = (this.u1G.userData.maxdim ?? 12) / 2;
    this.u1G.visible = !!live;
    if (live) place(this.u1G, onRNA(live.start), 0, r1 + 2);
    const bp = live ? (live.bp?.pos ?? live.end - 25) : null;
    const haveU2 = !!live && pos > bp;
    this.u2G.visible = haveU2;
    if (haveU2) place(this.u2G, onRNA(bp), 2, 9);
    const haveSpl = !!live && (com.has(live.id) || pos > live.end);
    this.splG.visible = haveSpl;
    if (haveSpl) place(this.splG, onRNA((live.start + live.end) / 2), 4, (this.splG.userData.maxdim ?? 30) / 2 + 4);
    const cut = this.lastCut && this.t - this.lastCut.t < 1.8 ? this.lastCut : null;
    this.lariat.visible = !!cut;
    if (cut) {
      const k = (this.t - cut.t) / 1.8;
      const circ = (cut.j.end - cut.j.start) * this.K.rna.value;      // nt -> nm of released intron
      const r = clamp(circ / (2 * Math.PI) * 0.12, 3, 18);
      this.lariat.geometry.dispose();
      this.lariat.geometry = new THREE.TorusGeometry(r, 0.5, 6, 26);
      this.lariat.position.copy(onRNA((cut.j.start + cut.j.end) / 2)).add(new THREE.Vector3(0, -12 - 40 * k, 10 * k));
      this.lariat.material.opacity = clamp(1 - k, 0, 1);
      this.lariat.rotation.set(this.t * .4, this.t * .3, 0);
    }
    this.controls.target.set(0, -16, 0);
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
    this.overlay(p, live, haveU2, haveSpl, cut);
  }
  overlay(p, live, haveU2, haveSpl, cut) {
    const w = this.host.clientWidth || 1, h = this.host.clientHeight || 1;
    const dist = this.camera.position.distanceTo(this.controls.target);
    const pxPerNm = h / (2 * dist * Math.tan(THREE.MathUtils.degToRad(this.camera.fov / 2)));
    this.bar.firstElementChild.style.width = clamp(10 * pxPerNm, 6, w * .6).toFixed(1) + 'px';
    const put = (key, text, obj, cls) => {
      const d = this.label(key, text, cls);
      if (!obj || !obj.visible) { d.style.display = 'none'; return; }
      const v = obj.getWorldPosition(new THREE.Vector3()).project(this.camera);
      if (v.z > 1 || Math.abs(v.x) > 1.3) { d.style.display = 'none'; return; }
      d.style.display = '';
      d.style.left = ((v.x * .5 + .5) * w).toFixed(0) + 'px';
      d.style.top = ((-v.y * .5 + .5) * h).toFixed(0) + 'px';
    };
    const pdb = (g) => (g.userData.meta?.pdb ? ` · PDB ${g.userData.meta.pdb}` : ' · approximate shape');
    put('pol', (p?.paused ? 'Pol II — paused' : 'Pol II') + pdb(this.polG), this.polG, 'pol');
    put('u1', 'U1 snRNP at the 5′ss' + pdb(this.u1G), this.u1G, 'u1');
    put('u2', 'U2 at the branch point · approximate shape', this.u2G, 'u2');
    put('spl', 'spliceosome' + pdb(this.splG), this.splG, 'spl');
    put('lar', 'intron lariat released', this.lariat, 'lar');
  }
  noteHTML() {
    const K = this.K;
    const nm = (k) => `${K[k].value} ${K[k].unit === 'nm' ? 'nm' : K[k].unit}`;
    const row = (k) => `<span title="${esc((K[k].src || 'no source given') + (K[k].note ? ' — ' + K[k].note : ''))}">${esc(K[k].label)}: <b>${nm(k)}</b>${K[k].found ? '' : '*'}</span>`;
    const asset = (g, what) => g?.userData?.meta?.pdb
      ? `${what} from PDB ${esc(g.userData.meta.pdb)}${g.userData.meta.title ? ` (${esc(String(g.userData.meta.title).slice(0, 60))})` : ''}`
      : `${what} is an approximate shape, not a structure`;
    const anyFb = ['rise', 'bpturn', 'diam', 'bubble', 'rna'].some(k => !K[k].found);
    return `<b>Real: the duplex geometry (${nm('rise')}/bp, ${K.bpturn.value} bp per turn, ${nm('diam')} across) and every size on screen, each measured from the cited structure.
      ${asset(this.polG, 'Pol II')}; ${asset(this.u1G, 'U1 snRNP')}; ${asset(this.splG, 'the spliceosome')}. U2 at the branch point is an approximate shape, not a structure.</b><br>
      Illustrative: the dynamics — bubble travel, strand folding, the order in which the snRNPs bind, lariat release — animate the kinetic model, and only the lead polymerase is drawn.
      True scale (bar = 10 nm).
      <span class="pk"><i style="background:#5ee0a0"></i>Pol II <i style="background:#4fc3f7"></i>U1 <i style="background:#c39bff"></i>U2 <i style="background:#ff8f8f"></i>spliceosome <i style="background:#e9e4d8"></i>exon <i style="background:#8a8f9c"></i>intron</span>
      <span class="kn">${['diam', 'rise', 'bpturn', 'bubble', 'rna', 'pol', 'spl', 'u1'].map(row).join(' · ')}${anyFb ? '<br>* approximate; no measured value available for this one' : ''}</span>`;
  }
  dispose() { this.renderer.dispose(); this.host.innerHTML = ''; }
}
