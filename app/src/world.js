import * as THREE from 'three';
import { CDS_STUB } from './sim.js';
import { motifsOf, motifStyle } from './motifs.js';

export const COLORS = {
  bg: 0x0b0d12,
  strand: 0xe9e4d8,
  pair: 0x4fc3f7,
  uorf: 0xffb547,
  overlap: 0xff6b4a,
  cds: 0xe85d9c,
  riboTC: 0xffffff,
  riboNoTC: 0x7c8494,
  riboUorf: 0xffb547,
  riboCds: 0xe85d9c,
};

const WIDTH = 14;
const MAX_RIBO = 400;

// One visual lane per isoform.
export class LaneView {
  constructor(scene, y, maxLen) {
    this.group = new THREE.Group();
    this.group.position.y = y;
    scene.add(this.group);
    this.scale = WIDTH / (maxLen + CDS_STUB);
    this.x0 = -WIDTH / 2;
    this.pickables = [];
    this.stallFlash = new Map();

    const rg = new THREE.SphereGeometry(1, 16, 12);
    rg.scale(0.256, 0.208, 0.208);
    // unlit HDR instance colours -> ribosomes are the brightest things in the scene and catch the bloom
    this.riboMesh = new THREE.InstancedMesh(rg, new THREE.MeshBasicMaterial({ color: 0xffffff }), MAX_RIBO);
    this.riboMesh.count = 0;
    this.riboMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.group.add(this.riboMesh);

    this.bursts = []; // protein particles
    const bg = new THREE.OctahedronGeometry(0.1);
    this.burstMesh = new THREE.InstancedMesh(bg, new THREE.MeshBasicMaterial({ color: new THREE.Color(COLORS.cds).multiplyScalar(2.2) }), 300);
    this.burstMesh.count = 0;
    this.group.add(this.burstMesh);
    this.content = new THREE.Group();
    this.group.add(this.content);
  }

  X(pos) { return this.x0 + pos * this.scale; }

  build(lane) {
    this.lane = lane;
    this.content.clear();
    this.pickables = [];
    const iso = lane.iso;
    const L = iso.len;

    // strand: one segment per nucleotide, colored by base
    const pts = [];
    for (let i = 0; i <= L; i++) pts.push(new THREE.Vector3(this.X(i), 0, 0));
    const strand = new THREE.Mesh(
      new THREE.TubeGeometry(new THREE.CatmullRomCurve3([pts[0], pts[pts.length - 1]]), 8, 0.025, 8),
      new THREE.MeshStandardMaterial({ color: COLORS.strand, roughness: 0.6 }));
    this.content.add(strand);

    // nucleotide ticks (instanced) if not too dense
    if (L <= 1500) {
      const g = new THREE.BoxGeometry(Math.max(this.scale * 0.6, 0.004), 0.09, 0.04);
      const m = new THREE.InstancedMesh(g, new THREE.MeshStandardMaterial({ roughness: 0.7 }), L);
      const base = { A: 0x9ad18b, C: 0x7fb8ff, G: 0xf2d06b, U: 0xff8f8f };
      const mat = new THREE.Matrix4(); const c = new THREE.Color();
      for (let i = 0; i < L; i++) {
        mat.makeTranslation(this.X(i + 0.5), 0, 0);
        m.setMatrixAt(i, mat);
        m.setColorAt(i, c.set(base[iso.utr5[i]] ?? 0x888888).multiplyScalar(0.55));
      }
      this.content.add(m);
    }

    // CDS stub
    const cds = new THREE.Mesh(new THREE.BoxGeometry(CDS_STUB * this.scale, 0.14, 0.14),
      new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.9 }));
    this.cds = cds; this.cdsFlash = 0;
    cds.position.x = this.X(L + CDS_STUB / 2);
    this.content.add(cds);
    // main AUG marker
    const aug = new THREE.Mesh(new THREE.ConeGeometry(0.08, 0.22, 12),
      new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.6 }));
    aug.position.set(this.X(L), -0.28, 0); aug.rotation.z = 0;
    this.content.add(aug);

    // base-pair arcs (arc diagram above the strand), grouped per hairpin for toggling
    this.arcGroups = [];
    const pairs = parsePairs(iso.fold?.dotbracket || '');
    const hpOf = (i, j) => lane.hairpins.findIndex(h => i >= h.start && j < h.end + 1);
    const buckets = new Map();
    for (const [i, j] of pairs) {
      const k = hpOf(i, j);
      if (!buckets.has(k)) buckets.set(k, []);
      buckets.get(k).push([i, j]);
    }
    for (const [k, list] of buckets) {
      const verts = [];
      for (const [i, j] of list) {
        const x1 = this.X(i + 0.5), x2 = this.X(j + 0.5);
        const r = (x2 - x1) / 2, cx = (x1 + x2) / 2;
        const h = Math.min(r, 1.0);
        const N = 14;
        for (let s = 0; s < N; s++) {
          const a1 = Math.PI * s / N, a2 = Math.PI * (s + 1) / N;
          verts.push(cx - r * Math.cos(a1), 0.05 + h * Math.sin(a1), 0, cx - r * Math.cos(a2), 0.05 + h * Math.sin(a2), 0);
        }
      }
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(verts, 3));
      // 28 entries in this data set have dG >= 0: not favourable structure, so they are not drawn
      // as structure (the model still stalls on them - see the hover tooltip).
      const unfav = k >= 0 && !(lane.hairpins[k]?.dG < 0);
      const mat = new THREE.LineBasicMaterial({ color: COLORS.pair, transparent: true, opacity: k >= 0 ? (unfav ? 0.08 : 0.35) : 0.14 });
      const ls = new THREE.LineSegments(g, mat);
      this.content.add(ls);
      if (k >= 0) {
        this.arcGroups[k] = ls;
        const h = lane.hairpins[k];
        const w = (h.end - h.start) * this.scale;
        const hit = new THREE.Mesh(new THREE.PlaneGeometry(Math.max(w, 0.15), Math.min(w / 2, 1.0) + 0.1),
          new THREE.MeshBasicMaterial({ visible: false }));
        hit.position.set(this.X((h.start + h.end) / 2), Math.min(w / 2, 1.0) / 2, 0.01);
        hit.userData = { kind: 'hairpin', idx: k, lane: this };
        this.content.add(hit); this.pickables.push(hit);
      }
    }

    // uORF bars below the strand
    this.uorfBars = [];
    lane.uorfs.forEach((u, idx) => {
      const over = u.type === 'overlapping' || u.end > L;
      const len = (Math.min(u.end, L + CDS_STUB) - u.start) * this.scale;
      const row = idx % 3;
      const bar = new THREE.Mesh(new THREE.BoxGeometry(Math.max(len, 0.04), 0.09, 0.09),
        new THREE.MeshStandardMaterial({ color: over ? COLORS.overlap : COLORS.uorf, emissive: over ? COLORS.overlap : COLORS.uorf, emissiveIntensity: 0.3 }));
      bar.position.set(this.X(u.start) + len / 2, -0.22 - row * 0.13, 0);
      bar.userData = { kind: 'uorf', idx, lane: this };
      this.content.add(bar); this.pickables.push(bar); this.uorfBars[idx] = bar;
    });
    // motif spans: thin coloured ticks just under the strand
    for (const mo of motifsOf(iso)) {
      if (mo.type === 'kozak_main') continue;
      const a = Math.max(0, mo.start), b = Math.min(mo.end, L);
      const w = Math.max((b - a) * this.scale, 0.03);
      const sp = new THREE.Mesh(new THREE.BoxGeometry(w, 0.035, 0.03),
        new THREE.MeshBasicMaterial({ color: new THREE.Color(motifStyle(mo.type).color) }));
      sp.position.set(this.X(a) + w / 2, -0.1, 0.02);
      this.content.add(sp);
    }
    this._rowEnd = [];
    this.refreshToggles();
  }

  get titlePos() { return [this.x0, 1.55, 0]; }
  get capPos() { return [this.x0 - 0.08, 0, 0]; }
  get cdsPos() { return [this.X(this.lane.iso.len + CDS_STUB / 2), 0.5, 0]; }
  // chips below the uORF rows; call in order of start; staggers into 3 rows to avoid overlap
  motifAnchor(m) {
    const x = this.X(Math.min(m.start, this.lane.iso.len)), W = 1.15;
    let row = this._rowEnd.findIndex(e => e < x);
    if (row < 0) row = this._rowEnd.length < 3 ? this._rowEnd.length : this._rowEnd.indexOf(Math.min(...this._rowEnd));
    this._rowEnd[row] = x + W;
    return { pos: [x, -0.62 - row * 0.27, 0], center: [0, 0.5], leader: true };
  }

  refreshToggles() {
    const lane = this.lane;
    lane.uorfs.forEach((u, i) => {
      const b = this.uorfBars[i]; if (!b) return;
      b.material.opacity = u.enabled ? 1 : 0.18; b.material.transparent = !u.enabled;
      b.material.emissiveIntensity = u.enabled ? 0.3 : 0;
    });
    lane.hairpins.forEach((h, i) => {
      const a = this.arcGroups[i]; if (!a) return;
      a.material.opacity = h.melted ? 0.06 : (h.dG < 0 ? 0.35 : 0.08);
    });
  }

  update(dt) {
    const lane = this.lane; if (!lane) return;
    const m = new THREE.Matrix4(); const c = new THREE.Color();
    let n = 0;
    // Simple mode: cap how many ribosomes are DRAWN (view only; the model runs all of them).
    // The stride only grows quickly and relaxes slowly, so beads do not flicker in and out.
    const cap = this.maxDraw ?? Infinity, cnt = (this.lane?.ribos ?? []).length;
    let stride = this._stride ?? 1;
    if (cnt / stride > cap) stride *= 2; else if (stride > 1 && cnt / (stride / 2) <= cap * 0.7) stride /= 2;
    if (!(cap < Infinity)) stride = 1;
    this._stride = stride;
    for (const r of lane.ribos) {
      if (n >= MAX_RIBO) break;
      if (stride > 1 && Math.floor(r.id * 1e6) % stride) continue;
      const wobble = r.mode === 'scan' && r.waitHp ? Math.sin(performance.now() * 0.03 + r.id * 50) * 0.02 : 0;
      m.makeTranslation(this.X(r.pos) + wobble, 0.0, 0.0);
      this.riboMesh.setMatrixAt(n, m);
      const col = r.mode === 'cds' ? COLORS.riboCds : r.mode === 'uorf' ? COLORS.riboUorf : r.hasTC ? COLORS.riboTC : COLORS.riboNoTC;
      this.riboMesh.setColorAt(n, c.set(col).multiplyScalar(r.mode === 'scan' ? (r.hasTC ? 1.25 : 0.9) : 1.9));
      n++;
    }
    this.riboMesh.count = n;
    this.riboMesh.instanceMatrix.needsUpdate = true;
    if (this.riboMesh.instanceColor) this.riboMesh.instanceColor.needsUpdate = true;

    // events -> visuals
    for (const e of lane.events) {
      if (e.kind === 'protein') {
        this.cdsFlash = 1;
        for (let k = 0; k < 7; k++) { const a = Math.random() * Math.PI * 2, v = 0.6 + Math.random() * 1.1; this.bursts.push({ x: this.X(lane.iso.len + CDS_STUB), y: 0, vx: Math.abs(Math.cos(a)) * v + 0.3, vy: Math.sin(a) * v, life: 2.6 }); }
      } else if (e.kind === 'stall') {
        this.stallFlash.set(e.h.idx, 0.5);
      }
    }
    lane.events.length = 0;
    for (const [k, t] of this.stallFlash) {
      const a = this.arcGroups[k];
      const left = t - dt;
      if (a && !lane.hairpins[k].melted) a.material.color.set(left > 0 ? 0xffffff : COLORS.pair);
      if (left <= 0) this.stallFlash.delete(k); else this.stallFlash.set(k, left);
    }
    this.cdsFlash = Math.max(0, this.cdsFlash - dt * 2.5);
    if (this.cds) this.cds.material.emissiveIntensity = 0.9 + 4 * this.cdsFlash;
    let bn = 0;
    this.bursts = this.bursts.filter(b => (b.life -= dt) > 0);
    for (const b of this.bursts) {
      if (bn >= 300) break;
      b.x += b.vx * dt; b.y += b.vy * dt;
      m.makeRotationZ(b.life * 3).setPosition(b.x, b.y, 0);
      m.scale(new THREE.Vector3(1, 1, 1).multiplyScalar(Math.min(1, b.life / 1.2)));
      this.burstMesh.setMatrixAt(bn++, m);
    }
    this.burstMesh.count = bn;
    this.burstMesh.instanceMatrix.needsUpdate = true;
  }

  dispose() { this.group.parent?.remove(this.group); }
}

export function parsePairs(db) {
  const st = [], out = [];
  for (let i = 0; i < db.length; i++) {
    if (db[i] === '(') st.push(i);
    else if (db[i] === ')') out.push([st.pop(), i]);
  }
  return out;
}
