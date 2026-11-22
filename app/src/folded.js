// Folded view: one isoform's 5' UTR drawn as its 2D secondary structure (fold.layout2d), lifted into 3D.
// Same interface as LaneView (group, pickables, build, update, refreshToggles, dispose) so main.js can swap them.
import * as THREE from 'three';
import { CDS_STUB } from './sim.js';
import { COLORS, parsePairs } from './world.js';

const MAX_RIBO = 400;
const BASE = { A: 0x9ad18b, C: 0x7fb8ff, G: 0xf2d06b, U: 0xff8f8f };

export class FoldView {
  // size = edge of the square box (world units) the structure is fitted into
  // box = [width, height] the structure is fitted into (aspect preserved)
  constructor(scene, box) {
    this.group = new THREE.Group();
    scene.add(this.group);
    this.box = box; this.size = Math.min(box[0], box[1]);
    this.pickables = [];
    this.stallFlash = new Map();
    this.content = new THREE.Group();
    this.group.add(this.content);
    this.bursts = [];
    const bg = new THREE.OctahedronGeometry(1);
    this.burstMesh = new THREE.InstancedMesh(bg, new THREE.MeshBasicMaterial({ color: new THREE.Color(COLORS.cds).multiplyScalar(2.2) }), 300);
    this.burstMesh.count = 0;
    this.group.add(this.burstMesh);
    this.riboMesh = new THREE.InstancedMesh(new THREE.SphereGeometry(1, 16, 12), new THREE.MeshBasicMaterial({ color: 0xffffff }), MAX_RIBO);
    this.riboMesh.count = 0;
    this.riboMesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    this.group.add(this.riboMesh);
  }

  // ---- geometry: nucleotide coordinates + start codon + CDS tail ----
  computeCoords(iso) {
    const L = iso.len;
    let xy = iso.fold?.layout2d;
    if (!Array.isArray(xy) || xy.length < L || !xy.every(p => Array.isArray(p) && Number.isFinite(p[0]))) {
      // no layout: fall back to a gentle spiral so the view still works
      xy = Array.from({ length: L }, (_, i) => { const a = Math.sqrt(i) * 1.6; return [Math.cos(a) * Math.sqrt(i), Math.sin(a) * Math.sqrt(i)]; });
    }
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < L; i++) { const [x, y] = xy[i]; x0 = Math.min(x0, x); x1 = Math.max(x1, x); y0 = Math.min(y0, y); y1 = Math.max(y1, y); }
    // leave room for the CDS tail
    const s = Math.min(this.box[0] * 0.82 / Math.max(x1 - x0, 1e-6), this.box[1] * 0.82 / Math.max(y1 - y0, 1e-6));
    const cx = (x0 + x1) / 2, cy = (y0 + y1) / 2;
    const pts = [];
    for (let i = 0; i < L; i++) {
      const x = (xy[i][0] - cx) * s, y = (xy[i][1] - cy) * s;
      pts.push(new THREE.Vector3(x, y, 0.18 * this.size / 6 * Math.sin(i * 0.045) ));
    }
    // typical nt spacing
    const d = []; for (let i = 1; i < L; i++) d.push(pts[i].distanceTo(pts[i - 1]));
    d.sort((a, b) => a - b);
    this.step = d.length ? d[d.length >> 1] : this.size / 30;
    // start codon sits one step beyond the last nt, pointing away from the structure's centre
    const last = pts[L - 1] ?? new THREE.Vector3();
    const cen = new THREE.Vector3(); pts.forEach(p => cen.add(p)); cen.multiplyScalar(1 / Math.max(L, 1));
    const out = last.clone().sub(cen); out.z = 0;
    if (out.lengthSq() < 1e-9) out.set(1, 0, 0);
    out.normalize();
    // blend with local backbone direction so the tail leaves smoothly
    if (L >= 2) { const loc = last.clone().sub(pts[L - 2]); loc.z = 0; if (loc.lengthSq() > 0) out.lerp(loc.normalize(), 0.3).normalize(); }
    this.dir = out; this.cen = cen;
    this.aug = last.clone().addScaledVector(out, Math.max(this.step * 2, this.size * 0.02));
    this.tailLen = this.size * 0.22;
    this.tailEnd = this.aug.clone().addScaledVector(out, this.tailLen);
    this.pts = pts;
  }

  // sim position (nt from cap, may run into the CDS stub) -> world point (in group space)
  P(pos, target = new THREE.Vector3()) {
    const L = this.pts.length;
    if (pos <= 0 || L === 0) return target.copy(this.pts[0] ?? this.aug);
    if (pos < L - 1) { const i = Math.floor(pos), f = pos - i; return target.copy(this.pts[i]).lerp(this.pts[i + 1], f); }
    if (pos < this.len) { const f = (pos - (L - 1)) / Math.max(this.len - (L - 1), 1e-6); return target.copy(this.pts[L - 1]).lerp(this.aug, f); }
    const f = Math.min((pos - this.len) / CDS_STUB, 1.2);
    return target.copy(this.aug).addScaledVector(this.dir, f * this.tailLen);
  }

  build(lane) {
    this.lane = lane;
    this.content.clear();
    this.pickables = [];
    const iso = lane.iso, L = iso.len;
    this.len = L;
    this.computeCoords(iso);
    const pts = this.pts, st = this.step;
    const nr = THREE.MathUtils.clamp(st * 0.42, 0.012, 0.09); // nucleotide radius
    this.nr = nr;
    this.riboR = Math.max(nr * 2.6, 0.085);

    // backbone tube
    const bb = [...pts, this.aug];
    if (bb.length >= 2) {
      const curve = new THREE.CatmullRomCurve3(bb, false, 'centripetal');
      const tube = new THREE.Mesh(new THREE.TubeGeometry(curve, Math.min(bb.length * 3, 6000), nr * 0.38, 6),
        new THREE.MeshStandardMaterial({ color: COLORS.strand, roughness: 0.55, metalness: 0.05 }));
      this.content.add(tube);
    }
    // nucleotides
    const sg = new THREE.SphereGeometry(nr, 10, 8);
    const nm = new THREE.InstancedMesh(sg, new THREE.MeshStandardMaterial({ roughness: 0.5 }), Math.max(L, 1));
    const m4 = new THREE.Matrix4(), c = new THREE.Color();
    for (let i = 0; i < L; i++) {
      m4.makeTranslation(pts[i].x, pts[i].y, pts[i].z); nm.setMatrixAt(i, m4);
      nm.setColorAt(i, c.set(BASE[iso.utr5?.[i]] ?? 0x888888).multiplyScalar(0.8));
    }
    this.content.add(nm);

    // base pairs as rungs, grouped per hairpin (for melt toggling / stall flash)
    this.arcGroups = [];
    const pairs = parsePairs(iso.fold?.dotbracket || '');
    const hpOf = (i, j) => lane.hairpins.findIndex(h => i >= h.start && j < h.end + 1);
    const buckets = new Map();
    for (const [i, j] of pairs) {
      if (i == null || j >= L) continue;
      const k = hpOf(i, j);
      if (!buckets.has(k)) buckets.set(k, []);
      buckets.get(k).push([i, j]);
    }
    for (const [k, list] of buckets) {
      const v = [];
      for (const [i, j] of list) v.push(pts[i].x, pts[i].y, pts[i].z, pts[j].x, pts[j].y, pts[j].z);
      const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(v, 3));
      const ls = new THREE.LineSegments(g, new THREE.LineBasicMaterial({ color: COLORS.pair, transparent: true, opacity: k >= 0 ? 0.75 : 0.35 }));
      this.content.add(ls);
      if (k >= 0) {
        this.arcGroups[k] = ls;
        // invisible pick sphere around the hairpin's nucleotides
        const h = lane.hairpins[k], cen = new THREE.Vector3(); let n = 0, r = 0;
        for (let i = h.start; i <= Math.min(h.end, L - 1); i++) { cen.add(pts[i]); n++; }
        cen.multiplyScalar(1 / Math.max(n, 1));
        for (let i = h.start; i <= Math.min(h.end, L - 1); i++) r = Math.max(r, pts[i].distanceTo(cen));
        const hit = new THREE.Mesh(new THREE.SphereGeometry(Math.max(r * 0.85, nr * 3), 12, 8), new THREE.MeshBasicMaterial({ visible: false }));
        hit.position.copy(cen);
        hit.userData = { kind: 'hairpin', idx: k, lane: this };
        this.content.add(hit); this.pickables.push(hit);
      }
    }

    // uORFs: amber sleeve over that stretch of backbone (pickable)
    this.uorfBars = [];
    lane.uorfs.forEach((u, idx) => {
      const over = u.type === 'overlapping' || u.end > L;
      const a = Math.max(0, u.start), b = Math.min(u.end, L + CDS_STUB);
      const n = Math.max(2, Math.min(200, Math.ceil(b - a)));
      const sp = []; for (let s = 0; s <= n; s++) sp.push(this.P(a + (b - a) * s / n));
      const curve = new THREE.CatmullRomCurve3(sp);
      const col = over ? COLORS.overlap : COLORS.uorf;
      const mesh = new THREE.Mesh(new THREE.TubeGeometry(curve, n * 2, nr * 0.75, 8),
        new THREE.MeshStandardMaterial({ color: col, emissive: col, emissiveIntensity: 0.35 }));
      mesh.userData = { kind: 'uorf', idx, lane: this };
      this.content.add(mesh); this.pickables.push(mesh); this.uorfBars[idx] = mesh;
    });

    // CDS tail (pink) + start codon cone + cap bead
    const tail = new THREE.Mesh(new THREE.TubeGeometry(new THREE.LineCurve3(this.aug, this.tailEnd), 8, Math.max(nr * 0.9, 0.035), 10),
      new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.9 }));
    this.cds = tail; this.cdsFlash = 0;
    this.content.add(tail);
    const cone = new THREE.Mesh(new THREE.ConeGeometry(Math.max(nr * 1.6, 0.06), Math.max(nr * 4, 0.16), 14),
      new THREE.MeshStandardMaterial({ color: COLORS.cds, emissive: COLORS.cds, emissiveIntensity: 0.7 }));
    cone.position.copy(this.aug);
    cone.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), this.dir);
    this.content.add(cone);
    const cap = new THREE.Mesh(new THREE.SphereGeometry(Math.max(nr * 2, 0.07), 16, 12),
      new THREE.MeshStandardMaterial({ color: 0xe9e4d8, emissive: 0xe9e4d8, emissiveIntensity: 0.4 }));
    cap.position.copy(pts[0] ?? this.aug);
    this.content.add(cap);
    this.refreshToggles();
  }

  // label anchors (group space)
  get titlePos() {
    let y = -Infinity; for (const p of this.pts) y = Math.max(y, p.y);
    return [0, Math.max(y, this.tailEnd.y, this.aug.y) + 0.45, 0];
  }
  get capPos() { return this.pts[0].toArray(); }
  get cdsPos() { return this.tailEnd.toArray(); }
  // motif chip anchor: pushed outward from the structure centre
  motifAnchor(m) {
    const p = this.P(Math.min(m.start, this.len - 0.01));
    const o = p.clone().sub(this.cen); o.z = 0;
    if (o.lengthSq() < 1e-9) o.set(0, 1, 0);
    o.normalize();
    return { pos: p.clone().addScaledVector(o, this.nr * 3), center: [o.x < 0 ? 1 : 0, o.y > 0 ? 1 : 0] };
  }

  refreshToggles() {
    const lane = this.lane;
    lane.uorfs.forEach((u, i) => {
      const b = this.uorfBars[i]; if (!b) return;
      b.material.opacity = u.enabled ? 1 : 0.15; b.material.transparent = !u.enabled;
      b.material.emissiveIntensity = u.enabled ? 0.35 : 0;
    });
    lane.hairpins.forEach((h, i) => {
      const a = this.arcGroups[i]; if (!a) return;
      a.material.opacity = h.melted ? 0.08 : 0.75;
    });
  }

  update(dt) {
    const lane = this.lane; if (!lane) return;
    const m = new THREE.Matrix4(), c = new THREE.Color(), v = new THREE.Vector3(), q = new THREE.Quaternion(), s = new THREE.Vector3();
    let n = 0;
    for (const r of lane.ribos) {
      if (n >= MAX_RIBO) break;
      this.P(r.pos, v);
      if (r.mode === 'scan' && r.waitHp) v.x += Math.sin(performance.now() * 0.03 + r.id * 50) * this.nr * 0.6;
      v.z += this.nr * 1.2;
      s.setScalar(this.riboR); m.compose(v, q, s);
      this.riboMesh.setMatrixAt(n, m);
      const col = r.mode === 'cds' ? COLORS.riboCds : r.mode === 'uorf' ? COLORS.riboUorf : r.hasTC ? COLORS.riboTC : COLORS.riboNoTC;
      this.riboMesh.setColorAt(n, c.set(col).multiplyScalar(r.mode === 'scan' ? (r.hasTC ? 1.25 : 0.9) : 1.9));
      n++;
    }
    this.riboMesh.count = n;
    this.riboMesh.instanceMatrix.needsUpdate = true;
    if (this.riboMesh.instanceColor) this.riboMesh.instanceColor.needsUpdate = true;

    for (const e of lane.events) {
      if (e.kind === 'protein') {
        this.cdsFlash = 1;
        for (let k = 0; k < 7; k++) {
          const a = (Math.random() - 0.5) * 2.2, sp = 0.5 + Math.random() * 0.9;
          const d = this.dir.clone().applyAxisAngle(new THREE.Vector3(0, 0, 1), a);
          this.bursts.push({ x: this.tailEnd.x, y: this.tailEnd.y, vx: d.x * sp, vy: d.y * sp, life: 2.4 });
        }
      } else if (e.kind === 'stall') this.stallFlash.set(e.h.idx, 0.5);
    }
    lane.events.length = 0;
    for (const [k, t] of this.stallFlash) {
      const a = this.arcGroups[k], left = t - dt;
      if (a && !lane.hairpins[k].melted) a.material.color.set(left > 0 ? 0xffffff : COLORS.pair);
      if (left <= 0) this.stallFlash.delete(k); else this.stallFlash.set(k, left);
    }
    this.cdsFlash = Math.max(0, this.cdsFlash - dt * 2.5);
    if (this.cds) this.cds.material.emissiveIntensity = 0.9 + 4 * this.cdsFlash;
    let bn = 0;
    this.bursts = this.bursts.filter(b => (b.life -= dt) > 0);
    const bs = Math.max(this.nr * 1.5, 0.06);
    for (const b of this.bursts) {
      if (bn >= 300) break;
      b.x += b.vx * dt; b.y += b.vy * dt;
      m.makeRotationZ(b.life * 3).setPosition(b.x, b.y, this.tailEnd.z);
      m.scale(s.setScalar(bs * Math.min(1, b.life / 1.2)));
      this.burstMesh.setMatrixAt(bn++, m);
    }
    this.burstMesh.count = bn;
    this.burstMesh.instanceMatrix.needsUpdate = true;
  }

  dispose() { this.group.parent?.remove(this.group); }
}
