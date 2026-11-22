// "One physics, all findings": scorecard + feasibility map over (λ = tcRate/scanSpeed, β = loadTCexp).
import { score } from './findings.js';
import { SOURCES } from './params.js';

const NX = 40, NY = 30, LMIN = 0.005, LMAX = 1, BMAX = 1.5;
const lamAt = (i) => LMIN * Math.pow(LMAX / LMIN, i / (NX - 1));
const betaAt = (j) => BMAX * j / (NY - 1);
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

// Robustness audit (n = 3000 Latin-hypercube draws over the declared range of the 16 then-exposed
// parameters, run outside the app): how often each claim survives re-sampling. Shown so that
// "5 / 5" is not read as "robust".
const ROBUST = { atf4: 3.5, all: 0.70, n: 3000, nparam: 16, nparamNow: 19 };

// Per-row display overrides where the model agrees with a paper for the wrong reason. The model's
// numbers are not changed — only what the row is allowed to claim.
const ROW_FIX = {
  atf4_u1: {
    // detail from findings.js is `${calm} → ${stress}` with uORF1 removed
    detail: (d) => {
      const m = String(d).match(/([\d.]+)\s*→\s*([\d.]+)/);
      if (!m) return `${d} (stress induction lost)`;
      const a = +m[1], b = +m[2];
      return `with uORF1 removed, output no longer rises under stress: ×${(b / Math.max(a, 1e-12)).toFixed(2)} induction (was ×3.5 with uORF1)`;
    },
    caveat: 'Only the loss of stress induction is claimed here. The model gets the basal direction wrong: removing uORF1 RAISES its basal output, while Vattem & Wek 2004 report that mutating uORF1 severely REDUCES ATF4 expression — uORF1 is a positive element this model does not have.',
  },
};

export class Scorecard {
  constructor(root, { onPick }) {
    this.root = root; this.onPick = onPick;
    root.innerHTML = `
      <h2>One physics, all findings</h2>
      <p class="sub2">Every published finding below must come out of the <b>same</b> scanning model with the <b>same</b> parameters — no per-gene tuning.</p>
      <div id="scHead"></div>
      <ul id="scList"></ul>
      <div class="maphead"><span>Feasibility map</span><span id="mapStat" class="mute"></span></div>
      <div class="mapwrap">
        <canvas id="fmap" width="560" height="330"></canvas>
        <div class="ylab">β · loading ∝ TC<sup>β</sup></div>
      </div>
      <div class="xlab"><span>0.005</span><span>λ = TC re-acquisition / scanning speed (per nt, log)</span><span>1</span></div>
      <div class="maplegend"><span><i class="sw all"></i>all findings hold</span><span><i class="sw part"></i>some fail</span><span><i class="sw dot"></i>current parameters</span></div>
      <p class="fine"><b class="cav">This map varies 2 of the ${ROBUST.nparam} parameters (λ and β); the other ${ROBUST.nparam - 2} are held at our chosen values, so the green area is not a measure of overall robustness.</b>
        Where all findings hold. Note: only the ratio λ is constrained — scanning speed itself (10–20 vs ~100 nt/s, disputed) is not identifiable from these data. Click the map to move there.</p>`;
    this.cv = root.querySelector('#fmap'); this.ctx = this.cv.getContext('2d');
    this.grid = null; this.job = 0;
    this.cv.addEventListener('click', (e) => {
      const r = this.cv.getBoundingClientRect();
      const fx = (e.clientX - r.left) / r.width, fy = (e.clientY - r.top) / r.height;
      const i = Math.round(fx * (NX - 1)), j = Math.round((1 - fy) * (NY - 1));
      this.onPick?.(lamAt(Math.max(0, Math.min(NX - 1, i))), betaAt(Math.max(0, Math.min(NY - 1, j))));
    });
    this.cv.addEventListener('mousemove', (e) => {
      const r = this.cv.getBoundingClientRect();
      const i = Math.round((e.clientX - r.left) / r.width * (NX - 1)), j = Math.round((1 - (e.clientY - r.top) / r.height) * (NY - 1));
      const g = this.grid?.[j * NX + i];
      this.cv.title = `λ = ${lamAt(i).toPrecision(2)} /nt, β = ${betaAt(j).toFixed(2)}${g != null ? ` · ${g} / ${this.total} findings` : ''}`;
    });
  }

  // cheap: list + dot
  update(data, P) {
    this.data = data; this.P = P;
    const s = score(data, P);
    this.total = s.total;
    // Denominator = every direction check that ran, including the ones the model fails. Counting
    // only the 'mechanism' rows would drop a red FAIL row out of the score sitting beside it.
    const tested = s.res.length, ok = s.res.filter(f => f.result.pass).length;
    const failed = s.res.filter(f => !f.result.pass);
    const all = ok === tested;
    this.root.querySelector('#scHead').innerHTML =
      `<div class="schead ${all ? 'ok' : 'bad'}"><b>${ok} of ${tested}</b> direction checks reproduced at our chosen parameter point${
        failed.length ? ` · <span class="scfail">${failed.map(f => esc(f.title.split(/[:(]/)[0].trim())).join('; ')} fail${failed.length > 1 ? '' : 's'}${failed.some(f => f.status === 'missing') ? ' (missing mechanism)' : ''}</span>` : ''}</div>
       <div class="scrob">Not a robustness claim. Re-sampling all ${ROBUST.nparam} UI-exposed parameters over their own declared ranges
         (n = ${ROBUST.n.toLocaleString()}): the ATF4 induction holds in <b>${ROBUST.atf4}%</b> of draws, and <b>all five mechanism findings at once in ${ROBUST.all}%</b> (CHOP is the sixth check and is a known failure, so it is not in that five). The audit sampled ${ROBUST.nparam} parameters; ${ROBUST.nparamNow} are exposed now, so ${ROBUST.nparamNow - ROBUST.nparam} added since are untested.</div>`;
    this.root.querySelector('#scList').innerHTML = s.res.map(f0 => {
      // display-only correction: never let a row claim more than the model actually supports
      const fix = ROW_FIX[f0.id];
      const f = fix ? { ...f0, result: { ...f0.result, detail: fix.detail ? fix.detail(f0.result.detail) : f0.result.detail } } : f0;
      const src = SOURCES[f.cite], ev = this.evidence?.findings?.[f.id];
      const url = ev?.url ?? src?.url, ttl = ev?.quote ? `What the paper says: “${ev.quote}” (${ev.cite ?? f.cite})` : src?.text ?? '';
      const cite = url ? `<a href="${esc(url)}" target="_blank" rel="noopener" title="${esc(ttl)}">${esc(f.cite)}</a>` : `<span title="${esc(ttl)}">${esc(f.cite)}</span>`;
      const cav = (f.status !== 'missing' && f.note ? `<div class="cav">${esc(f.note)}</div>` : '')
        + (fix?.caveat ? `<div class="cav wrong">⚠ ${esc(fix.caveat)}</div>` : '');
      if (f.status === 'missing') {
        return `<li class="missing"><span class="mk">${f.result.pass ? '✓' : '✗'}</span><div><div class="t">${esc(f.title)}</div>
          <div class="d">${esc(f.result.detail)} · ${cite}</div>
          ${f.result.pass ? '' : `<div class="note">✗ missing mechanism: ${esc(f.note)}</div>`}</div></li>`;
      }
      return `<li class="${f.result.pass ? 'pass' : 'fail'}" title="${esc(ttl)}"><span class="mk">${f.result.pass ? '✓' : '✗'}</span><div><div class="t">${esc(f.title)}</div>
        <div class="d">${esc(f.result.detail)} · ${cite}</div>${cav}</div></li>`;
    }).join('');
    this.draw();
    return s;
  }

  // expensive: grid sweep, chunked on idle time; key = params other than tcRate/loadTCexp/tc
  recompute(data, P) {
    const job = ++this.job;
    const base = { ...P };
    const next = new Array(NX * NY).fill(null);
    let k = 0;
    const idle = (f) => (window.requestIdleCallback ? requestIdleCallback(f, { timeout: 40 }) : setTimeout(() => f({ timeRemaining: () => 8 }), 0));
    const work = (dl) => {
      if (job !== this.job) return;
      const t0 = performance.now();
      while (k < NX * NY && (dl.timeRemaining() > 1 || performance.now() - t0 < 8)) {
        const i = k % NX, j = Math.floor(k / NX);
        next[k] = score(data, { ...base, tcRate: lamAt(i) * base.scanSpeed, loadTCexp: betaAt(j) }).passed;
        k++;
        if (performance.now() - t0 > 10) break;
      }
      this.grid = next; this.draw();
      if (k < NX * NY) idle(work); else {
        const nAll = next.filter(v => v === this.total).length;
        this.root.querySelector('#mapStat').textContent = `${Math.round(100 * nAll / next.length)}% of grid satisfies all`;
      }
    };
    this.root.querySelector('#mapStat').textContent = 'computing…';
    idle(work);
  }

  draw() {
    const { ctx, cv } = this, W = cv.width, H = cv.height, g = this.grid, tot = this.total || 1;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = '#12151d'; ctx.fillRect(0, 0, W, H);
    const cw = W / NX, ch = H / NY;
    if (g) for (let j = 0; j < NY; j++) for (let i = 0; i < NX; i++) {
      const v = g[j * NX + i]; if (v == null) continue;
      const x = i * cw, y = H - (j + 1) * ch;
      if (v === tot) ctx.fillStyle = '#5ee0a0';
      else { const f = v / tot; ctx.fillStyle = `rgba(79,195,247,${0.06 + 0.42 * f})`; }
      ctx.fillRect(x + 0.5, y + 0.5, cw - 1, ch - 1);
    }
    // grid ticks at decades of λ
    ctx.strokeStyle = 'rgba(233,228,216,0.12)'; ctx.lineWidth = 1;
    for (const l of [0.01, 0.1]) { const x = Math.log(l / LMIN) / Math.log(LMAX / LMIN) * W; ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke(); }
    // current params
    const P = this.P; if (!P) return;
    const lam = P.tcRate / P.scanSpeed, beta = P.loadTCexp;
    const x = Math.log(Math.max(lam, 1e-9) / LMIN) / Math.log(LMAX / LMIN) * W, y = H - beta / BMAX * H;
    ctx.strokeStyle = '#fff'; ctx.lineWidth = 1.5; ctx.setLineDash([4, 4]);
    ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = '#e85d9c'; ctx.strokeStyle = '#fff'; ctx.lineWidth = 2.5;
    ctx.beginPath(); ctx.arc(x, y, 8, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  }
}
