# UTR Race — Sundai Hack 150 (3D worlds), 2026-11-22

Three.js app comparing alternative 5' UTR isoforms: a ribosome-scanning sim run on real UTRs (folded with ViennaRNA), side by side.

- `pipeline/` Python precompute (env `utrrace`) -> `data/utr_data.json` -> copied to `app/public/utr_data.json`. Contract: `data/SCHEMA.md`.
- `app/` Vite vanilla + three. Node lives in `~/miniconda3/envs/web/bin` (add to PATH). `npm run dev` / `npm run build`.

Rules (adapted from the hack's starter system prompt):
- Vanilla JS, ES modules, `import * as THREE from 'three'`, addons from `three/addons/...`. No React/TS.
- Every change leaves the app running. Read the browser console error first when something breaks.
- InstancedMesh for >100 copies; target 60 fps on a laptop; cap pixel ratio at 2.
- Mouse + keyboard; a one-line controls hint overlay.
- Limited palette, one signature effect.
- Science: no invented numbers. Every default parameter shown in the UI is either sourced (citation in `src/params.js`) or labelled "assumed".
