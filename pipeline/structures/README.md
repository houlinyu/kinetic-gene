# pipeline/structures — real cryo-EM ribosome geometry for the 3D view

Turns deposited cryo-EM models of the **human 48S translation pre-initiation complex**
into small JSON assets that `app/` can draw at true physical scale (Angstroms).

Nothing here is modelled, idealised or invented. Every coordinate comes out of a
deposited mmCIF; every dimension is computed from those coordinates; every literature
number carries a citation and a URL that was actually opened.

## Outputs (owned by this directory)

| file | what |
| --- | --- |
| `app/public/struct/43s.json` | **primary** asset — PDB **8PJ1**, human 48S in the *open codon-scanning* state, 3.4 Å |
| `app/public/struct/48s_eif4f.json` | secondary asset — PDB **8OZ0**, human 48S with eIF4F/eIF4A at the mRNA **entry channel**, 3.5 Å |
| `app/public/struct/constants.json` | sourced physical constants for scaling (A-form helix, Ribo-seq footprint, 40S size, …) |
| `app/public/struct/README.json` | provenance: full citations, DOIs, URLs, licence |

Other files in `app/public/struct/` (`polii.json`, `spliceosome.json`, …) come from a
different pipeline and are **not** produced or documented here.

## Rerun from scratch

```bash
cd /home/unix/yuhoulin/kinetic-gene
PY=~/miniconda3/envs/utrrace/bin/python          # needs gemmi + numpy
# one-off: nice -n 19 $PY -m pip install gemmi

# 1. download the mmCIFs into the untrusted download dir
mkdir -p pipeline/structures/dl
nice -n 19 curl -sL -o pipeline/structures/dl/8pj1.cif.gz https://files.rcsb.org/download/8PJ1.cif.gz
nice -n 19 curl -sL -o pipeline/structures/dl/8oz0.cif.gz https://files.rcsb.org/download/8OZ0.cif.gz
nice -n 19 gunzip -kf pipeline/structures/dl/8pj1.cif.gz pipeline/structures/dl/8oz0.cif.gz

# 2. build the assets  (-I keeps python from importing anything out of the download dir)
nice -n 19 $PY -I pipeline/structures/build_struct.py
nice -n 19 $PY -I pipeline/structures/build_struct.py \
    --cif pipeline/structures/dl/8oz0.cif --out app/public/struct/48s_eif4f.json

# 3. check them
nice -n 19 $PY -I pipeline/structures/validate.py      # exit 0 = all good
```

`inspect.py` is a scratch helper that prints every chain of a cif with its entity
description, modelled-residue count and residue range:

```bash
nice -n 19 $PY -I pipeline/structures/inspect.py pipeline/structures/dl/8pj1.cif
```

`dl/` holds downloaded (untrusted) data only — never put scripts in it, and always
run Python with `-I` when pointing it at that directory.

## Asset format (`utr-race/struct-v1`)

Plain JSON, ~220-245 kB each (~50 kB gzipped, which is how Vite will serve it).
A binary format was not worth it at this size.

```jsonc
{
  "format": "utr-race/struct-v1",
  "meta": {
    "pdb_id", "title", "resolution_A", "resolution_method", "method", "organism",
    "citation", "doi", "pubmed", "url", "cif_url",
    "license": "PDB data are freely available",
    "origin":   "...",                      // what the coordinate origin is
    "origin_in_cif_frame_A": [x,y,z],       // so you can go back to the deposited frame
    "orientation": "...",                   // the convention, in words
    "rotation_matrix_row_major": [9 floats],// p_out = R @ (p_cif - origin)
    "parts_mapping": { "rule", "rules", "parts", "parts_absent" },
    "dimensions_nm": { ... },               // see below
    "n_mrna_nt_resolved": 47,
    "caveats": [ ... ]
  },
  "geometry": {
    "traces": [                             // one entry per continuous modelled stretch
      { "chain": "A", "part": "18S_rRNA", "label": "18S rRNA", "atom": "P",
        "step": 1, "res_first": 1, "res_last": 1988, "n": 1500,
        "xyz": [x0,y0,z0, x1,y1,z1, ...] }  // flat polyline, Angstroms, 0.1 A precision
    ],
    "surface_spheres": {                    // cheap blob for InstancedMesh
      "count": 681, "xyzr": [x,y,z,r, ...]  // one sphere per <=10 residues of 18S + 40S proteins
    },
    "mrna_channel": {
      "n_nt_resolved", "res_numbers", "sequence_resolved_5to3",
      "path_xyz_5to3_A", "contour_length_nm", "mean_PP_spacing_A",
      "end_5prime_exit_channel", "end_3prime_entry_channel",
      "p_site_geometric", "start_codon_in_resolved_path",
      "entry_edge_offset_from_p_site_nt"
    }
  },

  // ---- flat compat view for app/src/struct3d.js (redundant with `geometry`) ----
  "parts": [ { "name": "18S_rRNA", "kind": "spheres",      "n": 178,
               "points": [[x,y,z,r], ...] },
             { "name": "mRNA",     "kind": "trace_points", "n": 47,
               "points": [[x,y,z], ...] } ],
  "channel_path": [[x,y,z], ...],   // mRNA phosphates 5'->3'
  "pdb_id": "8PJ1", "citation": "...", "dimensions_nm": { ... }   // aliases
}
```

### Units and frame

* **All `geometry` coordinates and radii are Angstroms.** 1 Å = 0.1 nm. Divide by 10
  for nm. Everything under `meta.dimensions_nm` is already in nm.
* **Origin** = geometric centroid of all modelled non-hydrogen atoms of the 40S core
  (18S rRNA + 40S proteins). Unweighted, i.e. not mass-weighted.
* **Orientation**: `+x` = first principal component of the mRNA phosphate path, signed
  so that 5'→3' is `+x`; `+z` = the part of (tRNAi centroid − origin) perpendicular to
  `+x` (the intersubunit face); `+y = z × x`. The exact matrix is in
  `meta.rotation_matrix_row_major`, row-major, and the transform is
  `p_out = R @ (p_cif − origin)`. Each asset has its **own** frame — the two are not
  superimposed on each other.
* So in the view: **mRNA runs left→right, 5' on the left.** The 40S travels toward `+x`
  while the mRNA is drawn in at its `+x` (3', entry-channel) end.

### Parts

`traces[].part` is one of `18S_rRNA`, `40S_proteins`, `mRNA`, `tRNAi`, `eIF1`, `eIF1A`,
`eIF2`, `eIF3`, `eIF4A`, `eIF4G`, `eIF4E`, `eIF5`, `eIF5B`, `other` — assigned by
substring match on the mmCIF `_entity.pdbx_description` (first rule wins; the full rule
list and the resulting chain→part table are in `meta.parts_mapping`).

`40S_proteins` deliberately also absorbs RACK1, eL41 and the ubiquitin moiety of eS31,
because they are modelled as part of the small subunit in these entries. This is stated
in `meta.dimensions_nm["40S_definition"]`.

### Decimation

Backbone traces only: `P` for RNA, `CA` for protein. Step per part is in `STEP` in
`build_struct.py` (18S/mRNA/tRNAi/eIF1/eIF1A every residue, 40S proteins and eIF2 every
2nd, eIF3 every 3rd). A polyline is cut wherever consecutive kept atoms are further
apart than the gap threshold, so a trace never bridges an unmodelled loop. Each cut
produces a separate `traces[]` entry with the same `chain`.

Surface spheres: windows of up to 10 consecutive residues *within one continuous run*;
centre = mean of the window's non-H atoms, radius = RMS atom spread + 3 Å. The run
splitting matters — without it a window spanning a missing loop produced 65 Å spheres
(caught by `validate.py`).

## Feeding `app/src/struct3d.js` (don't break this)

`struct3d.js` already existed when these assets were built and reads a generic shape,
so the assets carry a redundant flat view rather than the renderer being changed:

* `parts` / `channel_path` / `pdb_id` / `citation` / `dimensions_nm` at the **top level**
  are what it reads. `geometry` is the richer, authoritative view — same frame, same
  Angstroms. Keep both in sync (the build script derives one from the other, and
  `validate.py` checks they agree).
* 40S-core compat points are `[x,y,z,radius]`, so `struct3d.js` picks up a bead radius
  for free. Everything else is `[x,y,z]`.
* **`constants.json` entry order is load-bearing.** `struct3d.js` fills six scale slots
  (`rise`, `bpturn`, `diam`, `ribo`, `chan`, `foot`) by regex over each entry's
  `name` + `label` and takes the **first** match. So:
  `mrna_nt_in_40S_channel` must come before any other entry whose name or label says
  "channel", and the word "footprint" must appear only in the 80S Ribo-seq entries
  (otherwise the "ribosome footprint" readout silently becomes the 40S number).
  `validate.py` asserts all six slots resolve to the intended entry. Current wiring:

  | slot | resolves to | shown as |
  | --- | --- | --- |
  | `rise` | `a_form_rna_rise_per_bp` | 2.8 Å |
  | `bpturn` | `a_form_rna_bp_per_turn` | 11 bp |
  | `diam` | `a_form_rna_helix_diameter` | 19.5 Å |
  | `ribo` | `s40_max_dimension` | 274 Å |
  | `chan` | `mrna_nt_in_40S_channel` | 24 nt |
  | `foot` | `riboseq_footprint_80S_yeast` | 28 nt |

  Note the `diam` slot gives **19.5 Å**, not the 23 Å fallback hard-coded in
  `struct3d.js` — see the rules at the bottom of this file.

## What the primary asset actually shows

PDB **8PJ1** is the *open codon-scanning* state (48S-1) of the human pathway in
Petrychenko et al. 2025. Things worth knowing before writing UI copy:

* 47 contiguous mRNA nucleotides are modelled, author residues **−28 … +19** (0 is
  skipped). The numbering is the depositors' **construct** numbering, **not** relative
  to the AUG.
* The construct's `ACCACC-AUG-G` start codon is at author residues **18–20**; only the
  A (18) and U (19) are modelled, right at the 3'/entry-channel edge. The start codon
  has **not** reached the P site — which is the point of this state.
* The mRNA residue nearest the tRNAi anticodon is **+3** (14.5 Å away; the
  codon–anticodon duplex is not formed in the open state, so this is a geometric proxy).
  The entry edge is therefore **+16 nt** from the P site — matching the published
  toeprint leading edge of +16 to +18 nt (Pisareva & Pisarev 2016).
* Independent geometric support for the channel labels, computed by the build script:
  the 5' end sits beside eIF3c (5.4 Å), eIF3a (7.2 Å) and eIF3d (14.1 Å) — the eIF3
  mRNA-**exit**-channel module; the 3' end sits beside eIF3g (5.2 Å), uS3/RPS3 (11.1 Å)
  and uS5/RPS2 (13.4 Å) — the canonical mRNA-**entry** channel.
* The mRNA is a synthetic in-vitro construct, not a cellular transcript.

Measured dimensions (see `meta.dimensions_nm`):

| quantity | value |
| --- | --- |
| 40S max diameter (backbone trace atoms) | **27.4 nm** |
| 40S extents along its own principal axes | 26.9 × 18.8 × 15.0 nm |
| 40S axis-aligned bbox, all atoms, output frame | 19.6 × 23.6 × 21.0 nm |
| whole 48S bbox, all atoms | 28.9 × 23.8 × 23.9 nm |
| mRNA path contour length (47 nt) | 29.5 nm |
| mRNA end-to-end distance | 12.6 nm |
| mean mRNA P–P spacing | 6.42 Å |

The 27.4 nm pair of atoms lies inside the 18S rRNA itself (human expansion segments),
not in a floppy protein tail. Figures of "~23 nm" quoted elsewhere usually describe the
compact core without expansion segments — don't mix the two.

## Rules this directory follows

* No invented coordinates, no invented numbers. The ~23 Å A-RNA helix diameter from
  textbooks is **omitted** from `constants.json` (listed under `omitted`) because no
  source for it could be opened; only the 19.5 Å phosphate-to-phosphate value is given.
* `validate.py` re-reads the published JSON and checks array lengths, coordinate ranges,
  rotation-matrix orthonormality, sphere radii, mRNA sequence/numbering consistency, and
  that `constants.json` agrees with the measured values in `43s.json`. It exits non-zero
  on any failure. Run it after any change.
