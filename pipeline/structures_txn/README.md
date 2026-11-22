# pipeline/structures_txn — transcription + splicing structures

Real experimental structures, at true physical scale, for the transcription+splicing tab.
Sibling of `pipeline/structures/` (translation side: `43s.json`, `48s_eif4f.json`,
`constants.json`). **This directory never writes those files**, and that pipeline never
writes ours.

Outputs (all in `app/public/struct/`):

| file | what | size |
| --- | --- | --- |
| `polii.json` | mammalian Pol II elongation complex, PDB **6GMH** | 91 KB |
| `spliceosome.json` | human spliceosome **C\*** complex, PDB **5XJC**, + a companion U1 snRNP, PDB **4PJO** | 157 KB |
| `txn_constants.json` | measured physical constants, each with a source | 15 KB |
| `txn_meta.json` | provenance, chosen vs. rejected structures, unverified list | 13 KB |

## Rerun

```bash
cd /home/unix/yuhoulin/kinetic-gene
PY=~/miniconda3/envs/utrrace/bin/python          # gemmi + numpy

# 1. mmCIFs -> pipeline/structures_txn/dl/ (skips files already there, ~24 MB total)
nice -n 19 bash pipeline/structures_txn/fetch.sh

# 2. the two 3D assets
nice -n 19 $PY -I pipeline/structures_txn/build_txn_structures.py

# 3. constants + provenance (reads the assets written in step 2, so run it second)
nice -n 19 $PY -I pipeline/structures_txn/build_txn_constants.py

# 4. check everything
nice -n 19 $PY -I pipeline/structures_txn/validate.py        # exit 0 = all checks passed
```

Helpers:

```bash
# chain / entity listing for any downloaded mmCIF
nice -n 19 $PY -I pipeline/structures/inspect.py pipeline/structures_txn/dl/6gmh.cif

# RCSB metadata for candidate entries (title, method, resolution, organism, entities)
nice -n 19 $PY -I pipeline/structures_txn/meta_probe.py 6GMH 5XJC 9EHZ 5FLM

# the GENCODE exon/intron measurement behind two of the constants
nice -n 19 $PY -I pipeline/structures_txn/gene_features.py \
    pipeline/discovery/dl/gencode.v47.basic.annotation.gtf.gz
```

## Format

`utr-race/struct-v1` — the same contract as `43s.json`, so a loader written for one
works for the other. Each asset is `{format, meta, geometry}`:

- `meta` — PDB provenance (`pdb_id`, `title`, `method`, `resolution_A`,
  `resolution_method`, `organism`, `organism_note`, `state`, `citation`,
  `citation_fields`, `doi`, `pubmed`, `url`, `cif_url`, `license`), plus `units`,
  `origin`, `origin_in_cif_frame_A`, `orientation`, `rotation_matrix_row_major`,
  `parts_mapping`, `dimensions_nm`, `caveats`.
- `geometry.traces[]` — one entry per **continuous modelled stretch**:
  `{chain, part, label, atom, step, res_first, res_last, n, xyz}`, where `xyz` is a flat
  `[x0,y0,z0,x1,…]` polyline in Å at 0.1 Å precision. `atom` is `CA` (protein) or `P`
  (nucleic acid); `step` is the decimation stride, so `step: 3` means every 3rd residue.
- `geometry.surface_spheres` — `{note, stride, count, xyzr}`, a flat `[x,y,z,r]`
  quadruple per sphere, one sphere per ≤10 consecutive residues, centre = mean of the
  window's non-H atoms, radius = RMS atom spread + 3.0 Å. Draw as an `InstancedMesh`.
- one structure-specific path block: `geometry.transcription_bubble` in `polii.json`,
  `geometry.splice_sites` in `spliceosome.json`.
- `spliceosome.json` additionally carries `companion_structures: [ … ]`, a full
  `utr-race/struct-v1` document for PDB 4PJO. **It is a different structure in its own
  frame — never superpose it on the C\* complex.**

### Orientation

Every asset is **translated and rotated only, never scaled**, so every distance in the
files is a true physical distance in Ångströms:

```
p_out = R @ (p_cif - origin)        R = meta.rotation_matrix_row_major (3×3, row-major)
```

- `polii.json` — `origin` = centroid of the 12 Pol II subunits' CA atoms.
  **+x = the direction DNA travels through the enzyme** (downstream duplex midpoint →
  upstream duplex midpoint), **+z** = the RNA exit direction orthogonalised against +x,
  **+y** = z × x. So: DNA enters at −x, the bubble is in the middle, DNA leaves at +x,
  and the transcript goes out towards +z.
- `spliceosome.json` — `origin` = centroid of all backbone trace atoms.
  **+x** = 1st principal component of the pre-mRNA phosphate path signed 5'→3',
  **+z** = towards the branch-point adenosine, **+y** = z × x.
- companion `4PJO` — **+x** = 5'→3' along the bound 5' splice site RNA.

## Parts

`polii.json`: `polii_core` (RPB1–RPB12), `clamp`, `elongation_factors` (DSIF, PAF1
complex, SPT6 — hide this part to show Pol II alone), `template_DNA`,
`nontemplate_DNA`, `nascent_RNA`.

`spliceosome.json`: `snRNA_U2`, `snRNA_U5`, `snRNA_U6`, `exon` (5' exon),
`intron` (lariat), `proteins`. `parts_absent` lists `snRNA_U1`, `snRNA_U4` and
`three_prime_exon` — U1 and U4 are gone **by biology** at the C\* stage, not by disorder.

Companion `4PJO`: `u1_proteins`, `snRNA_U1`, `five_prime_splice_site_RNA`.

## What is measured, and what is not

Every number in these files is either a deposited coordinate or something
`build_txn_structures.py` / `gene_features.py` measured from deposited coordinates. The
derivation is written next to each value in the JSON. Highlights:

- base-pair register of the scaffold: template *i* pairs with non-template *49 − i*,
  found as the antiparallel offset maximising C1′–C1′ contacts ≤ 12 Å. This independently
  reproduces the disordered gap: non-template residues 15–25 are the partners of template
  residues 24–34.
- transcription bubble = **11 nt** (template 24–34), measured, not quoted.
- DNA:RNA hybrid = **10 bp** (template 25–34 / RNA 37–46).
- nascent RNA modelled past the hybrid in the exit channel = **11 nt**, chord 6.3 nm.
- duplex geometry, measured on the 23-bp downstream duplex with a PCA-fitted axis:
  rise **3.47 Å/bp**, twist **36.1°** (**10.0 bp/turn**), P–P diameter **19.0 Å**.
- branch point: **A144**, found as the one nucleotide whose O2′ is 1.63 Å from the
  5′-terminal phosphorus of intron nucleotide +1 — i.e. the 2′–5′ lariat bond itself.
- 5′ splice site: intron +1/+2 read out of the model as **GU**.
- human intron/exon lengths: median **1556 bp** / **133 bp**, from GENCODE v47 basic,
  one Ensembl_canonical transcript per protein-coding gene.

**Not verified** (see `txn_meta.json.unverified` and `txn_constants.json.omitted`): no
*published* bubble size, hybrid length, canonical B-DNA rise/twist/diameter or
nuclease-protected RNA length appears anywhere in these files. Europe PMC's
`fullTextXML` endpoint returned HTTP 500 for every article tried, and the publisher pages
(Cell Press, Springer, PMC) returned 403 or a reCAPTCHA page, so those sentences could
not be read. They are listed as omitted rather than filled in from memory. The one
literature-derived item that *is* in an asset — the `clamp` residue range 1–346 — is
flagged `APPROXIMATE` in `polii.json` with a warning not to print its numbers.

## Why these entries

Full reasoning, plus the entries considered and rejected (9EHZ, 5FLM, 8XSO, 6ICZ, 5MQF)
and why, is in `app/public/struct/txn_meta.json`. In short:

- **6GMH** (Vos et al., *Nature* 560:607, 2018; cryo-EM 3.1 Å) — the only mammalian
  elongation complex combining high resolution with a 48-bp scaffold in which template
  DNA, non-template DNA **and** nascent RNA are all resolved. ⚠️ the polymerase is
  *Sus scrofa*, the elongation factors are human: label it **mammalian Pol II**, never
  "human Pol II". The genuinely human alternative, 9EHZ at 2.6 Å, models only 22 nt of
  template DNA and 11 nt of RNA, so there is almost no duplex to draw.
- **5XJC** (Zhang et al., *Cell* 169:918, 2017; cryo-EM 3.6 Å) — **human**, so no yeast
  substitution was needed. C\* state: 5' exon, intron lariat and branch point all
  modelled, with U2/U5/U6 snRNA present. 6ICZ is sharper (3.0 Å) but post-catalytic.
- **4PJO** (Kondo et al., *eLife* 4:e04986, 2015; X-ray 3.3 Å) — human U1 snRNP with a
  10-mer 5' splice site (5'-AGGUAAGUCC-3') paired to the 5' end of U1 snRNA: the one step
  the C\* complex cannot show. ⚠️ engineered minimal particle, so **do not quote its size
  as the size of U1 snRNP**.

PDB data are freely available (<https://www.rcsb.org/pages/usage-policy>). The short
attribution line for the UI is in `txn_meta.json.attribution_line_short`.
