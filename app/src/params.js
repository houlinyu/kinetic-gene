// Simulation parameters. Each entry: value, range, and provenance.
// provenance: "assumed" until replaced with a sourced value (cite = key into SOURCES).
export const PARAMS = {
  scanSpeed:   { value: 20,    min: 1,    max: 40,  step: 1,    unit: 'nt/s',  label: 'Scanning speed (contested: 10–20 vs ~100)', cite: 'wang2022' },
  loadRate:    { value: 0.5, min: 0.02, max: 1,   step: 0.01, unit: '/s',    label: 'Cap loading rate',      cite: null },
  helicase:    { value: 1.0,  min: 0,    max: 3,   step: 0.05, unit: '×',     label: 'eIF4A helicase',        cite: null },
  hpFree:      { value: 30,   min: 0,    max: 60,  step: 1,    unit: 'kcal/mol', label: 'Hairpin ΔG helicase-assisted scanners pass', cite: null },
  capZone:     { value: 15,   min: 0,    max: 60,  step: 1,    unit: 'nt',    label: 'Cap zone (hairpins here block loading)', cite: 'kozak1989' },
  capFree:     { value: 25,   min: 0,    max: 60,  step: 1,    unit: 'kcal/mol', label: 'Cap-zone hairpin ΔG tolerated', cite: null },
  capAUG:      { value: 0.3,  min: 0,    max: 1,   step: 0.05, unit: '×',     label: 'Initiation factor, AUG near cap', cite: null },
  capAUGnt:    { value: 12,   min: 0,    max: 40,  step: 1,    unit: 'nt',    label: 'Cap-proximal AUG zone', cite: null },
  hpKT:        { value: 4.0,  min: 1,    max: 10,  step: 0.1,  unit: 'kcal/mol', label: 'Hairpin barrier scale', cite: null },
  tc:          { value: 1.0,  min: 0.05, max: 1,   step: 0.01, unit: '',      label: 'eIF2·TC availability (1 = no stress)', cite: null },
  loadTCexp:   { value: 0.5,    min: 0,    max: 1.5, step: 0.05, unit: '',      label: 'Loading ∝ TC^β (β)', cite: null },
  tcRate:      { value: 1.5,  min: 0.05, max: 3,   step: 0.05, unit: '/s',    label: 'TC re-acquisition rate', cite: null },
  kozakStrong: { value: 0.97,  min: 0,    max: 1,   step: 0.01, unit: '',      label: 'P(init) strong Kozak',  cite: null },
  kozakAdeq:   { value: 0.6,  min: 0,    max: 1,   step: 0.01, unit: '',      label: 'P(init) adequate Kozak', cite: null },
  kozakWeak:   { value: 0.25, min: 0,    max: 1,   step: 0.01, unit: '',      label: 'P(init) weak Kozak',    cite: null },
  reinitMax:   { value: 0.8,  min: 0,    max: 1,   step: 0.01, unit: '',      label: 'P(reinitiate) short uORF', cite: null },
  reinitAA:    { value: 30,   min: 5,    max: 200, step: 1,    unit: 'aa',    label: 'Reinit decay length',   cite: null },
  elongSpeed:  { value: 15,   min: 1,    max: 60,  step: 1,    unit: 'nt/s',  label: 'Elongation speed',      cite: null },
  footprint:   { value: 30,   min: 10,   max: 60,  step: 1,    unit: 'nt',    label: 'Ribosome footprint',    cite: null },
};

export const SOURCES = {
  kozak1989: { text: 'Kozak 1989 Mol Cell Biol: −30 kcal/mol hairpin 12 nt from the cap blocked 40S loading; at 52 nt no effect; −61 kcal/mol stalled scanning', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC363665/' },

  wang2022:  { text: 'Wang et al. 2022 Cell: ~100 nt/s scanning (smFRET, yeast in vitro)', url: 'https://pubmed.ncbi.nlm.nih.gov/36334590/' },
  scan10:    { text: '10–20 nt/s scanning (single-molecule yeast scanning complexes)', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11446271/' },
  vattem2004:{ text: 'Vattem & Wek 2004 PNAS: ATF4 uORF1 (3 codons) promotes reinitiation; uORF2 overlaps CDS; eIF2α-P delays TC re-acquisition', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC509193/' },
  palam2011: { text: 'Palam, Baird & Wek 2011 J Biol Chem: eIF2α-P lowers polysomes generally; CHOP induction is attributed to ribosomes bypassing its uORF due to poor start-site context', url: 'https://pubmed.ncbi.nlm.nih.gov/21285359/' },
  noderer2014:{ text: 'Noderer et al. 2014 Mol Syst Biol: TIS efficiency of all 65,536 −6..+5 contexts', url: 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4299517/' },
};

export function values() {
  const v = {};
  for (const [k, p] of Object.entries(PARAMS)) v[k] = p.value;
  return v;
}
