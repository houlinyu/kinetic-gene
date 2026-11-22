// Parameters for the co-transcriptional splicing model (tab 2).
// Same shape as app/src/params.js: {value,min,max,step,unit,label,cite}.
// cite = key into SOURCES; cite: null means the value is ASSUMED and the UI must label it so.
// Every number with a cite was read out of the cited paper's own text (see SOURCES for the quote).
export const PARAMS = {
  // Initiation rate sets how many polymerases are on the gene, and so how much mRNA is made, but
  // polymerases in this model are independent (model.js draws each one's speed, pause and race
  // separately, and they cannot queue behind one another), so it cannot change WHICH isoform a
  // transcript becomes. It is animation-only, and the label says so rather than implying otherwise.
  initRate:    { value: 0.02,   min: 0.002, max: 0.5,  step: 0.002, unit: '/s',   label: 'Initiation rate (new Pol II per second) — animation only: polymerases do not interact, so this changes mRNA amount on screen, not the predicted isoform mix', cite: null },
  // Pausing enters the prediction exactly (model_exact.js pauseResolvent), but only through splice
  // sites that are already transcribed when the polymerase stalls: a pause delays every site
  // downstream of it by the same amount, and a uniform time shift cannot change a race.
  pauseProb:   { value: 0.95,   min: 0,     max: 1,    step: 0.01,  unit: '',     label: 'P(promoter-proximal pause) — affects the isoform mix only for 3′ splice sites within the pause site, which a uniform delay of everything downstream cannot', cite: 'jonkers2014' },
  pausePos:    { value: 50,     min: 0,     max: 200,  step: 1,     unit: 'nt',   label: 'Pause site, nt downstream of TSS (a 3′ss inside this window gets the whole pause to commit)', cite: 'jonkers2014' },
  pauseEscape: { value: 0.0017, min: 0.0002, max: 0.2, step: 0.0001, unit: '/s',  label: 'Pause escape rate (ln2 / 6.9 min half-life)', cite: 'jonkers2014' },
  elongRate:   { value: 25,     min: 8,     max: 67,   step: 1,     unit: 'nt/s', label: 'Pol II elongation rate (25 nt/s = 1.5 kb/min)', cite: 'veloso2014' },
  elongCV:     { value: 0.38,   min: 0,     max: 1,    step: 0.01,  unit: '',     label: 'Elongation-rate spread between polymerases (CV)', cite: 'jonkers2014' },
  commitScale: { value: 0.02,   min: 0.001, max: 0.5,  step: 0.001, unit: '/s',   label: 'Commitment rate of a median-strength site', cite: null },
  scoreSlope:  { value: 1.0,    min: 0,     max: 4,    step: 0.05,  unit: '',     label: 'How steeply commitment scales with splice-site score', cite: null },
  removalTime: { value: 480,    min: 30,    max: 1200, step: 10,    unit: 's',    label: 'Intron removal time after commitment (8 min)', cite: 'singh2009' },
  ctLag:       { value: 0,      min: 0,     max: 4000, step: 50,    unit: 'nt',   label: 'Extra nt Pol II must pass the 3′ss before the site can commit (0 = commit as soon as transcribed)', cite: null },
  exonWindow:  { value: 300,    min: 50,    max: 2000, step: 10,    unit: 'nt',   label: 'Exon-definition window (exons up to this size are favoured)', cite: null },
  exonDefGain: { value: 1.0,    min: 0,     max: 3,    step: 0.05,  unit: '×',    label: 'Exon-definition boost for a short exon', cite: null },
  postTx:      { value: 0.2,    min: 0,     max: 1,    step: 0.01,  unit: '',     label: 'Fraction of transcripts spliced post-transcriptionally (1 = no kinetic coupling)', cite: null },
};

export const SOURCES = {
  jonkers2014: {
    text: 'Jonkers, Kwak & Lis 2014 eLife 3:e02407 (GRO-seq + flavopiridol/triptolide, mouse ES cells): "the average half-life of paused Pol II is 7 min" (mean 6.9 min over 3181 genes); "Both inhibitors block transcription of more than 95% of genes"; "+20 to 100 nts downstream of the transcription start site (TSS), Pol II can be slowed down and paused"; "we measured elongation rates ranging from 0.5 to 4 kb/min at ~1000 genes"; stage averages 0.5 / 1.8 / 2.4 kb/min with SDs 0.43 / 0.69 / 0.79 kb/min (CV ~0.38 used here for between-polymerase spread).',
    url: 'https://pubmed.ncbi.nlm.nih.gov/24843027/',
  },
  veloso2014: {
    text: 'Veloso et al. 2014 Genome Res 24:896-905 (BruDRB-seq, >2000 human genes): "revealed gene-specific differences in elongation rates with a median rate of around 1.5 kb/min" (= 25 nt/s).',
    url: 'https://pubmed.ncbi.nlm.nih.gov/24714810/',
  },
  singh2009: {
    text: 'Singh & Padgett 2009 Nat Struct Mol Biol 16:1128-33 (DRB release + qRT-PCR, long human genes): "the rate of RNA polymerase II transcription over long genomic distances is about 3.8 kb min-1"; "co-transcriptional pre-mRNA splicing of U2-dependent introns occurs within 5-10 min of synthesis, irrespective of intron length between 1 kb and 240 kb".',
    url: 'https://pubmed.ncbi.nlm.nih.gov/19820712/',
  },
  drexler2020: {
    text: 'Drexler, Choquet & Churchman 2020 Mol Cell 77:985-998 (nano-COP, direct nascent-RNA nanopore sequencing): "the vast majority of human introns were not spliced until after Pol II transcribed 4 kb downstream, whereas in Drosophila S2 cells the majority of splicing occurred within 2 kb of introns"; consistent with "median splicing half-lives estimated from metabolic labeling data: 2 and 7-14 minutes in Drosophila and mammalian cells, respectively".',
    url: 'https://pubmed.ncbi.nlm.nih.gov/31839405/',
  },
  delamata2003: {
    text: 'de la Mata et al. 2003 Mol Cell 12:525-32: a human equivalent of the slow Drosophila C4 Pol II mutant "affects alternative splicing of the fibronectin EDI exon and adenovirus E1a pre-mRNA", and Ultrabithorax resplicing is stimulated in C4 mutant embryos — "a direct proof for the elongation control of alternative splicing in vivo".',
    url: 'https://pubmed.ncbi.nlm.nih.gov/14536091/',
  },
  dujardin2014: {
    text: 'Dujardin et al. 2014 Mol Cell 54:683-90 (Kornblihtt lab) states the direction of the kinetic-coupling effect: "In most cases, slow Pol II elongation allows weak splice sites to be recognized, leading to higher inclusion of alternative exons"; and for site strength, "Weakening of E9 5′ ss increases ETR-3 binding at the 3′ ss and subsequent E9 skipping, whereas strengthening of the 5′ ss usage has the opposite effect" (CFTR exon 9 is their counter-example, where slow Pol II instead causes skipping).',
    url: 'https://pubmed.ncbi.nlm.nih.gov/24793692/',
  },
  berget1995: {
    text: 'Berget 1995 J Biol Chem 270:2411-4, "Exon recognition in vertebrate splicing" — the exon-definition concept this model implements. The numeric exon-size window used here is NOT from this paper and is assumed.',
    url: 'https://doi.org/10.1074/jbc.270.6.2411',
  },
};

// Mechanism references for parameters whose numeric value is assumed (UI may show these as
// "mechanism from X, value assumed"). Keeps `cite` honest: cite = the number came from the paper.
export const MECHANISM = {
  exonWindow: 'berget1995',
  exonDefGain: 'berget1995',
  ctLag: 'drexler2020',
  postTx: 'drexler2020',
  commitScale: 'singh2009',
  scoreSlope: 'dujardin2014',
};

export function values() {
  const v = {};
  for (const [k, p] of Object.entries(PARAMS)) v[k] = p.value;
  return v;
}
