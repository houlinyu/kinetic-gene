#!/usr/bin/env python3
"""Hand-curated oncogene / therapeutic-lever table for the 'therapeutic window' view.

Every row below was read by the curator from the PubMed/PMC record named in `url`
(NCBI efetch, abstract text unless quote_source says otherwise) on 2026-11-22.
Quotes are verbatim and <= 15 words. Nothing here is computed; no effect sizes are
precomputed (SPEC_TARGET.md: the app scores live from the exact models).
Run: python3 pipeline/targets/build_targets.py  -> data/targets.json (+ app/public copy)
"""
import json, os, shutil, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

META = {
    "generated": "2026-11-22",
    "curator": "ONCOTARGET-DATA agent (manual curation; no automated text mining)",
    "method": ("Candidate genes and levers from SPEC_TARGET.md were each checked against a primary paper "
               "fetched with NCBI E-utilities (esearch + efetch, rettype=abstract) and read by the curator. "
               "A row ships only if that paper experimentally tested the gene's own 5' leader or its own splicing "
               "event, or (tier B) tested it only transcriptome-wide. Quotes are verbatim from the text fetched "
               "(abstract unless stated); every quote is <= 15 words. Candidates whose claim could not be verified "
               "were dropped and are listed in `rejected` with the reason."),
    "tiers": ("A: direct experiment on this gene's leader / splice event with isoform or uORF resolution. "
              "B: the gene's leader or event was measured only transcriptome-wide, or without resolving the "
              "isoforms the app compares. Same bar as data/evidence.json."),
    "honesty": ("Output of the app is PREDICTED differential sensitivity (model), never efficacy. Known failure "
                "modes, repeated from SPEC_TARGET.md: no validation against isoform-level translation efficiency; "
                "no IRES, RBP or ribosome-stalling mechanisms; kinetics assumed. The `param` field of each lever is "
                "OUR mapping of a drug class onto a model knob - the cited papers do not mention these parameters."),
    "sources_note": "cite keys are local to this file; url is the record the curator actually read.",
}

LEVERS = [
 {"id": "eif4a", "label": "eIF4A inhibition (rocaglates)",
  "param": "helicase", "also_params": ["hpFree"], "direction": "down", "model": "translation",
  "agent": "silvestrol / rocaglamide A (tool compounds); zotatifin (eFT226), first-in-class clinical eIF4A inhibitor",
  "mechanism": ("Rocaglates do not simply remove eIF4A unwinding activity: they clamp eIF4A onto polypurine "
                "sequences in the 5' leader, and the clamped protein blocks the scanning 43S subunit. In the model "
                "this is a loss of effective helicase support for a scanner meeting structure."),
  "status": ("silvestrol / RocA: preclinical tool compounds (leukaemia models in vitro and in vivo). zotatifin: "
             "described by its developers as a development candidate with 'significant antitumor activity'; "
             "IFN-stimulated-gene induction has been reported in breast-cancer patient biopsies after zotatifin "
             "treatment (the papers read here do not state a trial phase)."),
  "quote": "RocA specifically clamps eIF4A onto polypurine sequences in an ATP-independent manner",
  "quote_source": "abstract", "cite": "iwasaki2016", "url": "https://pubmed.ncbi.nlm.nih.gov/27309803/",
  "supporting": [
    {"cite": "wolfe2014", "url": "https://pubmed.ncbi.nlm.nih.gov/25079319/",
     "quote": "inhibition of eIF4A with silvestrol has powerful therapeutic effects against murine and human leukaemic cells"},
    {"cite": "ernst2020", "url": "https://pubmed.ncbi.nlm.nih.gov/32470302/",
     "quote": "eFT226 (Zotatifin), a compound with excellent physicochemical properties and significant antitumor activity"},
    {"cite": "zhao2023", "url": "https://pubmed.ncbi.nlm.nih.gov/37874652/",
     "quote": "induction of IFN-stimulated genes was observed in breast cancer patient biopsies following zotatifin treatment"}],
  "caveat": ("Our `helicase` knob is a scanning-support multiplier, i.e. the OLD 'structure needs unwinding' picture. "
             "Iwasaki 2016 showed rocaglates are gain-of-function clamps and that 5' UTR structure is only a minor "
             "determinant of their selectivity, so turning `helicase` down is at best a caricature of a rocaglate.")},

 {"id": "isr", "label": "ISR / eIF2 modulation",
  "param": "tc", "also_params": ["loadTCexp", "tcRate"], "direction": "down", "model": "translation",
  "direction_note": ("Stress / eIF2alpha phosphorylation moves `tc` DOWN; ISRIB moves it back UP. The lever is listed "
                     "as down because that is the direction a cancer cell experiences under ISR activation, and it is "
                     "the direction in which uORF mRNAs are de-repressed relative to bulk."),
  "agent": "ISRIB (tool compound; eIF2B activator) - and, in the opposite direction, any eIF2alpha-kinase stress",
  "mechanism": ("eIF2alpha phosphorylation lowers eIF2-GTP-Met-tRNAi (ternary complex) availability; scanning "
                "ribosomes then need longer to become reinitiation-competent, which changes how uORFs are used. "
                "ISRIB reverses the downstream effects of eIF2alpha phosphorylation, i.e. restores effective TC supply."),
  "status": "ISRIB: preclinical tool compound (cells and mice; cognition/ER-stress studies). Not an oncology drug.",
  "quote": "ISRIB, that potently (IC50 = 5 nM) reverses the effects of eIF2α phosphorylation",
  "quote_source": "abstract", "cite": "sidrauski2013", "url": "https://pubmed.ncbi.nlm.nih.gov/23741617/",
  "supporting": [
    {"cite": "vattem2004", "url": "https://pubmed.ncbi.nlm.nih.gov/15277680/",
     "quote": "increase the time required for the scanning ribosomes to become competent to reinitiate translation"}],
  "caveat": ("Direction is gene-specific, not global: low `tc` represses most mRNAs but DE-represses the "
             "reinitiation-dependent uORF mRNAs (ATF4 type). ISRIB is not a cancer agent; it is in the table as the "
             "cleanest published way to move this knob.")},

 {"id": "mtor", "label": "mTOR inhibition / cap-dependence",
  "param": "loadRate", "also_params": [], "direction": "down", "model": "translation",
  "agent": "Torin 1 (tool, full mTORC1 inhibition); rapamycin and rapalogs; INK128 (ATP-site mTOR inhibitor)",
  "mechanism": ("mTORC1 inhibition frees the 4E-BPs, which block the eIF4E-eIF4G1 interaction; eIF4E then loses much "
                "more of its capacity to bind TOP and TOP-like mRNAs than other mRNAs. In the model this is the rate "
                "at which a 43S subunit is loaded at the cap."),
  "status": ("rapamycin: 'a drug with many clinical uses' (Thoreen 2012). INK128: described by its authors as a "
             "clinically relevant ATP-site mTOR inhibitor, in prostate-cancer models. Torin 1: tool compound."),
  "quote": "consists almost entirely of transcripts with established 5' terminal oligopyrimidine (TOP) motifs",
  "quote_source": "abstract", "cite": "thoreen2012", "url": "https://pubmed.ncbi.nlm.nih.gov/22552098/",
  "supporting": [
    {"cite": "hsieh2012", "url": "https://pubmed.ncbi.nlm.nih.gov/22367541/",
     "quote": "we develop a clinically relevant ATP site inhibitor of mTOR, INK128"}],
  "caveat": ("Thoreen 2012 explicitly found NO support for the idea that mTORC1 preferentially hits long or complex "
             "5' UTRs ('We find no evidence to support proposals that mTORC1 preferentially regulates mRNAs' with "
             "increased 5' UTR length or complexity). So the mTOR lever must be presented as a TOP-motif lever only, "
             "and our structure features must not be used to rank mTOR sensitivity.")},

 {"id": "sso", "label": "Splice-switching ASO",
  "param": "site_score", "also_params": [], "direction": "down", "model": "splice",
  "param_note": ("The targeted site's own strength score in splice_data.json (donor.score / acceptor.score), not a global "
                 "knob: down at the masked site, which lets the competing site win."),
  "agent": "nusinersen (approved, SMA - proof of concept); MDM4 exon-6 ASO (preclinical oncology)",
  "mechanism": ("A steric-block oligonucleotide base-pairs over a splice site or regulatory element so the "
                "spliceosome cannot use it; the competing site wins instead. Nusinersen does this on SMN2 to "
                "include exon 7; an MDM4 exon-6 ASO does the reverse (forces skipping) to destabilise MDM4 mRNA."),
  "status": ("nusinersen: approved (December 2016) for spinal muscular atrophy. MDM4 exon-6 ASO: preclinical - "
             "melanoma cell lines and patient-derived xenografts, plus a DLBCL PDX."),
  "quote": "steric block ASO that binds the SMN2 messenger RNA and promotes exon 7 inclusion",
  "quote_source": "abstract", "cite": "wood2017", "url": "https://pubmed.ncbi.nlm.nih.gov/28977438/",
  "supporting": [
    {"cite": "dewaele2016", "url": "https://pubmed.ncbi.nlm.nih.gov/26595814/",
     "quote": "skipping of exon 6 decreased MDM4 abundance, inhibited melanoma growth"}],
  "caveat": ("No splice-switching ASO in an oncology clinical trial could be verified in the time available; the "
             "oncology example here is preclinical. Nusinersen is the approved proof of concept but is not an "
             "oncology drug. In the model an ASO is a per-site score change, which ignores delivery, "
             "pharmacokinetics and the ASO's own sequence constraints.")},

 {"id": "sf3b1", "label": "SF3B1 / SF3b complex modulation",
  "param": "commitScale", "also_params": ["removalTime"], "direction": "down", "model": "splice",
  "param_note": "Commitment to a transcribed site is slowed; the retained-intron arm of the drug would also raise removalTime.",
  "agent": "H3B-8800 (oral SF3b modulator)",
  "mechanism": ("H3B-8800 binds the SF3b complex, i.e. the machinery that recognises the branch point / 3' splice "
                "site region, so commitment to a 3' site is impaired and short GC-rich introns are retained. In the "
                "model: slower commitment to a site once transcribed."),
  "status": "Phase I first-in-human dose escalation in MDS, CMML and AML (84 patients enrolled).",
  "quote": "orally available modulator of the SF3b complex, H3B-8800",
  "quote_source": "abstract", "cite": "seiler2018", "url": "https://pubmed.ncbi.nlm.nih.gov/29457796/",
  "supporting": [
    {"cite": "steensma2021", "url": "https://pubmed.ncbi.nlm.nih.gov/34172893/",
     "quote": "We conducted a phase I clinical trial of H3B-8800, an oral small molecule"}],
  "caveat": ("The measured consequence is intron retention of short, GC-rich introns and selective killing of "
             "spliceosome-MUTANT cells; our model has no intron-retention output and no mutant-spliceosome state, so "
             "the mapping to commitScale is an approximation of one arm of the drug's effect.")},

 {"id": "clk_srpk", "label": "CLK / SRPK (SR-protein kinase) inhibition",
  "param": "scoreSlope", "also_params": ["commitScale"], "direction": "down", "model": "splice",
  "param_note": "Lower slope = the spliceosome discriminates less between strong and weak sites, so weak/alternative sites gain.",
  "agent": "TG003 (CLK1/Sty inhibitor, tool); SM08502 / cirtuvivint (CLK inhibitor, clinical); SRPIN340 (SRPK1/2 inhibitor, tool)",
  "mechanism": ("CLK and SRPK kinases phosphorylate SR proteins; SR phosphorylation is what makes a splice site "
                "efficiently recognised. Inhibiting them suppresses SR phosphorylation, disrupts spliceosome "
                "assembly and shifts alternative splicing. In the model: how steeply commitment scales with a site's "
                "score (site discrimination)."),
  "status": ("TG003 and SRPIN340: tool compounds (in vitro splicing, cells, Xenopus embryos / virus assays). "
             "SM08502: 'first-in-class CLK inhibitor being investigated in a Phase 1 clinical trial' (NCT03355066)."),
  "quote": "a benzothiazole compound, TG003, had a potent inhibitory effect on the activity of Clk1/Sty",
  "quote_source": "abstract", "cite": "muraki2004", "url": "https://pubmed.ncbi.nlm.nih.gov/15010457/",
  "supporting": [
    {"cite": "tam2020", "url": "https://pubmed.ncbi.nlm.nih.gov/31560935/",
     "quote": "SM08502 is a first-in-class CLK inhibitor being investigated in a Phase 1 clinical trial"},
    {"cite": "fukuhara2006", "url": "https://pubmed.ncbi.nlm.nih.gov/16840555/",
     "quote": "which preferentially inhibits SRPK1 and SRPK2 and down-regulates SRp75"}],
  "caveat": ("Mapping SR-protein phosphorylation onto a single 'score slope' knob is our modelling choice; none of "
             "these papers measured a per-site commitment rate. SRPIN340's published characterisation is antiviral, "
             "not oncological.")},

 {"id": "elong", "label": "Pol II elongation rate (kinetic coupling)",
  "param": "elongRate", "also_params": ["pauseEscape"], "direction": "down", "model": "splice",
  "param_note": "Slower Pol II (and slower pause escape) gives weak/alternative sites more time - the kinetic-coupling arm.",
  "agent": "DRB and flavopiridol (P-TEFb/CDK9 inhibitors, tool compounds); camptothecin (used as an elongation inhibitor)",
  "mechanism": ("P-TEFb (CDK9/cyclin T1) phosphorylation releases paused Pol II into productive elongation; "
                "flavopiridol and DRB inhibit P-TEFb, so Pol II spends longer upstream of downstream splice sites. "
                "With elongation inhibited, Pol II occupancy rises in the introns flanking alternative exons and "
                "those exons are generally included more."),
  "status": ("flavopiridol: 'a cyclin-dependent kinase (Cdk) inhibitor in clinical trials as a cancer therapy' "
             "(as stated in 2001); DRB and camptothecin: laboratory tools here."),
  "quote": "Flavopiridol is the most potent P-TEFb inhibitor reported",
  "quote_source": "abstract", "cite": "chao2001", "url": "https://pubmed.ncbi.nlm.nih.gov/11431468/",
  "supporting": [
    {"cite": "ip2011", "url": "https://pubmed.ncbi.nlm.nih.gov/21163941/",
     "quote": "increased Pol II occupancy was detected in the intronic regions flanking the alternative exons"},
    {"cite": "delamata2003", "url": "https://pubmed.ncbi.nlm.nih.gov/14536091/",
     "quote": "a human equivalent to Drosophila's C4 pol II in human cultured cells affects alternative splicing"}],
  "caveat": ("These agents inhibit transcription globally (flavopiridol blocks most Pol II transcription), so in a "
             "cell they are not a clean 'slow Pol II' knob; de la Mata 2003's clean version of the experiment is a "
             "genetic slow-Pol II mutant, not a drug. Ip 2011's direction (more inclusion) was a global tendency, "
             "with exceptions.")},
]

GENES = [
 {"gene": "MDM2", "role": "oncogene", "disease": "soft-tissue sarcoma / tumours with MDM2 overexpression",
  "mechanism": "alternative promoter (P1 vs P2) -> long leader with two uORFs vs short leader",
  "event": "uORF + alternative promoter", "isoform_resolved": True, "tier": "A",
  "measured": "polysome/ribosome loading of L-mdm2 vs S-mdm2 mRNA, and leader-reporter assays",
  "cells": "mouse/human tumour cells (as reported in the abstract: 'the tumor cells')",
  "quote": "this mRNA was loaded with ribosomes inefficiently in comparison with S-mdm2",
  "quote_source": "abstract", "cite": "brown1999", "url": "https://pubmed.ncbi.nlm.nih.gov/10523842/",
  "why_targetable": "Which promoter is used decides whether two uORFs sit in the leader, so MDM2 protein output is set by scanning/reinitiation, the part of the model the ISR lever moves.",
  "levers": ["isr", "eif4a"],
  "caveat": "Already tier A in data/evidence.json with a verified transcript mapping; the two Ensembl leaders must stay matched to L-/S-mdm2 before any isoform-level claim."},

 {"gene": "ERBB2", "role": "oncogene", "disease": "breast carcinoma (HER-2 overexpressed in 25-30%)",
  "mechanism": "short uORF in the transcript leader; leader represses, and repression is relieved in transformed cells",
  "event": "uORF", "isoform_resolved": False, "tier": "A",
  "measured": "polysome distribution of endogenous HER-2 mRNA; leader-reporter repression",
  "cells": "primary human fibroblasts and mammary epithelial cells vs BT474, MCF-7, COS-7",
  "quote": "The transcript leader in all cell types contains a short upstream open reading frame",
  "quote_source": "abstract", "cite": "child1999", "url": "https://pubmed.ncbi.nlm.nih.gov/10216954/",
  "why_targetable": "A uORF-controlled oncogene whose leader repression differs between normal and transformed cells - exactly the kind of difference a uORF-sensitive lever would exploit.",
  "levers": ["isr"],
  "caveat": "Not isoform-resolved: the same leader was found in both cell types (no mRNA size difference), so the effect is cell state, not an isoform switch."},

 {"gene": "PIM1", "role": "oncogene", "disease": "lymphoma / leukaemia (oncogenic serine-threonine kinase)",
  "mechanism": "long (400 nt), 76% GC 5' leader that represses translation; repression relieved by eIF4E overexpression",
  "event": "structured leader", "isoform_resolved": False, "tier": "A",
  "measured": "5'-UTR deletion reporters in vitro and in cells; eIF4E overexpression rescue",
  "cells": "Hut-78 lymphocytes (endogenous), rabbit reticulocyte lysate, NIH-3T3",
  "quote": "Deletion of the 5'-UTR of pim-1 increased translation of the Pim-1 protein approximately 10-fold",
  "quote_source": "abstract", "cite": "hoover1997", "url": "https://pubmed.ncbi.nlm.nih.gov/9419425/",
  "why_targetable": "A measured, structured, GC-rich leader whose output depends on cap-binding/unwinding capacity - the helicase and loading knobs.",
  "levers": ["eif4a", "mtor"],
  "caveat": "eIF4E dependence was shown by overexpression, not by an eIF4A inhibitor; no rocaglate was used on PIM1 in this paper."},

 {"gene": "ODC1", "role": "oncogene", "disease": "polyamine-driven proliferation; MYC-pathway tumours (gene role not characterised in the cited paper)",
  "mechanism": "GC-rich 5'-proximal leader segment with the potential to form a very stable hairpin; represses translation",
  "event": "structured leader", "isoform_resolved": False, "tier": "A",
  "measured": "ODC 5'UTR-luciferase / CAT reporters, transient expression and in vitro translation; 5'UTR truncations",
  "cells": "wild-type and ODC-deficient hamster cells; rabbit reticulocyte lysate (and rat ODC in vitro, Van Steeg 1991)",
  "quote": "shown to be due to the GC-rich 5' segment of the UTR",
  "quote_source": "abstract", "cite": "grens1990", "url": "https://pubmed.ncbi.nlm.nih.gov/2365700/",
  "why_targetable": "One of the few leaders where a structure-dependent block on scanning was mapped to a defined GC-rich segment, with a dG threshold measured in the follow-up (Van Steeg 1991).",
  "levers": ["eif4a"],
  "supporting": [{"cite": "vansteeg1991", "url": "https://pubmed.ncbi.nlm.nih.gov/2006916/",
                  "quote": "A truncated 5' UTR with a calculated free energy of less than -272 kJ"}],
  "caveat": "Rodent (hamster/rat) ODC leaders, not human ODC1; the human leader has not been tested in these papers. Van Steeg 1991 also reports the short uORF in the leader does NOT contribute to repression."},

 {"gene": "MYC", "role": "oncogene", "disease": "pancreatic adenocarcinoma, lymphoma, leukaemia",
  "mechanism": "long, highly structured leader (G-quadruplex-forming elements) -> eIF4A-dependent translation",
  "event": "structured leader / G-quadruplex", "isoform_resolved": False, "tier": "B",
  "measured": "transcriptome-scale ribosome footprinting under eIF4A inhibition (PDAC models)",
  "cells": "murine and human pancreatic adenocarcinoma models",
  "quote": "critical KRAS signaling molecules such as PI3K, RALA, RAC2, MET, MYC, and YAP1",
  "quote_source": "abstract", "cite": "singh2021", "url": "https://pubmed.ncbi.nlm.nih.gov/33632898/",
  "why_targetable": "MYC is named among the eIF4A-dependent transcripts in ribosome-footprinting data, so it is the canonical test case for the helicase lever.",
  "levers": ["eif4a", "mtor"],
  "caveat": ("Gene-level, not isoform-level: no MYC leader reporter in this paper. data/evidence.json keeps MYC at "
             "tier B and show_as_example=false because the MANE leader starts at a CUG (c-Myc1) and the published "
             "P0/P1/P2 leaders are not mapped to Ensembl transcripts. The older MYC claim (nanbru1997) is an IRES, "
             "which the scanning model does not implement.")},

 {"gene": "KRAS", "role": "oncogene", "disease": "pancreatic adenocarcinoma (PDAC)",
  "mechanism": "long, highly structured 5' leader with multiple G-quadruplex elements -> eIF4A-dependent",
  "event": "structured leader / G-quadruplex", "isoform_resolved": False, "tier": "B",
  "measured": "ribosome footprinting + computational 5'UTR feature analysis under eIF4A inhibition",
  "cells": "murine and human PDAC models",
  "quote": "Computational analyses identified these features in mRNAs encoding KRAS and key downstream molecules",
  "quote_source": "abstract", "cite": "singh2021", "url": "https://pubmed.ncbi.nlm.nih.gov/33632898/",
  "why_targetable": "If the leader features are real, the undruggable driver becomes reachable through its translation requirement rather than its protein.",
  "levers": ["eif4a"],
  "caveat": "The KRAS leader itself was identified computationally; the experimental read-out is transcriptome-wide footprinting, not a KRAS leader reporter. KRAS4A/4B alternative splicing was NOT verified here (see rejected)."},

 {"gene": "MET", "role": "oncogene", "disease": "stage IV lung adenocarcinoma (also gastric cancer line Hs746T)",
  "mechanism": "mutations in the exon 14 splice acceptor/donor sites cause exon 14 skipping, deleting the juxtamembrane CBL-binding site",
  "event": "alternative splice event (exon skipping)", "isoform_resolved": True, "tier": "A",
  "measured": "tumour sequencing + cDNA sequencing of the skipped transcript; clinical response to MET inhibitors",
  "cells": "patient tumours (MSKCC lung adenocarcinoma cohort); Hs746T gastric line (Yamamoto 2010)",
  "quote": "MET exon 14 RNA splice acceptor and donor sites, which lead to exon skipping",
  "quote_source": "abstract", "cite": "paik2015", "url": "https://pubmed.ncbi.nlm.nih.gov/25971939/",
  "why_targetable": "A splice event that is itself the actionable lesion: patients whose tumours skip exon 14 respond to MET inhibitors, so isoform fraction is the biomarker.",
  "levers": ["sso", "elong", "clk_srpk"],
  "supporting": [{"cite": "yamamoto2010", "url": "https://pubmed.ncbi.nlm.nih.gov/20331976/",
                  "quote": "Hs746T cells had a splice site mutation of exon 14"}],
  "caveat": "The skipping is driven by a cis mutation of the splice site, so a kinetic/trans lever would not restore exon 14; MET is here as a splice-event target, not as a drug-reversible one."},

 {"gene": "BCL2L1", "role": "oncogene", "disease": "lymphoid malignancy / apoptosis resistance (Bcl-xL)",
  "mechanism": "alternative 5' splice site inside exon 2: long mRNA -> anti-apoptotic Bcl-xL, short mRNA -> pro-apoptotic Bcl-xS",
  "event": "alternative 5' splice site", "isoform_resolved": True, "tier": "A",
  "measured": "cloning of both mRNAs; stable transfection and survival after growth-factor withdrawal; tissue expression",
  "cells": "IL-3-dependent murine cell line (FL5.12-type); mouse tissues",
  "quote": "Alternative splicing results in two distinct bcl-x mRNAs",
  "quote_source": "abstract", "cite": "boise1993", "url": "https://pubmed.ncbi.nlm.nih.gov/8358789/",
  "why_targetable": "The two isoforms have opposite functions from one 5' splice-site choice: shifting the ratio, not the level, is the therapeutic goal.",
  "levers": ["sso", "sf3b1", "clk_srpk", "elong"],
  "caveat": "Function was shown by overexpression of each isoform, not by shifting endogenous splicing with a drug; the paper is murine/rat-system work."},

 {"gene": "MCL1", "role": "oncogene", "disease": "myeloid leukaemia and many solid tumours (MCL-1L overexpressed)",
  "mechanism": "exon 2 skipping gives MCL-1S, a BH3-only pro-apoptotic protein lacking BH1, BH2 and the transmembrane domain",
  "event": "alternative splice event (exon skipping)", "isoform_resolved": True, "tier": "A",
  "measured": "cloning of the variant from human placenta; gene-structure analysis; yeast two-hybrid, co-precipitation, apoptosis after overexpression",
  "cells": "human placenta (mRNA source); CHO cells (apoptosis assays); transfected mammalian cells",
  "quote": "MCL-1S results from the splicing out of exon 2 during mRNA processing",
  "quote_source": "abstract", "cite": "bae2000", "url": "https://pubmed.ncbi.nlm.nih.gov/10837489/",
  "why_targetable": "Same logic as BCL2L1: one exon decides pro- vs anti-apoptotic product, and the long:short ratio is what cancer cells skew.",
  "levers": ["sso", "sf3b1", "clk_srpk", "elong"],
  "caveat": "Pro-apoptotic activity was shown by overexpression in CHO cells; no drug was used to shift the endogenous ratio in this paper."},

 {"gene": "AR", "role": "oncogene", "disease": "hormone-refractory / castration-resistant prostate cancer",
  "mechanism": "splicing of a cryptic 'intronic' exon onto the DNA-binding-domain exons -> AR-V7 lacking the ligand-binding domain, constitutively active",
  "event": "alternative splice event (cryptic exon inclusion)", "isoform_resolved": True, "tier": "A",
  "measured": "variant-specific qPCR in 25 HRPC vs 82 hormone-naive tumours; AR-V7-specific antibody; reporter and expression microarray activity",
  "cells": "human prostate tumour specimens; cultured prostate cancer cells under androgen depletion",
  "quote": "seven AR variant transcripts lacking the reading frames for the ligand-binding domain",
  "quote_source": "abstract", "cite": "hu2009", "url": "https://pubmed.ncbi.nlm.nih.gov/19117982/",
  "why_targetable": "The resistance mechanism IS an isoform: an agent that suppressed cryptic-exon inclusion would restore dependence on the ligand-binding domain that existing drugs target.",
  "levers": ["sso", "clk_srpk", "elong"],
  "caveat": "Isoform abundance and activity were measured; no splice-modulating drug was tested here. 20-fold higher expression in HRPC is a cohort comparison, not a per-cell measurement."},

 {"gene": "BRAF", "role": "oncogene", "disease": "melanoma, acquired resistance to vemurafenib",
  "mechanism": "aberrantly spliced BRAF(V600E) lacking exons 4-8 (p61BRAF(V600E)) removes the RAS-binding domain and dimerises constitutively",
  "event": "alternative splice event (multi-exon skipping)", "isoform_resolved": True, "tier": "A",
  "measured": "variant identification in resistant cells and patient tumours; dimerisation and ERK-signalling assays",
  "cells": "vemurafenib-resistant melanoma cells and patient-derived samples",
  "quote": "express a 61-kDa variant form of BRAF(V600E), p61BRAF(V600E), which lacks exons 4-8",
  "quote_source": "abstract", "cite": "poulikakos2011", "url": "https://pubmed.ncbi.nlm.nih.gov/22113612/",
  "why_targetable": "A drug-resistance isoform produced by splicing, so the isoform fraction is itself a resistance biomarker and a splicing target.",
  "levers": ["sso", "clk_srpk", "elong"],
  "supporting": [{"cite": "vido2018", "url": "https://pubmed.ncbi.nlm.nih.gov/30404005/",
                  "quote": "Expression of aberrantly spliced BRAF V600E isoforms (BRAF V600E ΔEx) mediates resistance"}],
  "caveat": "The papers do not establish which cis element or factor causes the skipping, so there is no verified handle for a specific splice lever."},

 {"gene": "CCND1", "role": "oncogene", "disease": "prostate adenocarcinoma",
  "mechanism": "use of the alternative 3' end at the exon 4 / intron 4 boundary gives cyclin D1b; the G/A870 polymorphism sits in that splice site",
  "event": "alternative splice event", "isoform_resolved": True, "tier": "A",
  "measured": "isoform-specific expression in prostate adenocarcinomas; minigenes testing G/A870; two case-control cohorts",
  "cells": "human prostate tumours and normal prostate; prostate cell lines (minigenes)",
  "quote": "A870 predisposed for transcript-b production in cells and in normal prostate",
  "quote_source": "abstract", "cite": "comstock2009", "url": "https://pubmed.ncbi.nlm.nih.gov/19706803/",
  "why_targetable": "A splice-site genotype changes the isoform ratio of a core cell-cycle oncogene, and tumours override that genotype - a splicing knob with a measured read-out.",
  "levers": ["sso", "clk_srpk", "elong"],
  "caveat": "The abstract states the polymorphism is 'within the splice site' and reports transcript-b production; the exact intron-4 retention mechanism was not re-verified here."},

 {"gene": "CDKN2A", "role": "tumour suppressor", "disease": "melanoma and many cancers (INK4a locus deleted/mutated frequently)",
  "mechanism": "one locus, two first exons and two reading frames: p16INK4a and p19ARF (human p14ARF) from alternative promoters/first exons",
  "event": "alternative promoter / alternative first exon (alternative reading frame)", "isoform_resolved": True, "tier": "A",
  "measured": "cDNA characterisation and ectopic expression of each product; cell-cycle arrest assays",
  "cells": "rodent fibroblasts (ectopic expression)",
  "quote": "An unrelated protein (p19ARF) arises in major part from an alternative reading frame",
  "quote_source": "abstract", "cite": "quelle1995", "url": "https://pubmed.ncbi.nlm.nih.gov/8521522/",
  "why_targetable": "A tumour suppressor where promoter/first-exon choice decides which protein is made - the cleanest example that a leader choice can change the protein, not just its amount.",
  "levers": [],
  "caveat": "Mouse INK4a locus and p19ARF; the human p14ARF equivalence is not established in this paper. No lever in our set moves promoter choice, so this gene is a mechanism illustration, not a drug target row."},

 {"gene": "MDM4", "role": "oncogene", "disease": "melanoma, diffuse large B-cell lymphoma (MDM4 upregulated in many cancers)",
  "mechanism": "exon 6 inclusion gives full-length MDM4; skipping gives an NMD-targeted isoform (MDM4-S). SRSF3 enhances inclusion",
  "event": "alternative splice event (cassette exon)", "isoform_resolved": True, "tier": "A",
  "measured": "isoform quantification in normal tissue vs tumours; ASO-mediated exon-6 skipping; tumour growth and drug sensitivity",
  "cells": "multiple human melanoma cell lines, melanoma PDX mice, DLBCL PDX",
  "quote": "skipping of exon 6 decreased MDM4 abundance, inhibited melanoma growth",
  "quote_source": "abstract", "cite": "dewaele2016", "url": "https://pubmed.ncbi.nlm.nih.gov/26595814/",
  "why_targetable": "The only gene in this table where a splice-switching oligonucleotide was actually given to tumour-bearing animals and changed the outcome - the proof of concept for the splice lever in oncology.",
  "levers": ["sso", "clk_srpk", "elong"],
  "caveat": "Preclinical only. The ASO works by steric block at exon 6, not by a global splicing knob; SRSF3 dependence means a CLK/SRPK lever could in principle move it, but that was not tested."},
]

REFERENCE_SET = [
 {"gene": "GAPDH", "role": "reference",
  "why": ("Measured comparator, not an assumption: in rat liver, when AMPK activation repressed translation of the "
          "TOP mRNAs rpS6, rpS8 and eEF1alpha, translation of GAPDH mRNA was unchanged, as was global protein "
          "synthesis. GAPDH therefore behaves as the insensitive bulk class for the mTOR/4E-BP axis our `loadRate` "
          "lever represents."),
  "quote": "neither global protein synthesis nor translation of mRNAs encoding GAPDH or beta-actin was changed",
  "quote_source": "abstract", "cite": "reiter2008", "url": "https://pubmed.ncbi.nlm.nih.gov/18638456/",
  "caveat": "Rat liver, AMPK activation by treadmill running (an upstream way of repressing mTORC1), not a rapalog and not human cells. No eIF4A-inhibitor data for GAPDH was verified here."},
 {"gene": "ACTB", "role": "reference",
  "why": "Same experiment, same result: beta-actin mRNA translation was unchanged while TOP mRNAs were repressed. Also a check on the common assumption that beta-actin is a TOP mRNA - in this experiment it did not behave like one.",
  "quote": "neither global protein synthesis nor translation of mRNAs encoding GAPDH or beta-actin was changed",
  "quote_source": "abstract", "cite": "reiter2008", "url": "https://pubmed.ncbi.nlm.nih.gov/18638456/",
  "caveat": "Same caveats as GAPDH (rat liver, AMPK). Being insensitive to mTORC1 repression does not make a leader unstructured; our pipeline must confirm the leader features once ACTB is built."},
 {"gene": "HBB", "role": "reference",
  "why": ("Beta-globin mRNA is the classic efficiently translated standard: in the ODC 5'UTR in vitro experiments it "
          "is the control message the repressed and de-repressed ODC constructs were compared against."),
  "quote": "or with that of control globin mRNA",
  "quote_source": "abstract", "cite": "vansteeg1991", "url": "https://pubmed.ncbi.nlm.nih.gov/2006916/",
  "caveat": "This documents globin mRNA's use as the efficiently translated comparator; it is not a measurement of the HBB leader's features. Weakest of the three reference justifications."},
]

DIFFERENTIAL = {
 "claim_under_test": ("The app's therapeutic-window view assumes that turning a knob hurts some transcripts more than "
                      "others, by a leader or splicing feature we can compute. Below is what the literature read here "
                      "actually shows, including where it contradicts that picture."),
 "eif4a": {
   "what_is_supported": ("Preferential repression is real and reproducible: eIF4A inhibition hits a minority of "
                         "transcripts. WHICH feature explains it is disputed. Wolfe 2014 (ribosome footprinting, "
                         "T-ALL cells and mouse models) reported a 12-nt (CGG)4 guanine-quartet motif and G-quadruplex "
                         "potential as hallmarks, enriched in oncogenes and super-enhancer transcription factors. "
                         "Rubio 2014 (ribosome profiling, MDA-MB-231 breast cancer cells, silvestrol) found 284 "
                         "eIF4A-dependent genes and reported 5'UTR complexity as the determinant. Singh 2021 (PDAC) "
                         "reported long, highly structured leaders with multiple G-quadruplexes and named MYC, MET, "
                         "YAP1, PI3K, RALA, RAC2 among the dependent transcripts. Iwasaki 2016 (in vitro + cells) "
                         "contradicts the structure explanation: RocA selectivity comes from polypurine sequence, with "
                         "structure only a minor determinant, and the clamped eIF4A blocks 43S scanning and causes "
                         "premature upstream initiation."),
   "quotes": [
     {"cite": "wolfe2014", "url": "https://pubmed.ncbi.nlm.nih.gov/25079319/",
      "quote": "the 12-nucleotide guanine quartet (CGG)4 motif that can form RNA G-quadruplex structures"},
     {"cite": "rubio2014", "url": "https://pubmed.ncbi.nlm.nih.gov/25273840/",
      "quote": "5′UTR complexity determines Silvestrol-sensitivity and altering 5′UTR structure modifies translational output"},
     {"cite": "singh2021", "url": "https://pubmed.ncbi.nlm.nih.gov/33632898/",
      "quote": "mRNAs with long and highly structured 5' untranslated regions, including those with multiple G-quadruplex elements"},
     {"cite": "iwasaki2016", "url": "https://pubmed.ncbi.nlm.nih.gov/27309803/",
      "quote": "secondary structure in 5' untranslated regions is only a minor determinant for RocA selectivity"},
     {"cite": "iwasaki2016", "url": "https://pubmed.ncbi.nlm.nih.gov/27309803/",
      "quote": "clamped eIF4A blocks 43S scanning, leading to premature, upstream translation initiation"}],
   "conflicts": [
     ("Method conflict, stated in the literature itself: Singh 2021 says its ribosome-footprinting result 'contrast[s] "
      "with a recent study that relied on an older method, polysome fractionation' (Rubio 2014, which implicated "
      "redox-related genes). The gene lists from different eIF4A studies do not fully agree."),
     ("Mechanism conflict: a structure-unwinding model (Wolfe/Rubio/Singh) vs a sequence-specific clamp model "
      "(Iwasaki). Our `helicase` knob implements the first; the second predicts that a polypurine tract, not a "
      "hairpin, marks the sensitive transcript, and predicts NEW upstream initiation events rather than slower scanning.")],
   "what_the_app_may_claim": ("That mRNAs whose leaders our pipeline measures as long, structured, G-quadruplex-bearing "
                              "lose more predicted output when `helicase` is lowered - and that this matches the feature "
                              "set reported for eIF4A-dependent transcripts in the three footprinting/profiling studies "
                              "above. It may NOT claim that a given gene will respond to a rocaglate in a patient, nor "
                              "that our structure score is the mechanism Iwasaki 2016 describes.")},
 "mtor": {
   "what_is_supported": ("Strongly supported but for a different feature: the mTORC1-sensitive class is almost entirely "
                         "TOP and TOP-like mRNAs (ribosome profiling in mouse cells, acute Torin 1), and the mechanism "
                         "is 4E-BP-mediated loss of eIF4E binding, which costs TOP mRNAs more than others. Crucially, "
                         "the same paper found NO support for leader length or complexity as an mTORC1 determinant."),
   "quotes": [
     {"cite": "thoreen2012", "url": "https://pubmed.ncbi.nlm.nih.gov/22552098/",
      "quote": "consists almost entirely of transcripts with established 5' terminal oligopyrimidine (TOP) motifs"},
     {"cite": "thoreen2012", "url": "https://pubmed.ncbi.nlm.nih.gov/22552098/",
      "quote": "We find no evidence to support proposals that mTORC1 preferentially regulates mRNAs"},
     {"cite": "rubio2014", "url": "https://pubmed.ncbi.nlm.nih.gov/25273840/",
      "quote": "eIF4A inhibition alters translation of an mRNA subset distinct from those affected by mTOR"}],
   "what_the_app_may_claim": ("That lowering `loadRate` should preferentially hurt transcripts our pipeline annotates "
                              "with a TOP / TOP-like motif, and that the eIF4A-sensitive and mTOR-sensitive sets are "
                              "reported to be DIFFERENT sets. It may not rank mTOR sensitivity by leader structure."),
   "positive_controls_not_in_the_reference_set": ("RPS6, RPS8 and EEF1A1 are the TOP mRNAs repressed in Reiter 2008 and "
                                                  "the TOP class of Thoreen 2012. They belong in the sensitive arm of "
                                                  "the mTOR comparison, NOT in the bulk reference set.")},
 "isr": {
   "what_is_supported": ("Preferential translation under low ternary-complex availability is mechanistically resolved "
                         "for the reinitiation-type uORF arrangement: uORF1 is positive-acting, uORF2 inhibitory, and "
                         "when eIF2-GTP falls the extra time needed to become reinitiation-competent lets ribosomes "
                         "scan through uORF2 into the main ORF. This is a uORF-architecture claim, not a structure claim."),
   "quotes": [
     {"cite": "vattem2004", "url": "https://pubmed.ncbi.nlm.nih.gov/15277680/",
      "quote": "increase the time required for the scanning ribosomes to become competent to reinitiate translation"},
     {"cite": "sidrauski2013", "url": "https://pubmed.ncbi.nlm.nih.gov/23741617/",
      "quote": "ISRIB, that potently (IC50 = 5 nM) reverses the effects of eIF2α phosphorylation"}],
   "what_the_app_may_claim": ("That the sign of the `tc` effect flips between uORF architectures (uORF1+inhibitory uORF2 "
                              "vs a single strong uORF), and that this is the measured mechanism for ATF4. Mouse ATF4 "
                              "leader; see data/evidence.json for the transcript-mapping caveats.")},
 "splicing": {
   "what_is_supported": ("Preferential effects are reported for both classes of splicing lever, but on features our "
                         "model does not yet compute. H3B-8800's selectivity for spliceosome-mutant cells is "
                         "attributed to retention of short, GC-rich introns. Elongation inhibition shifts alternative "
                         "exons towards inclusion as a global tendency, with Pol II piling up in the flanking introns, "
                         "and the affected genes were enriched for RNA-binding and RNA-processing functions."),
   "quotes": [
     {"cite": "seiler2018", "url": "https://pubmed.ncbi.nlm.nih.gov/29457796/",
      "quote": "preferential killing of spliceosome-mutant cells is due to retention of short, GC-rich introns"},
     {"cite": "ip2011", "url": "https://pubmed.ncbi.nlm.nih.gov/21163941/",
      "quote": "these exons generally became more included"}],
   "what_the_app_may_claim": ("That slowing `elongRate` shifts weak/alternative sites towards use (kinetic coupling), "
                              "citing Ip 2011 for the global direction and de la Mata 2003 for the clean genetic proof. "
                              "It may NOT claim an intron-retention or GC-content prediction, because the model has "
                              "neither output.")},
}

REJECTED = [
 {"item": "EGFRvIII (EGFR exons 2-7 deletion)", "reason": "The lesion is a genomic rearrangement/deletion, not a regulated splicing choice, so neither the leader model nor the splicing model represents it. No primary paper fetched in time; nothing shipped."},
 {"item": "KRAS4A vs KRAS4B alternative exon 4", "reason": "PubMed searches returned no primary paper the curator could read within the session; the KRAS row ships on its structured-leader/eIF4A evidence only, and the 4A/4B isoform claim is marked unverified and left out."},
 {"item": "start-codon-fidelity lever (kozakStrong/kozakAdeq/kozakWeak)", "reason": "No pharmacological agent with a verified published effect on start-codon selection / leaky scanning was found. SPEC_TARGET.md already flags this as '(tool compounds; verify)'. Dropped rather than invented."},
 {"item": "BCL2, CTNNB1, VEGFA, TP53 as target examples", "reason": "The verified mechanisms are IRES / cap-independent initiation (VEGFA huez1998, TP53 ray2006 in data/evidence.json) or unverified in this session (BCL2, CTNNB1). The scanning model has no IRES, so these cannot be scored honestly."},
 {"item": "MYC P0/P1/P2 leader-specific claims", "reason": "data/evidence.json already records that the published P1/P2 leaders are not mapped to Ensembl transcripts and that MANE MYC starts at a CUG; a leader-isoform claim would be unsupported. MYC ships gene-level, tier B."},
 {"item": "MCL1 / BCL2L1 'drug shifts the ratio' claims", "reason": "The primary papers establish the isoforms and their opposite functions by overexpression; neither shows a drug shifting the endogenous ratio. Kept as splice-event targets only."},
 {"item": "specific clinical phase for zotatifin", "reason": "The two zotatifin papers read here state 'development candidate', 'first-in-class' and report patient biopsies, but no phase. Status is recorded in those words instead of a phase number."},
 {"item": "DDIT3/CHOP as a second ISR differential source", "reason": "The PMID resolved to an unrelated paper and the correct record was not re-fetched in time; the ISR differential ships on Vattem 2004 alone. data/evidence.json's verified DDIT3 entry is unchanged and not re-used here."},
 {"item": "GAPDH leader-feature citation", "reason": "A PubMed search for a GAPDH 5'UTR reporter study returned an unrelated trypanosome paper. GAPDH ships as a reference gene on the Reiter 2008 insensitivity data only, not on a leader-feature claim."},
]

PIPELINE_REQUESTS = {
 "note": "Genes in this file that are NOT in pipeline/genes.txt; the lead must add them and rebuild utr_data.json (and splice_data.json for the splice rows) before the app can score them. pipeline/genes.txt was not edited by this agent.",
 "needed_for_leader_model": ["PIM1", "ODC1", "KRAS", "GAPDH", "ACTB", "HBB"],
 "needed_for_splice_model": ["MET", "BCL2L1", "MCL1", "AR", "BRAF", "CCND1", "MDM4", "CDKN2A"],
 "already_in_genes_txt": ["MDM2", "ERBB2", "MYC"],
}

doc = {"meta": META, "levers": LEVERS, "genes": GENES, "reference_set": REFERENCE_SET,
       "differential_evidence": DIFFERENTIAL, "rejected": REJECTED, "pipeline_requests": PIPELINE_REQUESTS}

out = os.path.join(ROOT, "data", "targets.json")
with open(out, "w") as f:
    json.dump(doc, f, indent=1, ensure_ascii=False)
    f.write("\n")
shutil.copyfile(out, os.path.join(ROOT, "app", "public", "targets.json"))
print("wrote", out, os.path.getsize(out), "bytes;", len(GENES), "genes,", len(LEVERS), "levers,", len(REFERENCE_SET), "reference genes")
