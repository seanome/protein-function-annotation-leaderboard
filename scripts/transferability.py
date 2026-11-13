"""Transferability T, scored once, by one person, with the six points written out.

The five blocks of rows were collected in parallel and each wrote its own reading of
the six points, so their scores are not comparable. Every row is scored again here,
from the molecule column and the human-target column only, before the drug column was
read for that row, against the six points as the prompt defines them:

  T1 one gene, not polygenic or life-history
  T2 secreted, not intracellular
  T3 catalytic or binding, not structural or network
  T4 target conserved across phyla, not a species-specific partner
  T5 peptide or single domain
  T6 works alone, no pathway or cofactor needed

Readings fixed before scoring, so the same question is asked of every row:
  - A row with no named molecule scores 0 and the points are written "no molecule".
  - A channel, pump or transporter is not catalytic and not a binder in the sense T3
    asks about, so T3 = 0 for those.
  - "Target conserved across phyla" means a protein or small molecule a human also has.
    Ice and a snake toxin are not human molecules, so those rows score T4 = 0.
  - "Single domain" is one structural domain. A P450, an acetylcholinesterase and a
    globin are single-domain; an immunoglobulin-like five-domain inhibitor, a sodium
    channel and a sodium pump are not.
  - T6 = 0 whenever the molecule needs a cofactor (NAD+, heme, CoA), a partner protein
    or a second gene product to do its job.

The key is the disease id and the position of the row within that disease.
"""

# (disease_id, row index within the disease): (six points, total)
POINTS = {
    ("D01", 0): "T1=1 one gene; T2=1 secreted serum glycoprotein; T3=1 binds the toxin; "
                "T4=0 the target is a snake toxin, not a human molecule; T5=0 five "
                "immunoglobulin-like domains; T6=1 the recombinant protein works alone",
    ("D02", 0): "T1=1 one gene; T2=0 a membrane channel; T3=0 a channel; T4=0 the change is "
                "in the animal's own channel, there is no target to act on; T5=0 four-domain "
                "channel; T6=0 needs the whole channel in a membrane",
    ("D03", 0): "T1=1 one gene, one residue; T2=0 a membrane pump; T3=0 a pump; T4=1 the "
                "sodium pump is conserved across animals; T5=0 a large multi-domain pump; "
                "T6=0 needs a membrane",
    ("D04", 0): "no molecule scored: no measured number and no accession read",
    ("D05", 0): "T1=1 one gene, one residue; T2=0 a membrane channel; T3=0 a channel; T4=1 "
                "voltage-gated sodium channels are conserved across animals; T5=0 a large "
                "multi-domain channel; T6=0 needs a membrane",
    ("D06", 0): "no molecule",
    ("D07", 0): "T1=1 one gene, CYP6B46; T2=0 a microsomal enzyme inside the cell; T3=1 "
                "catalytic; T4=0 it acts on nicotine, a plant metabolite, not a human "
                "molecule; T5=1 one P450 domain; T6=0 needs heme and a reductase partner",
    ("D08", 0): "T1=0 the paper's own analysis puts membrane fatty acids alongside the Adh "
                "pathway; T2=0 an enzyme inside the cell; T3=1 catalytic; T4=0 no human "
                "target to act on; T5=0 two domains per subunit; T6=0 needs NAD+",
    ("D09", 0): "no molecule",
    ("D10", 0): "no molecule",
    ("D11", 0): "no molecule",
    ("D12", 0): "no molecule",
    ("D13", 0): "T1=1 one gene, TST1; T2=0 inside the cell; T3=1 catalytic, a "
                "sulfurtransferase; T4=1 cyanide and thiosulfate are the same in a human; "
                "T5=1 a rhodanese domain; T6=0 needs a sulfur donor",
    ("D14", 0): "T1=0 cox1 and cox3 plus the nuclear SQR; T2=0 inner mitochondrial membrane; "
                "T3=1 catalytic; T4=1 humans have MT-CO1 and SQOR; T5=0 a multi-subunit "
                "complex; T6=0 needs the whole respiratory chain",
    ("D15", 0): "no molecule: the paper names no gene or protein variant",
    ("D16", 0): "T1=1 one gene, ace-1; T2=1 acetylcholinesterase is secreted into the "
                "synaptic cleft; T3=1 catalytic; T4=1 humans have the same enzyme and the "
                "same substrate; T5=1 one alpha/beta-hydrolase domain; T6=1 needs no cofactor",
    ("D17", 0): "T1=0 a V-ATPase complex plus Na+/H+ exchangers; T2=0 membrane proteins; "
                "T3=1 catalytic, an ATP-driven pump; T4=1 humans have both; T5=0 a "
                "multi-subunit complex; T6=0 needs ATP and the rest of the complex",
    ("D18", 0): "T1=1 one gene; T2=0 a membrane transporter; T3=0 a transporter; T4=1 humans "
                "have SLC14 urea transporters; T5=0 a multi-pass membrane protein; T6=0 "
                "needs a membrane",
    ("D19", 0): "T1=0 several transporters, GLUT1, GLUT2, GLUT3 and GLUT5; T2=0 membrane "
                "transporters; T3=0 transporters; T4=1 humans have all four; T5=0 "
                "multi-pass membrane proteins; T6=0 needs a membrane",
    ("D20", 0): "no molecule",
    ("D21", 0): "T1=1 one gene, UOX; T2=0 uricase works in the peroxisome; T3=1 catalytic; "
                "T4=1 urate is the same molecule in a human; T5=0 a tetramer of two-domain "
                "subunits; T6=1 needs no cofactor",
    ("D22", 0): "no molecule",
    ("D23", 0): "T1=1 one gene, HAMP; T2=1 hepcidin is secreted; T3=1 it binds ferroportin; "
                "T4=1 humans have ferroportin; T5=1 a 25-residue peptide; T6=1 works alone",
    ("D24", 0): "no molecule",
    ("D25", 0): "T1=1 one gene, PPARGC1A; T2=0 a nuclear coactivator; T3=0 a coactivator, "
                "which is network, not catalysis or binding to a target; T4=1 humans have "
                "PPARGC1A; T5=0 a large multi-region protein; T6=0 needs its partner "
                "transcription factors",
    ("D26", 0): "no molecule: the buffer is shell mineral",
    ("D27", 0): "T1=1 one gene, Tac1; T2=1 substance P is released outside the cell; T3=1 it "
                "binds the NK1 receptor; T4=1 humans have TACR1; T5=1 an 11-residue peptide; "
                "T6=1 works alone",
    ("D28", 0): "T1=1 one gene, ATG5; T2=0 autophagy machinery inside the cell; T3=0 a "
                "conjugation partner, which is network; T4=1 humans have ATG5; T5=0 more "
                "than one domain; T6=0 needs ATG12 and ATG16L1",
    ("D29", 0): "T1=0 two bacterial genes, frc and oxc; T2=0 both work inside the bacterium; "
                "T3=1 both are catalytic; T4=1 oxalate is the same molecule in a human; T5=0 "
                "multi-domain enzymes; T6=0 needs both enzymes and CoA",
    ("D30", 0): "T1=1 one gene, SDC4; T2=1 a cell-surface proteoglycan; T3=0 structural; "
                "T4=1 humans have SDC4; T5=0 a proteoglycan with a long glycan chain; T6=0 "
                "needs its glycan chains and membrane anchor",
    ("D31", 0): "T1=0 two genes, SLC2A5 and KHK; T2=0 a transporter and an intracellular "
                "kinase; T3=1 ketohexokinase is catalytic; T4=1 humans have both; T5=0 "
                "multi-domain; T6=0 needs the pair plus ATP",
    ("D31", 1): "no molecule: the paper measures lactate and shell buffering",
    ("D32", 0): "T1=1 one gene, SIRT3; T2=0 a mitochondrial enzyme; T3=1 catalytic; T4=1 "
                "humans have SIRT3; T5=0 more than one domain; T6=0 needs NAD+",
    ("D33", 0): "no molecule",
    ("D34", 0): "T1=1 one gene, GRIN2D; T2=0 a membrane receptor subunit; T3=0 a channel "
                "subunit; T4=1 humans have GRIN2D; T5=0 a large multi-domain subunit; T6=0 "
                "needs the rest of the receptor",
    ("D35", 0): "no molecule: the paper rules mechanisms out",
    ("D36", 0): "T1=1 one gene, COX3; T2=0 inner mitochondrial membrane; T3=0 a subunit of a "
                "complex, not catalytic on its own; T4=1 humans have MT-CO3; T5=0 a "
                "multi-pass membrane subunit; T6=0 needs the whole oxidase",
    ("D37", 0): "T1=1 one gene, ROCK2; T2=0 a cytoplasmic kinase; T3=1 catalytic; T4=1 humans "
                "have ROCK2; T5=0 a large multi-domain kinase; T6=0 needs ATP and its "
                "upstream Rho partner",
    ("D38", 0): "T1=1 one gene, MB; T2=0 myoglobin works inside the muscle cell; T3=1 it "
                "binds oxygen; T4=0 its partner is oxygen, not a human protein; T5=1 one "
                "globin domain; T6=0 needs heme",
    ("D39", 0): "T1=0 XIAP together with Akt1; T2=0 both inside the cell; T3=1 XIAP binds "
                "caspases; T4=1 humans have XIAP and the same caspases; T5=0 several BIR "
                "domains plus a RING domain; T6=0 needs the pair",
    ("D40", 0): "no molecule",
    ("D41", 0): "T1=1 one gene; T2=1 secreted into the hemolymph; T3=1 it binds ice crystal "
                "faces; T4=0 its target is ice, not a human molecule; T5=1 one small "
                "beta-helical domain; T6=1 works alone",
    ("D42", 0): "no molecule: the protective agent is glucose, a metabolite",
    ("D42", 1): "T1=1 one gene; T2=1 added to the medium outside the cell; T3=1 it binds ice; "
                "T4=0 its target is ice, not a human molecule; T5=1 one ice-binding domain; "
                "T6=1 the recombinant protein works alone",
    ("D43", 0): "no molecule: the paper reports whole-animal metabolism",
    ("D44", 0): "no molecule",
    ("D45", 0): "T1=1 one gene, p26; T2=0 a chaperone inside the cell; T3=1 it binds "
                "unfolding proteins; T4=0 its partners are the cell's own proteins; T5=1 one "
                "alpha-crystallin domain; T6=0 the human-cell result needed trehalose as well",
    ("D46", 0): "T1=1 one gene, Dsup; T2=0 it binds chromatin inside the nucleus; T3=1 it "
                "binds nucleosomes; T4=1 the nucleosome is the same in a human cell; T5=1 one "
                "disordered nucleosome-binding protein; T6=1 it works alone in human cells",
    ("D47", 0): "no molecule tied to the measured number",
    ("D48", 0): "no molecule: the paper measures bone structure and strength",
    ("D49", 0): "no molecule: the active agent is whole serum",
    ("D50", 0): "T1=1 one gene, SERPINH1; T2=0 an endoplasmic-reticulum chaperone; T3=1 it "
                "binds collagen; T4=1 humans have SERPINH1 and the same collagens; T5=0 a "
                "serpin fold with more than one region; T6=0 needs the collagen-folding "
                "machinery",
    ("D51", 0): "no molecule",
    ("D52", 0): "T1=1 one gene, ttk; T2=0 an intracellular kinase; T3=1 catalytic; T4=1 "
                "humans have TTK; T5=0 a large multi-domain kinase; T6=0 needs ATP and the "
                "rest of the checkpoint",
    ("D53", 0): "T1=1 one gene, ctgfa; T2=1 secreted; T3=1 it binds signalling partners; T4=1 "
                "humans have CCN2 and the four domains are 81% identical; T5=0 four domains, "
                "though the C-terminal part alone carries the effect; T6=1 the recombinant "
                "protein placed at the lesion worked alone",
    ("D54", 0): "no causal molecule: the proteins named are markers",
    ("D55", 0): "no molecule",
    ("D56", 0): "no molecule",
    ("D57", 0): "no molecule",
    ("D58", 0): "no molecule shown to cause replacement",
    ("D59", 0): "T1=1 one gene, BLMH; T2=0 inside the cell; T3=1 catalytic; T4=0 it acts on "
                "the drug bleomycin, not on a human molecule of fibrosis; T5=0 a papain-fold "
                "hexamer; T6=0 the paper itself doubts this is the cause",
    ("D60", 0): "T1=0 about 20 gene copies, not one gene; T2=0 an intracellular transcription "
                "factor; T3=0 the retrogene products are not shown to act as transcription "
                "factors; T4=1 humans have TP53; T5=0 a multi-domain protein; T6=0 needs the "
                "rest of the p53 response",
}

# Rows where the human-target column as collected starts with "yes" but does not mean
# "a paper shows this molecule acting in a human cell, on a human receptor or in a human
# assay" in the direction the step asks about.
HUMAN_OVERRIDE = {
    ("D22", 0): (
        "no. The human-cell result here runs the other way: spisulosine "
        "(1-deoxysphinganine), from the clam Spisula polynyma, kills human PC-3 and LNCaP "
        "cells at about 1 micromolar by raising ceramide. That is the lipid being toxic to "
        "human cells, not a tolerance molecule working in them."
    ),
}

# Rows where the drug column as collected does not mean "the molecule from the tolerant
# species became an approved drug". The text replaces the column; the fetched agency page
# is kept inside it, because the fact is real and only the reading changes.
DRUG_OVERRIDE = {
    ("D27", 0): (
        "no. Aprepitant (EMEND, US FDA NDA 021549, approved 26 March 2003, page fetched "
        "2026-10-19 at https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event="
        "overview.process&ApplNo=021549) blocks the human receptor for this molecule. It is "
        "not the naked mole-rat molecule made into a drug, which is what step 5 asks for, "
        "and it was developed without reference to the naked mole rat."
    ),
}
