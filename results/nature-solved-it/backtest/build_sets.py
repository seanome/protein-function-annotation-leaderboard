"""Build the transferability back-test sets from the records read on 2026-10-19.

Every row's facts come from a source opened on 2026-10-19: Drugs@FDA and the FDA
label through the openFDA API, ClinicalTrials.gov through its v2 API, and PubMed
through the PubMed tool. No year, DOI, PMID, application number or trial
identifier here was written from memory.

Outputs, all in this directory:
  sets.csv         one row per entry, with the set it belongs to
  excluded.csv     candidates checked and left out, with the reason
  blind_input.csv  opaque id and mechanism description only, shuffled
  key.csv          blind id to entry id

Run: python3 build_sets.py
"""

import csv
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHUFFLE_SEED = 20261019
OPENED = "2026-10-19"

# Citations used by more than one row, written once.
C_AVERIN = ("Averin, A. S., and Yu. N. Utkin. \"Cardiovascular Effects of Snake Toxins: "
            "Cardiotoxicity and Cardioprotection.\" Acta Naturae 13, no. 3 (2021): 4-14. "
            "https://doi.org/10.32607/actanaturae.11375. PMID 34707893.")
C_WHITE = ("White, Harvey D., and Derek P. Chew. \"Bivalirudin: An Anticoagulant for Acute "
           "Coronary Syndromes and Coronary Interventions.\" Expert Opinion on Pharmacotherapy "
           "3, no. 6 (2002): 777-88. https://doi.org/10.1517/14656566.3.6.777. PMID 12036417.")
C_BELL = ("Bell, W. R. \"Defibrinogenating Enzymes.\" Drugs 54, Suppl. 3 (1997): 18-30. "
          "https://doi.org/10.2165/00003495-199700543-00005. PMID 9360849.")

FDA = "Drugs@FDA record read through the openFDA drugsfda API on 2026-10-19"
LABEL = "FDA product label read through the openFDA label API on 2026-10-19"
CTG = "ClinicalTrials.gov study record read through the v2 API on 2026-10-19"

# set, id, name, species, molecule, human target, year, year_type, agency_or_registry,
# record_id, citation, source_type, mechanism (blind text), notes
ROWS = [
    # ---------------- Set A: approved ----------------
    dict(set="A", id="A01", name="captopril",
         species="Bothrops jararaca, a pit viper (the paper names a snake venom "
                 "bradykinin-potentiating peptide; PubMed indexes the record under "
                 "Bothrops jararaca)",
         molecule="bradykinin-potentiating peptide of pit viper venom; the marketed "
                  "molecule is a small synthetic inhibitor modelled on that peptide",
         target="angiotensin-converting enzyme",
         year=1981, year_type="FDA approval", agency="FDA", record="NDA018343 (1981-04-06)",
         citation=C_AVERIN + " Approval: " + FDA, source_type="review + agency record",
         mech="A short peptide found in the venom of a pit viper, secreted into another "
              "animal's bloodstream, that blocks a single human blood-pressure enzyme "
              "(angiotensin-converting enzyme) at its active site. A small synthetic "
              "molecule copying the peptide's contact points does the same. One gene "
              "makes the peptide; it needs no partner protein and no cofactor.",
         notes="The species is from the indexed record, not from the abstract text."),
    dict(set="A", id="A02", name="exenatide",
         species="Heloderma suspectum", molecule="exendin-4, a 39-residue peptide of the "
         "lizard's oral secretion", target="GLP-1 receptor",
         year=2005, year_type="FDA approval", agency="FDA", record="NDA021773 (2005-04-28)",
         citation="FDA label for BYETTA (exenatide): \"a synthetic peptide, GLP-1 receptor "
                  "agonist, that was originally identified in the lizard Heloderma "
                  "suspectum.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A 39-residue peptide from the oral secretion of a venomous lizard. It is "
              "released outside the animal's cells, and it switches on one human cell-surface "
              "receptor for a gut hormone (the GLP-1 receptor) by binding it. One gene, one "
              "chain, no cofactor, no partner protein.",
         notes=""),
    dict(set="A", id="A03", name="ziconotide",
         species="Conus magus", molecule="omega-conotoxin MVIIA, a 25-residue "
         "disulfide-bonded venom peptide", target="N-type (Cav2.2) voltage-gated calcium channel",
         year=2004, year_type="FDA approval", agency="FDA", record="NDA021060 (2004-12-28)",
         citation="FDA label for PRIALT: \"a synthetic equivalent of a naturally occurring "
                  "conopeptide found in the piscivorous marine snail, Conus magus.\" "
                  + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A 25-residue disulfide-bonded peptide from the venom of a fish-hunting marine "
              "snail, injected into prey. It blocks the pore of one human voltage-gated "
              "calcium channel in nerve endings. One gene, one small folded chain, works on "
              "its own.",
         notes=""),
    dict(set="A", id="A04", name="eptifibatide",
         species="southeastern pigmy rattlesnake (the source gives the common name only)",
         molecule="barbourin, a venom peptide; the marketed molecule is a cyclic "
                  "heptapeptide copying its binding loop",
         target="platelet integrin glycoprotein IIb/IIIa",
         year=1998, year_type="FDA approval", agency="FDA", record="NDA020718 (1998-05-18)",
         citation="Scarborough, R. M. \"Development of Eptifibatide.\" American Heart Journal "
                  "138, no. 6, pt. 1 (1999): 1093-1104. "
                  "https://doi.org/10.1016/s0002-8703(99)70075-x. PMID 10577440. Approval: " + FDA,
         source_type="review + agency record",
         mech="A venom peptide of a small rattlesnake that binds the fibrinogen-binding site "
              "of one human platelet surface receptor (an integrin) and so stops platelets "
              "sticking together. A seven-residue ring copying the "
              "peptide's binding loop has the same action. Secreted, one gene, no cofactor.",
         notes="The source names the snake by common name; no binomial was read today."),
    dict(set="A", id="A05", name="bivalirudin",
         species="Hirudo medicinalis", molecule="hirudin of leech saliva; the marketed "
         "molecule is a 20-residue synthetic peptide copying its two thrombin contacts",
         target="thrombin", year=2000, year_type="FDA approval", agency="FDA",
         record="NDA020873 (2000-12-15)",
         citation=C_WHITE + " Approval: " + FDA, source_type="review + agency record",
         mech="A peptide copying the two thrombin-contacting ends of a protein in the saliva "
              "of a blood-feeding leech. It is secreted into the host's blood and blocks one "
              "human clotting enzyme, thrombin, by binding both its active site and its "
              "substrate groove. Twenty residues, one chain, no partner needed.",
         notes=""),
    dict(set="A", id="A06", name="desirudin",
         species="Hirudo medicinalis", molecule="recombinant hirudin, a 65-residue leech "
         "salivary protein", target="thrombin",
         year=2003, year_type="FDA approval", agency="FDA", record="BLA021271 (2003-04-04)",
         citation=C_WHITE + " Approval: " + FDA, source_type="review + agency record",
         mech="The whole salivary anticoagulant protein of a blood-feeding leech, made in "
              "bacteria: a single 65-residue secreted chain that binds one human clotting "
              "enzyme, thrombin, and holds it inactive. One gene, no cofactor, no partner.",
         notes=""),
    dict(set="A", id="A07", name="defibrotide",
         species="pig (the label says porcine intestinal tissue; no binomial given)",
         molecule="mixture of single-stranded polydeoxyribonucleotides, mean weighted mass "
                  "14-19 kDa",
         target="vascular endothelium (the label describes endothelial protection; no single "
                "receptor is named)",
         year=2016, year_type="FDA approval", agency="FDA", record="NDA208114 (2016-03-30)",
         citation="FDA label for DEFITELIO: \"a polydisperse mixture of predominantly "
                  "single-stranded (ss) polydeoxyribonucleotide sodium salts derived from "
                  "porcine intestinal tissue.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A mixture of single-stranded DNA fragments, 14-19 kDa on average, extracted "
              "from pig gut tissue. It acts on the lining cells of blood vessels from the "
              "outside; no one receptor is named for it. It is not the product of a single "
              "gene, and its effect is on a cell layer rather than on one binding partner.",
         notes=""),
    dict(set="A", id="A08", name="heparin sodium",
         species="pig (the label says porcine intestinal mucosa; no binomial given)",
         molecule="heparin, a sulfated glycosaminoglycan chain",
         target="antithrombin III, which it activates to shut down thrombin and factor Xa",
         year=1972, year_type="FDA approval", agency="FDA", record="NDA017029 (1972-02-22)",
         citation="FDA label for heparin sodium injection: \"heparin sodium derived from "
                  "porcine intestinal mucosa.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A long sulfated sugar chain purified from pig gut lining, not a protein and not "
              "the product of one gene. Outside the cell it binds a human blood protein "
              "(antithrombin) and makes it trap clotting enzymes much faster. It works only "
              "through that partner protein.",
         notes=""),
    dict(set="A", id="A09", name="protamine sulfate",
         species="salmon (the label says the sperm of salmon and certain other fish; no "
                 "binomial given)",
         molecule="protamine, a short arginine-rich sperm nuclear peptide",
         target="heparin, which it binds and neutralises",
         year=1987, year_type="FDA approval (abbreviated application)", agency="FDA",
         record="ANDA089454 (1987-04-07)",
         citation="FDA label for protamine sulfate: \"They occur in the sperm of salmon and "
                  "certain other species of fish.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A short, strongly positive arginine-rich peptide from fish sperm nuclei, where "
              "its normal job is packing DNA. In blood it works "
              "by sticking to a negatively charged sugar chain and cancelling it. One gene, "
              "small peptide, no cofactor, but its action is charge binding rather than "
              "catalysis on a protein target.",
         notes="The record is an abbreviated (generic) application, not an original NDA."),
    dict(set="A", id="A10", name="calcitonin salmon",
         species="salmon (the label says the sequence found in calcitonin of salmon origin; "
                 "no binomial given)",
         molecule="salmon calcitonin, a 32-residue peptide hormone",
         target="calcitonin receptor on bone-resorbing cells",
         year=1986, year_type="FDA approval", agency="FDA", record="NDA017808 (1986-07-03)",
         citation="FDA label for calcitonin salmon injection: \"a synthetic polypeptide of 32 "
                  "amino acids in the same linear sequence that is found in calcitonin of "
                  "salmon origin.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A 32-residue hormone that a fish secretes into its own blood to control "
              "calcium. It binds the matching human cell-surface receptor on bone cells. One "
              "gene, one small secreted chain, no cofactor.",
         notes=""),
    dict(set="A", id="A11", name="hyaluronidase, bovine",
         species="cattle (the label says bovine testicular hyaluronidase; no binomial given)",
         molecule="testicular hyaluronidase, an enzyme",
         target="hyaluronan in human extracellular matrix",
         year=2004, year_type="FDA approval", agency="FDA", record="BLA021665 (2004-10-26)",
         citation="FDA label for AMPHADASE: \"a preparation of purified bovine testicular "
                  "hyaluronidase, a protein enzyme.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="An enzyme purified from bull testes that cuts the hyaluronan chains of the "
              "jelly between human cells, so fluid moves more freely through the tissue. It is secreted, acts "
              "outside cells, is a single catalytic protein, and needs no cofactor.",
         notes=""),
    dict(set="A", id="A12", name="hyaluronidase, ovine",
         species="sheep (Ovis aries in the cited paper's terms: sheep testicular "
                 "hyaluronidase)",
         molecule="sheep testicular hyaluronidase, two glycoforms of one enzyme",
         target="hyaluronan in human extracellular matrix",
         year=2004, year_type="FDA approval", agency="FDA", record="BLA021640 (2004-05-05)",
         citation="Yang, Subin, Mijung Jin, Chi Soo Park, Chulmin Moon, Mirae Kim, Jieun Kim, "
                  "Leeseul Jang, et al. \"Identification, Quantification, and Structural Role "
                  "of N-Glycans in Two Highly Purified Isoforms of Sheep Testicular "
                  "Hyaluronidase.\" International Journal of Biological Macromolecules 252 "
                  "(2023): 126437. https://doi.org/10.1016/j.ijbiomac.2023.126437. "
                  "PMID 37611686. Approval: " + FDA,
         source_type="primary paper + agency record",
         mech="An enzyme purified from ram testes, in two sugar-decorated forms of the same "
              "protein, that cuts the hyaluronan chains between human cells. Secreted, one "
              "catalytic chain, no cofactor; the sugar groups are part of the molecule.",
         notes=""),
    dict(set="A", id="A13", name="pancrelipase",
         species="pig (the label says porcine pancreatic glands; no binomial given)",
         molecule="mixture of lipases, proteases and amylases from pancreas extract",
         target="dietary fat, protein and starch in the human gut",
         year=2009, year_type="FDA approval", agency="FDA", record="BLA020725 (2009-04-30)",
         citation="FDA label for CREON: \"a pancreatic enzyme product consisting of a mixture "
                  "of enzymes including lipases, proteases, and amylases and is an extract "
                  "derived from porcine pancreatic glands.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A mixture of digestive enzymes extracted from pig pancreas. The enzymes are "
              "secreted ones and each is catalytic, but the product is many gene products "
              "together, and what they act on is food in the gut rather than a human protein.",
         notes=""),
    dict(set="A", id="A14", name="beractant",
         species="cattle (the label says natural bovine lung extract; no binomial given)",
         molecule="lung surfactant extract: phospholipids, neutral lipids, fatty acids and "
                  "surfactant proteins, with three lipids added",
         target="the air-liquid surface of human alveoli",
         year=1991, year_type="FDA approval", agency="FDA", record="BLA020032 (1991-07-01)",
         citation="FDA label for SURVANTA: \"a natural bovine lung extract containing "
                  "phospholipids, neutral lipids, fatty acids, and surfactant-associated "
                  "proteins.\" " + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A lipid-and-protein extract of cow lung, acting at the air-liquid surface of the lung. "
              "It works by lowering surface tension at the air-liquid boundary, a physical "
              "effect of the mixture, not by binding a receptor. Many components, several "
              "genes, and added lipids are needed for it to work.",
         notes=""),
    dict(set="A", id="A15", name="calfactant",
         species="cattle, calf (the label says natural bovine lung extract; no binomial given)",
         molecule="lung surfactant extract: 85% phospholipids, 8% neutral lipids, 7% "
                  "surfactant proteins B and C",
         target="the air-liquid surface of human alveoli",
         year=1998, year_type="FDA approval", agency="FDA", record="BLA020521 (1998-07-01)",
         citation="FDA label for INFASURF: \"an extract of natural bovine lung (pulmonary) "
                  "surfactant consisting of 85% phospholipids, 8% neutral lipids, and 7% "
                  "hydrophobic surfactant-associated proteins B and C.\" " + LABEL
                  + ". Approval: " + FDA,
         source_type="agency record",
         mech="A lipid-and-protein extract of calf lung, acting at the air-liquid surface of the lung. Two small "
              "hydrophobic proteins work together with the lipids to spread a film and lower "
              "surface tension. The effect is physical, needs the lipid partners, and comes "
              "from more than one gene.",
         notes=""),
    dict(set="A", id="A16", name="poractant alfa",
         species="pig (the label says natural porcine lung surfactant; no binomial given)",
         molecule="lung surfactant extract: 99% polar lipids, 1% hydrophobic surfactant "
                  "proteins",
         target="the air-liquid surface of human alveoli",
         year=1999, year_type="FDA approval", agency="FDA", record="BLA020744 (1999-11-18)",
         citation="FDA label for CUROSURF: \"an extract of natural porcine lung (pulmonary) "
                  "surfactant consisting of 99% polar lipids (mainly phospholipids) and 1% "
                  "hydrophobic low molecular weight surfactant associated proteins.\" "
                  + LABEL + ". Approval: " + FDA,
         source_type="agency record",
         mech="A lipid-and-protein extract of pig lung, acting at the air-liquid surface of the lung. Almost all of it "
              "is lipid; the small proteins only help the film spread. The action is physical "
              "surface tension lowering, needs the lipids, and is not one gene's product.",
         notes=""),
    dict(set="A", id="A17", name="eribulin mesylate",
         species="Halichondria okadai", molecule="synthetic analogue of halichondrin B, a "
         "macrocyclic polyketide of a marine sponge",
         target="tubulin, at the microtubule plus end",
         year=2010, year_type="FDA approval", agency="FDA", record="NDA201532 (2010-11-15)",
         citation="FDA label for eribulin mesylate: \"a synthetic analogue of halichondrin B, "
                  "a product isolated from the marine sponge Halichondria okadai.\" " + LABEL
                  + ". Approval: " + FDA,
         source_type="agency record",
         mech="A small molecule, not a protein, copying part of a large chemical made by a "
              "marine sponge. It binds tubulin inside human cells and stops microtubules "
              "growing. It is made by a multi-enzyme pathway in the animal, so no single gene "
              "carries it, but the molecule itself acts alone on one conserved target.",
         notes=""),
    dict(set="A", id="A18", name="trabectedin",
         species="a Caribbean tunicate (the paper read today gives no binomial; PubMed matched "
                 "the record to Ecteinascidia turbinata)",
         molecule="ET-743, a marine alkaloid with three fused tetrahydroisoquinoline rings",
         target="the DNA minor groove in human cells",
         year=2015, year_type="FDA approval", agency="FDA", record="NDA207953 (2015-10-23)",
         citation="Henry, Theophraste, Elizabeth Fabre, Laurent S. Baccar, and Michele "
                  "Lamuraglia. \"Longer Survival in Patients with Metastatic Uterine "
                  "Leiomyosarcoma Treated with Trabectedin: A Case Report.\" Molecular and "
                  "Clinical Oncology 10, no. 3 (2019): 387-90. "
                  "https://doi.org/10.3892/mco.2019.1797. PMID 30847179. Approval: " + FDA,
         source_type="primary paper + agency record",
         mech="A small ring-shaped chemical made by a sea squirt. It slips into the groove of "
              "DNA inside human cells and blocks transcription and repair. Not a protein and "
              "not one gene's product; the molecule itself needs no partner.",
         notes="The species binomial was not read in a source today."),
    dict(set="A", id="A19", name="omega-3-acid ethyl esters",
         species="fish (the label says omega-3 fatty acids sourced from fish oils; no species "
                 "named)",
         molecule="ethyl esters of eicosapentaenoic and docosahexaenoic acid",
         target="triglyceride handling in the liver (no single receptor named on the label)",
         year=2004, year_type="FDA approval", agency="FDA", record="NDA021654 (2004-11-10)",
         citation="FDA label for omega-3-acid ethyl esters: \"at least 900 mg of the ethyl "
                  "esters of omega-3 fatty acids sourced from fish oils.\" " + LABEL
                  + ". Approval: " + FDA,
         source_type="agency record",
         mech="A mixture of fatty acid esters from fish oil, absorbed from the gut. It changes how "
              "the liver handles fat. No single gene makes it, no receptor is named for it, "
              "and the effect is on a whole metabolic route rather than on one binding site.",
         notes=""),

    # ---------------- Set B(i): reached a registered trial, no approval ----------------
    dict(set="B_i", id="B101", name="dalazatide (ShK-186)",
         species="Stichodactyla helianthus", molecule="ShK, a 35-residue sea anemone venom "
         "peptide; the trial molecule is a synthetic analogue",
         target="Kv1.3 potassium channel on human effector memory T cells",
         year=2013, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT02446340 (start 2013-11, phase 1, completed)",
         citation="Chi, Victor, Michael W. Pennington, Raymond S. Norton, Eric J. Tarcha, Luz "
                  "M. Londono, Brian Sims-Fahey, Sanjeev K. Upadhyay, et al. \"Development of "
                  "a Sea Anemone Toxin as an Immunomodulator for Therapy of Autoimmune "
                  "Diseases.\" Toxicon 59, no. 4 (2011): 529-46. "
                  "https://doi.org/10.1016/j.toxicon.2011.07.016. PMID 21867724. Trial: " + CTG,
         source_type="review + registry record",
         mech="A small disulfide-bonded peptide from the venom of a sea anemone, released "
              "outside the animal, that blocks one human potassium channel on the surface of "
              "a subset of T cells. About 35 residues, a single folded domain, no cofactor and "
              "no partner protein.",
         notes=""),
    dict(set="B_i", id="B102", name="chlorotoxin (TM-601 and later conjugates)",
         species="scorpion; PubMed indexes the record under Leiurus quinquestriatus and the "
                 "abstract text read today drops the species name",
         molecule="chlorotoxin, a 36-residue disulfide-rich venom peptide",
         target="neuropilin-1 on human glioma cells (the paper read today), and matrix "
                "metalloproteinase 2 in earlier work",
         year=2002, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00040573 (start 2002-06, phase 1/2, completed)",
         citation="Sharma, Gagan, Carolyne B. Braga, Kai-En Chen, Xinying Jia, Venkatraman "
                  "Ramanujam, Brett M. Collins, Roberto Rittner, and Mehdi Mobli. \"Structural "
                  "Basis for the Binding of the Cancer Targeting Scorpion Toxin, ClTx, to the "
                  "Vascular Endothelia Growth Factor Receptor Neuropilin-1.\" Current Research "
                  "in Structural Biology 3 (2021): 179-86. "
                  "https://doi.org/10.1016/j.crstbi.2021.07.003. PMID 34401749. Trial: " + CTG,
         source_type="primary paper + registry record",
         mech="A 36-residue disulfide-rich peptide from scorpion venom that binds a receptor "
              "on the surface of human brain tumour cells and so marks them. One small folded "
              "domain, secreted by the animal, works without partners; it is often "
              "joined to a label or a cell.",
         notes=""),
    dict(set="B_i", id="B103", name="ancrod",
         species="Calloselasma rhodostoma", molecule="ancrod, a venom serine protease",
         target="human fibrinogen, from which it cleaves the A fibrinopeptides",
         year=2005, year_type="earliest matching trial start", agency="ClinicalTrials.gov",
         record="NCT00141011 (start 2005-09, phase 3, terminated)",
         citation=C_BELL + " Trial: " + CTG, source_type="review + registry record",
         mech="A snake venom enzyme that cuts one end off human fibrinogen, stripping clotting "
              "protein out of the blood. It is secreted, catalytic, a single chain of moderate "
              "size, and needs no cofactor.",
         notes="The registry text search also returned unrelated trials; the row uses the "
               "trials whose intervention is this molecule."),
    dict(set="B_i", id="B104", name="alfimeprase",
         species="Agkistrodon contortrix contortrix", molecule="truncated recombinant fibrolase, "
         "a zinc metalloproteinase of snake venom",
         target="the fibrin(ogen) A-alpha chain in human clots",
         year=2003, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00073554 (start 2003-06, phase 2, completed)",
         citation="Deitcher, Steven R., and Christopher F. Toombs. \"Non-Clinical and Clinical "
                  "Characterization of a Novel Acting Thrombolytic: Alfimeprase.\" "
                  "Pathophysiology of Haemostasis and Thrombosis 34, nos. 4-5 (2005): 215-20. "
                  "https://doi.org/10.1159/000092427. PMID 16707931. Trial: " + CTG,
         source_type="review + registry record",
         mech="A shortened form of a zinc-dependent enzyme from snake venom that cuts the "
              "fibrin of a human clot directly. Secreted, catalytic, one chain; it needs a "
              "zinc ion, and a common blood protein inhibitor switches it off outside the clot.",
         notes=""),
    dict(set="B_i", id="B105", name="desmoteplase",
         species="Desmodus rotundus", molecule="salivary plasminogen activator, a protease of "
         "vampire bat saliva", target="human plasminogen, activated only on fibrin",
         year=2001, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00638781 (start 2001-03, phase 2, completed)",
         citation="Patel, Rikesh, Sevasti Ispoglou, and Stavros Apostolakis. \"Desmoteplase as "
                  "a Potential Treatment for Cerebral Ischaemia.\" Expert Opinion on "
                  "Investigational Drugs 23, no. 6 (2014): 865-73. "
                  "https://doi.org/10.1517/13543784.2014.911285. PMID 24766516. Trial: " + CTG,
         source_type="review + registry record",
         mech="An enzyme in the saliva of a blood-feeding bat that switches on the human "
              "clot-dissolving protein plasminogen, and only does so where fibrin is present. "
              "Secreted, catalytic, one gene; it acts through a human partner protein, so it "
              "needs that protein to be there.",
         notes=""),
    dict(set="B_i", id="B106", name="recombinant nematode anticoagulant protein c2",
         species="hookworm (Ancylostoma species; the registry record names the protein but not "
                 "the species, and no species source was read today)",
         molecule="nematode anticoagulant protein c2, a small secreted inhibitor protein",
         target="the human tissue factor-factor VIIa complex",
         year=2005, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00116012 (start 2005-06, phase 2, completed)",
         citation="ClinicalTrials.gov record NCT00116012: \"Recombinant nematode anticoagulant "
                  "protein c2 (rNAPc2) is a novel and potent inhibitor of fVIIa/TF.\" " + CTG,
         source_type="registry record",
         mech="A small protein secreted by a blood-feeding worm into its host, which blocks "
              "the first step of human clotting by holding a two-protein starter complex "
              "inactive. One gene, small single domain, no cofactor; it binds a complex rather "
              "than a single protein.",
         notes="Species not confirmed in a source opened today."),
    dict(set="B_i", id="B107", name="pexiganan (MSI-78)",
         species="Xenopus laevis", molecule="magainin, a 21-27 residue skin peptide; the trial "
         "molecule is a 22-residue analogue",
         target="bacterial membranes on human skin (no human protein target)",
         year=1994, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00563394 (start 1994-08, phase 3, completed)",
         citation="Jacob, L., and M. Zasloff. \"Potential Therapeutic Applications of Magainins "
                  "and Other Antimicrobial Agents of Animal Origin.\" Ciba Foundation "
                  "Symposium 186 (1994): 197-216. https://doi.org/10.1002/9780470514658.ch12. "
                  "PMID 7768152. Trial: " + CTG,
         source_type="review + registry record",
         mech="A short positively charged peptide from frog skin secretion, about 22 residues "
              "in the version used, which kills bacteria by breaking their membranes. One "
              "gene, secreted, small, no cofactor; what it acts on is a microbial membrane, "
              "not a protein that humans and other animals share.",
         notes=""),
    dict(set="B_i", id="B108", name="crotoxin",
         species="Crotalus durissus terrificus", molecule="crotoxin, a two-subunit "
         "phospholipase A2 complex of rattlesnake venom",
         target="tumour cell membranes; the cited paper reports necrosis of human breast "
                "cancer cells and links activity to epidermal growth factor receptor levels",
         year=2026, year_type="trial start as registered", agency="ClinicalTrials.gov",
         record="NCT01481532 (start recorded as 2026-10, phase 1, not yet recruiting)",
         citation="De Andrade, Camila Marques, Cleni Mara Marzocchi Machado, Fabiola Attie de "
                  "Castro, Natercia Teixeira, Cristina Amaral, Felipe Campos, Nathalia "
                  "Cristina Canicoba, and Maria Regina Torqueti. \"Antitumor and "
                  "Antiproliferative Potential of Crotoxin in Triple Negative Breast Tumors.\" "
                  "Toxicon 258 (2025): 108322. https://doi.org/10.1016/j.toxicon.2025.108322. "
                  "PMID 40120860. Trial: " + CTG,
         source_type="primary paper + registry record",
         mech="A venom toxin of a rattlesnake made of two different subunits that must stay "
              "together: a membrane-cutting enzyme and a chaperone subunit that aims it. "
              "Secreted, catalytic, but it needs its partner subunit and calcium, and its "
              "effect is on cell membranes broadly.",
         notes="The registered start date is in the future as the record stands today."),
    dict(set="B_i", id="B109", name="SOR-C13",
         species="Blarina brevicauda", molecule="C-terminal 13-residue fragment of soricidin, a "
         "54-residue paralytic venom peptide",
         target="TRPV6 calcium channel on human tumour cells",
         year=2012, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT01578564 (start 2012-07, phase 1, completed)",
         citation="Bowen, Chris V., Drew DeBay, H. Stephen Ewart, Pamela Gallant, Sean Gormley, "
                  "T. Toney Ilenchuk, Umar Iqbal, et al. \"In Vivo Detection of Human "
                  "TRPV6-Rich Tumors with Anti-Cancer Peptides Derived from Soricidin.\" PLoS "
                  "ONE 8, no. 3 (2013): e58866. https://doi.org/10.1371/journal.pone.0058866. "
                  "PMID 23554944. Trial: " + CTG,
         source_type="primary paper + registry record",
         mech="A 13-residue piece of a peptide found in the paralytic saliva of a shrew. It "
              "blocks one human calcium channel that many tumours carry a lot of. One gene, "
              "secreted, very small, works alone.",
         notes=""),
    dict(set="B_i", id="B110", name="cenderitide (CD-NP)",
         species="Dendroaspis angusticeps (the cited paper names Dendroaspis natriuretic "
                 "peptide; PubMed indexes the record under Dendroaspis)",
         molecule="chimera: human C-type natriuretic peptide with the carboxyl tail of the "
                  "snake venom natriuretic peptide",
         target="the particulate guanylyl cyclase A and B receptors in human cells",
         year=2007, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00482937 (start 2007-01, phase 1, completed)",
         citation="Lee, Candace Y. W., Brenda K. Huntley, Daniel J. McCormick, Tomoko Ichiki, "
                  "S. Jeson Sangaralingham, Ondrej Lisy, and John C. Burnett. \"Cenderitide: "
                  "Structural Requirements for the Creation of a Novel Dual Particulate "
                  "Guanylyl Cyclase Receptor Agonist with Renal-Enhancing in Vivo and ex Vivo "
                  "Actions.\" European Heart Journal: Cardiovascular Pharmacotherapy 2, no. 2 "
                  "(2015): 98-105. https://doi.org/10.1093/ehjcvp/pvv040. PMID 27340557. "
                  "Trial: " + CTG,
         source_type="primary paper + registry record",
         mech="A short designed peptide: a human hormone peptide with the tail of a snake "
              "venom peptide stuck on the end. It switches on two human receptor enzymes at "
              "the cell surface. Small, secreted, no cofactor; only part of it comes from the "
              "animal.",
         notes="Half the molecule is a human peptide; the animal part is the venom tail."),
    dict(set="B_i", id="B111", name="batroxobin (hemocoagulase)",
         species="Bothrops atrox moojeni", molecule="batroxobin, a venom serine protease",
         target="human fibrinogen, from which it cleaves both fibrinopeptides",
         year=2020, year_type="earliest trial start found", agency="ClinicalTrials.gov",
         record="NCT04269954 (start 2020-03-01, phase 4, status unknown)",
         citation=C_BELL + " Trial: " + CTG, source_type="review + registry record",
         mech="A snake venom enzyme that cuts human fibrinogen into fibrin, used to strip "
              "clotting protein from the blood. Secreted, catalytic, a single chain, no "
              "cofactor.",
         notes="Marketed in some countries outside the United States and the European Union; "
               "no FDA or EMA approval was found today."),
    dict(set="B_i", id="B112", name="omiganan (MBI 226)",
         species="cattle (indolicidin comes from bovine neutrophils; the review read today "
                 "names bovine indolicidin, no binomial)",
         molecule="omiganan, a cationic peptide derived from indolicidin",
         target="microbial membranes on human skin (no human protein target)",
         year=2000, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00027248 (start 2000-09, phase 3, completed)",
         citation="Wiesner, Jochen, and Andreas Vilcinskas. \"Antimicrobial Peptides: The "
                  "Ancient Arm of the Human Immune System.\" Virulence 1, no. 5 (2010): "
                  "440-64. https://doi.org/10.4161/viru.1.5.12983. PMID 21178486. Trial: " + CTG,
         source_type="review + registry record",
         mech="A short positively charged peptide based on one that cattle immune cells keep "
              "in granules to kill bacteria. It breaks microbial membranes on the skin. One "
              "gene, small, no cofactor; it is stored inside the cell that makes it and "
              "released on demand, and its target is a microbe rather than a shared animal "
              "protein.",
         notes=""),
    dict(set="B_i", id="B113", name="iseganan (IB-367)",
         species="pig (protegrins come from porcine leukocytes; the review read today says "
                 "porcine protegrins, no binomial)",
         molecule="iseganan, a synthetic protegrin analogue",
         target="microbial membranes in the human mouth (no human protein target)",
         year=2000, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00022373 (start 2000-12, phase 3, status unknown)",
         citation="Lehrer, Robert I., and Tomas Ganz. \"Cathelicidins: A Family of Endogenous "
                  "Antimicrobial Peptides.\" Current Opinion in Hematology 9, no. 1 (2002): "
                  "18-22. https://doi.org/10.1097/00062752-200201000-00004. PMID 11753073. "
                  "Trial: " + CTG,
         source_type="review + registry record",
         mech="A small disulfide-bonded peptide copied from one that pig white blood cells use "
              "to kill microbes. It is stored as a precursor inside the cell and cut to the "
              "active form. It punches holes in microbial membranes, so its target is not a "
              "protein humans share.",
         notes=""),
    dict(set="B_i", id="B114", name="HBOC-201 (hemoglobin glutamer-250, bovine)",
         species="cattle (the registry record says isolated bovine red blood cells; no "
                 "binomial)",
         molecule="cross-linked, glutaraldehyde-polymerised bovine haemoglobin",
         target="oxygen, carried in human plasma (no receptor)",
         year=2003, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00317512 (start 2003-12, phase 2, completed)",
         citation="ClinicalTrials.gov record NCT00317512: \"a crosslinked and "
                  "glutaraldehyde-polymerized hemoglobin (Hb) extracted from isolated bovine "
                  "red blood cells.\" " + CTG,
         source_type="registry record",
         mech="Cow haemoglobin taken out of red cells and chemically cross-linked, infused to "
              "carry oxygen in plasma. It is a large four-chain protein with a haem iron "
              "group, works inside the circulation rather than on a receptor, and needs its "
              "haem cofactor and its multi-chain assembly.",
         notes=""),
    dict(set="B_i", id="B115", name="exendin(9-39)",
         species="Heloderma suspectum", molecule="exendin(9-39), a truncated form of the "
         "39-residue lizard oral secretion peptide",
         target="GLP-1 receptor, which it blocks",
         year=2004, year_type="earliest trial start", agency="ClinicalTrials.gov",
         record="NCT00980083 (start 2004-10, phase 1, completed)",
         citation="FDA label for BYETTA (exenatide), for the species and the parent peptide: "
                  "\"originally identified in the lizard Heloderma suspectum.\" " + LABEL
                  + ". Trial: " + CTG,
         source_type="agency record + registry record",
         mech="A shortened version of a 39-residue peptide from the oral secretion of a "
              "venomous lizard, missing its first eight residues. It sits on the same human "
              "gut-hormone receptor as the full peptide but blocks it instead of switching it "
              "on. One gene, small, secreted, no partner needed.",
         notes=""),

    # ---------------- Set B(ii): named gene, no candidate in trials ----------------
    dict(set="B_ii", id="B201", name="naked mole rat high-molecular-mass hyaluronan route",
         species="Heterocephalus glaber", molecule="very high mass hyaluronan, made by a "
         "distinctive hyaluronan synthase 2 (HAS2) and kept by low hyaluronidase activity",
         target="proposed for human cancer, through hyaluronan signalling",
         year=2013, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Tian, Xiao, Jorge Azpurua, Christopher Hine, Amita Vaidya, Max "
                  "Myakishev-Rempel, Julia Ablaeva, Zhiyong Mao, Eviatar Nevo, Vera Gorbunova, "
                  "and Andrei Seluanov. \"High-Molecular-Mass Hyaluronan Mediates the Cancer "
                  "Resistance of the Naked Mole Rat.\" Nature 499, no. 7458 (2013): 346-49. "
                  "https://doi.org/10.1038/nature12234. PMID 23783513. Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A sugar chain five times longer than the human one, secreted by the skin cells "
              "of a long-lived burrowing rodent and left intact because the animal's "
              "chain-cutting enzymes are slow. Cells become more sensitive to it, and "
              "removing it lets them turn cancerous. It is not one protein: it needs a "
              "specific synthase, low cutting activity, and the cell's response to the chain.",
         notes=""),
    dict(set="B_ii", id="B202", name="tardigrade Dsup",
         species="Ramazzottius varieornatus (named in the paper's keywords; the abstract text "
                 "read today drops the species name)",
         molecule="damage suppressor protein (Dsup), a nuclear nucleosome-binding protein",
         target="human chromatin in cultured cells, which it shields from hydroxyl radicals",
         year=2019, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Chavez, Carolina, Grisel Cruz-Becerra, Jia Fei, George A. Kassavetis, and "
                  "James T. Kadonaga. \"The Tardigrade Damage Suppressor Protein Binds to "
                  "Nucleosomes and Protects DNA from Hydroxyl Radicals.\" eLife 8 (2019): "
                  "e47682. https://doi.org/10.7554/eLife.47682. PMID 31571581. Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A nuclear protein of a tardigrade that sits on nucleosomes and protects the DNA "
              "wrapped on them from hydroxyl radicals. One gene, one chain, and it binds a "
              "structure conserved across animals, but it works inside the nucleus and must "
              "be expressed in the cell to do anything.",
         notes=""),
    dict(set="B_ii", id="B203", name="tardigrade CAHS D for dry storage of biologics",
         species="a tardigrade (the paper read today says the protein is derived from a "
                 "tardigrade and does not give the binomial in its abstract)",
         molecule="CAHS D, a cytoplasmic abundant heat-soluble protein that forms a hydrogel",
         target="human blood clotting factor VIII, kept working after drying and heating",
         year=2023, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Packebush, Maxwell H., Silvia Sanchez-Martinez, Sourav Biswas, Shraddha Kc, "
                  "Kenny H. Nguyen, John F. Ramirez, Vincent Nicholson, and Thomas C. Boothby. "
                  "\"Natural and Engineered Mediators of Desiccation Tolerance Stabilize Human "
                  "Blood Clotting Factor VIII in a Dry State.\" Scientific Reports 13, no. 1 "
                  "(2023): 4542. https://doi.org/10.1038/s41598-023-31586-9. PMID 36941331. "
                  "Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A disordered protein from the cytoplasm of a tardigrade that forms a gel as it "
              "dries and holds other proteins in shape. In a vial it keeps a human clotting "
              "protein working through drying and heat. One gene, one chain, but it acts as a "
              "bulk protectant rather than by binding a site, and it works inside the "
              "organism's own cells.",
         notes=""),
    dict(set="B_ii", id="B204", name="Botryllus histocompatibility factor",
         species="Botryllus schlosseri", molecule="BHF, a protein encoded in the fuhc locus, "
         "found at the inner face of the plasma membrane and in the nucleus",
         target="proposed for human transplant rejection, through self versus non-self "
                "recognition",
         year=2015, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Taketa, Daryl A., Marie L. Nydam, Adam D. Langenbacher, Delany Rodriguez, "
                  "Erin Sanders, and Anthony W. De Tomaso. \"Molecular Evolution and in Vitro "
                  "Characterization of Botryllus Histocompatibility Factor.\" Immunogenetics "
                  "67, no. 10 (2015): 605-23. https://doi.org/10.1007/s00251-015-0870-1. "
                  "PMID 26359175. Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A protein of a colonial sea squirt that decides whether two colonies fuse or "
              "scar. It is found on the inside of the cell membrane and in the nucleus, and "
              "the authors suggest it works as a scaffold holding other proteins of the same "
              "animal together. One gene, but intracellular, not catalytic, and dependent on "
              "partners that only that animal has.",
         notes=""),
    dict(set="B_ii", id="B205", name="opossum venom-neutralising serum protein",
         species="opossums, family Didelphidae (the review read today names no binomial)",
         molecule="serum proteins homologous to human alpha-1B-glycoprotein that inhibit "
                  "venom metalloproteinases or phospholipase A2 myotoxins",
         target="proposed for snakebite envenoming, by neutralising venom enzymes",
         year=2012, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found for this molecule)",
         record="no study found for this molecule on 2026-10-19",
         citation="Voss, Robert S., and Sharon A. Jansa. \"Snake-Venom Resistance as a "
                  "Mammalian Trophic Adaptation: Lessons from Didelphid Marsupials.\" "
                  "Biological Reviews of the Cambridge Philosophical Society 87, no. 4 (2012): "
                  "822-37. https://doi.org/10.1111/j.1469-185X.2012.00222.x. PMID 22404916. "
                  "Trials: " + CTG,
         source_type="review + registry search",
         mech="A blood protein of a snake-eating marsupial, circulating in its own plasma, "
              "that binds and switches off venom enzymes. It is secreted, binds a target, and "
              "comes from a single gene family, but several different proteins act together in "
              "the animal and each one blocks a different venom component.",
         notes="Review, not the primary biochemical papers it points to."),
    dict(set="B_ii", id="B206", name="salamander nAG for scar-free regrowth",
         species="Ambystoma mexicanum and newts (the paper compares axolotl and newt)",
         molecule="newt anterior gradient protein (nAG), made by nerve sheath Schwann cells "
                  "and gland cells and released by holocrine secretion",
         target="proposed for human failure to regenerate, through blastema growth",
         year=2010, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Kumar, Anoop, Graham Nevill, Jeremy P. Brockes, and Andrew Forge. \"A "
                  "Comparative Study of Gland Cells Implicated in the Nerve Dependence of "
                  "Salamander Limb Regeneration.\" Journal of Anatomy 217, no. 1 (2010): 16-25. "
                  "https://doi.org/10.1111/j.1469-7580.2010.01239.x. PMID 20456522. "
                  "Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A secreted protein that regenerating nerves switch on in a salamander's limb "
              "stump, released from gland cells under the wound skin, and enough on its own to "
              "let a nerve-cut stump regrow digits. One gene and secreted, but it only works "
              "inside a limb stump that already has the right cells and nerve signals.",
         notes=""),
    dict(set="B_ii", id="B207", name="naked mole rat fructose-driven glycolysis",
         species="Heterocephalus glaber", molecule="the GLUT5 fructose transporter and "
         "ketohexokinase, used to run glycolysis on fructose",
         target="proposed for human oxygen-loss injury, by rerouting energy metabolism",
         year=2017, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found for this route)",
         record="no study found for this route on 2026-10-19",
         citation="Park, Thomas J., Jane Reznick, Bethany L. Peterson, Gregory Blass, Damir "
                  "Omerbasic, Nigel C. Bennett, P. Henning J. L. Kuich, et al. "
                  "\"Fructose-Driven Glycolysis Supports Anoxia Resistance in the Naked "
                  "Mole-Rat.\" Science 356, no. 6335 (2017): 307-11. "
                  "https://doi.org/10.1126/science.aab3896. PMID 28428423. Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A burrowing rodent survives 18 minutes with no oxygen by running its energy "
              "metabolism on fructose instead of glucose, using a sugar transporter and a "
              "sugar kinase expressed throughout its tissues. This is a whole metabolic route "
              "inside cells, needing several proteins and the right substrate, not one "
              "molecule that could be handed over.",
         notes=""),
    dict(set="B_ii", id="B208", name="hibernating bear bone maintenance",
         species="Ursus americanus", molecule="no bear-specific molecule is named; the paper "
         "measures leptin and sclerostin and proposes suppression of sympathetic signalling",
         target="proposed for human disuse bone loss",
         year=2011, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (not applicable; no molecule to search)",
         record="no molecule to search on 2026-10-19",
         citation="Seger, Rita L., Randal A. Cross, Clifford J. Rosen, Robert C. Causey, Caren "
                  "M. Gundberg, Thomas O. Carpenter, Tai C. Chen, et al. \"Investigating the "
                  "Mechanism for Maintaining Eucalcemia despite Immobility and Anuria in the "
                  "Hibernating American Black Bear (Ursus americanus).\" Bone 49, no. 6 (2011): "
                  "1205-12. https://doi.org/10.1016/j.bone.2011.08.017. PMID 21893223.",
         source_type="primary paper",
         mech="A bear lies still for months and loses no bone. Bone turnover stays balanced "
              "and suppressed, and the measurements point to hormone signalling that the "
              "animal shares with humans, with no animal-specific molecule identified. The "
              "effect is a whole-body state, not a molecule that can be moved.",
         notes="No animal-specific molecule is named; the row records that."),
    dict(set="B_ii", id="B209", name="fish antifreeze protein for cold organ storage",
         species="arctic fish (the paper says antifreeze proteins derived from arctic fish and "
                 "gives no binomial)",
         molecule="type III antifreeze protein, used at 15 mg/mL in the storage solution",
         target="ice crystals in the stored organ, not a human protein",
         year=2004, year_type="year of the paper naming the route",
         agency="ClinicalTrials.gov (checked, no study found)",
         record="no study found for this molecule on 2026-10-19",
         citation="Amir, Gabriel, Boris Rubinsky, Liana Horowitz, Liron Miller, Jonathan Leor, "
                  "Yigal Kassif, David Mishaly, Aram K. Smolinsky, and Jacob Lavee. "
                  "\"Prolonged 24-Hour Subzero Preservation of Heterotopically Transplanted Rat "
                  "Hearts Using Antifreeze Proteins Derived from Arctic Fish.\" Annals of "
                  "Thoracic Surgery 77, no. 5 (2004): 1648-55. "
                  "https://doi.org/10.1016/j.athoracsur.2003.04.004. PMID 15111159. "
                  "Trials: " + CTG,
         source_type="primary paper + registry search",
         mech="A small protein that a cold-water fish secretes into its own body fluids, where "
              "it binds the surface of ice crystals and stops them growing. Added to a storage "
              "solution it let rat hearts be kept below zero for 24 hours. One gene, one small "
              "domain, works alone, but what it binds is ice rather than a protein.",
         notes=""),
]

EXCLUDED = [
    ("lixisenatide", "A", "Animal origin not confirmed in any source opened today. The FDA "
     "record gives the approval (BLA208471, 2016-07-27) but neither the label nor a paper read "
     "today names the source species for the peptide it is based on.", FDA),
    ("exenatide extended release", "A", "Same active molecule as an entry already in Set A "
     "(exenatide), so it would double-count one molecule.", FDA),
    ("insulin glargine and lixisenatide combination", "A", "Fixed combination of a human "
     "protein with a peptide whose animal origin was not confirmed today.", FDA),
    ("enalapril", "A", "A later synthetic inhibitor in the same class as an entry already in "
     "Set A; the review read today names it alongside captopril as a drug created from the "
     "venom peptide, but it is a follow-on chemical rather than a molecule taken from the "
     "animal.", C_AVERIN),
    ("tirofiban", "A", "Non-peptide molecule designed from a venom binding motif; no source "
     "opened today names an animal molecule as its starting point.", FDA),
    ("aprotinin", "A", "Approval confirmed (BLA020304, 1993-12-29) but no source opened today "
     "names the source tissue or species. The withdrawn product has no label in the openFDA "
     "label API.", FDA),
    ("pegloticase", "A", "The label says the enzyme is recombinant modified mammalian urate "
     "oxidase whose cDNA is based on mammalian sequences, and names no species, so Set A rule 5 "
     "cannot be met from a record opened today.", LABEL),
    ("secretin", "A", "The approved product read today (ChiRhoStim) is synthetic human "
     "secretin, so it is excluded as a human protein. The porcine product SecreFlo (NDA021136, "
     "2002-04-04) is on the Drugs@FDA record but has no label in the label API, so its species "
     "statement could not be read today.", LABEL),
    ("lepirudin", "A", "No agency record found today: the openFDA drugsfda API returns nothing "
     "for the product, so the approval could not be confirmed on an agency record.", FDA),
    ("cytarabine", "A", "A synthetic nucleoside analogue whose design came from a sponge "
     "nucleoside; no source opened today names an animal molecule as the active substance.",
     "PubMed search, 2026-10-19"),
    ("plitidepsin", "A", "No FDA or EMA approval found today.", "PubMed search, 2026-10-19"),
    ("crotalidae polyvalent immune Fab", "A", "Excluded by Set A rule 4: the active molecules "
     "are antibodies raised in sheep against venom, not an evolved animal molecule with its own "
     "function.", FDA),
    ("antithrombin alfa made in goats", "A", "Excluded by Set A rule 3: the active molecule is "
     "a human protein, produced in an animal.", FDA),
    ("thyroid, desiccated (porcine)", "A", "Pig thyroid extract supplies the same hormones "
     "humans make; no species-specific molecule, and no original approval record was read "
     "today.", FDA),
    ("corticotropin (porcine pituitary)", "A", "Could not confirm today whether the molecule "
     "differs from the human hormone, so the human-protein exclusion could not be ruled out.",
     FDA),
    ("tetrodotoxin", "B_i", "Trials exist (NCT00725114, start 2008-04), but the paper read "
     "today shows pufferfish acquire the toxin from flatworm prey, so it is not a molecule the "
     "fish evolved. Okabe, Taiki, Hikaru Oyama, Maho Kashitani, Yuta Ishimaru, Rei Suo, Haruo "
     "Sugita, and Shiro Itoi. \"Toxic Flatworm Egg Plates Serve as a Possible Source of "
     "Tetrodotoxin for Pufferfish.\" Toxins 11, no. 7 (2019): 402. "
     "https://doi.org/10.3390/toxins11070402. PMID 31373322.", CTG),
    ("Trichuris suis ova", "B_i", "Trials exist (NCT01070498, start 2010-02) but the "
     "intervention is live whole eggs of a pig whipworm, not a named molecule.", CTG),
    ("honey bee venom", "B_i", "Trials exist (NCT00253942, start 2005-11) but the intervention "
     "is whole venom, not a named molecule.", CTG),
    ("cinobufacini", "B_i", "Trials exist (NCT01236690, start 2010-11) but the intervention is "
     "a toad skin secretion extract, not a named molecule.", CTG),
    ("lumbrokinase (DLBS1033)", "B_i", "Trials exist (NCT01790997, start 2012-05) but the "
     "intervention is a protein fraction of an earthworm, not a named molecule.", CTG),
    ("contulakin-G (CGX-1160)", "B_i", "No study found on ClinicalTrials.gov today, so Set B(i) "
     "rule 1 is not met.", CTG),
    ("Xen2174", "B_i", "No study found on ClinicalTrials.gov today under that name or "
     "XEN-2174.", CTG),
    ("Vc1.1 (ACV1)", "B_i", "No study found on ClinicalTrials.gov today.", CTG),
    ("cobratoxin", "B_i", "No study found on ClinicalTrials.gov today.", CTG),
    ("melittin", "B_i", "No study found on ClinicalTrials.gov today.", CTG),
    ("tebanicline (ABT-594)", "B_i", "No study found on ClinicalTrials.gov today.", CTG),
    ("Arctica islandica longevity", "B_ii", "No molecule is named for the trait, so Set B(ii) "
     "rule 1 is not met.", "twenty-diseases.md, 2026-10-19"),
    ("arctic clam deoxysphingolipid tolerance", "B_ii", "The claim in twenty-diseases.md is "
     "marked unverified there and no paper naming a clam molecule was found today.",
     "twenty-diseases.md, 2026-10-19"),
    ("bar-headed goose haemoglobin", "B_ii", "No source naming the molecule was returned by "
     "PubMed today under the searches run.", "PubMed search, 2026-10-19"),
    ("hookworm anti-inflammatory protein AIP-2", "B_ii", "No source naming the molecule was "
     "returned by PubMed today under the searches run.", "PubMed search, 2026-10-19"),
]


def main():
    fields = ["set", "entry_id", "drug_or_programme", "species", "molecule", "human_target",
              "year", "year_type", "agency_or_registry", "record_id", "citation",
              "source_type", "source_opened", "mechanism_description", "notes"]
    with (HERE / "sets.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(fields)
        for r in ROWS:
            w.writerow([r["set"], r["id"], r["name"], r["species"], r["molecule"], r["target"],
                        r["year"], r["year_type"], r["agency"], r["record"], r["citation"],
                        r["source_type"], OPENED, r["mech"], r["notes"]])

    with (HERE / "excluded.csv").open("w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["candidate", "intended_set", "reason", "source_checked", "source_opened"])
        for cand, which, reason, src in EXCLUDED:
            w.writerow([cand, which, reason, src, OPENED])

    order = list(range(len(ROWS)))
    random.Random(SHUFFLE_SEED).shuffle(order)
    with (HERE / "blind_input.csv").open("w", newline="") as f, \
            (HERE / "key.csv").open("w", newline="") as k:
        wb = csv.writer(f, lineterminator="\n")
        wk = csv.writer(k, lineterminator="\n")
        wb.writerow(["blind_id", "mechanism_description"])
        wk.writerow(["blind_id", "entry_id", "set", "drug_or_programme", "shuffle_seed"])
        for pos, idx in enumerate(order, start=1):
            r = ROWS[idx]
            blind_id = f"E{pos:03d}"
            wb.writerow([blind_id, r["mech"]])
            wk.writerow([blind_id, r["id"], r["set"], r["name"], SHUFFLE_SEED])

    counts = {}
    for r in ROWS:
        counts[r["set"]] = counts.get(r["set"], 0) + 1
    print("entries per set:", counts, "excluded:", len(EXCLUDED), "seed:", SHUFFLE_SEED)


if __name__ == "__main__":
    main()
