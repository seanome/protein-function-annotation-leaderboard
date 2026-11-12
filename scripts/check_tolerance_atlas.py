"""Independent check of results/nature-solved-it/atlas/atlas.csv.

Written 2026-10-19 by the checking agent, which did not read the atlas README or
the build script's counting code.

Step 1 (fetch): pull every PMID and DOI out of columns 4, 5, 7 and 8 of atlas.csv,
fetch each PMID from NCBI E-utilities (title, DOI, abstract), fetch each UniProt
accession in column 6 from rest.uniprot.org, and cache the replies as JSON.
Step 2 (identifiers): mark an identifier FAILED when a PMID does not resolve, or
when the PMID's own DOI differs from the DOI written next to it.
Step 3 (judgments): merge the per-row reading judgments in check_judgments.csv
(written by hand after reading each abstract or full text) into check_rows.csv.
Step 4 (counts): recompute the five funnel counts and the negative-row count from
atlas.csv with the rules written in funnel() below, and only then compare with
counts.json.

Run: CHECK_CACHE=<dir> python scripts/check_tolerance_atlas.py fetch   (network, fills the cache)
     CHECK_CACHE=<dir> python scripts/check_tolerance_atlas.py report  (writes check_rows.csv)

CHECK_CACHE is where the fetched PubMed and UniProt replies are kept, outside the repo so
that the cache is never committed; it defaults to /tmp/atlas_check_cache. The report step
needs the cache only for the identifier checks; the counts come from atlas.csv alone.
"""

import csv
import os
import json
import re
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ATLAS = Path(__file__).resolve().parents[1] / "results" / "nature-solved-it" / "atlas"
CACHE = Path(os.environ.get("CHECK_CACHE", "/tmp/atlas_check_cache"))
COLS = {
    "c4": "citation_for_the_number",
    "c5": "independent_origins",
    "c6": "molecule_named",
    "c7": "moved_to_a_human_target",
    "c8": "approved_drug",
}
PMID_RE = re.compile(r"PMID[:\s]*(\d{6,9})")
DOI_RE = re.compile(r"(?:https?://doi\.org/)(10\.\d{4,9}/[^\s,]+)")
UNIPROT_RE = re.compile(
    r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b"
)


def rows():
    with open(ATLAS / "atlas.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def clean_doi(d):
    d = d.rstrip(".)")
    return d.lower()


def ids_in(text):
    return PMID_RE.findall(text), [clean_doi(d) for d in DOI_RE.findall(text)]


def efetch(pmids):
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id="
        + ",".join(pmids)
    )
    xml = urllib.request.urlopen(url, timeout=60).read()
    out = {}
    for art in ET.fromstring(xml).findall(".//PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID")
        title = "".join(art.find(".//ArticleTitle").itertext())
        abstract = " ".join("".join(a.itertext()) for a in art.findall(".//Abstract/AbstractText"))
        doi = None
        for aid in art.findall(".//PubmedData/ArticleIdList/ArticleId"):
            if aid.get("IdType") == "doi":
                doi = aid.text
        pmc = None
        for aid in art.findall(".//PubmedData/ArticleIdList/ArticleId"):
            if aid.get("IdType") == "pmc":
                pmc = aid.text
        out[pmid] = {"title": title, "abstract": abstract, "doi": doi, "pmc": pmc,
                     "year": art.findtext(".//PubDate/Year")}
    return out


def uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.json"
    try:
        j = json.loads(urllib.request.urlopen(url, timeout=60).read())
    except Exception as e:  # noqa: BLE001
        return {"error": str(e)}
    name = (j.get("proteinDescription", {}).get("recommendedName")
            or (j.get("proteinDescription", {}).get("submissionNames") or [{}])[0])
    return {
        "type": j.get("entryType"),
        "protein": name.get("fullName", {}).get("value") if name else None,
        "organism": j.get("organism", {}).get("scientificName"),
        "gene": ((j.get("genes") or [{}])[0].get("geneName") or {}).get("value"),
        "inactive": j.get("entryType") == "Inactive",
    }


def fetch():
    CACHE.mkdir(exist_ok=True)
    pmids, accs = set(), set()
    for r in rows():
        for key, col in COLS.items():
            p, _ = ids_in(r[col])
            pmids.update(p)
            if key == "c6":
                accs.update(UNIPROT_RE.findall(r[col]))
    pmids = sorted(pmids)
    pub = {}
    for i in range(0, len(pmids), 50):
        pub.update(efetch(pmids[i:i + 50]))
        time.sleep(0.5)
    (CACHE / "pubmed.json").write_text(json.dumps(pub, indent=1))
    uni = {}
    for a in sorted(accs):
        uni[a] = uniprot(a)
        time.sleep(0.2)
    (CACHE / "uniprot.json").write_text(json.dumps(uni, indent=1))
    print(len(pmids), "PMIDs asked,", len(pub), "returned;", len(uni), "UniProt accessions")


def identifier_checks():
    """One line per (row, column, identifier): does it resolve, and to the paper cited."""
    pub = json.loads((CACHE / "pubmed.json").read_text())
    uni = json.loads((CACHE / "uniprot.json").read_text())
    out = []
    for i, r in enumerate(rows(), 1):
        for key, col in COLS.items():
            text = r[col]
            pmids, dois = ids_in(text)
            # A PMID may be cited in one cell with its DOI written in another cell of
            # the same row, so compare against every DOI in the row.
            row_dois = set(dois)
            for c in COLS.values():
                row_dois.update(ids_in(r[c])[1])
            for p in pmids:
                rec = pub.get(p)
                if rec is None:
                    out.append((i, key, f"PMID {p}", "FAILED", "PMID did not resolve"))
                    continue
                if rec["doi"] and rec["doi"].lower() in row_dois:
                    out.append((i, key, f"PMID {p}", "HELD", f"DOI matches; {rec['title'][:90]}"))
                elif rec["doi"] is None:
                    out.append((i, key, f"PMID {p}", "HELD", f"PubMed has no DOI; {rec['title'][:90]}"))
                else:
                    out.append((i, key, f"PMID {p}", "HELD",
                                f"cited without a DOI; resolves to {rec['title'][:90]}"))
            pdois = {(v["doi"] or "").lower() for k, v in pub.items() if k in pmids}
            for d in dois:
                if d not in pdois:
                    out.append((i, key, f"DOI {d}", "CHECK", "DOI has no matching PMID in this cell"))
            if key == "c6":
                for a in UNIPROT_RE.findall(text):
                    u = uni.get(a, {})
                    if "error" in u or u.get("inactive"):
                        out.append((i, key, f"UniProt {a}", "FAILED", f"did not resolve: {u}"))
                    else:
                        out.append((i, key, f"UniProt {a}", "READ",
                                    f"{u['protein']} | {u['organism']} | gene {u['gene']} | {u['type']}"))
    return out


# ---------------------------------------------------------------- funnel counts
def says_yes(text):
    """Column 6/7/8 counts as yes only if the cell starts with 'yes'."""
    return text.strip().lower().startswith("yes")


def has_species_and_number(r):
    sp = r["species"].strip().lower()
    m = r["measured_ability"].strip().lower()
    if sp.startswith("none found"):
        return False
    return not (m.startswith("no number") or m.startswith("not found") or m.startswith("none found"))


def t_score(r):
    return int(re.match(r"T=(\d)", r["transferability_T"]).group(1))


def disease_id(r, scope):
    name = r["human_disease_and_insult_class"].rsplit(" [", 1)[0].strip()
    return scope[name]


def funnel(rs, scope, subset=None):
    """Diseases (not rows) clearing each step; a step needs all earlier steps in the same row."""
    if subset is not None:
        rs = [r for r in rs if subset(r)]
    step = [set(), set(), set(), set(), set()]
    for r in rs:
        d = disease_id(r, scope)
        step[0].add(d)
        if not has_species_and_number(r):
            continue
        step[1].add(d)
        if not says_yes(r["molecule_named"]):
            continue
        step[2].add(d)
        if not says_yes(r["moved_to_a_human_target"]):
            continue
        step[3].add(d)
        if says_yes(r["approved_drug"]):
            step[4].add(d)
    # Diseases with no row that has a species and a number.
    negatives = step[0] - step[1]
    return [len(s) for s in step], len(negatives), step


# ------------------------------------------------- what each paper actually says
# One entry per atlas row, filled in by hand after reading each paper: the abstract
# from the E-utilities reply in the cache, or the PubMed Central full text where the
# atlas says the number is in the body. "-" means the row cites nothing for that
# column, so there was nothing to open. Column keys: c3 the number and its
# condition, c4 the citation for the number, c5 independent origins, c6 the UniProt
# accessions, c7 the human target, c8 the approved drug.
JUDGMENTS = {
    1: {"c3": ("HELD", "abstract: cytotoxicity of myotoxin II inhibited by 92% at a 1:1 molar ratio"),
        "c4": ("HELD", "PMID 28759578 is the DM64 Pichia paper cited"),
        "c5": ("HELD", "abstract: 13 lineages plus 12 further motifs, 25 times in all"),
        "c6": ("HELD", "Q8MIS3 is Venom myotoxin inhibitor DM64, Didelphis marsupialis, as the row says"),
        "c7": ("HELD", "row claims nothing; the papers opened test DM64 on snake toxins"),
        "c8": ("HELD", "row claims nothing")},
    2: {"c3": ("HELD", "abstract gives no dose, time or fold-change, as the row says"),
        "c4": ("HELD", "PMID 24159039 is the grasshopper mouse Nav1.8 paper"),
        "c5": ("-", ""), "c6": ("-", "no accession claimed"), "c7": ("-", ""), "c8": ("-", "")},
    3: {"c3": ("HELD", "abstract: transfected HEK cells tolerated 50 mM ouabain; note the paper says "
                       "controls showed substantial mortality, the row says they died"),
        "c4": ("HELD", "PMID 24227116 is the Holzinger and Wink paper"),
        "c5": ("HELD", "abstract: N122H in four further species across four orders, 300 Myr of divergence"),
        "c6": ("HELD", "Q2XN70 is the Danaus plexippus Na+,K+-ATPase alpha subunit"),
        "c7": ("FAILED", "the pump put into human cells is the Drosophila ATPase carrying the monarch's "
                         "N122H, not the Danaus protein of column 6; no monarch molecule acts on a human target"),
        "c8": ("-", "")},
    4: {"c3": ("HELD", "neither abstract carries a dose, survival time or fold-change"),
        "c4": ("HELD", "both PMIDs resolve to the papers named"),
        "c5": ("HELD", "abstract: independent origins in Thamnophis couchii and T. sirtalis"),
        "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    5: {"c3": ("HELD", "abstract: 1,000-fold drop in affinity, whole-nerve assays"),
        "c4": ("HELD", "PMID 15815630 is the Nature clam paper"),
        "c5": ("HELD", "no origin count in the paper; RR = 70.0% confirmed in PMID 30075308"),
        "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    6: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
        "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    7: {"c3": ("HELD", "abstract: CYP6B46 moves 0.65% of ingested nicotine into the hemolymph"),
        "c4": ("HELD", "PMID 24379363 is the CYP6B46 paper"),
        "c5": ("-", ""), "c6": ("-", "no accession claimed"), "c7": ("-", ""), "c8": ("-", "")},
    8: {"c3": ("HELD", "abstract: seven isochromosomal lines on 4.5% v/v ethanol, four generations"),
        "c4": ("HELD", "PMID 28564174 is the Geer 1991 Evolution paper"),
        "c5": ("-", ""), "c6": ("HELD", "P00334 is Drosophila melanogaster alcohol dehydrogenase"),
        "c7": ("-", ""), "c8": ("-", "")},
    9: {"c3": ("HELD", "abstract: longer median lethal times for zinc only, not cadmium or copper"),
        "c4": ("HELD", "PMID 15712333 is the Nereis diversicolor paper"),
        "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    10: {"c3": ("HELD", "abstract: a review; the arsenic-tolerant earthworm lineage is named without a "
                        "binomial and without a number"),
         "c4": ("HELD", "PMID 24004878 is the ecotoxicological epigenetics review"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    11: {"c3": ("HELD", "both abstracts are concentration surveys with no tolerance number"),
         "c4": ("HELD", "PMIDs 35086557 and 31927288 resolve to the papers named"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    12: {"c3": ("HELD", "both Daphnia abstracts are about mine drainage, copper and cadmium, with no "
                        "lead dose"),
         "c4": ("HELD", "PMIDs 22278886 and 19379016 resolve to the papers named"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    13: {"c3": ("HELD", "neither abstract gives a cyanide dose the animal withstands"),
         "c4": ("HELD", "PMIDs 30018390 and 42572112 resolve; the second is Mol Biol Evol 43(8) 2026"),
         "c5": ("HELD", "abstract: CYP405A and CYP332A independently co-opted in Heliconiinae and Zygaena"),
         "c6": ("HELD", "A0A345BJE4 and A0A345BJE5 are Pieris rapae TST1 and TST2; PMID 28629163 is the "
                        "beta-cyanoalanine synthase paper"),
         "c7": ("-", ""), "c8": ("-", "")},
    14: {"c3": ("HELD", "full text: COX activity held across 0-0.9 uM H2S in the Puyacatengo and "
                        "Pichucalco populations; mitochondrial respiration measured across 0-90 uM"),
         "c4": ("HELD", "PMID 32586956 is the PNAS convergence paper"),
         "c5": ("HELD", "full text: three colonisations from 13,390,303 SNPs, 10 lineages across five "
                        "genera over 40 Myr, convergent substitutions at six COX1 and two COX3 codons"),
         "c6": ("HELD", "A0A023IMR0 is Poecilia mexicana cytochrome c oxidase subunit 1"),
         "c7": ("-", ""), "c8": ("-", "")},
    15: {"c3": ("HELD", "abstract: 8.7 +/- 0.3% mean, 10.4% maximum, 7.6 and 7.1% in juveniles and pups, "
                        "7% cut in total body oxygen stores"),
         "c4": ("HELD", "PMID 24829326 is the elephant seal carboxyhaemoglobin paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    16: {"c3": ("HELD", "abstract: 3.4- and 3.8-fold for MAmCq, 2.2- and 1.9-fold for HAmCq, G119S in ace-1"),
         "c4": ("HELD", "PMID 16363165 is the chlorpyrifos resistance paper"),
         "c5": ("HELD", "abstract: G119S arose independently at least four times in mosquitoes"),
         "c6": ("HELD", "Q867X2 is Culex quinquefasciatus acetylcholinesterase, gene ACE-1"),
         "c7": ("-", ""), "c8": ("-", "")},
    17: {"c3": ("HELD", "abstract: 10 mmol/l ambient ammonia built in 8 days, skin flux 0.0095 umol/cm2/min "
                        "across one pH unit, 30 mmol/l for 6 days given as the exposure, not as a stated limit"),
         "c4": ("HELD", "PMID 14747411 is the mudskipper NH4+ excretion paper"),
         "c5": ("HELD", "no origin count; the Boleophthalmus contrast is in PMID 12794674"),
         "c6": ("HELD", "A0A191B6Q3 is Periophthalmodon schlosseri vacuolar proton pump subunit B"),
         "c7": ("-", ""), "c8": ("-", "")},
    18: {"c3": ("HELD", "abstract: plasma urea fell from 340 to 225 mM at the lethal 685 ug/l silver dose; "
                        "the comparison with human kidney failure is the atlas's own, not the paper's"),
         "c4": ("HELD", "PMID 11489311 is the dogfish silver paper"),
         "c5": ("-", ""), "c6": ("HELD", "Q9I9F0 is a Squalus acanthias urea transporter; the three Rh "
                                         "accessions are Squalus ammonium transporters"),
         "c7": ("-", ""), "c8": ("-", "")},
    19: {"c3": ("HELD", "abstract: 17 mM fasted, up to 42 mM after feeding, glycated haemoglobin below "
                        "non-diabetic humans"),
         "c4": ("HELD", "PMID 9787825 is the hummingbird hyperglycaemia paper"),
         "c5": ("-", ""), "c6": ("HELD", "A0A7H1CRJ6 is Archilochus colubris GLUT5; the three Calypte anna "
                                         "accessions are GLUT6, GLUT10 and GLUT12, as the row says"),
         "c7": ("-", ""), "c8": ("-", "")},
    20: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    21: {"c3": ("HELD", "abstract: the paper measures resurrected enzyme activity, not what an animal "
                        "withstands, as the row says"),
         "c4": ("-", "no citation claimed for a number"),
         "c5": ("HELD", "abstract: multiple independent events, no count given"),
         "c6": ("HELD", "Q00511 is Aspergillus flavus uricase, alternative name Rasburicase; P16164 is "
                        "Sus scrofa uricase"),
         "c7": ("HELD", "abstract: ancient and modern uricases stably transfected into human HepG2 cells"),
         "c8": ("HELD", "FDA page fetched today: ELITEK, rasburicase, BLA 103946, Sanofi Synthelabo, "
                        "original approval 07/12/2002")},
    22: {"c3": ("HELD", "abstract: spisulosine IC50 1 uM in PC-3 and LNCaP, and no tolerance measurement"),
         "c4": ("HELD", "PMID 18343365 is the spisulosine paper"),
         "c5": ("-", ""), "c6": ("-", ""),
         "c7": ("HELD", "the row says the human-cell result runs the other way, which the abstract supports"),
         "c8": ("-", "")},
    23: {"c3": ("HELD", "abstract: 0.068 +/- 0.0032% absorbed, 6.1 mg Fe/day eaten, about 4.2 ug absorbed, "
                        "80 +/- 67 mg Fe/kg, turnover half-life 379 +/- 101 days"),
         "c4": ("HELD", "PMID 6128125 is the vampire bat iron balance paper"),
         "c5": ("HELD", "abstract: 27 bat species screened, REP15 loss on the vampire lineage, no origin count"),
         "c6": ("HELD", "P81172 is human hepcidin; the row states no Desmodus entry exists"),
         "c7": ("FAILED", "the cited paper measures bat hepcidin mRNA after iron challenge and compares it "
                          "with humans; no bat molecule acts in a human cell, on a human receptor or in a "
                          "human assay. The 'yes' rests on human hepcidin acting on human ferroportin"),
         "c8": ("-", "")},
    24: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    25: {"c3": ("HELD", "full text: 420 g per meal five times daily on days 18-24 in Landes geese, Sirius "
                        "Red with no fibrosis, and no liver weight or triglyceride figure, as the row says"),
         "c4": ("HELD", "PMID 41205434 is the goose PGC-1alpha paper"),
         "c5": ("-", ""), "c6": ("HELD", "A0A8B9ED99 and A0A8B9CCL4 are PGC-1alpha in Anser cygnoides and "
                                         "Anser brachyrhynchus, not Anser anser, as the row says"),
         "c7": ("HELD", "abstract of PMID 40651204: goose hepatocytes at 200 mmol/l glucose, not human cells"),
         "c8": ("-", "")},
    26: {"c3": ("HELD", "abstract: 9 days of anoxic submergence at 10 C against 2 days for Apalone "
                        "spinifera, with lactate in blood, bone and shell and rising plasma Ca and Mg"),
         "c4": ("HELD", "PMID 10893168 is the turtle shell buffering paper"),
         "c5": ("HELD", "abstract: a five-species ranking with no tree, as the row says"),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    27: {"c3": ("HELD", "abstract: 0-10% CO2 with ventilation up only at 10% and metabolism unchanged; "
                        "substance P 2 mg/kg and DiMe-C7 40-400 ug/kg restored the response at 2-10%"),
         "c4": ("HELD", "PMID 31967860 is the naked mole-rat hypercapnia paper"),
         "c5": ("HELD", "abstract of PMID 19219614: ammonia and capsaicin insensitivity, no origin count"),
         "c6": ("HELD", "A0A0P6JY17 is Heterocephalus glaber protachykinin-1; P25103 is the human substance-P "
                        "receptor"),
         "c7": ("FAILED", "the cited paper tests aprepitant, a synthetic antagonist, on human leukaemia "
                          "cells; naked mole-rat substance P was not put in a human cell or assay"),
         "c8": ("HELD", "FDA page fetched today: EMEND, aprepitant, NDA 021549, original approval "
                        "03/26/2003; the applicant is listed today as MSD, and the row claims no approval "
                        "for this animal molecule")},
    28: {"c3": ("HELD", "full text: 56% against 11% alive at 72 h, median survival 70 against 20 h, "
                        "hazard 1.9 against 7.2; the percentages are survival relative to the polyQ19 "
                        "control, not absolute shares, which the row does not say"),
         "c4": ("HELD", "PMID 32515471 is the aggresome paper"),
         "c5": ("-", ""), "c6": ("HELD", "G5BMN2 is Heterocephalus glaber autophagy protein 5"),
         "c7": ("-", ""), "c8": ("-", "")},
    29: {"c3": ("HELD", "full text: oxalate degradation exceeded 90% of the amount eaten at every dietary "
                        "concentration, up to the 3% oxalate diet"),
         "c4": ("HELD", "PMID 26896138 is the woodrat oxalate paper"),
         "c5": ("-", ""), "c6": ("HELD", "O06644 and P40149 are Oxalobacter formigenes frc and oxc"),
         "c7": ("HELD", "abstract of PMID 25645091: 80 patients, 700 million twice daily, hyperoxaluria "
                        "82.5% to 15%. The agent is the whole bacterium in people, and neither it nor its "
                        "enzymes are the woodrat's own molecule, which the row states"),
         "c8": ("-", "")},
    30: {"c3": ("FAILED", "the 1.5- to 2.0-fold syndecan-4 difference, the faster proliferation and the "
                          "lower cell-surface heparan sulfate all hold, but the added sentence that the two "
                          "breeds show no consistent difference in plasma cholesterol is in neither cited "
                          "paper's abstract"),
         "c4": ("HELD", "PMID 15596055 is the pigeon syndecan-4 paper"),
         "c5": ("HELD", "abstract: susceptibility behaves as one autosomal recessive trait"),
         "c6": ("HELD", "P31431 is human syndecan-4; the row states no Columba livia entry exists"),
         "c7": ("-", ""), "c8": ("-", "")},
    31: {"c3": ("HELD", "abstract: 18 minutes of anoxia survived without apparent injury"),
         "c4": ("HELD", "PMID 28428423 is the fructose glycolysis paper"),
         "c5": ("HELD", "abstract of PMID 30327830: independent invasions, different genes, no count"),
         "c6": ("HELD", "A0A0P6JEK5 and A0AAX6NXL4 are Heterocephalus glaber Slc2a5 and Khk; P22732 is "
                        "human SLC2A5"),
         "c7": ("-", ""), "c8": ("-", "")},
    32: {"c3": ("HELD", "abstract: 3 months of submergence anoxia at 3 C, shell lactate 133 mmol/kg, "
                        "plasma lactate 144 mmol/l"),
         "c4": ("HELD", "PMID 9320212 is the turtle shell lactate paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    33: {"c3": ("HELD", "full text: 45 min of deep hypothermic circulatory arrest at 18 C, troponin I "
                        "0.9 +/- 0.2 against 4.8 +/- 0.8 and 6.7 +/- 1.3 after 3 h; 6.7/0.9 is 7.4-fold, "
                        "and the row says about 8-fold"),
         "c4": ("HELD", "PMID 27187119 is the Anesthesiology proteomics paper"),
         "c5": ("HELD", "abstract of PMID 30327830: no origin count, as the row says"),
         "c6": ("HELD", "A0A8D2GQS9 is Urocitellus parryii SIRT3; Q9NTG7 is human SIRT3"),
         "c7": ("FAILED", "the cited paper raises human SIRT3 in human cardiomyocytes through HOTAIR and "
                          "FUS; the Arctic ground squirrel protein was not tested in a human cell"),
         "c8": ("-", "")},
    34: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    35: {"c3": ("HELD", "full text: no detectable calcium rise through the 10-minute hypoxia in P5-7 "
                        "slices, a rise only when extended to 15 minutes, and -3.8 +/- 0.75% against "
                        "-12.1 +/- 1.96% at P18-22, which is a 69% cut, against the row's 68%"),
         "c4": ("HELD", "PMID 22363676 is the calcium response paper"),
         "c5": ("HELD", "abstract of PMID 30327830: no origin count, as the row says"),
         "c6": ("HELD", "A0AAX6RX97 is Heterocephalus glaber Grin2d; O15399 is human GRIN2D; PMID 22155615 "
                        "is the GluN2D paper"),
         "c7": ("-", ""), "c8": ("-", "")},
    36: {"c3": ("HELD", "full text: 30 min of oxygen glucose deprivation at 36 C, bath PO2 0-2.9 mmHg "
                        "averaging 1.1 mmHg over 6 determinations, less lactate dehydrogenase release than "
                        "rat, and the same at pH 6.3 and in ionic-shift solution; the 2 h figure is in the "
                        "introduction, citing Christian et al. 2008, which the row says"),
         "c4": ("HELD", "PMID 28222226 is the Arctic ground squirrel OGD paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    37: {"c3": ("HELD", "abstract: 15 minutes at 7% oxygen at the normoxic maximum speed, heart rate "
                        "466 +/- 8 beats per minute, barnacle geese unable to finish"),
         "c4": ("HELD", "PMID 24710001 is the bar-headed goose running paper"),
         "c5": ("-", ""), "c6": ("HELD", "E5D5T8 is Anser indicus COX3; PMID 20685719 reports the Trp-116 "
                                         "to Arg substitution, the lower maximum activity and the higher "
                                         "affinity for reduced cytochrome c"),
         "c7": ("-", ""), "c8": ("-", "")},
    38: {"c3": ("HELD", "abstract: newborns at 3600 m had a lower pulmonary pressor response than those at "
                        "580 m, and the abstract carries no mmHg value, as the row says"),
         "c4": ("HELD", "PMID 29369354 is the llama nitric oxide paper"),
         "c5": ("-", ""), "c6": ("HELD", "O75116 is human ROCK2; the abstract reports lower ROCK2 and "
                                         "phosphorylated MLC20 and higher MYPT1 in highland lungs"),
         "c7": ("-", ""), "c8": ("-", "")},
    39: {"c3": ("HELD", "abstract: six apnoeas of 8-12 min, muscle blood flow to 31 +/- 8%, oxygen "
                        "depletion 1-2.3 ml O2/kg/min, high-energy phosphates unchanged"),
         "c4": ("HELD", "PMID 18840667 is the seal muscle apnoea paper"),
         "c5": ("HELD", "abstract of PMID 32324472: three lineages from separate recolonisations; it is a "
                        "review, which the row says"),
         "c6": ("HELD", "R9S002 is Mirounga angustirostris myoglobin; P02144 is human myoglobin"),
         "c7": ("-", ""), "c8": ("-", "")},
    40: {"c3": ("HELD", "abstract: 24 h of cold storage at 4 C then 24 h of rewarming to 37 C, with less "
                        "apoptosis than mouse cells"),
         "c4": ("HELD", "PMID 26457601 is the ground squirrel tubular cell paper"),
         "c5": ("-", ""), "c6": ("HELD", "I3M2E8 is Ictidomys tridecemlineatus XIAP; P98170 is human XIAP"),
         "c7": ("-", ""), "c8": ("-", "")},
    41: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    42: {"c3": ("HELD", "abstract: thermal hysteresis up to 5.5 C, over 4 C after 4 weeks at 4 C, with an "
                        "about 20-fold rise in thermal hysteresis protein"),
         "c4": ("HELD", "PMID 11029589 is the Tenebrio antifreeze protein paper"),
         "c5": ("HELD", "abstract of PMID 9108061: notothenioid and Arctic cod genes have different "
                        "progenitors, intron-exon layouts and spacers"),
         "c6": ("HELD", "O16119 is Tenebrio molitor thermal hysteresis protein YL-1; PMID 12372611 is the "
                        "ice-binding face paper"),
         "c7": ("HELD", "both papers cited are non-human: a transgenic mouse and porcine epidermis, which "
                        "is what the row's 'no' says"),
         "c8": ("-", "")},
    43: {"c3": ("HELD", "abstract: glucose-loaded frogs survived -5 C while saline controls died, against "
                        "a spring limit near -3 C, with lower ice content and red-cell injury"),
         "c4": ("HELD", "PMID 8409827 is the wood frog glucose paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    44: {"c3": ("HELD", "abstract: survival 75.0% to 89.0%, blastocysts per surviving oocyte 30.7% to "
                        "78.0%, double-strand-break-positive oocytes 65.2% to 30.8%"),
         "c4": ("HELD", "PMID 26202918 is the oocyte vitrification paper"),
         "c5": ("HELD", "abstract of PMID 9108061 supports more than one origin in fish"),
         "c6": ("-", "no accession claimed"),
         "c7": ("HELD", "the abstract does say these proteins should be tested on human oocytes before "
                        "clinical use, so the row's 'no' matches"),
         "c8": ("-", "")},
    45: {"c3": ("HELD", "abstract: steady-state torpor at -26 C through a 36-fold rise in torpid metabolic "
                        "rate over the minimum at 0 C, torpid down to -30 C when ramped, against -8 C and "
                        "13.5-fold in Callospermophilus lateralis"),
         "c4": ("HELD", "PMID 25590595 is the subzero thermogenesis paper"),
         "c5": ("-", ""), "c6": ("-", ""),
         "c7": ("HELD", "abstract of PMID 29576452: 13-lined ground squirrel cells, not Urocitellus "
                        "parryii, which is what the row says"),
         "c8": ("-", "")},
    46: {"c3": ("-", "no identifier to open"), "c4": ("-", ""), "c5": ("-", ""),
         "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    47: {"c3": ("HELD", "abstract: up to a 34-fold rise in immediate survival and colony formation after "
                        "drying to 0.2 g water per g dry weight, in human 293H cells"),
         "c4": ("HELD", "PMID 15963489 is the p26 and trehalose paper"),
         "c5": ("-", ""), "c6": ("HELD", "O44112 is the Artemia franciscana small heat shock protein p26"),
         "c7": ("HELD", "same paper: p26 expressed in human 293H cells; PMID 19228698 imports AfrLEA3m "
                        "into human hepatoma mitochondria"),
         "c8": ("-", "")},
    48: {"c3": ("HELD", "abstract: Dsup in human cultured cells suppressed X-ray DNA damage by about 40% "
                        "and improved radiotolerance"),
         "c4": ("HELD", "PMID 27649274 is the tardigrade genome paper"),
         "c5": ("HELD", "abstract of PMID 20817622: 5 kGy repaired without loss of viability, 25 kGy gives "
                        "frequent rearrangements, and no origin count, as the row says"),
         "c6": ("HELD", "P0DOW4 is the Ramazzottius varieornatus damage suppressor protein"),
         "c7": ("HELD", "same paper plus PMID 34681069 in human cells, and PMID 36858083 is the rat neuron "
                        "caution the row records"),
         "c8": ("-", "")},
    49: {"c3": ("HELD", "abstract: 112 mJ/cm2 of 222-nm light for a 1-log drop, the most resistant "
                        "vegetative cell tested"),
         "c4": ("HELD", "PMID 38670653 is the far-UVC paper"),
         "c5": ("-", ""), "c6": ("HELD", "Q9RXY7 is Deinococcus radiodurans R1 IrrE, tested against "
                                         "ionising radiation, not UV, as the row says"),
         "c7": ("HELD", "abstract of PMID 27974687: pprI in mice against cobalt-60 gamma rays, so not "
                        "human and not UV, which is what the row says"),
         "c8": ("-", "")},
    50: {"c3": ("HELD", "abstract: activation frequency of intracortical remodeling down 75%, mineral "
                        "apposition rate, geometry and whole-bone mechanics unchanged, hibernating femurs "
                        "less porous and more mineralised"),
         "c4": ("HELD", "PMID 18037367 is the grizzly bear bone paper"),
         "c5": ("-", ""), "c6": ("HELD", "abstract of PMID 21893223 reports the suppressed bone-specific "
                                         "alkaline phosphatase, the CTX range, the raised leptin and the "
                                         "sclerostin in a quarter of bears"),
         "c7": ("-", ""), "c8": ("-", "")},
    51: {"c3": ("HELD", "abstract: seven free-ranging female black bears, fibre number and cross-sectional "
                        "area unchanged, a moderate fall in protein concentration, and no percentage, as "
                        "the row says"),
         "c4": ("HELD", "PMID 9678502 is the black bear muscle paper"),
         "c5": ("-", ""), "c6": ("HELD", "P05019 is human IGF1, which the row marks as the human target "
                                         "rather than the bear agent"),
         "c7": ("HELD", "abstract of PMID 35077510: hibernating bear serum on human myotubes, raised total "
                        "protein, lower MuRF1, Akt/FOXO3a. The further claim that bear and human IGF-1 are "
                        "identical in sequence is not in the abstract and was not checked"),
         "c8": ("-", "")},
    52: {"c3": ("HELD", "abstract: HSP47 is the most substantially reduced platelet protein and the "
                        "abstract carries no number, as the row says"),
         "c4": ("-", "no citation claimed for a number"),
         "c5": ("-", ""), "c6": ("HELD", "P50454 is human Serpin H1"),
         "c7": ("HELD", "abstract: the same platelet signature in paralysed spinal cord injury patients"),
         "c8": ("-", "")},
    53: {"c3": ("HELD", "abstract: induced metamorphosis cut the regeneration rate twofold and caused "
                        "carpal and digit malformations, with no effect of body size"),
         "c4": ("HELD", "PMID 27499857 is the axolotl metamorphosis paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    54: {"c3": ("HELD", "abstract: hearts fully regenerated within 2 months of 20% ventricular resection, "
                        "by histology, and Mps1 mutants scarred"),
         "c4": ("HELD", "PMID 12481136 is the zebrafish heart paper"),
         "c5": ("-", ""), "c6": ("HELD", "Q8AYG3 is Danio rerio Ttk; P33981 is human TTK"),
         "c7": ("-", ""), "c8": ("-", "")},
    55: {"c3": ("HELD", "full text: ctgfa mutants showed about 71% less bridging than wild-type clutchmates "
                        "at 4 weeks"),
         "c4": ("HELD", "PMID 27811277 is the ctgfa paper"),
         "c5": ("-", ""), "c6": ("HELD", "Q5RI33 is Danio rerio ccn2a; P29279 is human CCN2"),
         "c7": ("HELD", "full text: human CTGF protein gave more than 10-fold more bridging at 1 week, in "
                        "the fish, which is the 'partial' the row records"),
         "c8": ("-", "")},
    56: {"c3": ("HELD", "full text: five of six 4 mm wounds fully re-epithelialized by day 3, and hair "
                        "follicles regrew in 1.5 cm wounds. The full text does name Acomys kempi and "
                        "Acomys percivali, so the two binomials were available"),
         "c4": ("HELD", "PMID 23018966 is the spiny mouse Nature paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    57: {"c3": ("HELD", "PMID 36629985 resolves to the IL11 essay in J Cardiovasc Transl Res 16(4) 2023, "
                        "an opinion piece with no number, as the row says"),
         "c4": ("-", "no citation claimed"), "c5": ("-", ""), "c6": ("-", ""),
         "c7": ("-", ""), "c8": ("-", "")},
    58: {"c3": ("HELD", "abstract: 10-day-old chicks, 1500 Hz at 120 dB SPL for 48 h, a patch with all hair "
                        "cells lost, new hair cells 48 h after exposure, bundles maturing over 10 days"),
         "c4": ("HELD", "PMID 3680064 is the chick cochlea paper"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    59: {"c3": ("HELD", "abstract: 7 days of constant 8000 lux, rods and cones dead by apoptosis within "
                        "24 h, all four cell types replaced within the 28-day recovery, mosaic not restored, "
                        "proliferation next to Mueller glia"),
         "c4": ("HELD", "PMID 10942883 resolves and its PubMed DOI matches the row's, including the ';2-h' "
                        "ending"),
         "c5": ("-", ""), "c6": ("-", ""), "c7": ("-", ""), "c8": ("-", "")},
    60: {"c3": ("HELD", "abstract: tooth families with up to five generations at hatching and 42 teeth in "
                        "function after hatching"),
         "c4": ("HELD", "PMID 2353674 is the alligator dentition paper"),
         "c5": ("-", ""), "c6": ("HELD", "PMID 23671090 reports the beta-catenin activation and sFRP1 loss; "
                                         "it does have a full text in PubMed Central, PMC3670376, although "
                                         "the row says the tool returned none"),
         "c7": ("-", ""), "c8": ("-", "")},
    61: {"c3": ("HELD", "full text: 15 U/kg in three animals, sub-lethal or lethal for C57BL/6 at 5 and "
                        "15 U/kg, collagen area about doubled in Mus by day 28 and did not rise in Acomys"),
         "c4": ("HELD", "PMID 40883787 is the Acomys bleomycin paper"),
         "c5": ("-", ""), "c6": ("HELD", "full text: about 2 times more bleomycin hydrolase by "
                                         "semi-quantitative analysis, with the authors' own BALB/c caveat; "
                                         "Q13867 is human BLMH"),
         "c7": ("-", ""), "c8": ("-", "")},
    62: {"c3": ("HELD", "abstract: cancer mortality 4.81% (95% CI 3.14-6.49) against 11-25% in humans, "
                        "from necropsies of 644 elephants"),
         "c4": ("HELD", "PMID 26447779 is the JAMA elephant paper"),
         "c5": ("HELD", "abstract of PMID 27642012: 20 TP53 copies, gained as body size grew, in the "
                        "Proboscidean lineage"),
         "c6": ("HELD", "G3UI57 is a Loxodonta africana p53 entry; P04637 is human TP53"),
         "c7": ("HELD", "abstract: 14.64% of elephant cells against 7.17% of human control and 2.71% of "
                        "Li-Fraumeni cells, which is the 'partial' the row records"),
         "c8": ("-", "")},
}
COLUMN_TITLES = {
    "c3": "the number, its unit and its condition",
    "c4": "citation for the number",
    "c5": "independent origins",
    "c6": "UniProt accessions",
    "c7": "human target",
    "c8": "approved drug",
}


def write_check_rows(rs, scope):
    """One line per row and column: held, failed, or nothing to open."""
    out = []
    for i, r in enumerate(rs, 1):
        j = JUDGMENTS[i]
        for col in ["c3", "c4", "c5", "c6", "c7", "c8"]:
            status, note = j[col]
            out.append([
                i, disease_id(r, scope), r["species"], col, COLUMN_TITLES[col],
                "NOTHING TO OPEN" if status == "-" else status, note,
            ])
    with open(ATLAS / "check_rows.csv", "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["row", "disease_id", "species", "column", "what_was_checked",
                    "status", "what_the_paper_says"])
        w.writerows(out)
    return out


def main():
    if sys.argv[1:] == ["fetch"]:
        fetch()
        return
    rs = rows()
    with open(ATLAS / "scope.csv", newline="") as fh:
        scope = {s["human_disease"]: s["disease_id"] for s in csv.DictReader(fh)}
    counts, neg, steps = funnel(rs, scope)
    # Step 0 for the whole atlas is every in-scope disease, rows or not.
    counts[0] = len(scope)
    neg_all = len(scope) - counts[1]
    t56, neg56, _ = funnel(rs, scope, subset=lambda r: t_score(r) >= 5)
    print("all rows:", counts, "negatives:", neg_all)
    print("T 5-6   :", t56, "diseases with a T5-6 row but no number:", neg56)
    for k, s in zip(["number", "molecule", "human", "drug"], steps[1:]):
        print(k, sorted(s))
    print("rows per disease >1:",
          [d for d in set(disease_id(r, scope) for r in rs)
           if sum(disease_id(r, scope) == d for r in rs) > 1])
    ids = identifier_checks()
    with open(CACHE / "identifier_checks.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["row", "column", "identifier", "status", "note"])
        w.writerows(ids)
    print(sum(1 for x in ids if x[3] == "FAILED"), "identifier failures")
    checks = write_check_rows(rs, scope)
    tally = {}
    for c in checks:
        tally[c[5]] = tally.get(c[5], 0) + 1
    print("checks:", tally)
    print("failed:", [(c[0], c[3]) for c in checks if c[5] == "FAILED"])
    # What the funnel becomes if every row with a failed human-target check is dropped
    # from step 4.
    failed_human = {c[1] for c in checks if c[5] == "FAILED" and c[3] == "c7"}
    print("diseases losing step 4:", sorted(failed_human))
    print("step 4 after dropping them:", sorted(steps[3] - failed_human))


if __name__ == "__main__":
    main()
