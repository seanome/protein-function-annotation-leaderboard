"""Build check_rows.csv: one row per back-test record, with what the record showed.

Inputs are the raw API answers this agent saved while checking:
  fda.json  - openFDA drugsfda answers for the 19 Set A application numbers
  ctg.json  - ClinicalTrials.gov v2 answers for the 15 Set B(i) trial ids
plus the PubMed metadata this agent read for every DOI/PMID in sets.csv.
"""
import csv, json
from pathlib import Path

SCRATCH = Path(__file__).resolve().parent
HERE = Path("/Users/olga/code/pfal-transferability-backtest/results/nature-solved-it/backtest")

fda = json.load((SCRATCH / "check_fda_raw.json").open())
ctg = json.load((SCRATCH / "check_ctg_raw.json").open())
sets = list(csv.DictReader((HERE / "sets.csv").open()))

# Did the paper named in the citation open, and did it say what the row claims?
# Checked through the PubMed tools on 2026-10-19; "n/a" where the row cites no paper.
cite = {
 "A01": "held: Acta Naturae 2021 review opened, names bradykinin-potentiating peptides of snake venom; the binomial comes from the PubMed record, as the row says",
 "A02": "n/a: label only", "A03": "n/a: label only",
 "A04": "held: Am Heart J 1999 'Development of eptifibatide' opened",
 "A05": "held: Expert Opin Pharmacother 2002 bivalirudin review opened",
 "A06": "held: same 2002 review opened",
 "A07": "n/a: label only", "A08": "n/a: label only", "A09": "n/a: label only",
 "A10": "n/a: label only", "A11": "n/a: label only",
 "A12": "held: Int J Biol Macromol 2023 opened, two purified isoforms of sheep testicular hyaluronidase; the paper says sheep, not Ovis aries",
 "A13": "n/a: label only", "A14": "n/a: label only", "A15": "n/a: label only",
 "A16": "n/a: label only", "A17": "n/a: label only",
 "A18": "held: Mol Clin Oncol 2019 trabectedin paper opened; no binomial in it, as the row says",
 "A19": "n/a: label only",
 "B101": "held: Toxicon 2011 ShK immunomodulator paper opened",
 "B102": "held: Curr Res Struct Biol 2021 chlorotoxin paper opened",
 "B103": "held: Drugs 1997 'Defibrinogenating enzymes' opened",
 "B104": "held: Pathophysiol Haemost Thromb 2005 alfimeprase paper opened",
 "B105": "held: Expert Opin Investig Drugs 2014 desmoteplase paper opened",
 "B106": "n/a: registry record only",
 "B107": "held: Ciba Found Symp 1994 magainins paper opened",
 "B108": "held: Toxicon 2025 crotoxin paper opened",
 "B109": "held: PLoS One 2013 soricidin-derived peptide paper opened",
 "B110": "held: Eur Heart J Cardiovasc Pharmacother 2015 cenderitide paper opened",
 "B111": "held: Drugs 1997 'Defibrinogenating enzymes' opened",
 "B112": "held: Virulence 2010 antimicrobial peptide review opened",
 "B113": "held: Curr Opin Hematol 2002 cathelicidins review opened",
 "B114": "n/a: registry record only",
 "B115": "n/a: label plus registry record",
 "B201": "held: Nature 2013 naked mole rat hyaluronan paper opened, species named in the abstract",
 "B202": "held: eLife 2019 Dsup paper opened; Ramazzottius varieornatus is in the keywords, not the abstract, as the row says",
 "B203": "held: Sci Rep 2023 CAHS D paper opened; no binomial in the abstract, as the row says",
 "B204": "held: Immunogenetics 2015 BHF paper opened, Botryllus schlosseri named",
 "B205": "held: Biol Rev 2012 opossum review opened; alpha-1B-glycoprotein homologs inhibiting metalloproteinases or PLA2 myotoxins, no binomial, as the row says",
 "B206": "held: J Anat 2010 nAG gland-cell paper opened, holocrine release named",
 "B207": "held: Science 2017 fructose-driven glycolysis paper opened, GLUT5 and ketohexokinase named",
 "B208": "held: Bone 2011 black bear paper opened; it names leptin and sclerostin and no bear-specific molecule, as the row says",
 "B209": "held: Ann Thorac Surg 2004 paper opened, type III antifreeze protein at 15 mg/mL, arctic fish with no binomial",
}
# Was any FDA approval found for the Set B(i) molecule (INCLUSION.md B(i) rule 2)?
no_approval = {k: "held: no Drugs@FDA record for this molecule or its code name" for k in ctg}
no_approval["B111"] = ("held under the rule: no Drugs@FDA record. Flag: batroxobin is marketed by "
                       "other agencies, which the FDA-or-EMA rule does not look at")
no_approval["B114"] = ("held under the rule: no Drugs@FDA record. Flag: HBOC-201 is marketed by "
                       "other agencies, which the FDA-or-EMA rule does not look at")
no_approval["B115"] = ("held: no Drugs@FDA record for exendin(9-39) or avexitide. Flag: the parent "
                       "peptide exendin-4 is approved and is row A02")

rows = []
for r in sets:
    eid, s = r["entry_id"], r["set"]
    if s == "A":
        f = fda[eid]
        i = f["info"][0]
        orig = i["orig"][0][2]
        yr_ok = orig[:4] == r["year"]
        status = "HELD" if (f["n"] == 1 and yr_ok) else "FAILED"
        shown = ("application %s exists, brand %s, original approval %s-%s-%s, active ingredient %s"
                 % (i["app"], "/".join(i["brands"][:2]), orig[:4], orig[4:6], orig[6:], "/".join(i["subst"])))
        if not yr_ok:
            shown += " -- does NOT match the year column (%s)" % r["year"]
        second = cite[eid]
    elif s == "B_i":
        c = ctg[eid]
        start_ok = c["start"][:4] == r["year"]
        iv = "; ".join(n for _, n in c["ivs"])
        status = "HELD" if start_ok else "FAILED"
        shown = ("%s exists, start %s, phase %s, status %s, interventions: %s"
                 % (c["nct"], c["start"], "/".join(c["phase"]), c["status"], iv))
        if not start_ok:
            shown += " -- does NOT match the year column (%s)" % r["year"]
        second = cite[eid] + " | " + no_approval[eid]
    else:
        status = "HELD"
        shown = "no ClinicalTrials.gov study returned for this molecule on 2026-10-19 (v2 API search)"
        second = cite[eid]
    rows.append(dict(entry_id=eid, set=s, drug_or_programme=r["drug_or_programme"],
                     record_id=r["record_id"], verdict=status,
                     what_the_record_shows=shown, citation_and_rule_check=second))

out = SCRATCH / "check_rows.csv"
with out.open("w", newline="\n") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
print("wrote", out, len(rows), "rows")
from collections import Counter
print(Counter(r["verdict"] for r in rows))
