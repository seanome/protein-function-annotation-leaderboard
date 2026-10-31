"""InterProScan XML to InterProScan's own TSV layout, one row per match location.

    python3 interproscan_xml_to_tsv.py out.tsv batch0.xml batch1.xml ...

The whole-proteome InterProScan 5.59 run (nf-core/proteinannotator, all member databases, in
$SCRATCH/botryllus-annotation-full/results/proteinannotator/interproscan on Sherlock) kept only
XML. Columns, as InterProScan writes them with -f tsv: protein, sequence MD5, length, member
database, signature, signature description, start, end, score, status, date, InterPro entry,
InterPro description. The score is the location's E-value where the database gives one, else
the match's E-value or score, else "-". The protein is the first part of the XML xref id
("FUN000432_DZANK1|FUN000432_DZANK1"), the unique name UniFIRE also uses.
"""
import csv
import sys
import xml.etree.ElementTree as ET


# The XML names member databases in capitals; the TSV, which build.py reads, uses these names
LIBRARY = {"PFAM": "Pfam", "GENE3D": "Gene3D", "NCBIFAM": "NCBIfam", "FUNFAM": "FunFam", "ANTIFAM": "AntiFam",
           "PROSITE_PATTERNS": "ProSitePatterns", "PROSITE_PROFILES": "ProSiteProfiles", "COILS": "Coils",
           "MOBIDB_LITE": "MobiDBLite", "PHOBIUS": "Phobius", "TIGRFAM": "TIGRFAM", "SIGNALP_EUK": "SignalP_EUK",
           "SIGNALP_GRAM_POSITIVE": "SignalP_GRAM_POSITIVE", "SIGNALP_GRAM_NEGATIVE": "SignalP_GRAM_NEGATIVE"}


def local(tag):
    return tag.rsplit("}", 1)[-1]


def main(out, *paths):
    n = 0
    with open(out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        for p in paths:
            for _, el in ET.iterparse(p, events=("end",)):
                if local(el.tag) != "protein":
                    continue
                seq = xref = None
                matches = []
                for ch in el:
                    t = local(ch.tag)
                    if t == "sequence":
                        seq = ch
                    elif t == "xref" and xref is None:
                        xref = ch.get("id").split("|")[0]
                    elif t == "matches":
                        matches = list(ch)
                md5, length = seq.get("md5"), len((seq.text or "").strip())
                for m in matches:
                    sig = next((c for c in m if local(c.tag) == "signature"), None)
                    if sig is None:
                        continue
                    lib = next((c for c in sig if local(c.tag) == "signature-library-release"), None)
                    entry = next((c for c in sig if local(c.tag) == "entry"), None)
                    locs = next((c for c in m if local(c.tag) == "locations"), [])
                    for loc in locs:
                        score = loc.get("evalue") or m.get("evalue") or loc.get("score") or m.get("score") or "-"
                        w.writerow([xref, md5, length, LIBRARY.get(lib.get("library"), lib.get("library")) if lib is not None else "-", sig.get("ac"),
                                    sig.get("desc") or sig.get("name") or "-", loc.get("start"), loc.get("end"), score,
                                    "T", "-", entry.get("ac") if entry is not None else "-",
                                    (entry.get("desc") or entry.get("name") or "-") if entry is not None else "-"])
                        n += 1
                el.clear()
    print(n, "rows written to", out)


if __name__ == "__main__":
    main(*sys.argv[1:])
