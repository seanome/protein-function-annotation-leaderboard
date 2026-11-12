"""Blind transferability scoring: T1..T6 from the description text alone."""

import csv
from pathlib import Path

BASE = Path(
    "/private/tmp/claude-501/-Users-olga-code-2024-kmerseek-analysis/"
    "9c1738be-f82c-4b99-b8ab-d5480f114584/scratchpad/blind"
)

# blind_id: ((T1, reason), (T2, reason), ... (T6, reason))
SCORES = {
    "E001": (
        (1, "one gene, stated"),
        (1, "secreted venom peptide"),
        (1, "binds the fibrinogen-binding site of a receptor"),
        (1, "an integrin, a receptor family not specific to one species"),
        (1, "small peptide; a seven-residue ring copies it"),
        (1, "no cofactor, stated"),
    ),
    "E002": (
        (1, "one gene, stated"),
        (1, "secreted into its own body fluids"),
        (1, "binds the surface of ice crystals"),
        (1, "what it binds is ice, which is not a species-specific partner"),
        (1, "one small domain"),
        (1, "works alone, stated"),
    ),
    "E003": (
        (1, "one gene, stated"),
        (1, "secreted into the host"),
        (1, "holds a starter complex inactive by binding"),
        (1, "a human clotting complex, not a species-specific partner"),
        (1, "small single domain"),
        (1, "no cofactor, stated"),
    ),
    "E004": (
        (0, "a designed hybrid; only part of it comes from the animal"),
        (1, "secreted, stated"),
        (1, "switches on two receptor enzymes by binding them"),
        (1, "human cell-surface receptor enzymes, not one species' partner"),
        (1, "a short peptide"),
        (1, "no cofactor, stated"),
    ),
    "E005": (
        (1, "one gene, one chain"),
        (0, "from the cytoplasm, and works inside the organism's own cells"),
        (0, "a bulk protectant, not binding a site"),
        (0, "no target named; it holds other proteins in shape generally"),
        (0, "size not stated"),
        (1, "in a vial it works with no partner added"),
    ),
    "E006": (
        (1, "two sugar-decorated forms of the same protein, so one gene"),
        (1, "secreted, stated"),
        (1, "cuts hyaluronan chains, catalytic"),
        (1, "hyaluronan between human cells, not a species-specific partner"),
        (0, "one catalytic chain, but size or domain count not stated"),
        (1, "no cofactor; the sugar groups are part of the molecule"),
    ),
    "E007": (
        (0, "several different proteins act together in the animal"),
        (1, "secreted, circulating in plasma"),
        (1, "binds and switches off venom enzymes"),
        (0, "its targets are snake venom components, not a conserved target"),
        (0, "size not stated"),
        (0, "several proteins act together, each on a different component"),
    ),
    "E008": (
        (0, "made by a multi-enzyme pathway, no single gene carries it"),
        (0, "acts on tubulin inside human cells; not stated secreted"),
        (1, "binds tubulin"),
        (1, "tubulin, one conserved target, stated"),
        (1, "a small molecule, smaller than a single domain"),
        (1, "the molecule itself acts alone"),
    ),
    "E009": (
        (1, "copies one salivary protein, so one gene"),
        (1, "secreted into the host's blood"),
        (1, "binds thrombin's active site and substrate groove"),
        (1, "thrombin, a human clotting enzyme, not one species' partner"),
        (1, "twenty residues, one chain"),
        (1, "no partner needed, stated"),
    ),
    "E010": (
        (0, "two different subunits that must stay together"),
        (1, "secreted, stated"),
        (1, "a membrane-cutting enzyme, catalytic"),
        (0, "no defined target; the effect is on cell membranes broadly"),
        (0, "two subunits, size not stated"),
        (0, "needs its partner subunit and calcium"),
    ),
    "E011": (
        (1, "one peptide from venom, one gene"),
        (1, "released outside the animal"),
        (1, "blocks a potassium channel by binding it"),
        (1, "a voltage-gated potassium channel, not one species' partner"),
        (1, "about 35 residues, a single folded domain"),
        (1, "no cofactor and no partner protein"),
    ),
    "E012": (
        (1, "copied from one peptide, so one gene"),
        (0, "stored as a precursor inside the cell"),
        (0, "punches holes in membranes, not catalysis or binding a target"),
        (0, "its target is not a protein humans share, stated"),
        (1, "a small disulfide-bonded peptide"),
        (0, "needs cutting from a precursor to become active"),
    ),
    "E013": (
        (0, "not one gene's product, stated"),
        (1, "acts at the air-liquid surface of the lung, outside cells"),
        (0, "physical surface tension lowering"),
        (0, "no molecular target named"),
        (0, "a lipid-and-protein extract, almost all lipid"),
        (0, "needs the lipids"),
    ),
    "E014": (
        (1, "one gene, stated"),
        (1, "secreted, released from gland cells"),
        (0, "no catalysis or binding named; it permits regrowth"),
        (0, "no target named"),
        (0, "size not stated"),
        (0, "needs a stump with the right cells and nerve signals"),
    ),
    "E015": (
        (0, "not the product of a single gene, stated"),
        (1, "acts on vessel lining cells from the outside"),
        (0, "no receptor named; the effect is on a cell layer"),
        (0, "no binding partner named"),
        (0, "a mixture of DNA fragments, 14-19 kDa, not a peptide"),
        (0, "not stated"),
    ),
    "E016": (
        (1, "one gene, stated"),
        (1, "a secreted salivary chain"),
        (1, "binds thrombin and holds it inactive"),
        (1, "thrombin, a human clotting enzyme, not one species' partner"),
        (1, "a single 65-residue chain"),
        (1, "no cofactor, no partner, stated"),
    ),
    "E017": (
        (1, "a single venom enzyme, one gene"),
        (1, "secreted, stated"),
        (1, "cuts fibrinogen, catalytic"),
        (1, "human fibrinogen, not a species-specific partner"),
        (0, "a single chain of moderate size, not a peptide or one domain"),
        (1, "needs no cofactor, stated"),
    ),
    "E018": (
        (1, "one gene, stated"),
        (1, "secreted in saliva"),
        (1, "blocks a calcium channel by binding it"),
        (1, "a human calcium channel, not one species' partner"),
        (1, "13 residues"),
        (1, "works alone, stated"),
    ),
    "E019": (
        (0, "a whole metabolic route needing several proteins"),
        (0, "inside cells, stated"),
        (0, "a metabolic route, not one catalytic or binding step"),
        (0, "needs the right substrate; no target"),
        (0, "several proteins, not one peptide"),
        (0, "needs several proteins and the right substrate"),
    ),
    "E020": (
        (1, "one gene, stated"),
        (1, "secreted in skin secretion"),
        (0, "kills bacteria by breaking membranes, not catalysis or binding"),
        (0, "acts on a microbial membrane, not a shared protein, stated"),
        (1, "about 22 residues"),
        (1, "no cofactor, stated"),
    ),
    "E021": (
        (0, "many gene products together, stated"),
        (1, "the enzymes are secreted ones"),
        (1, "each enzyme is catalytic"),
        (0, "acts on food in the gut, not a human protein"),
        (0, "a mixture of enzymes, size not stated"),
        (0, "many enzymes act together"),
    ),
    "E022": (
        (1, "one gene, stated"),
        (1, "secreted in oral secretion"),
        (1, "sits on a human receptor and blocks it"),
        (1, "a human gut-hormone receptor, not one species' partner"),
        (1, "a shortened 39-residue peptide"),
        (1, "no partner needed, stated"),
    ),
    "E023": (
        (1, "one gene, stated"),
        (1, "injected into prey, so released outside"),
        (1, "blocks the pore of a calcium channel by binding"),
        (1, "a human voltage-gated calcium channel, not one species' partner"),
        (1, "25 residues, one small folded chain"),
        (1, "works on its own, stated"),
    ),
    "E024": (
        (1, "a single catalytic protein, one gene"),
        (1, "secreted, acts outside cells"),
        (1, "cuts hyaluronan chains, catalytic"),
        (1, "hyaluronan between human cells, not one species' partner"),
        (0, "a single catalytic protein, but size or domain count not stated"),
        (1, "needs no cofactor, stated"),
    ),
    "E025": (
        (1, "a single venom enzyme, one gene"),
        (1, "secreted, stated"),
        (1, "cuts fibrinogen, catalytic"),
        (1, "human fibrinogen, not a species-specific partner"),
        (0, "a single chain, size not stated"),
        (1, "no cofactor, stated"),
    ),
    "E026": (
        (0, "a four-chain protein, so more than one gene"),
        (0, "normally inside red cells"),
        (0, "carries oxygen in plasma rather than acting on a receptor"),
        (0, "no target named"),
        (0, "a large four-chain protein"),
        (0, "needs its haem cofactor and multi-chain assembly"),
    ),
    "E027": (
        (0, "a whole-body state, no animal-specific molecule identified"),
        (0, "no molecule or location identified"),
        (0, "hormone signalling the animal shares, no molecule named"),
        (0, "no target named"),
        (0, "no molecule identified"),
        (0, "a whole-body state, not a molecule that can be moved"),
    ),
    "E028": (
        (0, "comes from more than one gene, stated"),
        (1, "acts at the air-liquid surface of the lung, outside cells"),
        (0, "a physical effect, lowering surface tension"),
        (0, "no molecular target named"),
        (0, "a lipid-and-protein extract"),
        (0, "needs the lipid partners, stated"),
    ),
    "E029": (
        (1, "one gene, stated"),
        (1, "secreted into its own blood"),
        (1, "binds the matching human receptor"),
        (1, "the matching human cell-surface receptor, not one species' partner"),
        (1, "a 32-residue chain"),
        (1, "no cofactor, stated"),
    ),
    "E030": (
        (1, "one gene, one chain"),
        (0, "a nuclear protein, works inside the nucleus"),
        (1, "sits on nucleosomes, binding"),
        (1, "binds a structure conserved across animals, stated"),
        (0, "one chain, but size or domain count not stated"),
        (0, "must be expressed in the cell to do anything"),
    ),
    "E031": (
        (1, "one gene, stated"),
        (1, "secreted in saliva"),
        (1, "catalytic, switches on plasminogen"),
        (1, "human plasminogen, not a species-specific partner"),
        (0, "size not stated"),
        (0, "needs the human partner protein and fibrin present"),
    ),
    "E032": (
        (1, "one gene makes the peptide, stated"),
        (1, "secreted into another animal's bloodstream"),
        (1, "blocks an enzyme at its active site by binding"),
        (1, "angiotensin-converting enzyme, not one species' partner"),
        (1, "a short peptide; a small molecule copies its contacts"),
        (1, "needs no partner protein and no cofactor"),
    ),
    "E033": (
        (1, "one gene, stated"),
        (0, "from sperm nuclei, where its job is packing DNA"),
        (1, "works by sticking to a sugar chain, binding"),
        (1, "a negatively charged sugar chain, not one species' partner"),
        (1, "a small peptide"),
        (1, "no cofactor, stated"),
    ),
    "E034": (
        (1, "one chain from venom, one gene"),
        (1, "secreted, stated"),
        (1, "cuts fibrin directly, catalytic"),
        (1, "the fibrin of a human clot, not one species' partner"),
        (0, "a shortened enzyme, size not stated"),
        (0, "needs a zinc ion and a blood inhibitor switches it off"),
    ),
    "E035": (
        (1, "one gene, stated"),
        (1, "released outside the animal's cells"),
        (1, "switches on a receptor by binding it"),
        (1, "the human GLP-1 receptor, not one species' partner"),
        (1, "39 residues, one chain"),
        (1, "no cofactor, no partner protein"),
    ),
    "E036": (
        (0, "not one gene's product, stated"),
        (0, "acts on DNA inside human cells, not secreted"),
        (1, "slips into the groove of DNA, binding"),
        (1, "DNA, the same in every phylum"),
        (1, "a small ring-shaped chemical"),
        (1, "needs no partner, stated"),
    ),
    "E037": (
        (0, "no single gene makes it, stated"),
        (0, "absorbed from the gut and acts on liver metabolism, not secreted"),
        (0, "changes a metabolic route, no receptor named"),
        (0, "no binding site named, stated"),
        (0, "a mixture of fatty acid esters"),
        (0, "the effect is a whole metabolic route"),
    ),
    "E038": (
        (0, "not one protein; needs a synthase and low cutting activity"),
        (1, "secreted by skin cells"),
        (0, "acts through the cell's response, not catalysis or binding"),
        (0, "no conserved target named; it needs the cell's own response"),
        (0, "a long sugar chain, not a peptide"),
        (0, "needs a synthase, low cutting activity and a cell response"),
    ),
    "E039": (
        (1, "one gene, stated"),
        (0, "stored inside the cell that makes it"),
        (0, "breaks microbial membranes, not catalysis or binding"),
        (0, "its target is a microbe, not a shared animal protein"),
        (1, "a short peptide"),
        (1, "no cofactor, stated"),
    ),
    "E040": (
        (0, "many components, several genes, stated"),
        (1, "acts at the air-liquid boundary of the lung, outside cells"),
        (0, "a physical effect of the mixture, not receptor binding"),
        (0, "no molecular target named"),
        (0, "a many-component lipid-and-protein extract"),
        (0, "added lipids are needed for it to work"),
    ),
    "E041": (
        (0, "not the product of one gene, stated"),
        (1, "works outside the cell"),
        (1, "binds antithrombin"),
        (1, "antithrombin, a human blood protein, not one species' partner"),
        (0, "a long sulfated sugar chain, not a peptide or domain"),
        (0, "works only through that partner protein"),
    ),
    "E042": (
        (1, "one peptide from venom, one gene"),
        (1, "secreted by the animal"),
        (1, "binds a receptor on the cell surface"),
        (1, "a human cell-surface receptor, not one species' partner"),
        (1, "36 residues, one small folded domain"),
        (1, "works without partners, stated"),
    ),
    "E043": (
        (1, "one gene, stated"),
        (0, "inside the cell membrane and in the nucleus"),
        (0, "a scaffold, not catalytic"),
        (0, "depends on partners that only that animal has"),
        (0, "size not stated"),
        (0, "dependent on partners that only that animal has"),
    ),
}


def main() -> None:
    rows_in = list(csv.DictReader((BASE / "blind_input.csv").open()))
    ids_in = [r["blind_id"] for r in rows_in]
    missing = [i for i in ids_in if i not in SCORES]
    extra = [i for i in SCORES if i not in ids_in]
    assert not missing, f"unscored rows: {missing}"
    assert not extra, f"scored rows not in input: {extra}"

    header = ["blind_id"]
    for n in range(1, 7):
        header += [f"T{n}", f"T{n}_reason"]
    header.append("T")

    out = BASE / "blind_scores.csv"
    dist = {k: 0 for k in range(7)}
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for bid in ids_in:
            items = SCORES[bid]
            assert len(items) == 6, bid
            row = [bid]
            total = 0
            for score, reason in items:
                assert score in (0, 1), (bid, score)
                assert reason.strip(), (bid, reason)
                row += [score, reason]
                total += score
            row.append(total)
            assert total == sum(s for s, _ in items)
            dist[total] += 1
            w.writerow(row)

    # read back and re-check the sum from the file itself
    check = list(csv.DictReader(out.open()))
    assert len(check) == len(ids_in), (len(check), len(ids_in))
    for r in check:
        assert int(r["T"]) == sum(int(r[f"T{n}"]) for n in range(1, 7)), r["blind_id"]

    print(f"wrote {out} with {len(check)} rows")
    for value in range(7):
        print(f"T={value}: {dist[value]}")


if __name__ == "__main__":
    main()
