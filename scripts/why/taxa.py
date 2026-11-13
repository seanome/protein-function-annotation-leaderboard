"""NCBI taxonomy from taxdump.tar.gz: parent, rank, scientific name, merged ids."""

import tarfile

# NCBI taxon ids used to group organisms.
BACTERIA, ARCHAEA, EUKARYOTA, VIRUSES = 2, 2157, 2759, 10239
METAZOA, FUNGI, VIRIDIPLANTAE, RHODOPHYTA, GLAUCOCYSTOPHYCEAE = 33208, 4751, 33090, 2763, 38254
SAR, HAPTISTA, CRYPTOPHYCEAE = 2698737, 2608109, 3027

# Catalogue of Life groups (Ruggiero et al. 2015 classification) as NCBI taxa. Protozoa is
# every eukaryote outside the other four, the way COL defines it by exclusion.
COL_GROUP_TAXA = {
    "Animalia": (METAZOA,),
    "Plantae": (VIRIDIPLANTAE, RHODOPHYTA, GLAUCOCYSTOPHYCEAE),
    "Fungi": (FUNGI,),
    "Chromista": (SAR, HAPTISTA, CRYPTOPHYCEAE),
    "Bacteria": (BACTERIA,),
    "Archaea": (ARCHAEA,),
}


class Taxonomy:
    def __init__(self, taxdump_path):
        self.parent, self.rank, self.name, self.merged = {}, {}, {}, {}
        with tarfile.open(taxdump_path) as tar:
            for line in tar.extractfile("nodes.dmp").read().decode().splitlines():
                f = line.split("\t|\t")
                self.parent[int(f[0])] = int(f[1])
                self.rank[int(f[0])] = f[2]
            for line in tar.extractfile("names.dmp").read().decode().splitlines():
                f = line.rstrip("\t|").split("\t|\t")
                if f[3] == "scientific name":
                    self.name[int(f[0])] = f[1]
            for line in tar.extractfile("merged.dmp").read().decode().splitlines():
                f = line.rstrip("\t|").split("\t|\t")
                self.merged[int(f[0])] = int(f[1])
        self._lineage = {}

    def resolve(self, taxid):
        return self.merged.get(taxid, taxid)

    def lineage(self, taxid):
        """Ancestors including the taxon itself, from the taxon up to the root."""
        taxid = self.resolve(taxid)
        if taxid in self._lineage:
            return self._lineage[taxid]
        out, t = [], taxid
        while t in self.parent and (not out or out[-1] != t):
            out.append(t)
            if self.parent[t] == t:
                break
            t = self.parent[t]
        self._lineage[taxid] = out
        return out

    def species(self, taxid):
        """The taxon of rank species at or above taxid; None when there is none."""
        for t in self.lineage(taxid):
            if self.rank.get(t) == "species":
                return t
        return None

    def superkingdom(self, taxid):
        lin = set(self.lineage(taxid))
        for t, label in ((BACTERIA, "Bacteria"), (ARCHAEA, "Archaea"), (EUKARYOTA, "Eukaryota"), (VIRUSES, "Viruses")):
            if t in lin:
                return label
        return "other"

    def col_group(self, taxid):
        lin = set(self.lineage(taxid))
        for group, taxa in COL_GROUP_TAXA.items():
            if lin & set(taxa):
                return group
        if EUKARYOTA in lin:
            return "Protozoa"
        if VIRUSES in lin:
            return "Viruses"
        return "other"

    def is_metazoan(self, taxid):
        return METAZOA in self.lineage(taxid)
