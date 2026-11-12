"""Small NCBI taxonomy helper built from the parquet caches (own code)."""
from pathlib import Path

import polars as pl

C = Path(__file__).resolve().parent / "cache"


class Tax:
    def __init__(self):
        n = pl.read_parquet(C / "nodes.parquet")
        self.parent = dict(zip(n["taxid"].to_list(), n["parent"].to_list()))
        self.rank = dict(zip(n["taxid"].to_list(), n["rank"].to_list()))
        nm = pl.read_parquet(C / "names.parquet")
        self.name = dict(zip(nm["taxid"].to_list(), nm["name"].to_list()))
        m = pl.read_parquet(C / "merged.parquet")
        self.merged = dict(zip(m["old"].to_list(), m["new"].to_list()))
        self.deleted = set(pl.read_parquet(C / "deleted.parquet")["taxid"].to_list())
        self._lin = {}

    def resolve(self, t):
        if t in self.parent:
            return t
        return self.merged.get(t)

    def lineage(self, t):
        t = self.resolve(t)
        if t is None:
            return None
        if t in self._lin:
            return self._lin[t]
        out = []
        cur = t
        while True:
            out.append(cur)
            p = self.parent[cur]
            if p == cur:
                break
            cur = p
        self._lin[t] = out
        return out

    def species_of(self, t):
        lin = self.lineage(t)
        if lin is None:
            return None
        for x in lin:
            if self.rank[x] == "species":
                return x
        return None

    def has(self, t, anc):
        lin = self.lineage(t)
        return lin is not None and anc in lin

    def col_group(self, t):
        lin = self.lineage(t)
        if lin is None:
            return None
        s = set(lin)
        if 33208 in s:
            return "Animalia"
        if s & {33090, 2763, 38254}:
            return "Plantae"
        if 4751 in s:
            return "Fungi"
        if s & {2698737, 2608109, 3027}:
            return "Chromista"
        if 2 in s:
            return "Bacteria"
        if 2157 in s:
            return "Archaea"
        if 2759 in s:
            return "Protozoa"
        return None

    def cellular(self, t):
        lin = self.lineage(t)
        return lin is not None and bool(set(lin) & {2, 2157, 2759})
