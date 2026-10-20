/* Organism pages. One file for the browser and for the build:
     computeAnswer(lookups, strata)  -> the answer for one species (cards, verdict, figures)
     renderAnswer(answer, strata)    -> that answer as HTML
   scripts/organism/build.py runs both in Node for every pre-built species, so a pre-built
   page and a page built here in the browser come from the same code.
   lookupSpecies(taxid, tables) does the four card lookups in the browser, for a species with
   no pre-built page: UniProt (taxonomy, proteomes, counts), AlphaFold DB, GoaT, and the
   tables this site ships (reviewed Swiss-Prot counts, TimeTree ages). */
(function (root) {
  "use strict";
  const BANDS = ["0-50", "50-70", "70-90", "90-100"];
  const PA_VERSION = "1.1.0";
  const MIN_MODELS = 100;

  // ------------------------------------------------------------------ small helpers
  const esc = (s) => String(s === null || s === undefined ? "" : s)
    .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  const pct = (x, d = 0) => (x === null || x === undefined || Number.isNaN(x) ? "—" : (100 * x).toFixed(d) + "%");
  const num = (x, d = 0) => (x === null || x === undefined || Number.isNaN(x) ? "—"
    : Number(x).toLocaleString("en-US", { minimumFractionDigits: d, maximumFractionDigits: d }));
  const bandLo = (b) => Number(b.split("-")[0]);
  const bandHi = (b) => Number(b.split("-")[1]);
  function bandOf(v) {
    if (v === null || v === undefined) return null;
    for (const b of BANDS) if ((v >= bandLo(b) && v < bandHi(b)) || (bandHi(b) === 100 && v === 100)) return b;
    return null;
  }
  const displayName = (sci, common) => (common ? `${common} (${sci})` : sci);

  // AlphaFold DB model counts in 5-point pLDDT bands -> mean pLDDT (estimate: each model is
  // counted at its band's middle, so the true mean is within 2.5 of this) and the share of
  // models in each of the four bands the results tables use.
  function plddtSummary(bands5) {
    if (!bands5) return null;
    let total = 0, sum = 0;
    const four = Object.fromEntries(BANDS.map((b) => [b, 0]));
    for (const [k, n] of Object.entries(bands5)) {
      const m = k.match(/plddt_(\d+)_(\d+)/);
      const lo = Number(m[1]), hi = Number(m[2]);
      total += n; sum += n * (lo + hi) / 2;
      four[bandOf(lo)] += n;
    }
    if (!total) return { models: 0, mean_plddt_est: null, share: null };
    return { models: total, mean_plddt_est: sum / total,
             share: Object.fromEntries(BANDS.map((b) => [b, four[b] / total])), counts: four };
  }

  // The nearest species with >= min reviewed entries, other than the species itself: walk up
  // the lineage (nearest first); the first ancestor with such a species under it wins, and
  // under it the species with most entries (ties: lower taxon id). Same rule as
  // scripts/organism/lookups.py.
  function nearestReviewed(lineage, selfTaxid, reviewed) {
    for (const a of lineage.slice(1)) {
      for (const [t, e] of reviewed.under[String(a.taxid)] || []) {
        if (t !== selfTaxid) return { taxid: t, entries: e, shared_ancestor: a };
      }
    }
    return null;
  }

  // Divergence time from human of a measured relative, as the page shows it (whole Mya).
  const relMya = (S, id) => { const r = S.relatives.find((x) => x.id === id); return r ? num(r.mya) : "—"; };
  // Two measured relatives with the same age, for the Figure 3 caption: "cow and dog, 89 Mya".
  function sharedAge(S) {
    const by = {};
    for (const r of S.relatives) if (r.mya !== null) (by[r.mya] = by[r.mya] || []).push(r.name);
    const m = Object.keys(by).map(Number).sort((a, b) => a - b).find((k) => by[k].length > 1);
    return m === undefined ? null : `${by[m].slice(0, 2).join(" and ")}, ${num(m)} Mya`;
  }

  function measuredRelative(mya, S) {
    if (mya === null || mya === undefined) return null;
    let best = null;
    for (const r of S.relatives) {
      if (r.mya === null) continue;
      const d = Math.abs(r.mya - mya);
      if (!best || d < best.d) best = { r, d };
    }
    return best && best.r;
  }

  function sharedDepth(lineage, otherLineage) {
    const other = new Set(otherLineage);
    for (let i = 0; i < lineage.length; i++) if (other.has(lineage[i].taxid)) return { i, taxid: lineage[i].taxid, name: lineage[i].name };
    return null;
  }

  // ------------------------------------------------------------------ the answer
  function computeAnswer(L, S) {
    const A = { taxid: L.taxid, name: L.name, common_name: L.common_name, slug: L.slug, built: L.built || null };
    const af = L.alphafold;
    const ps = af ? plddtSummary(af.bands5) : null;
    A.cards = {
      relative: L.relative ? {
        taxid: L.relative.taxid, name: L.relative.name, common_name: L.relative.common_name,
        entries: L.relative.entries, shared_ancestor: L.relative.shared_ancestor,
        mya: L.mya ? L.mya.value : null, mya_method: L.mya ? L.mya.method : null, mya_detail: L.mya || null,
      } : null,
      reviewed: { self: L.self_entries, order: L.ranks.order, class: L.ranks.class },
      alphafold: af ? {
        proteome: L.proteome ? L.proteome.upid : null, proteome_entries: af.proteome_entries,
        with_model: af.with_model, share_with_model: af.share_with_model,
        afdb_models_counted: ps ? ps.models : null, mean_plddt_est: ps ? ps.mean_plddt_est : null,
        band_share: ps ? ps.share : null, band_counts: ps ? ps.counts : null,
      } : null,
      ebp: L.ebp,
    };
    A.proteome = L.proteome;

    // Stratum: distance, then pLDDT band.
    const mya = A.cards.relative ? A.cards.relative.mya : null;
    const rel = measuredRelative(mya, S);
    // A mean over fewer than MIN_MODELS AlphaFold models does not place a proteome
    // (Priapulus caudatus: 1 model counted, of 18,404 proteins).
    const enough = ps && ps.models >= MIN_MODELS;
    const meanP = enough ? ps.mean_plddt_est : null;
    const asked = bandOf(meanP);
    const st = { relative: rel, band_asked: asked, mean_plddt_est: meanP, models_counted: ps ? ps.models : 0,
                 fallback: false, reason: null };
    if (!rel) {
      st.reason = mya === null ? "no_distance" : "no_relative";
    } else if (asked === null) {
      st.band = "all"; st.reason = ps && ps.models ? "few_models" : "no_plddt";
    } else {
      const c = S.cells[`${rel.id}|${asked}`];
      st.n_asked = c ? c.n : 0;
      if (c && c.n >= S.min_n) {
        st.band = asked;
      } else {
        st.fallback = true;
        const ok = BANDS.filter((b) => { const x = S.cells[`${rel.id}|${b}`]; return x && x.n >= S.min_n; });
        ok.sort((a, b) => Math.abs((bandLo(a) + bandHi(a)) / 2 - meanP) - Math.abs((bandLo(b) + bandHi(b)) / 2 - meanP) || bandLo(a) - bandLo(b));
        st.band = ok.length ? ok[0] : "all";
      }
    }
    A.stratum = st;

    // Verdict: the three tools with the highest Fmax in this stratum, each at the setting
    // that scored best on all domains at this distance.
    A.verdict = [];
    A.not_run = [];
    if (rel) {
      const cell = S.cells[`${rel.id}|${st.band}`];
      st.n = cell ? cell.n : 0;
      const ranked = Object.entries(cell ? cell.tools : {}).filter(([, v]) => v.fmax !== null)
        .sort((a, b) => b[1].fmax - a[1].fmax || (a[0] < b[0] ? -1 : 1));
      A.not_run = Object.keys(S.tools).filter((t) => !(cell && cell.tools[t]));
      for (const [tool, v] of ranked.slice(0, 3)) {
        const key = `${tool}|${v.setting}`;
        const sh = S.short[`${key}|${rel.id}`] || null;
        const pu = S.pileup && S.pileup.panels[key] ? S.pileup.panels[key].bands[st.band] || null : null;
        const cost = S.cost[key] || null;
        A.verdict.push({
          tool, name: S.tools[tool].name, kind: S.tools[tool].kind, shape: S.tools[tool].shape,
          setting: v.setting, setting_label: S.setting_label[key] || v.setting,
          fmax: v.fmax, recall: v.recall, precision: v.precision, overlap_median: v.iou, n: v.n,
          short_recall: sh ? sh.recall : null, short_precision: sh ? sh.precision : null, short_n: sh ? sh.n : null,
          pileup: pu ? { species: S.pileup.species, band: st.band, n_domains: pu.n, curve: pu.curve,
                         median_start_offset_aa: pu.median_start_offset_aa, n_tp_calls: pu.n_tp_calls } : null,
          cpu_h_per_1000: cost ? cost.cpu_h_per_1000 : null, gpu_h_per_1000: cost ? cost.gpu_h_per_1000 : null,
          cost_usd_est: cost ? cost.cost_usd_est : null, db_cpu_h_per_db: cost ? cost.db_cpu_h_per_db : null,
        });
      }
    }

    // Overlap of the three tools' correct finds (pile-up domains at zebrafish, same band).
    A.upset = null;
    if (S.pileup && A.verdict.length) {
      const sets = A.verdict.map((v) => {
        const p = S.pileup.panels[`${v.tool}|${v.setting}`];
        return p ? new Set(p.tp_rows) : null;
      });
      const inBand = Object.entries(S.pileup.region_band)
        .filter(([, b]) => st.band === "all" || b === st.band).map(([id]) => Number(id));
      if (sets.every((s) => s)) {
        const combos = {};
        for (const id of inBand) {
          const key = sets.map((s) => (s.has(id) ? "1" : "0")).join("");
          combos[key] = (combos[key] || 0) + 1;
        }
        A.upset = { species: S.pileup.species, band: st.band, n_domains: inBand.length,
                    tools: A.verdict.map((v) => v.tool), combos };
      } else {
        A.upset = { species: S.pileup.species, band: st.band, missing: A.verdict.filter((v, i) => !sets[i]).map((v) => v.tool) };
      }
    }

    // Distance curve: each verdict tool at its setting, Fmax on all domains, every relative.
    A.curve = A.verdict.map((v) => ({
      tool: v.tool, setting: v.setting,
      points: S.relatives.map((r) => ({ relative: r.id, mya: r.mya, fmax: (S.curve[`${v.tool}|${v.setting}`] || {})[r.id] ?? null })),
    }));

    // Dark proteins: the measured organism that shares the deepest ancestor with this species.
    const lin = L.lineage;
    let bestD = null;
    for (const d of S.dark) {
      const sd = sharedDepth(lin, d.lineage || [d.taxid]);
      if (sd && (!bestD || sd.i < bestD.sd.i)) bestD = { d, sd };
    }
    A.dark = bestD ? Object.assign({}, bestD.d, { shared_ancestor: { taxid: bestD.sd.taxid, name: bestD.sd.name } }) : null;
    A.dark_all = S.dark.map((d) => ({ organism: d.organism, name: d.name, dark_fraction: d.dark_fraction, proteins: d.proteins }));

    // expected_label_fraction_est: the first tool's recall in each pLDDT band at this distance,
    // weighted by the share of this species' AlphaFold models in that band.
    A.expected_label_fraction_est = null;
    if (rel && A.verdict.length && enough && ps.share) {
      const v0 = A.verdict[0];
      let w = 0, s = 0;
      const parts = [];
      for (const b of BANDS) {
        const c = S.cells[`${rel.id}|${b}`];
        const t = c && c.tools[v0.tool] && c.tools[v0.tool].setting === v0.setting ? c.tools[v0.tool] : null;
        parts.push({ band: b, share: ps.share[b], recall: t ? t.recall : null, n: c ? c.n : 0 });
        if (t && t.recall !== null) { w += ps.share[b]; s += ps.share[b] * t.recall; }
      }
      A.expected_label_fraction_est = { value: w ? s / w : null, share_of_models_covered: w, tool: v0.tool,
        setting: v0.setting, relative: rel.id, parts,
        rests_on: `human proteins (Pfam 38.2 domains) searched against ${rel.name} (${rel.scientific}), ${num(rel.mya)} Mya` };
    }
    A.command = command(A, S);
    return A;
  }

  // ------------------------------------------------------------------ the command
  function toolCommand(tool, S) {
    const cont = S.tools[tool].container;
    const c = {
      hmmer3_phmmer: [`# phmmer (${cont})`, "gunzip -kf proteome.fasta.gz swissprot.fasta.gz",
        "phmmer --domtblout phmmer.domtbl --noali -E 10.0 --cpu 8 proteome.fasta swissprot.fasta > /dev/null"],
      hmmer3_jackhmmer: [`# jackhmmer, three rounds (${cont})`, "gunzip -kf proteome.fasta.gz swissprot.fasta.gz",
        "jackhmmer -N 3 --domtblout jackhmmer.domtbl --noali -E 10.0 --cpu 8 proteome.fasta swissprot.fasta > /dev/null"],
      mmseqs2_seqseq: [`# MMseqs2, sensitivity 7, one round (${cont})`,
        "mmseqs createdb proteome.fasta.gz q_db", "mmseqs createdb swissprot.fasta.gz t_db",
        "mmseqs search q_db t_db result_s7 tmp_s7 --threads 16 -s 7 --max-seqs 1000 -e 10.0",
        "mmseqs convertalis q_db t_db result_s7 mmseqs2_s7.m8"],
      mmseqs2_iterative: [`# MMseqs2, sensitivity 7, three rounds (${cont})`,
        "mmseqs createdb proteome.fasta.gz q_db", "mmseqs createdb swissprot.fasta.gz t_db",
        "mmseqs search q_db t_db result_s7_n3 tmp_s7_n3 --threads 16 -s 7 --num-iterations 3 --max-seqs 1000 -e 10.0",
        "mmseqs convertalis q_db t_db result_s7_n3 mmseqs2_s7_3rounds.m8"],
    }[tool];
    if (c) return c;
    return [`# ${S.tools[tool].name}: its database needs build steps this page does not write out.`,
            `# The command the board ran is on its tool card: tools/${tool}.html`];
  }

  function command(A, S) {
    const tools = A.verdict.map((v) => v.tool);
    const sample = (A.slug || `taxon_${A.taxid}`).slice(0, 40);
    const lines = [`# ${displayName(A.name, A.common_name)}`];
    if (A.proteome) {
      lines.push(`# UniProt reference proteome ${A.proteome.upid}, ${num(A.proteome.protein_count)} proteins`,
        `curl -L -o proteome.fasta.gz 'https://rest.uniprot.org/uniprotkb/stream?query=proteome%3A${A.proteome.upid}&format=fasta&compressed=true'`);
    } else {
      lines.push("# UniProt has no reference proteome for this species: put your own protein FASTA at proteome.fasta.gz");
    }
    // Search all of Swiss-Prot, not only the nearest relative's few reviewed proteins: for the
    // oyster those are 158 cone-snail proteins. The relative sets the distance, nothing else.
    const r = A.cards.relative;
    lines.push(`# All reviewed Swiss-Prot proteins, the source of the labels${r ? ` (they include those of the nearest reviewed relative, ${r.name})` : ""}.`,
      "# Swiss-Prot has about 575,000 proteins; the board measured cost against one relative's proteome of about 20,000, so expect longer runs.",
      "curl -L -o swissprot.fasta.gz 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'");
    lines.push(`# nf-core/proteinannotator ${PA_VERSION}: Pfam, FunFam and InterProScan on your proteome.`,
      tools.length === 1 ? `# Version ${PA_VERSION} has no option to run ${S.tools[tools[0]].name}; its command follows.`
        : `# Version ${PA_VERSION} has no option to run ${tools.map((t) => S.tools[t].name).join("; ") || "the tools below"}; their commands follow.`,
      `printf 'id,fasta\\n${sample},proteome.fasta.gz\\n' > samplesheet.csv`,
      `nextflow run nf-core/proteinannotator -r ${PA_VERSION} -profile docker --input samplesheet.csv --outdir results`);
    // A step two tools share (unzipping, building the MMseqs2 databases) is written once.
    const done = new Set();
    if (tools.length) lines.push("", "# -E 10.0 keeps every hit, as on the board. The board's numbers count calls above the score cut-off",
      "# that gives each tool its best Fmax; filter the output before you take labels from it.");
    for (const t of tools) {
      lines.push("");
      for (const l of toolCommand(t, S)) {
        if (!l.startsWith("#") && done.has(l)) continue;
        done.add(l); lines.push(l);
      }
    }
    return lines.join("\n");
  }

  // ------------------------------------------------------------------ rendering
  const KIND_VAR = (k) => `var(--k-${k})`;
  function marker(shape, kind, x, y, r = 4.5) {
    const c = KIND_VAR(kind);
    if (shape === "square") return `<rect x="${x - r}" y="${y - r}" width="${2 * r}" height="${2 * r}" fill="${c}"/>`;
    if (shape === "diamond") return `<path d="M${x},${y - r - 1} L${x + r + 1},${y} L${x},${y + r + 1} L${x - r - 1},${y} Z" fill="${c}"/>`;
    if (shape === "triangle") return `<path d="M${x},${y - r - 1} L${x + r + 1},${y + r} L${x - r - 1},${y + r} Z" fill="${c}"/>`;
    return `<circle cx="${x}" cy="${y}" r="${r}" fill="${c}"/>`;
  }
  function legendMark(v) {
    return `<svg width="16" height="16" viewBox="0 0 16 16" aria-hidden="true">${marker(v.shape, v.kind, 8, 8, 5)}</svg>`;
  }

  // Pile-up thumbnail: share of domains covered by a call at each of 300 columns (100
  // residues before, the domain scaled to 100 columns, 100 residues after).
  function thumb(v, idx, S) {
    if (!v.pileup) return `<p class="muted small">Pile-up: not run for ${esc(v.name)} at this setting.</p>`;
    const W = 307, H = 70, P = 4;  // one column per x unit
    const pts = v.pileup.curve.map((y, j) => (y === null ? null : [P + j, H - 14 - y * (H - 22)]));
    let d = "", pen = false;
    for (const p of pts) { if (!p) { pen = false; continue; } d += (pen ? "L" : "M") + p[0] + "," + p[1].toFixed(1); pen = true; }
    const x0 = P + 100, x1 = P + 199;
    return `<figure class="thumb"><figcaption><b>Figure 1${"abc"[idx]}.</b> ${esc(v.name)}: at each residue, the share of ${num(v.pileup.n_domains)} Pfam domains${v.pileup.band === "all" ? "" : ` whose own stretch of the protein has pLDDT ${esc(v.pileup.band)} (grouped by the domain, not the whole protein, so the count differs from the one above)`} that a call covers. Human proteins searched against zebrafish (${relMya(S, "zebrafish")} Mya); the grey band is the domain, every domain stretched to one width.</figcaption>
<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Pile-up of ${esc(v.name)}" class="thumbsvg">
<rect x="${x0}" y="4" width="${x1 - x0}" height="${H - 18}" fill="var(--chip)"/>
<line x1="${P}" y1="${H - 14}" x2="${W - P}" y2="${H - 14}" stroke="var(--hair)"/>
<line x1="${P}" y1="8" x2="${W - P}" y2="8" stroke="var(--hair)" stroke-dasharray="2 3"/>
<path d="${d}" fill="none" stroke="${KIND_VAR(v.kind)}" stroke-width="2"/>
<text x="${P}" y="${H - 2}" class="tk">100 aa before</text><text x="${(x0 + x1) / 2}" y="${H - 2}" class="tk" text-anchor="middle">domain</text><text x="${W - P}" y="${H - 2}" class="tk" text-anchor="end">100 aa after</text>
<text x="${W - P}" y="20" class="tk" text-anchor="end">dashes: all covered</text>
</svg></figure>`;
  }

  function upsetSvg(A, S) {
    const U = A.upset;
    if (!U) return "";
    const head = `<figure><figcaption><b>Figure 2.</b> Which domains each of the three tools finds, and which only one of them finds. The tools searched human proteins against ${esc(U.species)} (${relMya(S, "zebrafish")} Mya)${U.band === "all" ? "" : `; domains whose own stretch of the protein has pLDDT ${esc(U.band)}`}. Per-domain calls were kept only for zebrafish, so this figure is at ${relMya(S, "zebrafish")} Mya whatever this species' distance. A filled dot under a bar means that tool found the domain; the bar is the number of domains found by that set of tools and no other.</figcaption>`;
    if (U.missing) return head + `<p class="muted">Not run: the pile-up tables hold no calls for ${esc(U.missing.join(", "))} at the setting chosen here, so the overlap cannot be counted.</p></figure>`;
    const order = Object.entries(U.combos).filter(([k]) => k !== "000").sort((a, b) => b[1] - a[1]);
    const none = U.combos["000"] || 0;
    const W = 560, barH = 120, rowH = 22, left = 150, top = 16;
    const n = order.length, colW = Math.min(52, (W - left - 10) / Math.max(n, 1));
    const max = Math.max(1, ...order.map(([, c]) => c));
    let s = `<svg viewBox="0 0 ${W} ${top + barH + 3 * rowH + 30}" role="img" aria-label="Overlap of correct finds" class="figsvg">`;
    order.forEach(([k, c], i) => {
      const x = left + i * colW + colW * 0.2, h = c / max * barH, w = colW * 0.6;
      s += `<rect x="${x}" y="${top + barH - h}" width="${w}" height="${h}" fill="var(--ink)"/>`;
      s += `<text x="${x + w / 2}" y="${top + barH - h - 3}" text-anchor="middle" class="tk">${c}</text>`;
      const ys = [];
      A.verdict.forEach((v, j) => {
        const y = top + barH + 14 + j * rowH;
        const on = k[j] === "1";
        if (on) ys.push(y);
        s += on ? `<circle cx="${x + w / 2}" cy="${y}" r="5" fill="var(--ink)"/>`
          : `<circle cx="${x + w / 2}" cy="${y}" r="5" fill="none" stroke="var(--hair)" stroke-width="1.5"/>`;
      });
      if (ys.length > 1) s += `<line x1="${x + w / 2}" y1="${ys[0]}" x2="${x + w / 2}" y2="${ys[ys.length - 1]}" stroke="var(--ink)" stroke-width="2"/>`;
    });
    A.verdict.forEach((v, j) => {
      const y = top + barH + 14 + j * rowH;
      s += `<line x1="${left}" y1="${y}" x2="${W - 6}" y2="${y}" stroke="var(--hair)" stroke-width="0.6"/>`;
      s += `<text x="${left - 8}" y="${y + 4}" text-anchor="end" class="lb">${esc(v.name)}</text>`;
    });
    s += `<text x="${left - 8}" y="${top + 10}" text-anchor="end" class="tk">domains</text></svg>`;
    return head + s + `<p class="small muted">${num(U.n_domains)} domains in this band; ${num(none)} of them found by none of the three. A domain counts as found when any call on it names the right Pfam family and overlaps it, at any score.</p></figure>`;
  }

  function curveSvg(A, S) {
    if (!A.curve.length) return "";
    const W = 600, H = 260, L0 = 48, R0 = 120, T0 = 14, B0 = 40;
    const xmin = 5;
    const X = (m) => L0 + (Math.log10(m) - Math.log10(xmin)) / (Math.log10(2500) - Math.log10(xmin)) * (W - L0 - R0);
    const ymax = Math.max(0.1, ...A.curve.flatMap((c) => c.points.map((p) => p.fmax || 0)));
    const Y = (f) => H - B0 - f / (Math.ceil(ymax * 10) / 10) * (H - T0 - B0);
    const ytop = Math.ceil(ymax * 10) / 10;
    let s = `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Fmax against distance" class="figsvg">`;
    for (let t = 0; t <= ytop + 1e-9; t += 0.1) {
      s += `<line x1="${L0}" y1="${Y(t)}" x2="${W - R0}" y2="${Y(t)}" stroke="var(--hair)" stroke-width="0.6"/><text x="${L0 - 6}" y="${Y(t) + 4}" text-anchor="end" class="tk">${t.toFixed(1)}</text>`;
    }
    for (const m of [10, 30, 100, 300, 1000, 2000]) {
      s += `<line x1="${X(m)}" y1="${T0}" x2="${X(m)}" y2="${H - B0}" stroke="var(--hair)" stroke-width="0.6"/><text x="${X(m)}" y="${H - B0 + 16}" text-anchor="middle" class="tk">${m}</text>`;
    }
    s += `<text x="${(L0 + W - R0) / 2}" y="${H - 6}" text-anchor="middle" class="lb">million years to the relative searched (Mya)</text>`;
    s += `<text x="12" y="${(T0 + H - B0) / 2}" transform="rotate(-90 12 ${(T0 + H - B0) / 2})" text-anchor="middle" class="lb">Fmax</text>`;
    const my = A.cards.relative ? A.cards.relative.mya : null;
    if (my) {
      const xm = X(Math.max(xmin, Math.min(2500, my)));
      s += `<line x1="${xm}" y1="${T0}" x2="${xm}" y2="${H - B0}" stroke="var(--ink)" stroke-width="1.4" stroke-dasharray="5 4"/>`;
      s += `<text x="${xm + 4}" y="${T0 + 10}" class="tk">${esc(A.common_name || A.name)}: ${num(my)} Mya</text>`;
    }
    A.curve.forEach((c, i) => {
      const v = A.verdict[i];
      const pts = c.points.filter((p) => p.fmax !== null && p.mya).sort((a, b) => a.mya - b.mya);
      // Two relatives at the same age (cow and dog, 89 Mya) are drawn as two marks; the
      // line goes through their mean.
      const byMya = {};
      for (const p of pts) (byMya[p.mya] = byMya[p.mya] || []).push(p.fmax);
      const line = Object.entries(byMya).map(([m, f]) => [Number(m), f.reduce((a, b) => a + b, 0) / f.length]).sort((a, b) => a[0] - b[0]);
      s += `<path d="${line.map(([m, f], j) => (j ? "L" : "M") + X(m).toFixed(1) + "," + Y(f).toFixed(1)).join("")}" fill="none" stroke="${KIND_VAR(v.kind)}" stroke-width="1.6"/>`;
      for (const p of pts) s += marker(v.shape, v.kind, X(p.mya), Y(p.fmax), 3.6);
      const last = line[line.length - 1];
      if (last) s += `<text x="${X(last[0]) + 8}" y="${Y(last[1]) + 4 + i * 12 - 12}" class="lb">${esc(v.name)}</text>`;
    });
    s += "</svg>";
    const legend = `<div class="legend">${A.verdict.map((v) => `<span class="item">${legendMark(v)}${esc(v.name)} (${esc(v.setting_label)})</span>`).join("")}<span class="item"><svg width="22" height="12" aria-hidden="true"><line x1="0" y1="6" x2="22" y2="6" stroke="var(--ink)" stroke-width="1.4" stroke-dasharray="5 4"/></svg>this species' distance to its nearest reviewed relative</span></div>`;
    return `<figure><figcaption><b>Figure 3.</b> How Fmax falls with distance, for the three tools above. Fmax is the best balance of finding domains and being right (1 is perfect). Each mark is one relative the human proteins were searched against, all Pfam domains${sharedAge(S) ? `; where two relatives share an age (${sharedAge(S)}) the line goes through their mean` : ""}. Divergence times are TimeTree 5's${S.relatives.some((r) => r.mya === null) ? `; ${S.relatives.filter((r) => r.mya === null).map((r) => r.name).join(" and ")}, with no TimeTree age, is not drawn` : ""}. The step from 10 to 100 Mya takes as much room as the step from 100 to 1,000.</figcaption>${legend}${s}</figure>`;
  }

  function card(title, value, sub, source) {
    return `<div class="mcard"><div class="mtitle">${title}</div><div class="mvalue">${value}</div><div class="msub">${sub}</div><div class="msrc">${source}</div></div>`;
  }

  function renderAnswer(A, S, src) {
    const c = A.cards, st = A.stratum, rel = c.relative;
    const nm = displayName(A.name, A.common_name);
    const mya = rel && rel.mya !== null ? num(rel.mya) : "—";
    const v = A.verdict;
    const dark = A.dark;
    // The three sentences.
    // The three sentences. Reviewed = checked by a UniProt curator (Swiss-Prot).
    const s1 = rel ? `The closest species with at least 100 curator-checked (Swiss-Prot) proteins is ${esc(displayName(rel.name, rel.common_name))}; ${rel.mya !== null ? `the two lineages split ${mya} million years ago (Mya)` : "TimeTree has no age for when the two lineages split"}.`
      : "No species with at least 100 curator-checked (Swiss-Prot) proteins could be placed on this species' lineage.";
    const s2 = v.length ? `In human proteins searched against a species that far away (${esc(st.relative.name)}, ${num(st.relative.mya)} Mya), the three tools with the best Fmax (the balance of finding domains and being right) were ${v.map((x) => esc(x.name)).join("; ")}, and ${esc(v[0].name)} found ${pct(v[0].recall)} of the known Pfam domains.`
      : "No tool on this board was measured at this distance.";
    const s3 = dark ? `In the measured organism closest to this one, ${esc(dark.name)}, ${pct(dark.dark_fraction)} of proteins had no hit from any sequence tool, so expect a share like that to stay without a label.` : "";
    let h = `<h1>${esc(nm)}</h1>
<p class="lede three">${s1} ${s2} ${s3}</p>
<div class="cmdrow"><button type="button" class="copy" data-copy="cmd-all">Copy command</button> <span class="small muted">downloads this species' proteins and Swiss-Prot, runs nf-core/proteinannotator ${PA_VERSION}, then the three tools; needs a Unix computer with curl, Nextflow and Docker</span></div>
<pre class="cmd" id="cmd-all"><code>${esc(A.command)}</code></pre>`;

    // Four cards.
    const relSub = rel ? `${esc(displayName(rel.name, rel.common_name))}, ${num(rel.entries)} reviewed proteins; lineages meet at ${esc(rel.shared_ancestor.name)} (${esc(rel.shared_ancestor.rank)})` : "none found";
    const mm = rel && rel.mya_detail ? rel.mya_detail : null;
    const myaSrc = !mm ? "TimeTree: no age for this pair" : mm.method === "pairwise"
      ? `TimeTree 5 (${esc(S.timetree_release)}), pairwise age of the two species`
      : mm.method === "node" ? `TimeTree 5 (${esc(S.timetree_release)}), age of their shared ancestor ${esc(rel.shared_ancestor.name)}`
        : `TimeTree 5 (${esc(S.timetree_release)}), pairwise age of two reviewed species that split at ${esc(rel.shared_ancestor.name)}`;
    const ord = c.reviewed.order, cls = c.reviewed.class;
    const afc = c.alphafold;
    h += `<h2 id="cards">Four facts about this species</h2><div class="mcards">
${card("Nearest reviewed relative", rel ? `${mya} Mya` : "—", relSub, `${myaSrc}; Swiss-Prot counts from ${esc(S.uniprot_release)}`)}
${card("Reviewed proteins nearby", `${ord ? num(ord.entries) : "—"} <span class="unit">in the order</span>`,
    `${ord ? esc(ord.name) : "NCBI gives no order"}; class ${cls ? esc(cls.name) + ": " + num(cls.entries) : "not given by NCBI"}; this species: ${num(c.reviewed.self)}`,
    `Swiss-Prot reviewed entries, ${esc(S.uniprot_release)}, summed over NCBI Taxonomy`)}
${card("Proteins with an AlphaFold model", afc ? pct(afc.share_with_model) : "—",
    afc ? `${num(afc.with_model)} of ${num(afc.proteome_entries)} proteins in ${esc(A.proteome.upid)}; mean pLDDT ${afc.mean_plddt_est === null ? "—" : num(afc.mean_plddt_est, 1)} <span class="est">estimate</span>, from ${num(afc.afdb_models_counted)} AlphaFold DB models of this organism's reference-proteome entries, counted in 5-point bands (AlphaFold DB keeps entries UniProt has since dropped, so this count can differ from the protein count)`
      : A.proteome ? `not counted for ${esc(A.proteome.upid)}` : "no UniProt reference proteome", afc ? `${esc(S.uniprot_release)} cross-references; AlphaFold DB search, models version 6` : "UniProt proteomes")}
${card("Earth BioGenome Project target", !c.ebp ? "—" : c.ebp.on_list ? "on the list" : "not on the list",
    !c.ebp ? "GoaT did not answer" : c.ebp.on_list ? `listed by ${esc(String(c.ebp.projects || "").split(",").join(", "))}` : "no project lists it",
    `GoaT, attribute long_list (${esc(S.goat_release)})`)}
</div>`;

    // Verdict.
    h += `<h2 id="verdict">Three tools, in the order to run them</h2>`;
    if (!st.relative) {
      h += `<p class="notice">No verdict: ${st.reason === "no_distance" ? "TimeTree gives no age between this species and its nearest reviewed relative" : "no measured relative"}, so the page cannot place it at a distance.</p>`;
    } else {
      if (st.fallback) {
        h += `<p class="notice">Fewer than ${S.min_n} measured domains in this species' stratum (${esc(st.relative.name)}, pLDDT ${esc(st.band_asked)}: ${num(st.n_asked)} domains). The tools below are ranked in the nearest stratum that has ${S.min_n} or more: pLDDT ${esc(st.band)}.</p>`;
      } else if (st.reason === "no_plddt") {
        h += `<p class="notice">No AlphaFold models counted for this species, so the ranking uses all domains at this distance.</p>`;
      } else if (st.reason === "few_models") {
        h += `<p class="notice">Only ${num(st.models_counted)} AlphaFold models counted for this species, fewer than ${MIN_MODELS}, so the ranking uses all domains at this distance.</p>`;
      }
      const gap = v.length > 1 ? v[0].fmax - v[v.length - 1].fmax : null;
      h += `<p class="small">How this was measured: ${num(S.n_queries)} human proteins were searched against the proteome of ${esc(st.relative.name)} (<i>${esc(st.relative.scientific)}</i>), whose lineage split from human ${num(st.relative.mya)} Mya. Of the ${S.relatives.filter((r) => r.mya !== null).length} relatives measured that have a TimeTree age, its split age is the closest to this species' ${mya} Mya. Split age is a stand-in for how much proteins have changed; lineages change at different speeds. Counted: ${num(st.n)} Pfam domains in human proteins ${st.band === "all" ? "of any AlphaFold confidence" : `whose mean AlphaFold confidence (pLDDT, 0–100; AlphaFold DB calls 70 to 90 confident and above 90 very high) is ${esc(st.band)}`}. Tools are ranked by Fmax, the best balance of finding domains (recall) and being right (precision), 1 being perfect${gap !== null ? `; the three differ by ${num(gap, 3)} in Fmax, a gap this one set of proteins cannot settle, so if cost matters pick the cheapest of them` : ""}.</p><ol class="verdict">`;
      v.forEach((x, i) => {
        const place = x.pileup && x.pileup.median_start_offset_aa !== null ? `; starts a found domain a median ${num(x.pileup.median_start_offset_aa)} residues from where the domain starts (measured at zebrafish, ${relMya(S, "zebrafish")} Mya)` : "; placement: not run";
        const sh = x.short_recall !== null ? `; finds ${pct(x.short_recall)} of short Swiss-Prot features 2–15 residues long, such as motifs and binding sites (${pct(x.short_precision)} of those calls right)` : "; short features: not run";
        const sameName = x.name.toLowerCase().includes(x.setting_label.split(",")[0].toLowerCase());
        h += `<li class="vtool"><h3>${legendMark(x)} ${i + 1}. ${esc(x.name)}${sameName ? "" : ` <span class="muted small">${esc(x.setting_label)}</span>`}</h3>
<p>Finds ${pct(x.recall)} of Pfam domains at this distance, and ${pct(x.precision)} of the domains it reports (its calls) are right${place}${sh}.</p>
<p class="small">Counted at the score cut-off that gives this tool its best Fmax: ${num(x.fmax, 3)} on ${num(x.n)} domains${x.short_n ? `; short features from ${num(x.short_n)} Swiss-Prot features` : ""}. A found domain and its call overlap over a median ${pct(x.overlap_median)} of the span they cover together.</p>
<p class="small cost">Computer time per 1,000 proteins searched against one relative's proteome: ${x.cpu_h_per_1000 === null ? "— processor-hours (not run)" : num(x.cpu_h_per_1000, 2) + " processor-hours (CPU-h)"}${x.gpu_h_per_1000 ? ", " + num(x.gpu_h_per_1000, 2) + " GPU-hours" : ""}; cloud price ${x.cost_usd_est === null ? "—" : "$" + num(x.cost_usd_est, 3)} <span class="est">estimate</span> (cost_usd_est). Search step only; ${x.db_cpu_h_per_db === null ? "no database-building step was recorded for it" : "building the relative's database took " + num(x.db_cpu_h_per_db, 2) + " CPU-h more, once per relative"}.</p>
${thumb(x, i, S)}
<button type="button" class="copy" data-copy="cmd-${i}">Copy command</button><pre class="cmd" id="cmd-${i}"><code>${esc(commandFor(A, S, i))}</code></pre></li>`;
      });
      h += "</ol>";
      if (A.not_run.length) h += `<p class="small muted">Not run at this distance and pLDDT: ${A.not_run.map((t) => esc(S.tools[t].name)).join(", ")} (shown as — in the tables). Pfam's own models (hmmscan) are left out of the ranking, because the answer key is made from those models.</p>`;
    }

    // Dark box.
    h += `<section class="darkbox"><h2 id="dark">What no tool on this board will find for you</h2>`;
    if (dark) {
      h += `<p><b>${pct(dark.dark_fraction, 1)}</b> of ${num(dark.proteins)} proteins of ${esc(dark.name)} had no hit from any sequence tool, even at a loose E-value below 1 (the E-value is the number of hits that size expected by chance). They were searched against Swiss-Prot with the organism's own group removed; we call these proteins dark. Of the five organisms measured, this one is closest to this species: their lineages meet at ${esc(dark.shared_ancestor.name)}.</p>
<table class="mini"><thead><tr><th>Part of the dark share</th><th>Measured</th></tr></thead><tbody>
<tr><td>no relative in Swiss-Prot (no hit at all)</td><td>${pct(dark.no_relative_fraction, 1)}</td></tr>
<tr><td>disordered (no stable fold)</td><td>— not run</td></tr>
<tr><td>short motif only</td><td>— not run</td></tr></tbody></table>
<p class="small muted">The tables hold the dark share but do not split it into disordered proteins and short motifs; those two parts are not measured. Measured on: ${A.dark_all.map((d) => `${esc(d.name)} ${pct(d.dark_fraction, 1)} of ${num(d.proteins)}`).join("; ")}.</p>`;
    } else {
      h += `<p>— not run: no measured organism shares an ancestor with this species.</p>`;
    }
    h += "</section>";

    // Estimate.
    const E = A.expected_label_fraction_est;
    h += `<h2 id="estimate">Expected to get a label</h2>`;
    if (E && E.value !== null) {
      h += `<p><span class="bignum">${pct(E.value)}</span> <span class="est">estimate</span> (expected_label_fraction_est): the share of this species' Pfam domains, not proteins, that ${esc(S.tools[E.tool].name)} would be expected to find, if it finds them as often as it did in ${esc(E.rests_on)}.</p>
<p class="small">How: ${esc(S.tools[E.tool].name)}'s measured recall in each pLDDT band at that distance, weighted by the share of this species' AlphaFold models in the band. The bands with a measured recall hold ${pct(E.share_of_models_covered)} of its models; the rest (pLDDT 0–50, no measured domains) are left out, so for the whole proteome the share is likely lower.</p>
<table class="mini"><thead><tr><th>pLDDT band</th><th>share of this species' models</th><th>recall measured</th><th>domains measured</th></tr></thead><tbody>${E.parts.map((p) => `<tr><td>${p.band}</td><td>${pct(p.share, 1)}</td><td>${p.recall === null ? "— not run" : pct(p.recall, 1)}</td><td>${num(p.n)}</td></tr>`).join("")}</tbody></table>`;
    } else {
      h += `<p>— Cannot be estimated: ${!afc ? "no AlphaFold models counted for this species" : "no verdict"}.</p>`;
    }

    // Figures.
    h += `<h2 id="figures">Figures</h2>${upsetSvg(A, S)}${curveSvg(A, S)}`;
    return h;
  }

  function commandFor(A, S, i) {
    const B = Object.assign({}, A, { verdict: [A.verdict[i]] });
    return command(B, S);
  }

  // ------------------------------------------------------------------ browser lookups
  async function getJSON(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return { json: await r.json(), headers: r.headers };
  }
  async function count(url) {
    const r = await fetch(url);
    if (!r.ok) throw new Error(`${r.status} ${url}`);
    return Number(r.headers.get("X-Total-Results"));
  }
  const UNI = "https://rest.uniprot.org";
  const AFDB = "https://alphafold.ebi.ac.uk/api/search";
  const GOAT = "https://goat.genomehubs.org/api/v2";
  // encodeURIComponent leaves ( and ) as they are; GoaT refuses them unencoded (HTTP 400).
  const q = (s) => encodeURIComponent(s).replace(/\(/g, "%28").replace(/\)/g, "%29");

  async function lookupSpecies(taxid, T, progress = () => {}) {
    // NCBI's own record, the same source as the build's taxdump. UniProt's taxonomy lags NCBI
    // (Phaeodactylum tricornutum: order Naviculales in UniProt 2026_03, Surirellales at NCBI).
    progress("NCBI Taxonomy lineage");
    const r = await fetch(`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=${taxid}&retmode=xml`);
    if (!r.ok) throw new Error(`${r.status} NCBI Taxonomy`);
    const doc = new DOMParser().parseFromString(await r.text(), "application/xml");
    const top = doc.querySelector("TaxaSet > Taxon");
    if (!top) throw new Error(`NCBI Taxonomy has no taxon ${taxid}`);
    const kid = (el, tag) => { const c = Array.from(el.children).find((x) => x.tagName === tag); return c ? c.textContent : null; };
    const self = { taxid: Number(kid(top, "TaxId")), rank: kid(top, "Rank"), name: kid(top, "ScientificName") };
    const lex = top.querySelector("LineageEx");
    const anc = lex ? Array.from(lex.children).map((t) => ({ taxid: Number(kid(t, "TaxId")), rank: kid(t, "Rank"), name: kid(t, "ScientificName") })) : [];
    const lineage = [self].concat(anc.reverse());
    const other = top.querySelector("OtherNames");
    const common = other ? (kid(other, "GenbankCommonName") || kid(other, "CommonName")) : null;
    const L = { taxid: self.taxid, name: self.name, common_name: common, lineage, built: "browser" };
    L.slug = self.name.toLowerCase().replace(/[^a-z0-9]+/g, "_").replace(/^_|_$/g, "");
    L.self_entries = T.reviewed.self_entries[String(L.taxid)] || 0;
    L.ranks = {};
    for (const r of ["order", "class"]) {
      const a = lineage.find((x) => x.rank === r);
      L.ranks[r] = a ? { taxid: a.taxid, name: a.name, entries: T.reviewed.rank_entries[String(a.taxid)] || 0 } : null;
    }
    const nr = nearestReviewed(lineage, L.taxid, T.reviewed);
    if (nr) {
      const nmc = T.reviewed.names[String(nr.taxid)] || [null, null];
      L.relative = { taxid: nr.taxid, name: nmc[0], common_name: nmc[1], entries: nr.entries,
                     shared_ancestor: { taxid: nr.shared_ancestor.taxid, name: nr.shared_ancestor.name, rank: nr.shared_ancestor.rank } };
      const ag = T.mrca.ages[String(nr.shared_ancestor.taxid)];
      L.mya = ag ? { value: ag[0], method: ag[1] === "node" ? "node" : "node_pair", node: nr.shared_ancestor.taxid,
                     pair: ag[1] === "pair" ? [ag[2], ag[3]] : null } : null;
    } else { L.relative = null; L.mya = null; }

    progress("UniProt reference proteome");
    const pr = (await getJSON(`${UNI}/proteomes/search?query=${q(`taxonomy_id:${taxid} AND reference:true`)}&fields=upid,organism_id,protein_count&size=50`)).json.results || [];
    let best = null;
    for (const p of pr) {
      const oid = p.taxonomy.taxonId;
      const rev = await count(`${UNI}/uniprotkb/search?query=${q(`reviewed:true AND organism_id:${oid}`)}&size=0`);
      const cand = { upid: p.id, taxid: oid, protein_count: p.proteinCount, rev };
      if (!best || rev > best.rev || (rev === best.rev && (cand.protein_count > best.protein_count
          || (cand.protein_count === best.protein_count && cand.upid < best.upid)))) best = cand;
    }
    L.proteome = best ? { upid: best.upid, taxid: best.taxid, protein_count: best.protein_count } : null;
    if (best) {
      progress("AlphaFold models (UniProt and AlphaFold DB)");
      const withM = await count(`${UNI}/uniprotkb/search?query=${q(`proteome:${best.upid} AND database:alphafolddb`)}&size=0`);
      const all = await count(`${UNI}/uniprotkb/search?query=${q(`proteome:${best.upid}`)}&size=0`);
      const base = `taxId:${best.taxid} AND isUniProtReferenceProteome:true`;
      const bands5 = {};
      const jobs = [];
      for (let lo = 0; lo < 100; lo += 5) {
        const hi = lo + 5;
        const rng = hi === 100 ? `[${lo} TO ${hi}]` : `[${lo} TO ${hi}}`;
        const key = `plddt_${String(lo).padStart(2, "0")}_${String(hi).padStart(3, "0")}`;
        // AlphaFold DB answers 404 when nothing matches: that is a count of 0.
        jobs.push(getJSON(`${AFDB}?q=${q(`${base} AND globalMetricValue:${rng}`)}&type=main&rows=0`)
          .then((r) => { bands5[key] = r.json.numFound; })
          .catch((e) => { if (String(e.message).startsWith("404")) bands5[key] = 0; else throw e; }));
      }
      await Promise.all(jobs);
      const ordered = Object.fromEntries(Object.keys(bands5).sort().map((k) => [k, bands5[k]]));
      L.alphafold = { proteome_entries: all, with_model: withM, share_with_model: all ? withM / all : null, bands5: ordered };
    } else L.alphafold = null;

    progress("Earth BioGenome Project list (GoaT)");
    try {
      const g = (await getJSON(`${GOAT}/search?query=${q(`tax_eq(${taxid})`)}&result=taxon&taxonomy=ncbi&fields=long_list&size=1`)).json;
      const f = g.results && g.results[0] && g.results[0].result.fields && g.results[0].result.fields.long_list;
      const val = f ? (Array.isArray(f.value) ? f.value.join(",") : f.value) : null;
      L.ebp = { on_list: Boolean(val), projects: val };
    } catch (e) { L.ebp = null; }
    return L;
  }

  const api = { BANDS, bandOf, plddtSummary, nearestReviewed, measuredRelative, computeAnswer, renderAnswer, lookupSpecies, command };
  root.Organism = api;
  if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof window !== "undefined" ? window : globalThis);
