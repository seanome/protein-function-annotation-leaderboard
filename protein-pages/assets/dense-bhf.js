// One protein page as one grid: program rows that unfold into hit rows, and the residues under an open hit.
// Built from assets/protein.js (its helpers, alignment and dot-plot code) plus the grid layout at the end.
// PAGE is set by the page: the protein sequence (q), the drawing scale (px, w), the hits
// (k: track, l: label, q: range on this protein, s: score values, t: target, r: range on the
// target, a: alignment if made at build time, m: kmerseek alphabet), score names per track,
// kmerseek's class tables (alpha) and the number of Swiss-Prot sequence files (tseq_files).
// This is the BHF page's dense.js (claude.ai artifact JgFXM7JTxrdpgiw9B5dEXF, 2026-11-03), shipped as
// dense-bhf.js; changed here: the "not run" line gives each program's own reason (window.NOTRUN), and a
// ruler tick is left out where its label would touch the end label; PAGE.noseq (protein.html, no sequence on
// the site) gives hit rows without composition and opens no residues; the domain lane draws InterProScan's
// Pfam matches where hmmscan did not run.
(function () {
  const tip = document.getElementById("tip");
  const panel = document.getElementById("panel");
  const pbody = document.getElementById("pbody");
  const SVGNS = "http://www.w3.org/2000/svg";
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  // Links to other sites open in a new tab. accHtml escapes s and makes every Pfam (PF00069,
  // PF00069.32) and InterPro (IPR000719) accession in it a link; the accessions are letters and
  // digits, so escaping leaves them as they were. uniHtml links a UniProt accession (P12345,
  // A0A0B4J2F0, P12345-2) and only escapes anything else, such as a PDB chain from hhsearch.
  const ext = (url, text) => `<a class="mono" href="${esc(url)}" target="_blank" rel="noopener">${esc(text)}</a>`;
  const accHtml = (s) => esc(s).replace(/\b(PF\d{5}|IPR\d{6})(\.\d+)?\b/g, (m, a) =>
    ext(`https://www.ebi.ac.uk/interpro/entry/${a.startsWith("PF") ? "pfam" : "InterPro"}/${a}/`, m));
  const UNIPROT = /^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$/;
  const uniHtml = (acc) => UNIPROT.test(acc || "") ? ext(`https://www.uniprot.org/uniprotkb/${encodeURIComponent(acc)}/entry`, acc) : esc(acc);
  const H = PAGE.hits;
  // A target's Pfam domains are stored once per page (PAGE.td); give each hit its own reference
  if (PAGE.td) for (const h of H) if (h.t && PAGE.td[h.t]) h.td = PAGE.td[h.t];
  // and where they come from (interproscan, hmmscan or interpro), also once per page
  if (PAGE.ds) for (const h of H) if (h.t && PAGE.ds[h.t]) h.ds = PAGE.ds[h.t];

  // the programs' own column names, in plain words
  const GLOSS = { "i-Evalue": "E-value of this domain match alone (hmmer's i-Evalue)", "score": "score (bits)", "target env from-to": "envelope: hmmer's wider bound on the matched part of the other protein",
    "region_evalue": "E-value of this matched stretch (kmerseek's region_evalue)", "evalue": "E-value", "fident": "identical fraction", "its score or E-value": "its score or E-value (the member database's own figure)", "bits": "score (bits)", "alntmscore": "TM-score of the alignment",
    "Pfam": "Pfam accession", "member database": "member database", "signature": "signature", "InterPro entry": "InterPro entry", "GO terms": "GO terms" };
  function tipHtml(h) {
    let rows = `<tr><td>this protein</td><td>${h.q[0]}-${h.q[1]}</td></tr>`;
    if (h.r) rows += `<tr><td>${h.t ? esc(h.t) : h.k === "hmmscan" ? "Pfam model" : "matched protein"}</td><td>${h.r[0]}-${h.r[1]}</td></tr>`;
    const names = PAGE.names[h.k] || [];
    h.s.forEach((v, i) => { if (/^no .* column|^not in the file$/.test(String(v))) return; const vv = /e-?value/i.test(names[i] || "") && isFinite(parseFloat(v)) ? fmt2(parseFloat(v)) : v; rows += `<tr><td>${esc(GLOSS[names[i]] || names[i] || "")}</td><td>${esc(vv)}</td></tr>`; });
    if (h.id) rows += `<tr><td>identical</td><td>${h.id[0]} of ${h.id[1]} aligned columns (${Math.round((100 * h.id[0]) / h.id[1])}%)</td></tr>`;
    if (PAGE.plddt) { let s = 0, n = 0; for (let p = h.q[0]; p <= h.q[1]; p++) { const x = PAGE.plddt[p - 1]; if (x != null) { s += x; n++; } } if (n) rows += `<tr><td>pLDDT here</td><td>${Math.round(s / n)} (mean over ${h.q[0]}-${h.q[1]})</td></tr>`; }
    return `<div>${esc(PAGE.tracks[h.k].name)}</div><b>${esc(h.l)}</b>${about(h, esc)}<table>${rows}</table>` +
      "<div>Click to show or hide the residues.</div>";
  }
  // Full name, gene and species of the matched Swiss-Prot protein, and its Pfam families over
  // the matched part (for hmmscan and InterProScan: the family itself). txt turns text into HTML:
  // esc in the tooltip, which cannot be clicked; accHtml in the panel, so accessions are links.
  function about(h, txt) {
    let s = "";
    if (h.n) s += `<div class="full">${esc(h.n)}</div>`;
    const gs = [h.g ? "gene " + h.g : "", h.o || ""].filter(Boolean).join(" · ");
    if (gs) s += `<div>${esc(gs)}</div>`;
    if (h.k === "hmmscan" || h.k === "interproscan") {
      if (h.f && h.f.length) s += `<div>${h.k === "hmmscan" ? "Pfam family" : "Family"}: ${txt(h.f.join("; "))}</div>`;
    } else if (h.f === null || h.f === undefined) {
      s += "<div>Pfam families on the matched part: not known (no Pfam search output for this protein)</div>";
    } else {
      const shown = h.f.slice(0, 5).join("; ") + (h.f.length > 5 ? `; and ${h.f.length - 5} more` : "");
      const from = h.ds && SRC[h.ds] ? ` (${SRC[h.ds](h.t)})` : "";
      s += `<div>Pfam families on the matched part: ${h.f.length ? txt(shown) : "none"}${esc(from)}</div>`;
    }
    return s;
  }
  // the panel's own line: only what the row above does not show (the Pfam families of the matched part)
  function aboutShort(h) {
    if (h.k === "hmmscan" || h.k === "interproscan") return h.f && h.f.length ? `<div>${h.k === "hmmscan" ? "Pfam family" : "Family"}: ${accHtml(h.f.join("; "))}</div>` : "";
    if (h.f === null || h.f === undefined) return "<div>Pfam families on the matched part: not known (no Pfam search output for this protein)</div>";
    const shown = h.f.slice(0, 5).join("; ") + (h.f.length > 5 ? `; and ${h.f.length - 5} more` : "");
    return `<div>Pfam families on the matched part: ${h.f.length ? accHtml(shown) : "none"}</div>`;
  }
  function place(e) {
    const r = tip.getBoundingClientRect();
    let x = e.clientX + 14, y = e.clientY + 14;
    if (x + r.width > innerWidth - 8) x = Math.max(8, e.clientX - r.width - 14);
    if (y + r.height > innerHeight - 8) y = Math.max(8, e.clientY - r.height - 14);
    tip.style.left = x + "px"; tip.style.top = y + "px";
  }
  function wire(root) {
    root.querySelectorAll(".hit").forEach((g) => {
      const h = H[+g.dataset.i];
      g.addEventListener("mouseenter", (e) => { tip.innerHTML = tipHtml(h); tip.hidden = false; place(e); if (window.STRUCT) STRUCT.show(h.q[0], h.q[1]); });
      g.addEventListener("mousemove", place);
      g.addEventListener("mouseleave", () => { tip.hidden = true; if (window.STRUCT) STRUCT.show(null); });
      g.addEventListener("focus", () => {
        const b = g.getBoundingClientRect();
        tip.innerHTML = tipHtml(h); tip.hidden = false; place({ clientX: b.left, clientY: b.bottom });
      });
      g.addEventListener("blur", () => { tip.hidden = true; });
      g.addEventListener("click", () => { if (window.STRUCT) STRUCT.pin(h.q[0], h.q[1]); open(h); });
      g.addEventListener("keydown", (e) => { if (e.key === "Enter") open(h); });
    });
  }

  // "Show all": redraw each tool's row with all of its hits, with the same lanes, boxes and
  // labels build.py draws for the best few, then move the rows below it down.
  const LANE_H = 18, BOX_H = 14, ROW_PAD = 8;
  // The same box as build.py's shape_svg: the kind's colour, the program's shape.
  // E-value as opacity, as build.py's e_opacity: 0.25 at 10, full at 1e-30 or below
  const eOpacity = (e) => (e === null || e === undefined || e <= 0 ? 1 : Math.round(100 * Math.max(0.25, Math.min(1, 0.25 + (0.75 * (1 - Math.log10(e))) / 31))) / 100);
  function shapeSvg(x, y, w, hh, tr, op) {
    const f = (v) => v.toFixed(1), o = op < 1 ? ` opacity="${op.toFixed(2)}"` : "";
    if (tr.shape === "open") return `<rect x="${f(x + 1)}" y="${f(y + 1)}" width="${f(Math.max(1, w - 2))}" height="${f(hh - 2)}" rx="2" fill="#fff" stroke="${tr.color}" stroke-width="2"${o}/>`;
    if (tr.shape === "hatched") return `<rect x="${f(x + 0.5)}" y="${f(y + 0.5)}" width="${f(Math.max(1, w - 1))}" height="${f(hh - 1)}" rx="2" fill="url(#hatch-${tr.cat})" stroke="${tr.color}" stroke-width="1"${o}/>`;
    if (tr.shape === "half") return `<rect x="${f(x + 0.5)}" y="${f(y + 0.5)}" width="${f(Math.max(1, w - 1))}" height="${f(hh - 1)}" rx="2" fill="#fff" stroke="${tr.color}"${o}/><rect x="${f(x)}" y="${f(y + hh / 2)}" width="${f(w)}" height="${f(hh / 2)}" fill="${tr.color}"${o}/>`;
    return `<rect x="${f(x)}" y="${f(y)}" width="${f(w)}" height="${f(hh)}" rx="${f(tr.shape === "rounded" ? hh / 2 : 2)}" fill="${tr.color}"${o}/>`;
  }
  // Label placing, as build.py's place_labels: inside the box when it fits, else beside it
  // (right, or left at the drawing's right edge) with a leader line; never cut or shrunk.
  const textW = (t) => t.length * 11 * 0.58, GAP = 8;
  function placeLabels(list) {
    const lanes = [], out = [];
    for (const h of list) {
      const x0 = (h.q[0] - 1) * PAGE.px, x1 = h.q[1] * PAGE.px, w = Math.max(2, x1 - x0), lw = textW(h.l);
      let side, lx, a, b;
      if (lw + 6 <= w) { side = "in"; lx = x0 + 3; a = x0; b = x0 + w; }
      else if (x1 + GAP + lw <= PAGE.w || x0 - GAP - lw < 0) { side = "right"; lx = Math.min(x1 + GAP, PAGE.w - lw); a = x0; b = Math.max(x0 + w, lx + lw); }
      else { side = "left"; lx = x0 - GAP - lw; a = lx; b = x0 + w; }
      let li = lanes.findIndex((ends) => ends.every(([e0, e1]) => a > e1 + 4 || b < e0 - 4));
      if (li < 0) { lanes.push([]); li = lanes.length - 1; }
      lanes[li].push([a, b]);
      out.push([li, side, lx]);
    }
    return [out, Math.max(1, lanes.length)];
  }
  // An InterProScan feature prediction (h.x), as build.py's feature_svg: a dashed box, never counted
  const featureSvg = (x, y, w, hh, col) => `<rect x="${(x + 1).toFixed(1)}" y="${(y + 1).toFixed(1)}" width="${Math.max(1, w - 2).toFixed(1)}" height="${(hh - 2).toFixed(1)}" rx="2" fill="transparent" stroke="${col}" stroke-width="1.5" stroke-dasharray="3 2"/>`;
  function boxSvg(h, i, lane, side, lx, dim) {
    const tr = PAGE.tracks[h.k];
    const x0 = (h.q[0] - 1) * PAGE.px, w = Math.max(2, h.q[1] * PAGE.px - x0), y = ROW_PAD + lane * LANE_H;
    const op = h.x ? 1 : eOpacity(h.e);
    const rect = h.x ? featureSvg(x0, y, w, BOX_H, tr.color) : shapeSvg(x0, y, w, BOX_H, tr, op), ty = y + BOX_H - 3;
    let lab;
    if (side === "in") {
      const solid = (tr.shape === "filled" || tr.shape === "rounded") && op >= 0.6 && !h.x;
      lab = `<text x="${(x0 + 3).toFixed(1)}" y="${ty}" fill="${solid ? "#ffffff" : "#1a1a1a"}" class="bl${solid ? "" : " halo"}">${esc(h.l)}</text>`;
    } else {
      const lw = textW(h.l), [l0, l1] = side === "right" ? [x0 + w, lx - 1] : [lx + lw + 1, x0];
      const lead = l1 > l0 ? `<line x1="${l0.toFixed(1)}" y1="${y + BOX_H / 2}" x2="${l1.toFixed(1)}" y2="${y + BOX_H / 2}" stroke="${tr.color}" stroke-width="1"/>` : "";
      lab = lead + `<text x="${lx.toFixed(1)}" y="${ty}" fill="#1a1a1a" class="bl out">${esc(h.l)}</text>`;
    }
    return `<g class="hit${dim ? " dim" : ""}" tabindex="0" data-i="${i}">${rect}${lab}</g>`;
  }
  // Search: each hit's text is what its tooltip shows (entry name, accession, full name, gene,
  // species, Pfam families of the matched part and of the target, GO terms when the input has
  // them). Lower case; every word typed must appear in it, as on the index page of the full site.
  const SEARCHED_SCORES = /GO|Pfam|signature|member database|InterPro/i;
  const hay = H.map((h) => {
    const names = PAGE.names[h.k] || [];
    const parts = [h.l, h.t, h.n, h.g, h.o].concat(h.f || [], (h.td || []).map((d) => d[2]));
    h.s.forEach((v, i) => { if (SEARCHED_SCORES.test(names[i] || "") && !/^(no |none$|not in )/.test(String(v))) parts.push(v); });
    return parts.filter(Boolean).join(" ").toLowerCase();
  });
  let words = [];
  const matches = (i) => words.every((w) => hay[i].includes(w));
  const fmtE = (v) => (v >= -3 ? String(+Math.pow(10, v).toFixed(3)) : "1e" + v);
  const num = (v) => { const x = parseFloat(v); return isNaN(x) ? null : x; };
  function scoreOf(h) {  // [value for sorting, text]: bits, or hhsearch's Prob, or InterProScan's score
    const k = h.k;
    if (k === "mmseqs2" || k === "mmseqs2_3rounds" || k === "prostt5" || k === "foldseek") return [num(h.s[2]), h.s[2] + " bits"];
    if (k === "phmmer" || k === "jackhmmer" || k === "hmmscan") return [num(h.s[1]), h.s[1] + " bits"];
    if (k === "hhsearch") return [num(h.s[0]), "Prob " + h.s[0]];
    return [null, ""];
  }
  function identOf(h) {
    const k = h.k;
    if (k === "mmseqs2" || k === "mmseqs2_3rounds" || k === "prostt5" || k === "foldseek") { const f = num(h.s[1]); return f === null ? null : 100 * f; }
    return null;
  }
  function k_e(h) {
    const names = PAGE.names[h.k] || [];
    const n = names.findIndex((x) => /e-?value/i.test(x));
    return n < 0 ? -1 : n;
  }
  // Swiss-Prot sequences live in ../targets/NNN.js, one file per group of accessions.
  function shardOf(acc) {
    let h = 0;
    for (const ch of acc) h = (Math.imul(h, 31) + ch.charCodeAt(0)) >>> 0;
    return h % PAGE.tseq_files;
  }
  const loading = {};
  function targetSeq(acc) {
    if (PAGE.noseq) return Promise.resolve(null);  // protein.html: no sequences on the site for this protein
    window.TSEQ = window.TSEQ || {};
    if (acc in TSEQ) return Promise.resolve(TSEQ[acc]);
    const n = shardOf(acc);
    if (!loading[n]) loading[n] = new Promise((res) => {
      const s = document.createElement("script");
      s.src = `${PAGE.tbase ?? "../"}targets/${String(n).padStart(3, "0")}.js`;  // protein.html sits at the site's root: tbase ""
      s.onload = res; s.onerror = res;
      document.head.appendChild(s);
    });
    return loading[n].then(() => TSEQ[acc]);
  }

  // Local alignment, BLOSUM62, gap open -11 (first gap residue), extend -1, the same scores
  // build.py gives Biopython's PairwiseAligner. Used only for hits beyond the best 5 per tool.
  const AA = "ARNDCQEGHILKMFPSTWYVBZX*";
  const B62 = `4 -1 -2 -2 0 -1 -1 0 -2 -1 -1 -1 -1 -2 -1 1 0 -3 -2 0 -2 -1 0 -4
-1 5 0 -2 -3 1 0 -2 0 -3 -2 2 -1 -3 -2 -1 -1 -3 -2 -3 -1 0 -1 -4
-2 0 6 1 -3 0 0 0 1 -3 -3 0 -2 -3 -2 1 0 -4 -2 -3 3 0 -1 -4
-2 -2 1 6 -3 0 2 -1 -1 -3 -4 -1 -3 -3 -1 0 -1 -4 -3 -3 4 1 -1 -4
0 -3 -3 -3 9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
-1 1 0 0 -3 5 2 -2 0 -3 -2 1 0 -3 -1 0 -1 -2 -1 -2 0 3 -1 -4
-1 0 0 2 -4 2 5 -2 0 -3 -3 1 -2 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4
0 -2 0 -1 -3 -2 -2 6 -2 -4 -4 -2 -3 -3 -2 0 -2 -2 -3 -3 -1 -2 -1 -4
-2 0 1 -1 -3 0 0 -2 8 -3 -3 -1 -2 -1 -2 -1 -2 -2 2 -3 0 0 -1 -4
-1 -3 -3 -3 -1 -3 -3 -4 -3 4 2 -3 1 0 -3 -2 -1 -3 -1 3 -3 -3 -1 -4
-1 -2 -3 -4 -1 -2 -3 -4 -3 2 4 -2 2 0 -3 -2 -1 -2 -1 1 -4 -3 -1 -4
-1 2 0 -1 -3 1 1 -2 -1 -3 -2 5 -1 -3 -1 0 -1 -3 -2 -2 0 1 -1 -4
-1 -1 -2 -3 -1 0 -2 -3 -2 1 2 -1 5 0 -2 -1 -1 -1 -1 1 -3 -1 -1 -4
-2 -3 -3 -3 -2 -3 -3 -3 -1 0 0 -3 0 6 -4 -2 -2 1 3 -1 -3 -3 -1 -4
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4 7 -1 -1 -4 -3 -2 -2 -1 -2 -4
1 -1 1 0 -1 0 0 0 -1 -2 -2 0 -1 -2 -1 4 1 -3 -2 -2 0 0 0 -4
0 -1 0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1 1 5 -2 -2 0 -1 -1 0 -4
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1 1 -4 -3 -2 11 2 -3 -4 -3 -2 -4
-2 -2 -2 -3 -2 -1 -2 -3 2 -1 -1 -2 -1 3 -3 -2 -2 2 7 -1 -3 -2 -1 -4
0 -3 -3 -3 -1 -2 -2 -3 -3 3 1 -2 1 -1 -2 -2 0 -3 -1 4 -3 -2 -1 -4
-2 -1 3 4 -3 0 1 -1 0 -3 -4 0 -3 -3 -2 0 -1 -4 -3 -3 4 1 -1 -4
-1 0 0 1 -3 3 4 -2 0 -3 -3 1 -1 -3 -1 0 -1 -3 -2 -2 1 4 -1 -4
0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2 0 0 -2 -1 -1 -1 -1 -1 -4
-4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 1`.split("\n").map((r) => r.trim().split(/\s+/).map(Number));
  const idx = (c) => { const i = AA.indexOf(c); return i < 0 ? AA.indexOf("X") : i; };
  function localAlign(a, b) {
    const n = a.length, m = b.length, OPEN = -11, EXT = -1, NEG = -1e9;
    const W = m + 1;
    const H = new Float64Array((n + 1) * W), E = new Float64Array((n + 1) * W).fill(NEG), F = new Float64Array((n + 1) * W).fill(NEG);
    const ai = Array.from(a, idx), bi = Array.from(b, idx);
    let best = 0, bi_ = 0, bj = 0;
    for (let i = 1; i <= n; i++) for (let j = 1; j <= m; j++) {
      const k = i * W + j;
      E[k] = Math.max(H[k - 1] + OPEN, E[k - 1] + EXT);
      F[k] = Math.max(H[k - W] + OPEN, F[k - W] + EXT);
      const h = Math.max(0, H[k - W - 1] + B62[ai[i - 1]][bi[j - 1]], E[k], F[k]);
      H[k] = h;
      if (h > best) { best = h; bi_ = i; bj = j; }
    }
    let i = bi_, j = bj, qa = [], ta = [], state = "H";
    while (i > 0 && j > 0) {
      const k = i * W + j;
      if (state === "H") {
        if (H[k] === 0) break;
        if (H[k] === H[k - W - 1] + B62[ai[i - 1]][bi[j - 1]]) { qa.push(a[i - 1]); ta.push(b[j - 1]); i--; j--; continue; }
        state = H[k] === F[k] ? "F" : "E";
      }
      if (state === "F") { qa.push(a[i - 1]); ta.push("-"); state = F[k] === F[k - W] + EXT ? "F" : "H"; i--; continue; }
      if (state === "E") { qa.push("-"); ta.push(b[j - 1]); state = E[k] === E[k - 1] + EXT ? "E" : "H"; j--; continue; }
    }
    return { q: qa.reverse().join(""), t: ta.reverse().join(""), qs: i + 1, ts: j + 1 };
  }

  // ---------------------------------------------------------------- residue panel ----
  // The tools' alignments are in <page>.aln.js next to the page, loaded on the first click.
  let alnLoad = null;
  function alignments() {
    if (window.ALN) return Promise.resolve(window.ALN);  // the standalone page loads aln.js up front
    if (!alnLoad) alnLoad = new Promise((res) => {
      const s = document.createElement("script");
      s.src = location.pathname.split("/").pop().replace(/(\.dense)?\.html$/, ".aln.js");
      s.onload = () => res(window.ALN || {}); s.onerror = () => res({});
      document.head.appendChild(s);
    });
    return alnLoad;
  }
  // Residue classes, only for a kmerseek hit and only in that hit's own alphabet (h.m), the one
  // kmerseek compared the two proteins in. Other programs compare residues, so their panels
  // show no classes. Fill colours mark hydrophobic and polar only in the two-class hp_ alphabets.
  const CLASS_WORD = { h: "hydrophobic", p: "polar" };
  let ALPHA_NAME = null, ALPHA = null, HP = false, groups = {};
  function setAlphabet(name) {
    ALPHA_NAME = name && PAGE.alpha && PAGE.alpha[name] ? name : null;
    ALPHA = ALPHA_NAME ? PAGE.alpha[ALPHA_NAME] : null;
    HP = !!ALPHA_NAME && /^hp_/.test(ALPHA_NAME);
    groups = {};
    if (ALPHA) for (const [r, c] of Object.entries(ALPHA)) (groups[c] = groups[c] || []).push(r);
  }
  const classOf = (c) => (ALPHA && c && c !== "-" ? ALPHA[c.toUpperCase()] || null : null);
  const cssClass = (c) => {
    if (!ALPHA) return "n";
    const k = classOf(c);
    return HP && (k === "h" || k === "p") ? k : (k ? "o" : "x");
  };
  const similar = (a, b) => a !== "-" && b !== "-" && B62[idx(a.toUpperCase())][idx(b.toUpperCase())] > 0;

  // Counts over the shown alignment: identical columns, more similar (BLOSUM62 above 0; for
  // kmerseek, the same class in its alphabet), and columns with a gap.
  function stats(a, kmer) {
    let ident = 0, sim = 0, gaps = 0;
    for (let i = 0; i < a.q.length; i++) {
      const x = a.q[i], y = a.t[i] || "-";
      if (x === "-" || y === "-") { gaps++; continue; }
      if (x.toUpperCase() === y.toUpperCase()) ident++;
      else if (kmer ? classOf(x) && classOf(x) === classOf(y) : similar(x, y)) sim++;
    }
    return { cols: a.q.length, ident, sim, gaps };
  }
  const pct = (a, b) => (b ? Math.round((100 * a) / b) : 0);

  // One line of W columns: [seed bar], your protein, [its classes], a match line, [target
  // classes], the target, [seed bar]. inSeed(i) says whether column i is in a kmerseek seed.
  function chunkHtml(a, start, n, qpos, tpos, tname, kmer, flank, W, inSeed) {
    const rowsQ = [], rowsT = [], mid = [], cq = [], ct = [], sd = [];
    let qi = qpos, ti = tpos;
    const qFirst = qi, tFirst = ti;
    for (let i = start; i < Math.min(start + W, n); i++) {
      const x = a.q[i], y = a.t[i] || "-";
      const fl = flank && (i < flank[0] || i >= n - flank[1]);
      const same = x !== "-" && y !== "-" && x.toUpperCase() === y.toUpperCase();
      const cell = (c, pos, who) => c === "-" ? `<i class="gap">-</i>` :
        `<i class="${cssClass(c)}${same && !fl ? " same" : ""}${fl ? " flank" : ""}" data-c="${c}" data-p="${pos}" data-w="${esc(who)}">${c}</i>`;
      rowsQ.push(cell(x, x === "-" ? 0 : qi, "your protein"));
      rowsT.push(cell(y, y === "-" ? 0 : ti, tname));
      if (inSeed) sd.push(inSeed(i) ? `<i class="sd"></i>` : "<i></i>");
      const kx = classOf(x), ky = classOf(y);
      if (kmer) {
        cq.push(`<i class="${fl ? "flank" : ""}">${kx || (x === "-" ? "-" : "x")}</i>`);
        ct.push(`<i class="${fl ? "flank" : ""}">${ky || (y === "-" ? "-" : "x")}</i>`);
        mid.push(fl ? " " : same ? "|" : kx && kx === ky ? "+" : " ");
      } else {
        mid.push(same ? "|" : similar(x, y) ? "+" : " ");
      }
      if (x !== "-") qi++;
      if (y !== "-") ti++;
    }
    const line = (lab, body, end, cls) => `<div class="lab">${lab}</div><div class="cells ${cls || ""}">${body}</div><div class="lab r">${end}</div>`;
    const an = esc(ALPHA_NAME || "classes");
    const seedRow = inSeed && sd.some((c) => c.length > 7) ? line("seed", sd.join(""), "", "seedrow") : "";
    let h = seedRow + line(`your protein aa ${qFirst}`, rowsQ.join(""), qi - 1);
    if (kmer) h += line(`your protein ${an}`, cq.join(""), "", "enc");
    h += line("", mid.map((m) => `<i>${m === " " ? "&nbsp;" : m}</i>`).join(""), "", "marks");
    if (kmer) h += line(`${esc(tname)} ${an}`, ct.join(""), "", "enc");
    h += line(`${esc(tname)} aa ${tFirst}`, rowsT.join(""), ti - 1, "under");
    h += seedRow;
    return { html: `<div class="chunk">${h}</div>`, qi, ti };
  }

  // The protein's name for an axis: its ID once (FUN008084 for FUN008084_FUN008084), and its
  // name from protein_names.tsv up to the first ";" or "(" when there is one.
  function shortId(id) {
    const m = /^([^_]+)_(.*)$/.exec(id);
    return m && m[1] === m[2].replace(/\.\d+$/, "") ? m[1] : id;
  }
  // Text cut to fit w pixels at the map's 11 px font, with "..." at the end when cut.
  const MAP_CH = 5.6;
  const fit = (t, w) => (t.length * MAP_CH <= w ? t : t.slice(0, Math.max(1, Math.floor(w / MAP_CH) - 1)) + "…");
  // Greedy: best first, keep each item that overlaps none kept.
  function apart(items) {
    const kept = [];
    for (const it of items) if (kept.every((k) => it.a > k.b || it.b < k.a)) kept.push(it);
    return kept;
  }
  const SRC = {
    interproscan: (t) => `from InterProScan run on ${t}`,
    hmmscan: (t) => `from hmmscan with Pfam run on ${t}`,
    interpro: (t) => `from InterPro's Pfam matches for ${t}`,
  };

  // A small map: your protein across, the other sequence down, every hit on that sequence as a
  // line in its tool's colour, the one shown drawn thick. Along the top: what hmmscan and
  // InterProScan found on your protein, in the marks of the drawing above. Along the right: the
  // Pfam domains of the other protein. Each axis is scaled to its own length: x fills PW pixels,
  // y is PW times the length ratio, kept between Y_MIN and Y_MAX, so a short and a long protein
  // both get a readable plot.
  const PW = 220, Y_MIN = 90, Y_MAX = 170;
  function dotPlot(h, tlen, tname) {
    const qL = PAGE.q.length, sl = document.getElementById("emax");
    const maxE = sl ? Math.pow(10, +sl.value) : Infinity;
    const shown = (x) => x.e === null || x.e === undefined || x.e <= maxE;
    const ph = Math.round(Math.max(Y_MIN, Math.min(Y_MAX, (PW * tlen) / qL)));
    const scx = PW / qL, scy = ph / tlen;
    // top edge: one lane of family matches (hmmscan, InterProScan), one lane of predictions
    const fam = apart(H.filter((x) => (x.k === "hmmscan" || x.k === "interproscan") && !x.x && shown(x))
      .sort((a, b) => (a.e ?? 0) - (b.e ?? 0)).map((x) => ({ a: x.q[0], b: x.q[1], l: x.l, x, e: x.e })));
    const pred = [];  // predictions of one kind that overlap are drawn as one box
    H.filter((x) => x.x).sort((a, b) => a.q[0] - b.q[0]).forEach((x) => {
      const p = pred.find((p) => p.l === x.l && x.q[0] <= p.b + 1 && x.q[1] >= p.a - 1);
      if (p) { p.a = Math.min(p.a, x.q[0]); p.b = Math.max(p.b, x.q[1]); } else pred.push({ a: x.q[0], b: x.q[1], l: x.l, x });
    });
    const LH = 16, BH = 12, lanes = [fam, pred].filter((l) => l.length);
    const ml = 40, mt = 8 + Math.max(1, lanes.length) * LH + 6, mr = 112;
    const W = ml + PW + mr, Hh = mt + ph + 40;
    const sx = (p) => ml + (p - 1) * scx, sy = (p) => mt + (p - 1) * scy;
    let g = "";
    lanes.forEach((lane, li) => {
      const y = 8 + li * LH;
      g += `<line x1="${ml}" x2="${ml + PW}" y1="${y + BH / 2}" y2="${y + BH / 2}" stroke="#C8C8C8" stroke-width="1.5"/>`;
      for (const it of lane) {
        const x0 = sx(it.a), w = Math.max(2, (it.b - it.a + 1) * scx), tr = PAGE.tracks[it.x.k];
        const rect = it.x.x ? featureSvg(x0, y, w, BH, tr.color) : shapeSvg(x0, y, w, BH, tr, eOpacity(it.e));
        const what = it.x.x ? `${it.l} (InterProScan prediction)` : `${it.l} (${tr.name}${it.e !== null && it.e !== undefined ? ", E-value " + it.e : ""})`;
        const solid = !it.x.x && (tr.shape === "filled" || tr.shape === "rounded") && eOpacity(it.e) >= 0.6;
        const lab = it.l.length * MAP_CH + 4 <= w ? `<text x="${(x0 + 2).toFixed(1)}" y="${y + BH - 2}" class="dl${solid ? " on" : ""}">${esc(it.l)}</text>` : "";
        g += `<g><title>${esc(what)}, residues ${it.a}-${it.b}</title>${rect}${lab}</g>`;
      }
    });
    // right edge: the other protein's Pfam domains, as open boxes in the Pfam colour
    const tx = ml + PW + 8;
    g += `<line x1="${tx + 5}" x2="${tx + 5}" y1="${mt}" y2="${mt + ph}" stroke="#C8C8C8" stroke-width="1.5"/>`;
    for (const [a, b, l] of h.td || []) {
      const y0 = sy(a), hh = Math.max(2, (b - a + 1) * scy);
      g += `<g><title>${esc(l)}, residues ${a}-${b} of ${esc(tname)}</title><rect x="${tx}" y="${y0.toFixed(1)}" width="10" height="${hh.toFixed(1)}" rx="1.5" fill="#fff" stroke="${PAGE.pfam_color}" stroke-width="1.5"/>`;
      if (hh >= 11) g += `<text x="${tx + 14}" y="${(y0 + Math.min(hh, 22) / 2 + 4).toFixed(1)}" class="dl">${esc(fit(l.replace(/ \(PF\d+\)$/, ""), mr - 18))}</text>`;
      g += "</g>";
    }
    // the plot: a frame, every hit on this protein as a line, the shown one thick
    // the two axes are the two proteins, drawn as protein lines: yours along the bottom, the other on the left
    g += `<line x1="${ml}" x2="${ml + PW}" y1="${mt + ph}" y2="${mt + ph}" stroke="#C8C8C8" stroke-width="1.5"/>` +
      `<line x1="${ml}" x2="${ml}" y1="${mt}" y2="${mt + ph}" stroke="#C8C8C8" stroke-width="1.5"/>`;
    const same = H.filter((x) => (h.t ? x.t === h.t : x.k === h.k && x.l === h.l) && x.r);
    const seg = (x, w, op) => `<line x1="${sx(x.q[0]).toFixed(1)}" y1="${sy(Math.min(x.r[0], x.r[1])).toFixed(1)}" x2="${sx(x.q[1]).toFixed(1)}" y2="${sy(Math.max(x.r[0], x.r[1])).toFixed(1)}" stroke="${PAGE.tracks[x.k].color}" stroke-width="${w}" stroke-linecap="round" opacity="${op}"/>`;
    for (const x of same) if (x !== h) g += seg(x, 1.5, 0.55);
    g += seg(h, 4, 1);
    g += `<text x="${ml}" y="${mt + ph + 12}" class="ax">1</text><text x="${ml + PW}" y="${mt + ph + 12}" class="ax" text-anchor="end">${qL}</text>`;
    g += `<text x="${ml - 4}" y="${mt + 9}" class="ax" text-anchor="end">1</text><text x="${ml - 4}" y="${mt + ph}" class="ax" text-anchor="end">${tlen}</text>`;
    // the protein's short name: a short word in brackets in its name (BHF), else the name up to ";" or "("
    const ab = /\(([^()]{1,12})\)/.exec(PAGE.pn || ""), pn = ab ? ab[1] : (PAGE.pn || "").split(/[;(]/)[0].trim();
    const xl = `your protein, ${shortId(PAGE.id)}${pn ? " (" + pn + ")" : ""}, aa`;
    g += `<text x="${ml}" y="${mt + ph + 28}" class="ax"><title>${esc(xl)}</title>${esc(fit(xl, PW))}</text>`;
    const yl = `${tname}, aa`, yc = mt + ph / 2;
    g += `<text x="12" y="${yc}" class="ax" text-anchor="middle" transform="rotate(-90 12 ${yc})"><title>${esc(yl)}</title>${esc(fit(yl, ph))}</text>`;
    // the key: what each mark is, in plain words, read before the plot
    const kinds = {};
    for (const x of same) { const tr = PAGE.tracks[x.k]; (kinds[tr.color] = kinds[tr.color] || new Set()).add(tr.name); }
    const tools = Object.entries(kinds).map(([c, ns]) => `<span><i style="background:${c}"></i>${esc([...ns].join(", "))}</span>`).join("");
    const top = [];
    if (fam.length) top.push(`solid boxes: Pfam families and other InterPro signatures that hmmscan or InterProScan found on your protein${sl ? ` (E-value at most ${fmtE(+sl.value)}, the chosen cutoff)` : ""}; where two overlap, the one with the better E-value`);
    if (pred.length) top.push(`dashed boxes: ${[...new Set(pred.map((p) => p.l))].join(", ")}, predicted by InterProScan`);
    const topLine = top.length ? `Along the top, ${top.join("; ")}.` : "Along the top: nothing; hmmscan and InterProScan found no family or prediction on your protein.";
    let rightLine;
    if (h.k === "hmmscan") rightLine = "The y axis is the Pfam model itself, so nothing is drawn on its edge.";
    else if (h.td && h.td.length) rightLine = `Along the right, open boxes: Pfam domains of ${uniHtml(h.t)}, ${esc((SRC[h.ds] || (() => "source not recorded"))(h.t))}; where two overlap, only the first is drawn.`;
    else if (h.ds) rightLine = `Along the right: nothing; no Pfam domain is known on ${uniHtml(h.t)} (${esc(SRC[h.ds](h.t))}).`;
    else rightLine = `Along the right: nothing; the Pfam domains of ${uniHtml(h.t)} are not known (no Pfam search output for it).`;
    return `<div class="map" style="max-width:${W}px"><p class="mapkey">Across: your protein; down: ${esc(tname)}. Each line is one hit between the two (${same.length} here), the thick one is this hit. Top edge: hmmscan and InterProScan on your protein; right edge: the other protein's Pfam domains.</p>` +
      `<div class="mapkey tools">${tools}</div><details class="mapmore"><summary>more about the edges</summary><p class="mapkey">${esc(topLine)}</p><p class="mapkey">${rightLine}</p></details>` +
      `<svg viewBox="0 0 ${W} ${Hh}" width="${W}" height="${Hh}" role="img" aria-label="Map of the hits between your protein and ${esc(tname)}">${g}</svg></div>`;
  }

  function open(h, box) {
    if (PAGE.noseq) {  // protein.html: the hit's numbers only; its residues are on the 51 drawn pages
      box.innerHTML = `<h3>${esc(PAGE.tracks[h.k].name)}: ${esc(h.l)}</h3><p class="note">${(PAGE.names[h.k] || []).map((n, j) => `${esc(GLOSS[n] || n)} ${esc(h.s[j] ?? "")}`).join(" · ")}. ` +
        "This protein's sequence is not on the site, so its residues cannot be shown; the 51 proteins with a drawn page show them.</p>";
      return;
    }
    Promise.all([alignments(), h.t ? targetSeq(h.t) : Promise.resolve(null)]).then(([alns, tseq]) => show(h, alns[H.indexOf(h)], tseq, box));
  }

  function show(h, a, tseq, box) {
    const q = PAGE.q, tool = PAGE.tracks[h.k].name;
    setAlphabet(h.k === "kmerseek" ? h.m : null);
    const tname = h.k === "hmmscan" ? `${h.l} model` : h.l;
    const names = PAGE.names[h.k] || [];
    let html = `<h3>${esc(tool)}: ${esc(h.l)}` + (h.t && h.t !== h.l ? ` <span class="w">Swiss-Prot ${uniHtml(h.t)}</span>` : "") + "</h3>" + aboutShort(h);
    const fmtS = (n, v) => (/e-?value/i.test(n || "") && isFinite(parseFloat(v)) ? String(+parseFloat(v).toPrecision(3)) : v);
    const scoreBits = h.s.map((v, i) => [names[i], v]).filter(([, v]) => v !== "not in the file" && !/^no .* column/.test(v)).map(([n, v]) => `${esc(GLOSS[n] || n || "")} <b>${accHtml(fmtS(n, v))}</b>`).join(" · ");
    let src = "", kmer = false, flank = null, note = "";
    if (h.k === "kmerseek" && tseq) {
      kmer = true;
      const F = 10, q0 = Math.max(1, h.q[0] - F), t0 = Math.max(1, h.r[0] - F);
      const lf = Math.min(h.q[0] - q0, h.r[0] - t0);
      const rf = Math.min(q.length - h.q[1], tseq.length - h.r[1], F);
      a = { q: q.slice(h.q[0] - 1 - lf, h.q[1] + rf), t: tseq.slice(h.r[0] - 1 - lf, h.r[1] + rf), qs: h.q[0] - lf, ts: h.r[0] - lf };
      flank = [lf, rf];
      src = `kmerseek matches the two segments side by side without gaps; ${lf} residues before and ${rf} after the match are shown faded.`;
    } else if (h.k === "kmerseek" && h.tsub) {  // protein.html: the Swiss-Prot residues of the region, kept with the hit
      kmer = true; flank = [0, 0];
      a = { q: q.slice(h.q[0] - 1, h.q[1]), t: h.tsub, qs: h.q[0], ts: h.r[0] };
      src = "kmerseek matches the two segments side by side without gaps.";
    } else if (!a && h.r && tseq && h.k !== "hhsearch") {
      const lo = Math.min(h.r[0], h.r[1]), hi = Math.max(h.r[0], h.r[1]);
      a = localAlign(q.slice(h.q[0] - 1, h.q[1]), tseq.slice(lo - 1, hi));
      if (a && a.q) { a.qs += h.q[0] - 1; a.ts += lo - 1; }
      src = "recomputed in your browser, because the tool's own alignment was not given: the two ranges aligned with BLOSUM62, gap open -11, gap extend -1";
    } else if (a) {
      src = a.src || "";
    }
    const tlen = h.tl || (tseq ? tseq.length : null);
    let statLine = `<p class="stats">${scoreBits}`;
    statLine += ` · your protein ${h.q[0]}-${h.q[1]} of ${q.length} aa (${pct(h.q[1] - h.q[0] + 1, q.length)}%)`;
    if (h.r) statLine += ` · ${esc(tname)} ${Math.min(...h.r)}-${Math.max(...h.r)}` + (tlen ? ` of ${tlen} aa (${pct(Math.abs(h.r[1] - h.r[0]) + 1, tlen)}%)` : "");
    html += statLine + "</p>";
    if (h.k !== "hhsearch" && h.t && !tseq && !(a && a.q)) html += `<p class="note">${uniHtml(h.t)} is not in swissprot.fasta, so its residues cannot be shown.</p>`;
    let right = "";
    if (a && a.q) {
      const leg = [], kk = h.kk ? `, k = ${h.kk}` : "";
      if (kmer && ALPHA) {
        leg.push(`<span><b>${esc(ALPHA_NAME)}${kk}</b>: the alphabet and k-mer size kmerseek compared the two proteins in</span>`);
        if (HP) for (const c of Object.keys(groups).sort()) leg.push(`<span><i class="${c === "h" || c === "p" ? c : "o"}">${c}</i>${esc(CLASS_WORD[c] || "class " + c)}: <code>${groups[c].sort().join("")}</code></span>`);
        else if (ALPHA_NAME === "protein20") leg.push("<span>protein20 does not group residues: each class is one amino acid</span>");
        else leg.push(`<span>classes:${Object.keys(groups).sort().map((c) => `<code>${groups[c].sort().join("")}</code> → <b class="mono">${esc(c)}</b>`).join(", ")}</span>`);
      } else if (kmer) {
        leg.push(`<span>The classes of ${esc(h.m || "this hit's alphabet")}${kk} are not in kmerseek_alphabets.rs, so no classes are shown.</span>`);
      }
      // the seed: columns where k residues in a row are in the same class in both (build.py seed_runs)
      const seeds = kmer && Array.isArray(h.sd) ? h.sd : null;
      const inSeed = seeds && seeds.length ? (i) => seeds.some(([s0, s1]) => a.qs + i >= s0 && a.qs + i <= s1) : null;
      if (inSeed) leg.push(`<span><i class="sd"></i>seed (see the Seed line above)</span>`);
      leg.push(`<span><b class="mono">A</b> bold: the same residue in both</span>`);
      leg.push(kmer ? `<span><b class="mono">|</b> same residue, <b class="mono">+</b> same class</span>` : `<span><b class="mono">|</b> same residue, <b class="mono">+</b> similar (BLOSUM62 above 0), <b class="mono">-</b> gap</span>`);
      const st = stats(flank ? { q: a.q.slice(flank[0], a.q.length - flank[1]), t: a.t.slice(flank[0], a.t.length - flank[1]) } : a, kmer);
      right += `<p class="asum"><b>${st.ident} of ${st.cols}</b> aligned columns identical (${pct(st.ident, st.cols)}%), ` +
        `${st.sim} more ${kmer ? `in the same ${esc(ALPHA_NAME)} class` : "similar (BLOSUM62 score above 0)"}, ` +
        `${st.gaps} gap column${st.gaps === 1 ? "" : "s"}${flank ? " (the faded flanks not counted)" : ""}</p>`;
      if (seeds) {
        const dt = h.r[0] - h.q[0];
        right += seeds.length ? `<p class="asum">Seed: ${h.sn} overlapping window${h.sn === 1 ? "" : "s"} of ${h.kk} residues in the same ${esc(ALPHA_NAME)} class in both, at ` +
          seeds.map(([s0, s1]) => `your protein aa ${s0}-${s1}, ${esc(tname)} aa ${s0 + dt}-${s1 + dt}`).join("; ") + "</p>" :
          `<p class="asum">No seed: nowhere in this hit are ${h.kk} residues in a row in the same ${esc(ALPHA_NAME)} class in both. kmerseek starts from such a run, so its positions or its alphabet here do not match what kmerseek used.</p>`;
      }
      right += `<p class="asrc">Source: ${esc(src.replace(/\.$/, ""))}.</p>`;
      right += `<div class="alegend">${leg.join("")}</div>`;
      let qi = a.qs, ti = a.ts, body = "";
      // columns per line: as many tens as fit beside the map (15 px a column, about 190 px of labels), 30 to 60
      const linked = H.filter((x) => (h.t ? x.t === h.t : x.k === h.k && x.l === h.l) && x.r).length;
      const avail = (box.clientWidth || innerWidth) - 60 - (tlen && h.r && linked > 1 && innerWidth > 760 ? 400 : 0);
      const CW = Math.max(30, Math.min(60, 10 * Math.floor((avail - 190) / 150)));
      for (let st = 0; st < a.q.length; st += CW) {
        const c = chunkHtml(a, st, a.q.length, qi, ti, tname, kmer, flank, CW, inSeed);
        body += c.html; qi = c.qi; ti = c.ti;
      }
      right += `<div class="seqs">${body}</div>`;
      if (!kmer && h.k !== "hhsearch" && h.k !== "hmmscan") {
        const qe = a.qs + a.q.replace(/-/g, "").length - 1;
        if (a.qs > h.q[0] || qe < h.q[1]) right += `<p class="note">This alignment covers residues ${a.qs}-${qe} of your protein; the tool reported ${h.q[0]}-${h.q[1]}.</p>`;
      }
    } else if (h.k === "hmmscan") {
      right += `<p class="note">hmmscan's alignment of ${esc(h.l)} to your protein was not given (hmmscan.out is not in the input folder). Your protein's residues ${h.q[0]}-${h.q[1]}:</p><div class="seqs"><div class="chunk"><div class="lab">your protein aa ${h.q[0]}</div><div class="cells">${Array.from(q.slice(h.q[0] - 1, h.q[1]), (c, i) => `<i class="${cssClass(c)}" data-c="${c}" data-p="${h.q[0] + i}" data-w="your protein">${c}</i>`).join("")}</div><div class="lab r">${h.q[1]}</div></div></div>`;
    } else {
      if (h.k !== "unifire" && h.k !== "interproscan") right += "<p>No alignment could be made from these residues.</p>";  // a rule or a signature has no other protein to align
    }
    const linkedN = H.filter((x) => (h.t ? x.t === h.t : x.k === h.k && x.l === h.l) && x.r).length;
    const map = tlen && h.r && linkedN > 1 ? dotPlot(h, tlen, tname) : "";
    html += `<div class="pgrid${map ? "" : " one"}">${map}<div class="pright">${right}</div></div>`;
    box.innerHTML = html;
    box.querySelectorAll(".cells i[data-c]").forEach((c) => {
      c.addEventListener("mouseenter", (e) => {
        const k = classOf(c.dataset.c);
        tip.innerHTML = `<b>${esc(c.dataset.c)} ${c.dataset.p} of ${esc(c.dataset.w)}</b>` +
          (k ? `<div>${esc((HP && CLASS_WORD[k]) || "class")} (${esc(k)}) in ${esc(ALPHA_NAME)}: <span class="mono">${groups[k].sort().join("")}</span></div>` : "");
        tip.hidden = false; place(e);
      });
      c.addEventListener("mouseleave", () => { tip.hidden = true; });
    });
  }


  // ============================ the one-grid layout ============================
  // One CSS grid (--cols) for the ruler, the domain lane, the column names, every program row,
  // every hit row and the residues under an open hit, so the position track shares one x axis.
  const L = PAGE.q.length;
  const PROG = Object.keys(PAGE.tracks);
  const WHOLE = window.WHOLE || [];      // programs that name a hit with no position (Biomni)
  const NOTRUN = window.NOTRUN || {};    // programs not run, and why
  const gridEl = document.getElementById("grid");
  const emax = document.getElementById("emax"), hq = document.getElementById("hq");
  const countEl = document.getElementById("count"), sentEl = document.getElementById("sent"), sentNote = document.getElementById("sentnote");
  const maxE = () => Math.pow(10, +emax.value);
  const shownE = (h) => h.e === null || h.e === undefined || h.e <= maxE();
  const openProg = {}, openHit = {}, shownHits = {}, FIRST_HITS = 200, MORE_HITS = 500;
  let region = null;                     // [a, b, label] from a click on the domain lane: keep hits overlapping it
  let sortKey = "e", sortUp = true;      // order of the hit rows inside an unfolded program
  let TW = 400; const TH = 34, BH = 10, BY = 3;   // track width (measured), row height, box height and top

  function measure() {
    const t = document.querySelector(".colhead .track");
    if (t && t.clientWidth) TW = Math.max(120, t.clientWidth);
  }
  const sx = (p) => ((p - 1) / L) * TW, wx = (a, b) => Math.max(2, ((b - a + 1) / L) * TW);
  const fmtNum = (x) => (x === null || x === undefined ? "" : x.toLocaleString("en-US"));
  // the E-value to 2 significant figures (the program's own text is in the tooltip and the residues)
  const fmt2 = (e) => (e === 0 ? "0" : e >= 0.001 ? String(+e.toPrecision(2)) : e.toExponential(1).replace(".0e", "e"));
  const eText = (h) => (h.e === null || h.e === undefined ? "" : fmt2(h.e));
  // identity: the program's own figure where it gives one (ProstT5, MMseqs2), else counted on the
  // program's own alignment (phmmer, jackhmmer, hmmscan), else on kmerseek's gapless pair (h.id, build time)
  const idPct = (h) => { const p = identOf(h); if (p !== null) return p; return h.id ? (100 * h.id[0]) / h.id[1] : null; };
  const idText = (h) => { const p = idPct(h); return p === null ? "" : Math.round(p) + "%"; };
  // where a hit's identical count comes from: the program's own alignment, kmerseek's gapless pair, or a
  // realignment made when the site was built (ProstT5 beyond its own fident figure)
  const idSource = (h) => (h.k === "kmerseek" ? "kmerseek's gapless pair" : h.k === "diamond" ? "realigned (BLOSUM62, local, gap open -11, extend -1); DIAMOND's table kept no positions, and its own identity is in the row's numbers" : h.k === "prostt5" || h.k === "foldseek" ? "realigned (BLOSUM62, local, gap open -11, extend -1); the program's own figure is the percentage" : "the alignment the program wrote");
  const idFrom = (h) => (identOf(h) !== null ? "the program's own figure" : h.id ? idSource(h).split(";")[0] : "");
  // fainter = fewer identical residues: under 20%, 20-40, 40-60, 60 and over (the hits page's four steps)
  const idOpacity = (h) => { const p = idPct(h); return p === null ? 0.55 : p < 20 ? 0.3 : p < 40 ? 0.5 : p < 60 ? 0.75 : 1; };
  const idxOf = (k) => H.map((h, i) => i).filter((i) => H[i].k === k);
  const inRegion = (h) => !region || (h.q[0] <= region[1] && h.q[1] >= region[0]);
  const passes = (i) => shownE(H[i]) && inRegion(H[i]) && (!words.length || matches(i));
  function keptOf(k) { return idxOf(k).filter((i) => !H[i].x && shownE(H[i])); }   // at the E-value, before search and region
  function matchOf(k) { return keptOf(k).filter(passes).sort(cmp); }
  const featsOf = (k) => idxOf(k).filter((i) => H[i].x);
  const SORT = {
    e: (h) => (h.e === null || h.e === undefined ? Infinity : h.e), score: (h) => { const s = scoreOf(h)[0]; return s === null ? -Infinity : -s; },
    aa: (h) => h.q[0], id: (h) => { const p = idPct(h); return p === null ? Infinity : -p; }, tl: (h) => h.tl || 0,
  };
  function cmp(a, b) { const f = SORT[sortKey] || SORT.e; const d = f(H[a]) - f(H[b]); return (sortUp ? d : -d) || (H[a].e ?? 0) - (H[b].e ?? 0) || H[a].q[0] - H[b].q[0]; }
  const svgOpen = (hh, cls) => `<svg width="${TW}" height="${hh}" viewBox="0 0 ${TW} ${hh}" aria-hidden="true"${cls ? ` class="${cls}"` : ""}>`;
  const base = (y) => `<line x1="0" x2="${TW}" y1="${y}" y2="${y}" stroke="#C8C8C8" stroke-width="1.5"/>`;
  const markSvg = (tr) => `<svg width="12" height="10" viewBox="0 0 12 10" class="mk">${tr.shape === "half" ? `<rect x=".5" y=".5" width="11" height="9" rx="2" fill="#fff" stroke="${tr.color}"/><rect y="5" width="12" height="5" fill="${tr.color}"/>` : tr.shape === "open" ? `<rect x="1" y="1" width="10" height="8" rx="2" fill="#fff" stroke="${tr.color}" stroke-width="2"/>` : tr.shape === "hatched" ? `<rect x=".5" y=".5" width="11" height="9" rx="2" fill="url(#hatch-${tr.cat})" stroke="${tr.color}"/>` : `<rect width="12" height="10" rx="${tr.shape === "rounded" ? 5 : 2}" fill="${tr.color}"/>`}</svg>`;
  // the alphabet lane: the query rewritten in the alphabet of the kmerseek hit under the pointer
  // (default: the best kmerseek hit's). Class colours are used for nothing else on the page.
  const CLASS_COLS = ["#9b8bb4", "#9aa15c", "#c08a9e", "#7d8ea3", "#a07f5f", "#6fa39a", "#b59a3f", "#8f8f8f"];  // alphabet classes only
  let laneAlpha = null;
  function alphaOf() {
    if (laneAlpha && PAGE.alpha[laneAlpha]) return laneAlpha;
    const b = H.filter((h) => h.k === "kmerseek" && h.m && PAGE.alpha[h.m]).sort((a, b) => a.e - b.e)[0];
    return b ? b.m : null;
  }
  function classColours(name) {
    const A = PAGE.alpha[name], groups = {};
    for (const [r, c] of Object.entries(A)) (groups[c] = groups[c] || []).push(r);
    const names = Object.keys(groups).sort((a, b) => groups[b].length - groups[a].length);
    const col = {}; let extra = 0;
    names.forEach((c, i) => { col[c] = /^hp_/.test(name) && c === "h" ? "#f3b08f" : /^hp_/.test(name) && c === "p" ? "#9dbdeb" : CLASS_COLS[/^hp_/.test(name) ? extra++ : i % CLASS_COLS.length]; });
    return { groups, names, col };
  }
  const PL = (p) => (p >= 90 ? "#0053D6" : p >= 70 ? "#65CBF3" : p >= 50 ? "#FFDB13" : "#FF7D45");  // AlphaFold's pLDDT colours
  const meanPl = (a, b) => { if (!PAGE.plddt) return null; let s = 0, n = 0; for (let p = a; p <= b; p++) { const v = PAGE.plddt[p - 1]; if (v !== null && v !== undefined) { s += v; n++; } } return n ? s / n : null; };

  // the ruler, with the pLDDT strip under it (same x), so structure confidence sits with the hits
  function rulerRow() {
    const step = L <= 400 ? 50 : L <= 1000 ? 100 : L <= 2500 ? 250 : 500;
    let g = svgOpen(34) + `<line x1="0" x2="${TW}" y1="14" y2="14" stroke="#c3c2b7"/>`;
    // a tick label is drawn only where it clears the "1" at the start and the "L aa" at the end (about 6 px a character)
    const endX = TW - (String(L).length + 3) * 6.2 - 6;
    for (let p = step; p < L; p += step) if (sx(p + 1) + String(p).length * 3.2 < endX && sx(p + 1) - String(p).length * 3.2 > 12) g += `<line x1="${sx(p + 1).toFixed(1)}" x2="${sx(p + 1).toFixed(1)}" y1="10" y2="14" stroke="#c3c2b7"/><text x="${sx(p + 1).toFixed(1)}" y="9" class="ax" text-anchor="middle">${p}</text>`;
    g += `<text x="0" y="9" class="ax">1</text><text x="${TW}" y="9" class="ax" text-anchor="end">${L} aa</text>`;
    if (PAGE.plddt) {
      for (let a = 1; a <= L;) {
        let b = a; const c = PL(PAGE.plddt[a - 1]);
        while (b < L && PL(PAGE.plddt[b]) === c) b++;
        g += `<rect x="${sx(a).toFixed(1)}" y="19" width="${wx(a, b).toFixed(1)}" height="6" fill="${c}"><title>pLDDT ${Math.round(PAGE.plddt[a - 1])}: residues ${a}-${b}</title></rect>`;
        a = b + 1;
      }
    }
    const an = alphaOf();
    let alphaNote = "";
    if (an) {
      const { groups, names, col } = classColours(an), A = PAGE.alpha[an];
      for (let a = 1; a <= L;) {
        let b = a; const c = A[PAGE.q[a - 1]];
        while (b < L && A[PAGE.q[b]] === c) b++;
        if (c) g += `<rect x="${sx(a).toFixed(1)}" y="27" width="${wx(a, b).toFixed(1)}" height="5" fill="${col[c]}"><title>${esc(PAGE.q.slice(a - 1, b))} at ${a}-${b}: class ${esc(c)} (${groups[c].sort().join("")}) in ${esc(an)}</title></rect>`;
        a = b + 1;
      }
      const ak = document.getElementById("alphakey");
      if (ak) ak.innerHTML = `<b class="mono">${esc(an)}</b> classes: ` + names.map((c) => `<i class="pl" style="background:${col[c]}"></i><span class="mono">${esc(groups[c].sort().join(""))}</span>`).join(" ") + ` <span class="dim">(the alphabet of the kmerseek hit under the pointer)</span>`;
    }
    g += "</svg>";
    const pl = PAGE.plddt ? ` Under it: pLDDT, then your protein in a kmerseek alphabet (key above).` : "";
    return `<div class="g ruler"><div class="lbl">residue 1 to ${L} of your protein, the same axis in every row.${pl}</div><div class="track">${g}</div></div>`;
  }
  // the domain lane: hmmscan's Pfam domains, InterProScan's predictions (dashed) and the families the
  // headline sentence rests on (brackets with the count of programs). Click one to keep only hits in it.
  function laneRow() {
    // Pfam domains: hmmscan's; where hmmscan did not run (protein.html), InterProScan's Pfam matches, as on the first page
    let famSrc = "hmmscan", fam = H.map((h, i) => i).filter((i) => H[i].k === "hmmscan" && shownE(H[i]));
    if (!fam.length && !H.some((h) => h.k === "hmmscan")) { famSrc = "InterProScan"; fam = H.map((h, i) => i).filter((i) => H[i].k === "interproscan" && H[i].s[0] === "Pfam" && shownE(H[i])); }
    fam.sort((a, b) => (H[a].e ?? 0) - (H[b].e ?? 0));
    const kept = []; for (const i of fam) if (kept.every((j) => H[i].q[0] > H[j].q[1] || H[i].q[1] < H[j].q[0])) kept.push(i);
    const pred = [];
    H.filter((h) => h.x).sort((a, b) => a.q[0] - b.q[0]).forEach((h) => {
      const p = pred.find((p) => p.l === h.l && h.q[0] <= p.b + 1 && h.q[1] >= p.a - 1);
      if (p) { p.a = Math.min(p.a, h.q[0]); p.b = Math.max(p.b, h.q[1]); } else pred.push({ a: h.q[0], b: h.q[1], l: h.l });
    });
    const agree = (PAGE.agree || []).slice().sort((a, b) => a.a - b.a);
    if (!kept.length && !pred.length && !agree.length) return "";  // nothing to draw (protein.html: no hmmscan row)
    const lanes = [];  // brackets that overlap go on separate lines
    const laneOf = agree.map((a) => { let li = lanes.findIndex((end) => a.a > end + 2); if (li < 0) { lanes.push(a.b); li = lanes.length - 1; } else lanes[li] = a.b; return li; });
    const hh = 16 + (pred.length ? 14 : 0) + lanes.length * 15;
    let g = svgOpen(hh, "lane") + base(9);
    const on = (a, b) => (region && region[0] === a && region[1] === b ? " on" : "");
    for (const i of kept) {
      const h = H[i], x0 = sx(h.q[0]), w = wx(h.q[0], h.q[1]);
      g += `<g class="dom${on(h.q[0], h.q[1])}" data-a="${h.q[0]}" data-b="${h.q[1]}" data-l="${esc(h.l)}"><title>${esc(h.l)} (${esc((h.f || [])[0] || "Pfam")}), ${famSrc}, residues ${h.q[0]}-${h.q[1]}, E-value ${eText(h)}. Click to keep only hits in it.</title>` +
        `<rect x="${x0.toFixed(1)}" y="3" width="${w.toFixed(1)}" height="12" rx="2" fill="#efeee9" stroke="#3D3D3D" stroke-width=".8"/>` +
        (h.l.length * 6.4 + 4 <= w ? `<text x="${(x0 + w / 2).toFixed(1)}" y="12" class="dl" text-anchor="middle">${esc(h.l)}</text>` : "") + `</g>`;
    }
    let y = 16;
    if (pred.length) {
      for (const p of pred) {
        const x0 = sx(p.a), w = wx(p.a, p.b);
        g += `<g class="dom${on(p.a, p.b)}" data-a="${p.a}" data-b="${p.b}" data-l="${esc(p.l)}"><title>${esc(p.l)}, predicted by InterProScan, residues ${p.a}-${p.b}; not counted as a hit. Click to keep only hits in it.</title>` +
          `<rect x="${(x0 + 1).toFixed(1)}" y="${y + 2}" width="${Math.max(1, w - 2).toFixed(1)}" height="10" rx="2" fill="transparent" stroke="#3D3D3D" stroke-width="1.2" stroke-dasharray="3 2"/>` +
          `<text x="${(x0 + w / 2).toFixed(1)}" y="${y + 10}" class="dl" text-anchor="middle">${esc(p.l)} ${p.a}–${p.b}</text></g>`;
      }
      y += 14;
    }
    agree.forEach((a, n) => {
      const x0 = sx(a.a), x1 = sx(a.b) + wx(a.b, a.b), yy = y + laneOf[n] * 15;
      const lab = `${a.name}: ${a.n} program${a.n === 1 ? "" : "s"}`;
      g += `<g class="dom${on(a.a, a.b)}" data-a="${a.a}" data-b="${a.b}" data-l="${esc(a.name)}"><title>${esc(a.name)} (${a.acc}): ${a.hits} hits from ${esc(a.programs.join(", "))} whose matched proteins carry this family; residues ${a.a}-${a.b}. Click to keep only hits in it.</title>` +
        `<rect x="${x0.toFixed(1)}" y="${yy}" width="${(x1 - x0).toFixed(1)}" height="15" fill="transparent"/>` +
        `<path d="M${x0.toFixed(1)},${yy + 1}v5h${(x1 - x0).toFixed(1)}v-5" fill="none" stroke="#0b0b0b" stroke-width="1"/>` +
        `<text x="${((x0 + x1) / 2).toFixed(1)}" y="${yy + 14}" class="dl halo" text-anchor="middle">${esc(lab)}</text></g>`;
    });
    g += "</svg>";
    const note = agree.length ? `Brackets: the Pfam families of the matched proteins that the sentence at the top rests on, with how many programs' hits carry them.` : "";
    return `<div class="g row lane"><div class="name"><b>domains on your protein</b><div class="s" title="Pfam domains from ${famSrc === "hmmscan" ? "hmmscan" : "InterProScan's Pfam matches (hmmscan was run only on the 51 drawn proteins)"} (grey boxes); InterProScan predictions (dashed, never counted). ${note} Click any to keep only the hits landing in it.">Pfam (${famSrc}), predictions (dashed)${agree.length ? ", and the families the sentence rests on (brackets)" : ""}. Click one to filter.</div></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt"></div><div class="track">${g}</div></div>`;
  }
  // programs with a hit at each residue, each program once (the same steps build.py draws)
  function agreeRow() {
    const per = {};
    H.forEach((h) => { if (!h.x && shownE(h) && inRegion(h) && (!words.length || matches(H.indexOf(h)))) (per[PAGE.tracks[h.k].p] = per[PAGE.tracks[h.k].p] || []).push(h.q); });
    const cnt = new Array(L + 2).fill(0);
    for (const ivs of Object.values(per)) {
      const cover = new Uint8Array(L + 2);
      for (const [a, b] of ivs) for (let i = Math.max(1, a); i <= Math.min(L, b); i++) cover[i] = 1;
      for (let i = 1; i <= L; i++) cnt[i] += cover[i];
    }
    const hh = 40, top = 6, bot = hh - 8, yv = (c) => (bot - ((bot - top) * c) / Math.max(1, PAGE.m)).toFixed(1);
    let d = "", prev = null, mx = 0;
    for (let i = 0; i < L; i++) {
      const c = cnt[i + 1]; if (c > mx) mx = c;
      if (c !== prev) { d += prev === null ? `M0,${yv(c)}` : `L${sx(i + 1).toFixed(1)},${yv(prev)}L${sx(i + 1).toFixed(1)},${yv(c)}`; prev = c; }
    }
    d += `L${TW},${yv(prev)}`;
    const g = svgOpen(hh) + `<line x1="0" x2="${TW}" y1="${bot}" y2="${bot}" stroke="#e1e0d9"/>` +
      (mx ? `<line x1="0" x2="${TW}" y1="${yv(mx)}" y2="${yv(mx)}" stroke="#e1e0d9" stroke-dasharray="2 3"/>` : "") +
      `<path d="${d}" fill="none" stroke="#0b0b0b" stroke-width="1.4"/>` +
      `<text x="${TW + 4}" y="${+yv(mx) + 4}" class="ax">${mx} programs</text><text x="${TW + 4}" y="${bot + 4}" class="ax">0</text></svg>`;
    return `<div class="g row agree"><div class="name"><b>programs with a hit</b><div class="s">at each residue, each program once; the most anywhere: <b>${mx}</b> of ${PAGE.m} that place hits</div></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt"></div><div class="track">${g}</div></div>`;
  }
  // the folded track: a histogram of how many hits cover each residue, with the best hit's box on top
  let YMAX = 1;  // the tallest coverage on the page, so every program's histogram shares one scale
  function coverMax(list) { const diff = new Int32Array(L + 2); for (const i of list) { const [a, b] = H[i].q; diff[Math.max(1, a)]++; diff[Math.min(L, b) + 1]--; } let c = 0, m = 0; for (let p = 1; p <= L; p++) { c += diff[p]; if (c > m) m = c; } return m; }
  function coverSvg(list, tr, best) {
    const diff = new Int32Array(L + 2), dop = new Float64Array(L + 2);
    for (const i of list) { const [a, b] = H[i].q, o = idOpacity(H[i]); diff[Math.max(1, a)]++; diff[Math.min(L, b) + 1]--; dop[Math.max(1, a)] += o; dop[Math.min(L, b) + 1] -= o; }
    let c = 0, max = 0, so = 0; const cnt = new Int32Array(L + 1), op = new Float64Array(L + 1);
    for (let p = 1; p <= L; p++) { c += diff[p]; so += dop[p]; cnt[p] = c; op[p] = c ? so / c : 0; if (c > max) max = c; }
    const base_y = TH - 4, hmax = TH - 6;
    let g = svgOpen(TH) + `<line x1="0" x2="${TW}" y1="${base_y}" y2="${base_y}" stroke="#C8C8C8" stroke-width="1"/>`;
    for (let a = 1; a <= L;) {
      let b = a; while (b < L && cnt[b + 1] === cnt[a] && Math.abs(op[b + 1] - op[a]) < 0.01) b++;
      const n = cnt[a];
      if (n > 0) { const hgt = Math.max(1.5, (hmax * n) / YMAX); g += `<rect x="${sx(a).toFixed(1)}" y="${(base_y - hgt).toFixed(1)}" width="${wx(a, b).toFixed(1)}" height="${hgt.toFixed(1)}" fill="${tr.color}" opacity="${(0.15 + 0.6 * op[a]).toFixed(2)}"><title>residues ${a}-${b}: ${n} hit${n === 1 ? "" : "s"} (this program's peak: ${max}; the scale is the same in every row, top = ${YMAX}); shade = identity of the hits here</title></rect>`; }
      a = b + 1;
    }
    if (best !== undefined) {
      const h = H[best], x0 = sx(h.q[0]), w = wx(h.q[0], h.q[1]);
      g += `<g class="hit" data-i="${best}">${shapeSvg(x0, BY, w, BH, tr, idOpacity(h))}${seedSvg(h)}</g>`;
    }
    return g + "</svg>";
  }
  // a hit row's track: the box on your protein, and under it the matched protein as a thin line to its
  // own scale, the matched part dark, its Pfam domains as ticks (so the eye sees domain vs tail)
  // kmerseek's seed: the k residues in a row in the same class in both proteins, where its search started
  const seedSvg = (h) => (h.k === "kmerseek" && Array.isArray(h.sd) ? h.sd.map(([a, b]) => `<rect x="${sx(a).toFixed(1)}" y="${BY + 3}" width="${wx(a, b).toFixed(1)}" height="${BH - 6}" fill="#b8860b"><title>seed: residues ${a}-${b}, ${b - a + 1} in a row whose ${esc(h.m)} classes are the same in both proteins at every position (k = ${h.kk}); where kmerseek's search started</title></rect>`).join("") : "");
  // a row of predictions only (InterProScan on a protein with no signature hit): dashed boxes, labelled
  function featTrack(feats, tr) {
    let g = svgOpen(TH) + base(BY + BH / 2);
    const merged = [];  // predictions of one kind that overlap are drawn as one box, with their count
    feats.map((i) => H[i]).sort((a, b) => a.q[0] - b.q[0]).forEach((h) => {
      const p = merged.find((p) => p.l === h.l && h.q[0] <= p.b + 1 && h.q[1] >= p.a - 1);
      if (p) { p.a = Math.min(p.a, h.q[0]); p.b = Math.max(p.b, h.q[1]); p.n++; p.i.push(H.indexOf(h)); } else merged.push({ a: h.q[0], b: h.q[1], l: h.l, n: 1, i: [H.indexOf(h)] });
    });
    for (const p of merged) { const x0 = sx(p.a), w = wx(p.a, p.b); g += `<g class="hit" data-i="${p.i[0]}"><title>${esc(p.l)}: ${p.n} overlapping prediction${p.n === 1 ? "" : "s"} merged, residues ${p.a}-${p.b}</title>${featureSvg(x0, BY, w, BH, tr.color)}<text x="${(x0 + w / 2).toFixed(1)}" y="${BY + BH + 9}" class="ax halo" text-anchor="middle">${esc(p.l)} ${p.a}–${p.b}${p.n > 1 ? ` (${p.n} predictions)` : ""}</text></g>`; }
    return g + "</svg>";
  }
  function boxTrack(i) {
    const h = H[i], tr = PAGE.tracks[h.k], x0 = sx(h.q[0]), w = wx(h.q[0], h.q[1]);
    let g = svgOpen(TH) + base(BY + BH / 2) + `<g class="hit" data-i="${i}">${h.x ? featureSvg(x0, BY, w, BH, tr.color) : shapeSvg(x0, BY, w, BH, tr, idOpacity(h))}${seedSvg(h)}</g>`;
    if (h.r && h.tl) {
      const tl = h.tl, ty = 21, lo = Math.min(...h.r), hi = Math.max(...h.r), tsx = (p) => ((p - 1) / tl) * TW, twx = (a, b) => Math.max(2, ((b - a + 1) / tl) * TW);
      const who = h.k === "hmmscan" ? `${h.l} model` : h.l;
      g += `<g class="tline"><title>${esc(who)}, ${tl} aa, drawn to its own length; dark: the matched part ${lo}-${hi}${(h.td || []).length ? "; ticks: its Pfam domains" : ""}</title>` +
        `<line x1="0" x2="${TW}" y1="${ty}" y2="${ty}" stroke="#d6d5cf" stroke-width="2"/>`;
      for (const [a, b, l] of h.td || []) g += `<rect x="${tsx(a).toFixed(1)}" y="${ty - 3}" width="${twx(a, b).toFixed(1)}" height="6" fill="none" stroke="#8a8a8a" stroke-width=".8"><title>${esc(l)}, ${a}-${b} of ${esc(who)}</title></rect>`;
      g += `<line x1="${tsx(lo).toFixed(1)}" x2="${(tsx(lo) + twx(lo, hi)).toFixed(1)}" y1="${ty}" y2="${ty}" stroke="${tr.color}" stroke-width="3"/>` +
        ``;
      const td = h.td || [];
      if (td.length) {
        const nm = (l) => l.replace(/ \(PF\d+\)$/, "");
        const inside = td.find(([a, b]) => a <= hi && b >= lo);
        let lab;
        if (inside) lab = `in ${nm(inside[2])} ${inside[0]}–${inside[1]}`;
        else { const near = td.map((d) => [Math.min(Math.abs(d[0] - hi), Math.abs(lo - d[1])), d]).sort((x, y) => x[0] - y[0])[0]; lab = `${nm(near[1][2])} ${near[1][0]}–${near[1][1]}, ${near[0]} aa away`; }
        const lw = lab.length * 4.9, lx = Math.max(0, Math.min(tsx(lo), TW - lw));
        g += `<text x="${lx.toFixed(1)}" y="${ty + 9}" class="ax tl halo" text-anchor="start"><title>the matched part of the other protein lies ${inside ? "inside" : "outside"} its Pfam domains; "away" is the distance in residues to the nearest one</title>${esc(lab)}</text>`;
      }
      g += `</g>`;
    }
    return g + "</svg>";
  }
  // one tick per hit on a log E-value axis, 0.001 to 10, so a cloud of weak hits reads as a cloud
  function eStrip(list) {
    const W = 60, hh = 8;
    const x = (e) => (e <= 0 ? 0 : Math.max(0, Math.min(W, ((Math.log10(e) + 3) / 4) * W)));
    let g = `<svg width="${W}" height="${hh}" viewBox="0 0 ${W} ${hh}" class="estrip" role="img"><title>one tick per hit on a log axis: E-value 0.001 at the left edge, 10 at the right; ${list.length} hits. Ticks bunched at the right are many weak hits; one tick far left is a single strong one.</title><line x1="0" x2="${W}" y1="${hh - 1}" y2="${hh - 1}" stroke="#e1e0d9"/>`;
    for (const i of list) { const e = H[i].e; if (e === null || e === undefined) continue; g += `<line x1="${x(e).toFixed(1)}" x2="${x(e).toFixed(1)}" y1="0" y2="${hh - 2}" stroke="currentColor" stroke-width="1" opacity=".7"/>`; }
    return g + `</svg><div class="eax"><span>0.001</span><span>10</span></div>`;
  }
  const howOf = (k) => (window.HOW && window.HOW[k]) || "";
  // two programs with the same hit list (phmmer and jackhmmer here): say so, do not merge them
  function sameAs(k, match) {
    const key = (i) => `${H[i].t}|${H[i].q}|${H[i].e}`;
    const mine = new Set(match.map(key));
    for (const k2 of PROG) {
      if (k2 === k || PROG.indexOf(k2) > PROG.indexOf(k)) continue;
      const other = matchOf(k2);
      if (other.length && other.length === mine.size && other.every((i) => mine.has(key(i)))) return PAGE.tracks[k2].name;
    }
    return null;
  }
  function progRow(k) {
    const tr = PAGE.tracks[k], kept = keptOf(k), match = matchOf(k), feats = featsOf(k);
    const filtered = words.length || region;
    const best = match.slice().sort((a, b) => (H[a].e ?? 0) - (H[b].e ?? 0) || H[a].q[0] - H[b].q[0])[0];  // always the best-E-value hit, whatever the sort
    const open = !!openProg[k] && match.length > 0;
    const name = `<div class="name"><span class="tri">${match.length ? (open ? "&#9662;" : "&#9656;") : ""}</span>${markSvg(tr)}<b>${esc(tr.name)}</b><div class="s" title="${esc(howOf(k))}">${esc(howOf(k))}</div></div>`;
    let cells, track;
    if (NOTRUN[k]) {
      cells = `<div class="num dim">—</div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt dim" title="${esc(NOTRUN[k])}">${esc(NOTRUN[k])}</div>`;
      track = `<div class="track"></div>`;
    } else if (!kept.length) {
      const names = [...new Set(feats.map((i) => H[i].l))].join(", ");
      cells = feats.length
        ? `<div class="num"><span class="dim">0</span><div class="s">${feats.length} prediction${feats.length === 1 ? "" : "s"}</div></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt">${esc(names)}: predicted, drawn dashed, never counted as a hit<div class="s">no signature with an E-value at most ${fmtE(+emax.value)}</div></div>`
        : `<div class="num dim">0</div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt dim">no hit with E-value at most ${fmtE(+emax.value)}</div>`;
      track = `<div class="track">${feats.length ? featTrack(feats, tr) : ""}</div>`;
    } else if (!match.length) {
      cells = `<div class="num"><span class="dim">0/${kept.length}</span></div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt dim">none of its ${kept.length} hits ${region ? "lands in " + esc(region[2]) : "matches the search"}</div>`;
      track = `<div class="track"></div>`;
    } else {
      const b = H[best], sc = scoreOf(b);
      const nProt = new Set(match.map((i) => H[i].t || H[i].l)).size;
      const nSet = k === "kmerseek" ? new Set(match.map((i) => `${H[i].m}|${H[i].kk}`)).size : 0;
      const same = sameAs(k, match);
      const hitsTxt = filtered ? `<b>${match.length}</b><span class="dim">/${kept.length}</span>` : `<b>${match.length}</b>`;
      const peak = coverMax(match);
      const sub = `${nProt} protein${nProt === 1 ? "" : "s"}${nSet ? `<br>${nSet} alphabet+k setting${nSet === 1 ? "" : "s"}` : ""}${peak > 1 ? `<br>peak ${peak} on one residue` : ""}${same ? `<br><span class="same" title="its hit list is identical to ${esc(same)}'s: same proteins, residues and E-values">same hits as ${esc(same)}</span>` : ""}`;
      const gen = {}; for (const i of match) { const o = H[i].o; if (o) { const gname = o.split(/[\s(]/)[0]; gen[gname] = (gen[gname] || 0) + 1; } }
      const nSp = new Set(match.map((i) => H[i].o).filter(Boolean)).size;
      const topG = Object.entries(gen).sort((x, y) => y[1] - x[1]).slice(0, 3).map(([g, n]) => `<i>${esc(g)}</i> ${n}`).join(", ");
      const spTxt = nSp ? `<div class="s" title="distinct species among the ${match.length} hits, and the genera with most hits">${nSp} species${topG ? ": " + topG : ""}</div>` : "";
      cells = `<div class="num">${hitsTxt}<div class="s">${sub}</div></div>` +
        `<div class="num" style="color:${tr.color}"><span class="ink">${esc(eText(b))}</span>${match.length > 1 ? eStrip(match) : ""}</div>` +
        `<div class="num">${esc(sc[1].replace(" bits", "").replace("Prob ", ""))}</div><div class="num">${idText(b)}</div>` +
        `<div class="txt"><span class="mono">${esc(b.l)}</span>${b.n ? `<span class="s"> ${esc(b.n)}</span>` : ""}${spTxt}</div>`;
      track = `<div class="track">${coverSvg(match, tr, best)}</div>`;
    }
    const cls = `g row prog${open ? " open" : ""}${match.length ? "" : " nohit"}`;
    let out = `<div class="${cls}" data-k="${k}" role="button" tabindex="0" aria-expanded="${open}">${name}${cells}${track}</div>`;
    if (open) {  // a program with thousands of hits opens its best ones first; drawing them all at once takes seconds
      const lim = shownHits[k] || FIRST_HITS, left = match.length - lim;
      match.slice(0, lim).forEach((i, n) => { out += hitRow(i, n + 1); });
      if (left > 0) out += `<div class="g row morehits"><button type="button" data-more="${k}">Show ${Math.min(MORE_HITS, left).toLocaleString("en-US")} more of the ${left.toLocaleString("en-US")} not shown</button></div>`;
    }
    return out;
  }
  // the two most common residues of a stretch of your protein and their share; over 40% is a biased stretch
  function topTwo(a, b) {
    const cnt = {}; for (let p = a; p <= b; p++) { const c = PAGE.q[p - 1]; cnt[c] = (cnt[c] || 0) + 1; }
    const top = Object.entries(cnt).sort((x, y) => y[1] - x[1]).slice(0, 2);
    const share = top.reduce((s, x) => s + x[1], 0) / (b - a + 1);
    return { res: top.map((x) => x[0]).join("+"), pct: Math.round(100 * share) };
  }
  function hitRow(i, n) {
    const h = H[i], tr = PAGE.tracks[h.k], sc = scoreOf(h), isOpen = !!openHit[i];
    const fams = (h.f || []).map((f) => f.replace(/\s*\(PF\d+\)(:.*)?$/, "")).slice(0, 3).join(", ") + ((h.f || []).length > 3 ? ` +${h.f.length - 3}` : "");
    const sub = [h.n, h.o && `<i>${esc(h.o)}</i>`, fams && `Pfam: ${esc(fams)}`].filter(Boolean).map((s, j) => (j === 0 ? esc(s) : s)).join(" · ");
    const name = `<div class="name"><span class="rank">${n}</span><span class="nm"><span class="mono">${esc(h.l)}</span>${h.t && h.t !== h.l ? ` ${uniHtml(h.t)}` : ""}</span>` +
      `<div class="s" title="${esc([h.n, h.o, fams].filter(Boolean).join(" · "))}">${sub}</div></div>`;
    const setting = h.k === "kmerseek" && h.m ? `<span class="chip" title="the alphabet and k-mer size kmerseek compared the two proteins in; ${h.sn || 0} seed window${h.sn === 1 ? "" : "s"}">${esc(h.m)} k${h.kk}</span>` : esc(sc[1].replace(" bits", "").replace("Prob ", ""));
    const re = h.k === "prostt5" || h.k === "foldseek" || h.k === "diamond";  // the count below is from a realignment, not from Foldseek's own alignment
    const idc = h.id ? `<span title="${esc(idFrom(h))}">${idText(h)}</span><div class="s" title="${esc(idSource(h))}">${h.id[0]} of ${h.id[1]}${re ? `<br><span class='re'>= ${Math.round((100 * h.id[0]) / h.id[1])}% realigned</span>` : ""}</div>` : idText(h);
    const tt = topTwo(h.q[0], h.q[1]);
    const cells = `<div class="num">${h.q[0]}–${h.q[1]}` + (PAGE.noseq ? `<div class="s">${h.q[1] - h.q[0] + 1} aa</div></div>` :
      `<div class="s${tt.pct >= 40 ? " bias" : ""}" title="${h.q[1] - h.q[0] + 1} residues; the two most common, ${esc(tt.res)}, make up ${tt.pct}% of them${tt.pct >= 40 ? ": a composition-biased stretch" : ""}">${h.q[1] - h.q[0] + 1} aa · ${esc(tt.res)} ${tt.pct}%</div></div>`) + `<div class="num">${esc(eText(h))}</div><div class="num">${setting}</div><div class="num">${idc}</div>` +
      `<div class="txt">${h.r ? `<span class="s">${esc(h.k === "hmmscan" ? "model" : h.t || "")} ${Math.min(...h.r)}–${Math.max(...h.r)}${h.tl ? ` of ${h.tl}` : ""}</span>` : ""}${PAGE.plddt ? `<div class="s">pLDDT here ${Math.round(meanPl(h.q[0], h.q[1]))}</div>` : ""}</div>`;
    let out = `<div class="g row hit${isOpen ? " open" : ""}" data-i="${i}" role="button" tabindex="0" aria-expanded="${isOpen}" style="--c:${tr.color}">${name}${cells}<div class="track">${boxTrack(i)}</div></div>`;
    if (isOpen) out += `<div class="panel" data-i="${i}"><p class="note">Loading the residues...</p></div>`;
    return out;
  }
  // every program with nothing to show, on one line, so absence costs one row and not seven
  function silentRow(silent, notRun, wholes) {
    const chip = (tr, text, cut) => `<span class="sil" title="${esc(text)}">${esc(tr.name)}<span class="dim"> (${esc(cut)})</span></span>`;
    const cutOf = (k) => (PAGE.nohit && PAGE.nohit[k]) || (k === "hhsearch" ? "Prob ≥ 50" : `E ≤ ${fmtE(+emax.value)}`);
    const parts = [];
    const none = silent.map((k) => chip(PAGE.tracks[k], howOf(k), cutOf(k))).concat(wholes.map((w) => `<span class="sil" title="${esc(w[2])}">${esc(w[1])}<span class="dim"> (${esc(w[3].replace(/^no (hit|prediction) ?/, "").replace(/[()]/g, "") || "no prediction")})</span></span>`));
    if (none.length) parts.push(`<b>no hit:</b> ` + none.join(", "));
    if (notRun.length) parts.push(`<b>not run:</b> ` + notRun.map((k) => `<span class="sil" title="${esc(NOTRUN[k])}">${esc(PAGE.tracks[k].name)}<span class="dim"> (${esc(NOTRUN[k].replace(/^not run:?\s*/, ""))})</span></span>`).join(", "));
    return parts.length ? `<div class="g row silent"><div class="sil-body">${parts.join(" · ")}</div></div>` : "";
  }
  function wholeRow(w) {
    const [k, name, how, text, hit, trackText] = w;
    return `<div class="g row whole" data-k="${k}"><div class="name"><span class="tri"></span>${markSvg(PAGE.tracks[k.replace(/^w_/, "")] || { shape: "open", color: "#8a8a8a", cat: "seq" })}<b>${esc(name)}</b><div class="s" title="${esc(how)}">${esc(how)}</div></div>` +
      `<div class="num">${hit ? "1" : '<span class="dim">0</span>'}${hit === 2 ? '<div class="s">weak: see the distance</div>' : ""}</div><div class="num"></div><div class="num"></div><div class="num"></div><div class="txt${hit ? "" : " dim"}" title="${esc(text.replace(/<[^>]+>/g, ""))}">${text}</div><div class="track"><span class="s${hit === 2 ? " bias" : ""}">${trackText || (hit ? "whole protein, no position" : "")}</span></div></div>`;
  }
  function render() {
    placeColhead();
    document.documentElement.style.setProperty("--hh", headerEl.offsetHeight + "px");
    measure();
    YMAX = Math.max(1, ...PROG.map((k) => coverMax(matchOf(k))));
    let html = rulerRow() + laneRow() + agreeRow();
    const silent = [], notRun = [];
    for (const k of PROG) {
      if (NOTRUN[k]) notRun.push(k); else if (!keptOf(k).length && !featsOf(k).length) silent.push(k); else html += progRow(k);
    }
    for (const w of WHOLE) if (w[4]) html += wholeRow(w);           // 1: a hit; 2: a protein named but not counted (FANTASIA)
    html += silentRow(silent, notRun, WHOLE.filter((w) => !w[4]));
    gridEl.innerHTML = html;
    const v = +emax.value;
    let st = PAGE.sent[0]; PAGE.sent.forEach((x) => { if (x[0] <= v) st = x; });
    sentEl.title = PAGE.twelve ? "The programs counted: " + PAGE.twelve : ""; sentEl.innerHTML = st[1]; sentNote.innerHTML = (st[2] || "") + (region ? `<p class="note filt">Rows are filtered to ${esc(region[2])} ${region[0]}–${region[1]}; the sentence and the ladder count all hits.</p>` : "");
    const all = H.map((h, i) => i).filter((i) => !H[i].x), shown = all.filter((i) => shownE(H[i])), m = shown.filter(passes);
    const nProt = new Set(m.map((i) => H[i].t || H[i].l)).size;
    countEl.innerHTML = `<b>${m.length}</b> of ${all.length} hits · ${nProt} protein${nProt === 1 ? "" : "s"}` + (region ? ` · in <b>${esc(region[2])}</b> ${region[0]}–${region[1]} <button class="btn x" id="unregion" title="Show hits anywhere on the protein">×</button>` : "");
    document.querySelectorAll(".colhead [data-sort]").forEach((b) => { b.classList.toggle("on", b.dataset.sort === sortKey); b.querySelector("i").textContent = b.dataset.sort === sortKey ? (sortUp ? "▲" : "▼") : ""; });
    gridEl.querySelectorAll(".panel[data-i]").forEach((p) => open(H[+p.dataset.i], p));
    const un = document.getElementById("unregion"); if (un) un.addEventListener("click", () => { region = null; render(); });
  }
  // one click: a program row folds or unfolds; a hit row shows or hides its residues; a domain filters
  gridEl.addEventListener("click", (e) => {
    const more = e.target.closest("[data-more]");
    if (more) { const k = more.dataset.more; shownHits[k] = (shownHits[k] || FIRST_HITS) + MORE_HITS; render(); return; }
    if (e.target.closest("a") || e.target.closest("button")) return;
    const d = e.target.closest(".dom[data-a]");
    if (d) { const a = +d.dataset.a, b = +d.dataset.b; region = region && region[0] === a && region[1] === b ? null : [a, b, d.dataset.l]; render(); return; }
    const pr = e.target.closest(".row.prog");
    if (pr) { tip.hidden = true; const k = pr.dataset.k; if (!pr.classList.contains("nohit")) { openProg[k] = !openProg[k]; render(); focusRow(`.row.prog[data-k="${k}"]`); } return; }
    const hr = e.target.closest(".row.hit");
    if (hr) { tip.hidden = true; const i = +hr.dataset.i; openHit[i] = !openHit[i]; if (window.STRUCT) STRUCT.pin(openHit[i] ? H[i].q[0] : null, H[i].q[1]); render(); focusRow(`.row.hit[data-i="${i}"]`); }
  });
  gridEl.addEventListener("keydown", (e) => {
    const r = e.target.closest(".row[role=button]");
    if (!r) return;
    if (e.key === "Enter" || e.key === " ") { e.preventDefault(); r.click(); }
    if (e.key === "j" || e.key === "k" || e.key === "ArrowDown" || e.key === "ArrowUp") {
      const rows = [...gridEl.querySelectorAll(".row[role=button]")], n = rows.indexOf(r) + (e.key === "j" || e.key === "ArrowDown" ? 1 : -1);
      if (rows[n]) { e.preventDefault(); rows[n].focus(); }
    }
  });
  const focusRow = (sel) => { const r = gridEl.querySelector(sel); if (r) r.focus({ preventScroll: true }); };
  gridEl.addEventListener("mouseover", (e) => {
    const t = e.target.closest(".hit[data-i], .row.hit[data-i]");
    if (!t) return;
    const h = H[+t.dataset.i];
    tip.innerHTML = tipHtml(h); tip.hidden = false; place(e);
    if (window.STRUCT) STRUCT.show(h.q[0], h.q[1]);
    overlapShow(+t.dataset.i);
  });
  // hits in other programs that overlap the pointed one stay full, the rest fade; the count says how many
  const overlapNote = document.getElementById("overlap");
  function overlapShow(i) {
    const h = H[i], me = PAGE.tracks[h.k].p;
    if (h.k === "kmerseek" && h.m) { if (laneAlpha !== h.m) { laneAlpha = h.m; const r = gridEl.querySelector(".ruler"); if (r) r.outerHTML = rulerRow(); } }
    let n = 0; const progs = new Set();
    gridEl.querySelectorAll(".row.hit[data-i], .hit[data-i]").forEach((el) => {
      const j = +el.dataset.i, x = H[j], ov = x.q[0] <= h.q[1] && x.q[1] >= h.q[0];
      el.classList.toggle("faded", !ov);
      if (ov && j !== i && !x.x && el.classList.contains("row") && PAGE.tracks[x.k].p !== me) { n++; progs.add(PAGE.tracks[x.k].name); }
    });
    if (overlapNote) overlapNote.textContent = n ? `${h.l} ${h.q[0]}–${h.q[1]} overlaps ${n} hit${n === 1 ? "" : "s"} from ${progs.size} other unfolded program${progs.size === 1 ? "" : "s"} (folded programs are not counted)` : `${h.l} ${h.q[0]}–${h.q[1]}: no hit from another unfolded program overlaps it (folded programs are not counted)`;
  }
  function overlapClear() { gridEl.querySelectorAll(".faded").forEach((el) => el.classList.remove("faded")); if (overlapNote) overlapNote.textContent = ""; }
  gridEl.addEventListener("mousemove", (e) => { if (!tip.hidden) place(e); });
  gridEl.addEventListener("mouseout", (e) => {
    const t = e.target.closest(".hit[data-i], .row.hit[data-i]");
    if (t && !t.contains(e.relatedTarget)) { tip.hidden = true; if (window.STRUCT) STRUCT.show(null); overlapClear(); }
  });
  document.querySelectorAll(".colhead [data-sort]").forEach((b) => b.addEventListener("click", () => {
    const k = b.dataset.sort; if (sortKey === k) sortUp = !sortUp; else { sortKey = k; sortUp = true; } render();
  }));
  document.getElementById("expand-all").addEventListener("click", () => { for (const k of PROG) if (matchOf(k).length) openProg[k] = true; render(); });
  document.getElementById("collapse-all").addEventListener("click", () => { for (const k in openProg) delete openProg[k]; for (const i in openHit) delete openHit[i]; render(); });
  document.getElementById("clear").addEventListener("click", () => { emax.value = emax.dataset.start; hq.value = ""; words = []; region = null; sortKey = "e"; sortUp = true; laneAlpha = null; render(); if (ladder) ladder.innerHTML = ladderHtml(); });
  emax.addEventListener("change", render);
  // the ladder: one chip per E-value setting with what the sentence says there; the current one is dark
  const ladder = document.getElementById("ladder");
  function ladderHtml() {
    return [1, 0, -1, -2, -3].map((v) => {
      let st = PAGE.sent[0]; PAGE.sent.forEach((x) => { if (x[0] <= v) st = x; });
      const txt = st[1].replace(/<[^>]+>/g, ""), m = /^(\d+) of (?:the )?(\d+) programs?/.exec(txt);
      let what = /agree on (.+?)(?: \(|\.$)/.exec(txt) ? "agree on " + /agree on (.+?)(?: \(|\.$)/.exec(txt)[1] : /different Pfam families/.test(txt) ? "hits name different families" : /different parts/.test(txt) ? "hits on different parts" : "";  // no hit: "0 hits, 0 programs" says it
      const n = m ? m[1] : "0";
      const all = H.filter((h) => !h.x && (h.e === null || h.e === undefined || h.e <= Math.pow(10, v))).length;
      return `<button class="chip${+emax.value === v ? " on" : ""}" data-v="${v}" title="${esc(txt)}">E ≤ ${fmtE(v)}: ${fmtNum(all)} hit${all === 1 ? "" : "s"}, ${n} program${n === "1" ? "" : "s"}${what ? ", " + what : ""}</button>`;
    }).join("");
  }
  if (ladder) {
    ladder.innerHTML = ladderHtml();
    ladder.addEventListener("click", (e) => { const b = e.target.closest("button[data-v]"); if (!b) return; emax.value = b.dataset.v; render(); ladder.innerHTML = ladderHtml(); });
  }
  hq.addEventListener("input", () => { words = hq.value.toLowerCase().split(/\s+/).filter(Boolean); render(); });
  // download: the rows shown, one line per hit, as a tab-separated file
  document.getElementById("dl").addEventListener("click", () => {
    const rows = [["program", "hit", "accession", "protein", "species", "your_protein_start", "your_protein_end", "evalue", "score", "identity_pct", "identity_pct_source", "identical", "aligned_columns", "identical_count_source", "target_start", "target_end", "target_length", "kmerseek_alphabet", "k", "pfam_families_of_hit"]];
    for (const k of PROG) for (const i of matchOf(k)) { const h = H[i], sc = scoreOf(h), p = idPct(h);
      rows.push([PAGE.tracks[k].name, h.l, h.t || "", h.n || "", h.o || "", h.q[0], h.q[1], h.e ?? "", sc[0] ?? "", p === null ? "" : p.toFixed(1), p === null ? "" : idFrom(h), h.id ? h.id[0] : "", h.id ? h.id[1] : "", h.id ? idSource(h).split(";")[0] : "", h.r ? Math.min(...h.r) : "", h.r ? Math.max(...h.r) : "", h.tl || "", h.m || "", h.kk || "", (h.f || []).join("; ")]); }
    const blob = new Blob([rows.map((r) => r.map((c) => String(c).replace(/\t/g, " ")).join("\t")).join("\n") + "\n"], { type: "text/tab-separated-values" });
    const a = document.createElement("a"); a.href = URL.createObjectURL(blob); a.download = `${PAGE.id}_hits_E${fmtE(+emax.value)}.tsv`; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  });
  const help = document.getElementById("help"), helpBtn = document.getElementById("help-btn");
  helpBtn.addEventListener("click", () => { const on = help.classList.toggle("open"); helpBtn.setAttribute("aria-expanded", on); });
  document.addEventListener("keydown", (e) => {
    if (e.target.matches("input, select, textarea")) return;
    if (e.key === "Escape") { help.classList.remove("open"); let any = false; for (const i in openHit) { delete openHit[i]; any = true; } if (window.STRUCT) STRUCT.pin(null); if (any) render(); }
    if (e.key === "?") helpBtn.click();
    if ((e.key === "j" || e.key === "k") && !e.target.closest(".row")) { const r = gridEl.querySelector(".row[role=button]"); if (r) { e.preventDefault(); r.focus(); } }
    if (e.key === "s") document.getElementById("structure").scrollIntoView({ behavior: "smooth", block: "start" });
    if (e.key === "t") scrollTo({ top: 0, behavior: "smooth" });
  });
  // the pTM of the predicted structure in its heading: below 0.5 the arrangement of the parts is a guess
  const sh = document.getElementById("struct-h");
  if (sh && PAGE.ptm) sh.innerHTML += ` <span class="how">pTM ${esc(PAGE.ptm)}${+PAGE.ptm < 0.5 ? ": the arrangement of the parts is a guess, even where each part is confident" : ""}; mean pLDDT ${Math.round(meanPl(1, L))}</span>`;
  // on a phone the column names sit at the top of the grid, so they scroll sideways with it
  const colhead = document.querySelector(".colhead"), headerEl = document.querySelector("header.top"), mainEl = document.querySelector("main");
  function placeColhead() { if (innerWidth <= 900) { if (colhead.parentNode !== mainEl) mainEl.prepend(colhead); } else if (colhead.parentNode !== headerEl) headerEl.appendChild(colhead); }
  let rt; addEventListener("resize", () => { clearTimeout(rt); rt = setTimeout(render, 150); });
  if (hq.value) words = hq.value.toLowerCase().split(/\s+/).filter(Boolean);
  render();
  // the sticky ruler sits under the header: keep --hh equal to the header's height as the ladder, the
  // fonts and the window width change it (set here, at the end, once headerEl and the grid exist)
  if (window.ResizeObserver) new ResizeObserver(() => document.documentElement.style.setProperty("--hh", headerEl.offsetHeight + "px")).observe(headerEl);
})();
