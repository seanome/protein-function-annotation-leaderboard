// Tooltips, the "show all" redraw and the residue panel for one protein page.
// PAGE is set by the page: the protein sequence (q), the drawing scale (px, w), the hits
// (k: track, l: label, q: range on this protein, s: score values, t: target, r: range on the
// target, a: alignment if made at build time, m: kmerseek alphabet), score names per track,
// kmerseek's class tables (alpha) and the number of Swiss-Prot sequence files (tseq_files).
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

  function tipHtml(h) {
    let rows = `<tr><td>this protein</td><td>${h.q[0]}-${h.q[1]}</td></tr>`;
    if (h.r) rows += `<tr><td>${h.t ? esc(h.t) : h.k === "hmmscan" ? "Pfam model" : "matched protein"}</td><td>${h.r[0]}-${h.r[1]}</td></tr>`;
    const names = PAGE.names[h.k] || [];
    h.s.forEach((v, i) => { rows += `<tr><td>${esc(names[i] || "")}</td><td>${esc(v)}</td></tr>`; });
    return `<div>${esc(PAGE.tracks[h.k].name)}</div><b>${esc(h.l)}</b>${about(h, esc)}<table>${rows}</table>` +
      "<div>Click to see the residues.</div>";
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
  wire(document);

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
  // Folded programs: one row shaded by how many hits cover each residue. Remembered per program.
  // A row with more than FOLD_OVER hits at the slider's E-value starts folded; a click on its name
  // folds or unfolds it for the rest of this visit (folded[k] true or false). foldNow: as drawn now.
  const FOLD_OVER = 100;
  const folded = {}, foldNow = {};
  // A program can be folded when it has at least one hit (feature predictions do not count).
  const FOLDABLE = {};
  H.forEach((h) => { if (!h.x) FOLDABLE[h.k] = true; });
  // The folded row: the count of hits over each residue, darker for more, in the program's colour.
  function densitySvg(list, tr) {
    const L = PAGE.q.length, diff = new Int32Array(L + 2);
    for (const i of list) { const [a, b] = H[i].q; diff[Math.max(1, a)]++; diff[Math.min(L, b) + 1]--; }
    let c = 0, max = 0;
    const cnt = new Int32Array(L + 1);
    for (let p = 1; p <= L; p++) { c += diff[p]; cnt[p] = c; if (c > max) max = c; }
    const op = (n) => (max <= 1 ? 1 : 0.15 + (0.85 * (n - 1)) / (max - 1));
    // a see-through area over the whole row first, so a click anywhere on it unfolds the row
    let g = `<rect class="dens" x="0" y="0" width="${PAGE.w}" height="${ROW_PAD + BOX_H + 18}" fill="transparent"><title>Click to show each hit</title></rect>` +
      `<line x1="0" y1="${ROW_PAD + BOX_H / 2}" x2="${PAGE.w}" y2="${ROW_PAD + BOX_H / 2}" stroke="#C8C8C8" stroke-width="1.5"/>`;
    for (let a = 1; a <= L;) {
      let b = a;
      while (b < L && cnt[b + 1] === cnt[a]) b++;
      const n = cnt[a];
      if (n > 0) {
        g += `<rect class="dens" x="${((a - 1) * PAGE.px).toFixed(2)}" y="${ROW_PAD}" width="${((b - a + 1) * PAGE.px).toFixed(2)}" height="${BOX_H}" fill="${tr.color}" opacity="${op(n).toFixed(2)}">` +
          `<title>residues ${a}-${b}: ${n} hit${n === 1 ? "" : "s"}</title></rect>`;
      }
      a = b + 1;
    }
    const what = max ? `up to ${max} hit${max === 1 ? "" : "s"} on one residue; darker means more` : "no hit matches the search";
    g += `<text x="${PAGE.w - 2}" y="${ROW_PAD + BOX_H + 13}" class="tm" text-anchor="end">${what}</text>`;
    return g;
  }
  const tracksDiv = document.getElementById("tracks");
  const original = tracksDiv.innerHTML;
  // Redraw the rows with only the hits below maxE, and all of them or the best PAGE.top per tool.
  // slider position v is log10 of the E-value limit: 1 -> 10, -3 -> 0.001, -7 -> 1e-7
  const fmtE = (v) => (v >= -3 ? String(+Math.pow(10, v).toFixed(3)) : "1e" + v);
  function draw(maxE, showAll) {
    const tmp = document.createElement("div");
    tmp.innerHTML = original;
    const [ls, rs] = tmp.querySelectorAll("svg");
    const byTrack = {};
    H.forEach((h, i) => { (byTrack[h.k] = byTrack[h.k] || []).push(i); });
    const lr = ls.querySelectorAll("g.row"), rr = rs.querySelectorAll("g.row");
    let y = 0;
    lr.forEach((lg, n) => {
      const rg = rr[n], k = rg.dataset.k;
      let hgt = +rg.dataset.h;
      const idx = byTrack[k];
      const nh = !idx && rg.querySelector("text.nh");  // a row with no hit even at E-value 10
      if (nh) nh.textContent = `no hit with E-value at most ${fmtE(Math.round(Math.log10(maxE)))}`;
      if (idx) {
        const all = idx.filter((i) => H[i].e === null || H[i].e === undefined || H[i].e <= maxE);
        // feature predictions are always drawn and never count among the best PAGE.top hits
        const kept = all.filter((i) => !H[i].x), feats = all.filter((i) => H[i].x);
        // While searching, the best PAGE.top matching hits are drawn too, after the best PAGE.top hits.
        const best = showAll ? kept : kept.slice(0, PAGE.top);
        const shown = best.concat(words.length && !showAll ? kept.slice(PAGE.top).filter(matches).slice(0, PAGE.top - best.filter(matches).length) : [], feats);
        const sep = rg.querySelector("line.sep");
        lg.querySelectorAll(".more").forEach((o) => o.remove());
        const nMatch = words.length ? kept.filter(matches).length : 0;
        const leftLines = [];  // lines under the tool's name: [text, action or ""]; an action line is a link
        const hits = (n) => `${n.toLocaleString()} hit${n === 1 ? "" : "s"}`;
        if (words.length && kept.length) leftLines.push([`${nMatch} of ${hits(kept.length)} match`, ""]);
        if (!all.length) {
          rg.innerHTML = (sep ? sep.outerHTML : "") + `<text x="4" y="${ROW_PAD + 11}" class="nh">no hit with E-value at most ${fmtE(Math.round(Math.log10(maxE)))}</text>`;
          hgt = ROW_PAD + 22;
        } else if ((foldNow[k] = k in folded ? folded[k] : kept.length > FOLD_OVER) && kept.length) {
          rg.innerHTML = (sep ? sep.outerHTML : "") + densitySvg(words.length ? kept.filter(matches) : kept, PAGE.tracks[k]);
          hgt = ROW_PAD + BOX_H + 22;
          leftLines.push([`${hits(kept.length)}, folded`, ""]);
        } else {
          const [place, nLanes] = placeLabels(shown.map((i) => H[i]));
          const parts = shown.map((i, n) => boxSvg(H[i], i, ...place[n], words.length > 0 && !matches(i)));
          hgt = ROW_PAD + nLanes * LANE_H + 4;
          if (!kept.length) {  // only feature predictions: the row still says it has no hit
            parts.push(`<text x="4" y="${ROW_PAD + nLanes * LANE_H + 11}" class="nh">no hit with E-value at most ${fmtE(Math.round(Math.log10(maxE)))}; dashed: predictions, not counted</text>`);
            hgt += LANE_H;
          }
          const line = `<line x1="0" y1="${ROW_PAD + BOX_H / 2}" x2="${PAGE.w}" y2="${ROW_PAD + BOX_H / 2}" stroke="#C8C8C8" stroke-width="1.5"/>`;
          rg.innerHTML = (sep ? sep.outerHTML : "") + line + parts.join("");
          if (kept.length && !words.length) leftLines.push([hits(kept.length), ""]);
        }
        leftLines.forEach(([label, act], n) => {
          const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
          t.setAttribute("x", 26); t.setAttribute("y", ROW_PAD + 26 + 15 * n);
          t.setAttribute("class", act ? "tm more act" : "tm more");
          if (act) { t.dataset.k = k; t.dataset.act = act; t.setAttribute("role", "button"); t.setAttribute("tabindex", "0"); }
          t.textContent = label; lg.appendChild(t);
          hgt = Math.max(hgt, ROW_PAD + 34 + 15 * n);
        });
      }
      lg.setAttribute("transform", `translate(0,${y})`);
      rg.setAttribute("transform", `translate(0,${y})`);
      y += hgt;
    });
    // "programs with a hit" line, recounted from the hits below maxE (the same steps build.py draws)
    const ag = rs.querySelector('g.row[data-k="agree"] path');
    if (ag) {
      const L = PAGE.q.length, per = {};
      H.forEach((h) => { if (!h.x && (h.e === null || h.e === undefined || h.e <= maxE)) (per[PAGE.tracks[h.k].p] = per[PAGE.tracks[h.k].p] || []).push(h.q); });
      const diff = new Array(L + 2).fill(0);
      for (const ivs of Object.values(per)) {
        const cover = new Uint8Array(L + 2);
        for (const [a, b] of ivs) for (let i = Math.max(1, a); i <= Math.min(L, b); i++) cover[i] = 1;
        for (let i = 1; i <= L; i++) diff[i] += cover[i];
      }
      const yv = (c) => (10 + 48 - (c * 48) / Math.max(1, PAGE.m)).toFixed(1);
      let d = "", prev = null;
      for (let i = 0; i < L; i++) {
        const c = diff[i + 1];
        if (c !== prev) { d += prev === null ? `M0,${yv(c)}` : `L${(i * PAGE.px).toFixed(2)},${yv(prev)}L${(i * PAGE.px).toFixed(2)},${yv(c)}`; prev = c; }
      }
      d += `L${(L * PAGE.px).toFixed(2)},${yv(prev)}`;
      ag.setAttribute("d", d);
    }
    ls.setAttribute("height", y); ls.removeAttribute("viewBox");
    rs.setAttribute("height", y); rs.setAttribute("viewBox", `0 0 ${PAGE.w} ${y}`);
    return tmp.innerHTML;
  }
  const tg = document.getElementById("toggle"), slider = document.getElementById("emax"), ev = document.getElementById("ev");
  const showAll = true;  // every hit up to the slider's E-value is drawn
  // The drawing as built shows the hits at the counting threshold (data-start); the slider may start
  // elsewhere (build.py --show-at), and then the page redraws once on load.
  const START = slider ? +slider.dataset.start : 1;
  function redraw() {
    const v = slider ? +slider.value : START;
    if (ev) ev.textContent = fmtE(v);
    // The sentence at the top counts the hits at the slider's E-value (build.py sentence_steps).
    if (PAGE.sent && PAGE.sent.length) {
      let st = PAGE.sent[0];
      PAGE.sent.forEach((x) => { if (x[0] <= v) st = x; });
      document.getElementById("sent").innerHTML = st[1];  // HTML from build.py: escaped, accessions as links
      document.getElementById("sentnote").innerHTML = st[2];
    }
    tracksDiv.innerHTML = draw(Math.pow(10, v), showAll);
    wire(tracksDiv);
    decorate();
    hitTable(Math.pow(10, v));
    if (hqn) {
      const maxE = Math.pow(10, v);
      const shownHits = H.map((h, i) => i).filter((i) => !H[i].x && (H[i].e === null || H[i].e === undefined || H[i].e <= maxE));
      if (!words.length) hqn.textContent = "";
      else {
        const m = shownHits.filter(matches), ks = new Set(m.map((i) => H[i].k));
        const tools = Object.keys(PAGE.tracks).filter((k) => ks.has(k)).map((k) => PAGE.tracks[k].name);
        hqn.textContent = `${m.length.toLocaleString()} of ${shownHits.length.toLocaleString()} hits match` +
          (tools.length ? `: ${tools.join(", ")}` : "") + ".";
      }
    }
  }
  // Every hit the drawing shows (the slider's E-value and the search), in a table sorted by a
  // column the reader clicks: E-value first. Feature predictions (dashed boxes) are left out.
  const htDiv = document.getElementById("ht");
  let htSort = { c: 3, up: true };
  const num = (v) => { const x = parseFloat(v); return isNaN(x) ? null : x; };
  function scoreOf(h) {  // [value for sorting, text]: bits, or hhsearch's Prob, or InterProScan's score
    const k = h.k;
    if (k === "mmseqs2" || k === "mmseqs2_3rounds" || k === "prostt5") return [num(h.s[2]), h.s[2] + " bits"];
    if (k === "phmmer" || k === "jackhmmer" || k === "hmmscan") return [num(h.s[1]), h.s[1] + " bits"];
    if (k === "hhsearch") return [num(h.s[0]), "Prob " + h.s[0]];
    return [null, ""];
  }
  function identOf(h) {
    const k = h.k;
    if (k === "mmseqs2" || k === "mmseqs2_3rounds" || k === "prostt5") { const f = num(h.s[1]); return f === null ? null : 100 * f; }
    return null;
  }
  const HT_COLS = [
    ["Program", (h) => PAGE.tracks[h.k].name],
    ["Hit", (h) => h.l || ""],
    ["Protein, species", (h) => [h.n, h.o].filter(Boolean).join(", ")],
    ["E-value", (h) => (h.e === null || h.e === undefined ? null : h.e)],
    ["Score", (h) => scoreOf(h)[0]],
    ["Identity (%)", identOf],
    ["Your protein (aa)", (h) => h.q[0]],
    ["Pfam families on the hit", (h) => (h.f || []).join("; ")],
  ];
  function hitTable(maxE) {
    if (!htDiv) return;
    let rows = H.map((h, i) => i).filter((i) => !H[i].x && (H[i].e === null || H[i].e === undefined || H[i].e <= maxE));
    if (words.length) rows = rows.filter(matches);
    const get = HT_COLS[htSort.c][1];
    rows.sort((a, b) => {
      const x = get(H[a]), y = get(H[b]);
      if (x === null || x === "") return (y === null || y === "") ? 0 : 1;  // empty values last either way
      if (y === null || y === "") return -1;
      const d = typeof x === "number" ? x - y : String(x).localeCompare(String(y));
      return htSort.up ? d : -d;
    });
    const head = HT_COLS.map((c, n) => `<th data-c="${n}" class="${n === htSort.c ? (htSort.up ? "up" : "down") : ""}" ` +
      `aria-sort="${n === htSort.c ? (htSort.up ? "ascending" : "descending") : "none"}">${esc(c[0])}</th>`).join("");
    const body = rows.map((i) => {
      const h = H[i], sc = scoreOf(h), idt = identOf(h);
      return `<tr data-i="${i}"><td>${esc(PAGE.tracks[h.k].name)}</td><td class="mono">${esc(h.l || "")}` +
        (h.t && h.t !== h.l ? ` ${uniHtml(h.t)}` : "") + `</td><td>${esc(h.n || "")}${h.o ? ", <i>" + esc(h.o) + "</i>" : ""}</td>` +
        `<td class="r">${h.e === null || h.e === undefined ? "" : esc(h.s[k_e(h)] || h.e)}</td><td class="r">${esc(sc[1])}</td>` +
        `<td class="r">${idt === null ? "" : idt.toFixed(1)}</td><td class="r">${h.q[0]}-${h.q[1]}</td>` +
        `<td>${accHtml((h.f || []).slice(0, 3).join("; ") + ((h.f || []).length > 3 ? "; and " + (h.f.length - 3) + " more" : ""))}</td></tr>`;
    }).join("");
    htDiv.innerHTML = `<p class="note">${rows.length.toLocaleString()} hit${rows.length === 1 ? "" : "s"}. Click a column name to sort by it; click a row to see the residues.</p>` +
      `<div class="tablewrap"><table class="small hits"><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`;
  }
  // the E-value as the program printed it: the score column that holds it
  function k_e(h) {
    const names = PAGE.names[h.k] || [];
    const n = names.findIndex((x) => /e-?value/i.test(x));
    return n < 0 ? -1 : n;
  }
  if (htDiv) htDiv.addEventListener("click", (e) => {
    const th = e.target.closest("th[data-c]");
    if (th) { const c = +th.dataset.c; htSort = { c, up: htSort.c === c ? !htSort.up : c !== 4 }; hitTable(currentMaxE()); return; }
    if (e.target.closest("a")) return;
    const tr = e.target.closest("tr[data-i]");
    if (tr) { if (window.STRUCT) STRUCT.pin(H[+tr.dataset.i].q[0], H[+tr.dataset.i].q[1]); open(H[+tr.dataset.i]); }
  });
  // hovering a row shows that hit's residues on the 3-D structure, when the page has one
  if (htDiv) htDiv.addEventListener("mouseover", (e) => {
    const tr = e.target.closest("tr[data-i]");
    if (window.STRUCT) STRUCT.show(tr ? H[+tr.dataset.i].q[0] : null, tr ? H[+tr.dataset.i].q[1] : null);
  });
  if (htDiv) htDiv.addEventListener("mouseleave", () => { if (window.STRUCT) STRUCT.show(null); });
  const currentMaxE = () => Math.pow(10, slider ? +slider.value : START);
  // Each program's name (mark, name and a small triangle) becomes a button that folds its row.
  function decorate() {
    tracksDiv.querySelectorAll(".names g.row").forEach((lg) => {
      const k = lg.dataset.k, t = lg.querySelector("text.tn");
      if (!FOLDABLE[k] || !t || lg.querySelector(".fold")) return;
      const g = document.createElementNS(SVGNS, "g");
      g.setAttribute("class", "fold"); g.setAttribute("role", "button"); g.setAttribute("tabindex", "0");
      g.dataset.k = k;
      g.setAttribute("aria-expanded", foldNow[k] ? "false" : "true");
      const name = PAGE.tracks[k].name;
      g.setAttribute("aria-label", foldNow[k] ? `${name}: show each hit again` : `${name}: fold its hits into one row`);
      const tt = document.createElementNS(SVGNS, "title");
      tt.textContent = foldNow[k] ? "Click to show each hit again" : "Click to fold these hits into one row";
      const mark = t.previousElementSibling;
      g.appendChild(tt);
      if (mark && mark.tagName === "g") g.appendChild(mark);
      g.appendChild(t);
      lg.appendChild(g);
      // a small triangle after the name: pointing down while the hits are shown, right when folded
      const x = +t.getAttribute("x") + t.getComputedTextLength() + 6, y = +t.getAttribute("y") - 5;
      const tri = document.createElementNS(SVGNS, "path");
      tri.setAttribute("class", "tri");
      tri.setAttribute("d", foldNow[k] ? `M${x},${y - 5}l7,5l-7,5z` : `M${x},${y - 3}l9,0l-4.5,7z`);
      g.appendChild(tri);
      const area = document.createElementNS(SVGNS, "rect");  // one click target over mark, name and triangle
      area.setAttribute("x", 0); area.setAttribute("y", y - 9); area.setAttribute("width", x + 12); area.setAttribute("height", 18);
      area.setAttribute("fill", "transparent");
      g.insertBefore(area, g.firstChild.nextSibling);
    });
  }
  function toggleFold(k) {
    folded[k] = !foldNow[k];
    redraw();
  }
  // Click a program's name (or its triangle) to fold or unfold its row; click a folded row to unfold it.
  tracksDiv.addEventListener("click", (e) => {
    const f = e.target.closest(".fold");
    if (f) { toggleFold(f.dataset.k); const g = tracksDiv.querySelector(`.fold[data-k="${f.dataset.k}"]`); if (g) g.focus(); return; }
    const d = e.target.closest(".dens"), row = d && d.closest("g.row");
    if (row && foldNow[row.dataset.k]) toggleFold(row.dataset.k);
  });
  tracksDiv.addEventListener("keydown", (e) => {
    const f = e.target.closest(".fold");
    if (f && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); toggleFold(f.dataset.k); }
  });
  const hq = document.getElementById("hq"), hqn = document.getElementById("hqn");
  if (hq) hq.addEventListener("input", () => { words = hq.value.toLowerCase().split(/\s+/).filter(Boolean); redraw(); });
  if (slider) slider.addEventListener("input", redraw);
  if (hq && hq.value) words = hq.value.toLowerCase().split(/\s+/).filter(Boolean);  // text the browser kept on reload
  redraw();  // the drawing as built has no links; this adds them and applies the slider, folds and search
  document.getElementById("close").addEventListener("click", () => { panel.hidden = true; });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") panel.hidden = true; });

  // Swiss-Prot sequences live in ../targets/NNN.js, one file per group of accessions.
  function shardOf(acc) {
    let h = 0;
    for (const ch of acc) h = (Math.imul(h, 31) + ch.charCodeAt(0)) >>> 0;
    return h % PAGE.tseq_files;
  }
  const loading = {};
  function targetSeq(acc) {
    window.TSEQ = window.TSEQ || {};
    if (acc in TSEQ) return Promise.resolve(TSEQ[acc]);
    const n = shardOf(acc);
    if (!loading[n]) loading[n] = new Promise((res) => {
      const s = document.createElement("script");
      s.src = `../targets/${String(n).padStart(3, "0")}.js`;
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
    if (!alnLoad) alnLoad = new Promise((res) => {
      const s = document.createElement("script");
      s.src = location.pathname.split("/").pop().replace(/\.html$/, ".aln.js");
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
  const PW = 220, Y_MIN = 110, Y_MAX = 320;
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
    if (fam.length) top.push(`solid boxes: Pfam families and other InterPro signatures that hmmscan or InterProScan found on your protein${sl ? ` (E-value at most ${fmtE(+sl.value)}, the slider's value)` : ""}; where two overlap, the one with the better E-value`);
    if (pred.length) top.push(`dashed boxes: ${[...new Set(pred.map((p) => p.l))].join(", ")}, predicted by InterProScan`);
    const topLine = top.length ? `Along the top, ${top.join("; ")}.` : "Along the top: nothing; hmmscan and InterProScan found no family or prediction on your protein.";
    let rightLine;
    if (h.k === "hmmscan") rightLine = "The y axis is the Pfam model itself, so nothing is drawn on its edge.";
    else if (h.td && h.td.length) rightLine = `Along the right, open boxes: Pfam domains of ${uniHtml(h.t)}, ${esc((SRC[h.ds] || (() => "source not recorded"))(h.t))}; where two overlap, only the first is drawn.`;
    else if (h.ds) rightLine = `Along the right: nothing; no Pfam domain is known on ${uniHtml(h.t)} (${esc(SRC[h.ds](h.t))}).`;
    else rightLine = `Along the right: nothing; the Pfam domains of ${uniHtml(h.t)} are not known (no Pfam search output for it).`;
    return `<div class="map" style="max-width:${W}px"><p class="mapkey">Across: your protein. Down: ${esc(tname)}. Each line is one hit between the two, in the colour of its kind of program; the thick one is the hit shown on the right.</p>` +
      `<div class="mapkey tools">${tools}</div><p class="mapkey">${esc(topLine)}</p><p class="mapkey">${rightLine}</p>` +
      `<svg viewBox="0 0 ${W} ${Hh}" width="${W}" height="${Hh}" role="img" aria-label="Map of the hits between your protein and ${esc(tname)}">${g}</svg></div>`;
  }

  function open(h) {
    tip.hidden = true;
    pbody.innerHTML = "<p>Loading the residues...</p>"; panel.hidden = false;
    Promise.all([alignments(), h.t ? targetSeq(h.t) : Promise.resolve(null)]).then(([alns, tseq]) => show(h, alns[H.indexOf(h)], tseq));
  }

  function show(h, a, tseq) {
    const q = PAGE.q, tool = PAGE.tracks[h.k].name;
    setAlphabet(h.k === "kmerseek" ? h.m : null);
    const tname = h.k === "hmmscan" ? `${h.l} model` : h.l;
    const names = PAGE.names[h.k] || [];
    let html = `<h3>${esc(tool)}: ${esc(h.l)}` + (h.t && h.t !== h.l ? ` <span class="w">Swiss-Prot ${uniHtml(h.t)}</span>` : "") + "</h3>" + about(h, accHtml);
    const scoreBits = h.s.map((v, i) => [names[i], v]).filter(([, v]) => v !== "not in the file" && !/^no .* column/.test(v)).map(([n, v]) => `${esc(n || "")} <b>${accHtml(v)}</b>`).join(" · ");
    let src = "", kmer = false, flank = null, note = "";
    if (h.k === "kmerseek" && tseq) {
      kmer = true;
      const F = 10, q0 = Math.max(1, h.q[0] - F), t0 = Math.max(1, h.r[0] - F);
      const lf = Math.min(h.q[0] - q0, h.r[0] - t0);
      const rf = Math.min(q.length - h.q[1], tseq.length - h.r[1], F);
      a = { q: q.slice(h.q[0] - 1 - lf, h.q[1] + rf), t: tseq.slice(h.r[0] - 1 - lf, h.r[1] + rf), qs: h.q[0] - lf, ts: h.r[0] - lf };
      flank = [lf, rf];
      src = `kmerseek matches the two segments side by side without gaps; ${lf} residues before and ${rf} after the match are shown faded.`;
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
    if (h.k !== "hhsearch" && h.t && !tseq) html += `<p class="note">${uniHtml(h.t)} is not in swissprot.fasta, so its residues cannot be shown.</p>`;
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
      if (inSeed) leg.push(`<span><i class="sd"></i>seed: ${h.kk} residues in a row in the same ${esc(ALPHA_NAME)} class in both, where kmerseek's search started</span>`);
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
      const avail = (panel.clientWidth || innerWidth) - 40 - (tlen && h.r && innerWidth > 760 ? 400 : 0);
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
      right += "<p>No alignment could be made from these residues.</p>";
    }
    const map = tlen && h.r ? dotPlot(h, tlen, tname) : "";
    html += `<div class="pgrid">${map}<div class="pright">${right}</div></div>`;
    pbody.innerHTML = html;
    panel.hidden = false;
    pbody.querySelectorAll(".cells i[data-c]").forEach((c) => {
      c.addEventListener("mouseenter", (e) => {
        const k = classOf(c.dataset.c);
        tip.innerHTML = `<b>${esc(c.dataset.c)} ${c.dataset.p} of ${esc(c.dataset.w)}</b>` +
          (k ? `<div>${esc((HP && CLASS_WORD[k]) || "class")} (${esc(k)}) in ${esc(ALPHA_NAME)}: <span class="mono">${groups[k].sort().join("")}</span></div>` : "");
        tip.hidden = false; place(e);
      });
      c.addEventListener("mouseleave", () => { tip.hidden = true; });
    });
  }
})();
