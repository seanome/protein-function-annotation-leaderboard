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
  const H = PAGE.hits;

  function tipHtml(h) {
    let rows = `<tr><td>this protein</td><td>${h.q[0]}-${h.q[1]}</td></tr>`;
    if (h.r) rows += `<tr><td>${esc(h.t)}</td><td>${h.r[0]}-${h.r[1]}</td></tr>`;
    const names = PAGE.names[h.k] || [];
    h.s.forEach((v, i) => { rows += `<tr><td>${esc(names[i] || "")}</td><td>${esc(v)}</td></tr>`; });
    return `<div>${esc(PAGE.tracks[h.k].name)}</div><b>${esc(h.l)}</b>${about(h)}<table>${rows}</table>` +
      "<div>Click to see the residues.</div>";
  }
  // Full name, gene and species of the matched Swiss-Prot protein, and its Pfam families over
  // the matched part (for hmmscan and InterProScan: the family itself).
  function about(h) {
    let s = "";
    if (h.n) s += `<div class="full">${esc(h.n)}</div>`;
    const gs = [h.g ? "gene " + h.g : "", h.o || ""].filter(Boolean).join(" · ");
    if (gs) s += `<div>${esc(gs)}</div>`;
    if (h.k === "hmmscan" || h.k === "interproscan") {
      if (h.f && h.f.length) s += `<div>${h.k === "hmmscan" ? "Pfam family" : "Family"}: ${esc(h.f.join("; "))}</div>`;
    } else if (h.f === null || h.f === undefined) {
      s += "<div>Pfam families on the matched part: not known (no hmmscan or InterProScan output for this protein)</div>";
    } else {
      const shown = h.f.slice(0, 5).join("; ") + (h.f.length > 5 ? `; and ${h.f.length - 5} more` : "");
      s += `<div>Pfam families on the matched part: ${h.f.length ? esc(shown) : "none"}</div>`;
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
      g.addEventListener("mouseenter", (e) => { tip.innerHTML = tipHtml(h); tip.hidden = false; place(e); });
      g.addEventListener("mousemove", place);
      g.addEventListener("mouseleave", () => { tip.hidden = true; });
      g.addEventListener("focus", () => {
        const b = g.getBoundingClientRect();
        tip.innerHTML = tipHtml(h); tip.hidden = false; place({ clientX: b.left, clientY: b.bottom });
      });
      g.addEventListener("blur", () => { tip.hidden = true; });
      g.addEventListener("click", () => open(h));
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
  function boxSvg(h, i, lane, side, lx) {
    const tr = PAGE.tracks[h.k];
    const x0 = (h.q[0] - 1) * PAGE.px, w = Math.max(2, h.q[1] * PAGE.px - x0), y = ROW_PAD + lane * LANE_H;
    const op = eOpacity(h.e);
    const rect = shapeSvg(x0, y, w, BOX_H, tr, op), ty = y + BOX_H - 3;
    let lab;
    if (side === "in") {
      const solid = (tr.shape === "filled" || tr.shape === "rounded") && op >= 0.6;
      lab = `<text x="${(x0 + 3).toFixed(1)}" y="${ty}" fill="${solid ? "#ffffff" : "#1a1a1a"}" class="bl${solid ? "" : " halo"}">${esc(h.l)}</text>`;
    } else {
      const lw = textW(h.l), [l0, l1] = side === "right" ? [x0 + w, lx - 1] : [lx + lw + 1, x0];
      const lead = l1 > l0 ? `<line x1="${l0.toFixed(1)}" y1="${y + BOX_H / 2}" x2="${l1.toFixed(1)}" y2="${y + BOX_H / 2}" stroke="${tr.color}" stroke-width="1"/>` : "";
      lab = lead + `<text x="${lx.toFixed(1)}" y="${ty}" fill="#1a1a1a" class="bl out">${esc(h.l)}</text>`;
    }
    return `<g class="hit" tabindex="0" data-i="${i}">${rect}${lab}</g>`;
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
        const kept = idx.filter((i) => H[i].e === null || H[i].e === undefined || H[i].e <= maxE);
        const shown = showAll ? kept : kept.slice(0, PAGE.top);
        const sep = rg.querySelector("line.sep"), old = lg.querySelector(".more");
        if (old) old.remove();
        if (!kept.length) {
          rg.innerHTML = (sep ? sep.outerHTML : "") + `<text x="4" y="${ROW_PAD + 11}" class="nh">no hit with E-value at most ${fmtE(Math.round(Math.log10(maxE)))}</text>`;
          hgt = ROW_PAD + 22;
        } else {
          const [place, nLanes] = placeLabels(shown.map((i) => H[i]));
          const parts = shown.map((i, n) => boxSvg(H[i], i, ...place[n]));
          hgt = ROW_PAD + nLanes * LANE_H + 4;
          const line = `<line x1="0" y1="${ROW_PAD + BOX_H / 2}" x2="${PAGE.w}" y2="${ROW_PAD + BOX_H / 2}" stroke="#C8C8C8" stroke-width="1.5"/>`;
          rg.innerHTML = (sep ? sep.outerHTML : "") + line + parts.join("");
          const more = kept.length - shown.length;
          const label = more > 0 ? `+${more} more hit${more === 1 ? "" : "s"}` : (showAll && kept.length > PAGE.top ? `all ${kept.length} hits` : "");
          if (label) {
            const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
            t.setAttribute("x", 26); t.setAttribute("y", ROW_PAD + 26); t.setAttribute("class", "tm more");
            t.textContent = label; lg.appendChild(t);
            hgt = Math.max(hgt, ROW_PAD + 34);
          }
        }
      }
      lg.setAttribute("transform", `translate(0,${y})`);
      rg.setAttribute("transform", `translate(0,${y})`);
      y += hgt;
    });
    // "programs with a hit" line, recounted from the hits below maxE (the same steps build.py draws)
    const ag = rs.querySelector('g.row[data-k="agree"] path');
    if (ag) {
      const L = PAGE.q.length, per = {};
      H.forEach((h) => { if (h.e === null || h.e === undefined || h.e <= maxE) (per[PAGE.tracks[h.k].p] = per[PAGE.tracks[h.k].p] || []).push(h.q); });
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
  let showAll = false;
  // The drawing as built shows the hits at the slider's starting value (the counting threshold).
  const START = slider ? +slider.dataset.start : 1;
  function redraw() {
    const v = slider ? +slider.value : START;
    if (ev) ev.textContent = fmtE(v);
    tracksDiv.innerHTML = v === START && !showAll ? original : draw(Math.pow(10, v), showAll);

    wire(tracksDiv);
    if (tg) {
      const maxE = Math.pow(10, v), n = H.filter((h) => h.e === null || h.e === undefined || h.e <= maxE).length;
      const byK = {};
      H.forEach((h) => { if (h.e === null || h.e === undefined || h.e <= maxE) byK[h.k] = (byK[h.k] || 0) + 1; });
      const more = Object.values(byK).reduce((s, c) => s + Math.max(0, c - PAGE.top), 0);
      tg.hidden = !showAll && more === 0;
      tg.textContent = showAll ? `Show the best ${PAGE.top} hits per tool` : `Show all ${n.toLocaleString()} hits (${more.toLocaleString()} more)`;
    }
  }
  if (tg) tg.addEventListener("click", () => { showAll = !showAll; redraw(); });
  if (slider) slider.addEventListener("input", redraw);
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
  // Residue classes: the kmerseek alphabet of this page when there is one (hp_lehninger2 here),
  // used to colour residues in every tool's alignment, so one colour means one class everywhere.
  const ALPHA_NAME = Object.keys(PAGE.alpha || {})[0] || null;
  const ALPHA = ALPHA_NAME ? PAGE.alpha[ALPHA_NAME] : null;
  const CLASS_WORD = { h: "hydrophobic", p: "polar" };
  const groups = {};
  if (ALPHA) for (const [r, c] of Object.entries(ALPHA)) (groups[c] = groups[c] || []).push(r);
  const classOf = (c) => (ALPHA && c && c !== "-" ? ALPHA[c.toUpperCase()] || null : null);
  const cssClass = (c) => { const k = classOf(c); return k === "h" || k === "p" ? k : (k ? "o" : "x"); };
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

  // One line of 60 columns: your protein, [its classes], a match line, [target classes], the target.
  function chunkHtml(a, start, n, qpos, tpos, tname, kmer, flank) {
    const W = 60, rowsQ = [], rowsT = [], mid = [], cq = [], ct = [];
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
    let h = line(`your protein aa ${qFirst}`, rowsQ.join(""), qi - 1);
    if (kmer) h += line(`your protein ${esc(ALPHA_NAME)}`, cq.join(""), "", "enc");
    h += line("", mid.map((m) => `<i>${m === " " ? "&nbsp;" : m}</i>`).join(""), "", "marks");
    if (kmer) h += line(`${esc(tname)} ${esc(ALPHA_NAME)}`, ct.join(""), "", "enc");
    h += line(`${esc(tname)} aa ${tFirst}`, rowsT.join(""), ti - 1, "under");
    return { html: `<div class="chunk">${h}</div>`, qi, ti };
  }

  // A small map: your protein across, the other sequence down, every hit on that sequence as a
  // line in its tool's colour, the one shown drawn thick; Pfam domains of both along the edges.
  function dotPlot(h, tlen, tname) {
    const qL = PAGE.q.length, S = 200, sc = S / Math.max(qL, tlen);
    const pw = qL * sc, ph = tlen * sc, ml = 34, mt = 26, tw = 22;
    const W = ml + pw + tw + 70, Hh = mt + ph + 34;
    const sx = (p) => ml + (p - 1) * sc, sy = (p) => mt + (p - 1) * sc;
    let g = "";
    // Pfam domains of your protein (hmmscan, best E-value first, not overlapping) along the top
    const qd = [];
    H.filter((x) => x.k === "hmmscan").sort((a, b) => (a.e ?? 1) - (b.e ?? 1)).forEach((x) => {
      if (qd.every((d) => x.q[1] < d[0] || x.q[0] > d[1])) qd.push([x.q[0], x.q[1], x.l]);
    });
    g += `<line x1="${ml}" x2="${(ml + pw).toFixed(1)}" y1="12" y2="12" stroke="#C8C8C8" stroke-width="1.5"/>`;
    for (const [a, b, l] of qd) g += `<rect x="${sx(a).toFixed(1)}" y="7" width="${((b - a + 1) * sc).toFixed(1)}" height="10" rx="1.5" fill="#fff" stroke="${PAGE.pfam_color}"><title>${esc(l)} ${a}-${b}</title></rect>`;
    const tx = ml + pw + 6;
    g += `<line x1="${tx + 6}" x2="${tx + 6}" y1="${mt}" y2="${(mt + ph).toFixed(1)}" stroke="#C8C8C8" stroke-width="1.5"/>`;
    for (const [a, b, l] of h.td || []) {
      g += `<rect x="${tx + 1}" y="${sy(a).toFixed(1)}" width="10" height="${((b - a + 1) * sc).toFixed(1)}" rx="1.5" fill="#fff" stroke="${PAGE.pfam_color}"><title>${esc(l)} ${a}-${b}</title></rect>`;
      if ((b - a + 1) * sc > 9) g += `<text x="${tx + 14}" y="${(sy((a + b) / 2) + 3).toFixed(1)}" class="dl">${esc(l.split(" ")[0])}</text>`;
    }
    const same = H.filter((x) => (h.t ? x.t === h.t : x.k === h.k && x.l === h.l) && x.r);
    const seg = (x, w, op) => `<line x1="${sx(x.q[0]).toFixed(1)}" y1="${sy(Math.min(x.r[0], x.r[1])).toFixed(1)}" x2="${sx(x.q[1]).toFixed(1)}" y2="${sy(Math.max(x.r[0], x.r[1])).toFixed(1)}" stroke="${PAGE.tracks[x.k].color}" stroke-width="${w}" stroke-linecap="round" opacity="${op}"/>`;
    for (const x of same) if (x !== h) g += seg(x, 1.5, 0.55);
    g += seg(h, 4, 1);
    g += `<text x="${ml}" y="${(mt + ph + 12).toFixed(1)}" class="ax">1</text><text x="${(ml + pw).toFixed(1)}" y="${(mt + ph + 12).toFixed(1)}" class="ax" text-anchor="end">${qL}</text>`;
    g += `<text x="${ml - 4}" y="${mt + 8}" class="ax" text-anchor="end">1</text><text x="${ml - 4}" y="${(mt + ph).toFixed(1)}" class="ax" text-anchor="end">${tlen}</text>`;
    g += `<text x="${(ml + pw / 2).toFixed(1)}" y="${(mt + ph + 26).toFixed(1)}" class="ax" text-anchor="middle">${esc(PAGE.id)} (aa)</text>`;
    g += `<text x="10" y="${(mt + ph / 2).toFixed(1)}" class="ax" text-anchor="middle" transform="rotate(-90 10 ${(mt + ph / 2).toFixed(1)})">${esc(tname)} (aa)</text>`;
    // the lines take the colour of their kind of program; the key names the programs of each colour
    const kinds = {};
    for (const x of same) { const tr = PAGE.tracks[x.k]; (kinds[tr.color] = kinds[tr.color] || new Set()).add(tr.name); }
    const key = Object.entries(kinds).map(([c, ns]) => `<span><i style="background:${c}"></i>${esc([...ns].join(", "))}</span>`).join("");
    return `<div class="map"><p class="mapkey">Each line is one hit on ${esc(tname)}, in the colour of its kind of program; the thick one is shown on the right. Boxes on the edges: Pfam domains.</p>` +
      `<div class="mapkey tools">${key}</div><svg width="${W.toFixed(0)}" height="${Hh.toFixed(0)}">${g}</svg></div>`;
  }

  function open(h) {
    tip.hidden = true;
    pbody.innerHTML = "<p>Loading the residues...</p>"; panel.hidden = false;
    Promise.all([alignments(), h.t ? targetSeq(h.t) : Promise.resolve(null)]).then(([alns, tseq]) => show(h, alns[H.indexOf(h)], tseq));
  }

  function show(h, a, tseq) {
    const q = PAGE.q, tool = PAGE.tracks[h.k].name;
    const tname = h.k === "hmmscan" ? `${h.l} model` : h.l;
    const names = PAGE.names[h.k] || [];
    let html = `<h3>${esc(tool)}: ${esc(h.l)}` + (h.t && h.t !== h.l ? ` <span class="w">Swiss-Prot ${esc(h.t)}</span>` : "") + "</h3>" + about(h);
    const scoreBits = h.s.map((v, i) => [names[i], v]).filter(([, v]) => v !== "not in the file" && !/^no .* column/.test(v)).map(([n, v]) => `${esc(n || "")} <b>${esc(v)}</b>`).join(" · ");
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
    if (h.k !== "hhsearch" && h.t && !tseq) html += `<p class="note">${esc(h.t)} is not in swissprot.fasta, so its residues cannot be shown.</p>`;
    let right = "";
    if (a && a.q) {
      const leg = [];
      if (ALPHA) {
        for (const c of Object.keys(groups).sort()) leg.push(`<span><i class="${c === "h" || c === "p" ? c : "o"}">${c}</i>${esc(CLASS_WORD[c] || "class " + c)} in ${esc(ALPHA_NAME)}: <code>${groups[c].sort().join("")}</code></span>`);
      }
      leg.push(`<span><b class="mono">A</b> bold: the same residue in both</span>`);
      leg.push(kmer ? `<span><b class="mono">|</b> same residue, <b class="mono">+</b> same class</span>` : `<span><b class="mono">|</b> same residue, <b class="mono">+</b> similar (BLOSUM62 above 0), <b class="mono">-</b> gap</span>`);
      const st = stats(flank ? { q: a.q.slice(flank[0], a.q.length - flank[1]), t: a.t.slice(flank[0], a.t.length - flank[1]) } : a, kmer);
      right += `<p class="asum"><b>${st.ident} of ${st.cols}</b> aligned columns identical (${pct(st.ident, st.cols)}%), ` +
        `${st.sim} more ${kmer ? `in the same ${esc(ALPHA_NAME)} class` : "similar (BLOSUM62 score above 0)"}, ` +
        `${st.gaps} gap column${st.gaps === 1 ? "" : "s"}${flank ? " (the faded flanks not counted)" : ""}</p>`;
      right += `<p class="asrc">Source: ${esc(src.replace(/\.$/, ""))}.</p>`;
      right += `<div class="alegend">${leg.join("")}</div>`;
      let qi = a.qs, ti = a.ts, body = "";
      for (let st = 0; st < a.q.length; st += 60) {
        const c = chunkHtml(a, st, a.q.length, qi, ti, tname, kmer, flank);
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
          (k ? `<div>${esc(CLASS_WORD[k] || "class")} (${esc(k)}) in ${esc(ALPHA_NAME)}: <span class="mono">${groups[k].sort().join("")}</span></div>` : "");
        tip.hidden = false; place(e);
      });
      c.addEventListener("mouseleave", () => { tip.hidden = true; });
    });
  }
})();
