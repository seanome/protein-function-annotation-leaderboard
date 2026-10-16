/* Region pile-up. Rebuilds each panel's matrix from the runs in data/pileup.json (built by
   scripts/pileup.py), draws it, and computes every readout and the caption from that same
   matrix with one function, readouts(). Cell codes: 0 not covered, 1 covered, 255 past the
   protein's end (null). */
(function () {
  "use strict";
  const NCOL = 300, FL = 100;
  const NULL = 255;
  const META = JSON.parse(document.getElementById("pu-meta").textContent);
  const css = getComputedStyle(document.documentElement);
  const tok = (n) => css.getPropertyValue(n).trim();
  // A kind's colour from the CSS tokens at draw time, so the canvases, curves and legends
  // agree in light and dark mode.
  const kindColor = (kind) => tok("--k-" + kind) || META.kinds[kind];
  let D, R, P, N, ALL_ROWS, FAMILY_ROWS, COV = [], state = { sort: "coverage", focus: 0, color: "coverage", all: false };

  // ---------- matrix ----------
  function isNull(r, j) { return r.c0 === null || j < r.c0 || j > r.c1; }
  function buildCov(p) {
    const m = new Uint8Array(N * NCOL);
    for (let i = 0; i < N; i++) {
      const r = R[i], base = i * NCOL;
      for (let j = 0; j < NCOL; j++) if (isNull(r, j)) m[base + j] = NULL;
      for (const [a, b] of p.rows[i].r) for (let j = a; j <= b; j++) m[base + j] = 1;
    }
    return m;
  }
  function buildScore(p) {
    const m = new Float32Array(N * NCOL).fill(NaN);
    for (let i = 0; i < N; i++) for (const [a, b, v] of p.rows[i].r) for (let j = a; j <= b; j++) m[i * NCOL + j] = v;
    return m;
  }

  // ---------- readouts: the one function every number on the page comes from ----------
  // rows: the row indices to use (all rows by default, or FAMILY_ROWS for the profile
  // capped at one row per Pfam family, since rows from one protein or family are not
  // independent).
  function readouts(m, rows) {
    const curve = new Float64Array(NCOL), nn = new Float64Array(NCOL);
    let rowsScored = 0, missed = 0, bodySum = 0, flankCov = 0, flankN = 0;
    for (const i of rows || ALL_ROWS) {
      const base = i * NCOL;
      let bodyCov = 0, bodyN = 0;
      for (let j = 0; j < NCOL; j++) {
        const v = m[base + j];
        if (v === NULL) continue;
        nn[j]++; curve[j] += v;
        if (j >= FL && j < FL + 100) { bodyN++; bodyCov += v; } else { flankN++; flankCov += v; }
      }
      if (bodyN) { rowsScored++; bodySum += bodyCov / bodyN; if (bodyCov === 0) missed++; }
    }
    for (let j = 0; j < NCOL; j++) curve[j] = nn[j] ? curve[j] / nn[j] : NaN;
    const bodyCurve = Array.from(curve.slice(FL + 10, FL + 90)).sort((a, b) => a - b);
    const plateau = bodyCurve[Math.floor(bodyCurve.length / 2)];
    const colWidth = (j) => (j >= FL && j < FL + 100 ? META.median_len / 100 : 1);
    function width(side, hi, lo) {
      // Residues between the outermost column at or above `hi` and the first column beyond it
      // at or below `lo`, moving out of the region from its centre: leftwards at the start,
      // rightwards at the end. Fewer residues means a sharper edge at the region's boundary.
      const cols = side === "start" ? range(FL + 50, -1, -1) : range(FL + 50, NCOL, 1);
      let j80 = null, j20 = null;
      for (const j of cols) {
        if (isNaN(curve[j])) continue;
        if (curve[j] >= hi) { j80 = j; continue; }
        if (curve[j] <= lo) { if (j80 !== null) j20 = j; break; }
      }
      if (j80 === null) return { text: "not reached", value: null };
      if (j20 === null) return { text: "does not fall that far", value: null };
      let w = 0;
      for (let j = Math.min(j80, j20); j < Math.max(j80, j20); j++) w += colWidth(j);
      return { text: w.toFixed(0) + " aa", value: w };
    }
    return {
      curve, rowsScored, missed, missedFrac: missed / rowsScored,
      bodyMean: bodySum / rowsScored, flankMean: flankCov / flankN, plateau,
      startAbs: width("start", 0.8, 0.2), endAbs: width("end", 0.8, 0.2),
      startRel: width("start", 0.8 * plateau, 0.2 * plateau), endRel: width("end", 0.8 * plateau, 0.2 * plateau),
    };
  }
  function range(a, b, s) { const o = []; for (let i = a; s > 0 ? i < b : i > b; i += s) o.push(i); return o; }
  function readoutText(ro) {
    return `${ro.missed.toLocaleString("en-US")} of ${ro.rowsScored.toLocaleString("en-US")} domains missed (${(100 * ro.missedFrac).toFixed(1)}%). ` +
      `Mean coverage ${ro.bodyMean.toFixed(3)} inside the domain and ${ro.flankMean.toFixed(3)} in the flanks. ` +
      `Edge sharpness, residues for the curve to fall from 0.8 to 0.2: start ${ro.startAbs.text}, end ${ro.endAbs.text}. ` +
      `From 80% to 20% of the curve's own height inside the domain (${ro.plateau.toFixed(2)}): start ${ro.startRel.text}, end ${ro.endRel.text}`;
  }
  // The edge in a small panel: a number of residues, or "none" (the panel below says why).
  const edgeShort = (w) => (w.value == null ? "none" : w.text);

  // ---------- order ----------
  function order() {
    const idx = range(0, N, 1);
    const key = {
      coverage: (i) => { let c = 0, n = 0; const m = COV[state.focus]; for (let j = FL; j < FL + 100; j++) { const v = m[i * NCOL + j]; if (v !== NULL) { n++; c += v; } } return n ? c / n : -1; },
      identity: (i) => (R[i].pid === null ? -1 : R[i].pid),
      length: (i) => R[i].len,
      plddt: (i) => (R[i].plddt === null ? -1 : R[i].plddt),
    }[state.sort];
    const k = idx.map(key);
    return idx.sort((a, b) => k[b] - k[a] || (R[a].acc < R[b].acc ? -1 : R[a].acc > R[b].acc ? 1 : R[a].s - R[b].s));
  }

  // ---------- drawing ----------
  function hexToRgb(h) { h = h.replace("#", ""); return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16)); }
  function mix(a, b, t) { return a.map((x, i) => Math.round(x + (b[i] - x) * t)); }
  function colors(p) {
    return {
      on: hexToRgb(kindColor(p.kind)), off: hexToRgb(tok("--surface") || "#ffffff"),
      nul: hexToRgb(tok("--null") || "#b9bdbb"), stripA: hexToRgb(tok("--strip-a") || "#8c6d46"), stripB: hexToRgb(tok("--strip-b") || "#dccbab"),
    };
  }
  function drawOverview(canvas, m, p, ord) {
    const H = canvas.height, W = NCOL, ctx = canvas.getContext("2d"), img = ctx.createImageData(W, H), c = colors(p);
    for (let y = 0; y < H; y++) {
      const r0 = Math.floor((y * N) / H), r1 = Math.max(r0 + 1, Math.floor(((y + 1) * N) / H));
      for (let j = 0; j < W; j++) {
        let cov = 0, n = 0;
        for (let k = r0; k < r1; k++) { const v = m[ord[k] * NCOL + j]; if (v !== NULL) { n++; cov += v; } }
        const rgb = n === 0 ? c.nul : mix(c.off, c.on, cov / n);
        const o = (y * W + j) * 4; img.data[o] = rgb[0]; img.data[o + 1] = rgb[1]; img.data[o + 2] = rgb[2]; img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
  }
  function drawFocus(canvas, m, s, p, ord, scoreRange) {
    // One pixel row per region, two pixel columns per matrix column, 8 px protein strip.
    const SW = STRIP - 2, W = STRIP + NCOL * 2, ctx = canvas.getContext("2d"), img = ctx.createImageData(W, N), c = colors(p);
    let tone = 0;
    for (let y = 0; y < N; y++) {
      const i = ord[y];
      if (y > 0 && R[i].acc !== R[ord[y - 1]].acc) tone ^= 1;
      const set = (x, rgb) => { const o = (y * W + x) * 4; img.data[o] = rgb[0]; img.data[o + 1] = rgb[1]; img.data[o + 2] = rgb[2]; img.data[o + 3] = 255; };
      for (let x = 0; x < SW; x++) set(x, tone ? c.stripB : c.stripA);
      set(SW, c.off); set(SW + 1, c.off);
      for (let j = 0; j < NCOL; j++) {
        const v = m[i * NCOL + j];
        let rgb;
        if (v === NULL) rgb = c.nul;
        else if (!v) rgb = c.off;
        else if (state.color === "score") {
          const t = Math.min(1, Math.max(0, (s[i * NCOL + j] - scoreRange[0]) / (scoreRange[1] - scoreRange[0] || 1)));
          rgb = mix(c.off, c.on, 0.4 + 0.6 * t);
        } else rgb = c.on;
        set(SW + 2 + 2 * j, rgb); set(SW + 3 + 2 * j, rgb);
      }
    }
    ctx.putImageData(img, 0, 0);
  }
  const SVGNS = "http://www.w3.org/2000/svg";
  function el(n, a, p) { const e = document.createElementNS(SVGNS, n); for (const k in a) e.setAttribute(k, a[k]); if (p) p.appendChild(e); return e; }
  // Geometry shared by the full-size curve and the heatmap under it: the heatmap canvas is
  // STRIP + 2 px per column wide and sits FOCUS_L units in from the left, so column j has
  // the same x in both.
  const STRIP = 10, FOCUS_L = 46, FOCUS_W = FOCUS_L + STRIP + 2 * NCOL;
  function drawCurve(svg, curve, color, big, capped) {
    svg.innerHTML = "";
    // A small panel is stretched to its box (90 px tall) and carries no text: its numbers
    // are printed under it, and its lines keep their width (non-scaling strokes).
    const H = big ? 150 : 90, T = big ? 6 : 3, B = big ? 22 : 3;
    const L = big ? FOCUS_L + STRIP : 0, CW = big ? 2 : 1, W = big ? FOCUS_W : NCOL;
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    if (!big) svg.setAttribute("preserveAspectRatio", "none");
    const ns = big ? {} : { "vector-effect": "non-scaling-stroke" };
    const sx = (j) => L + (j + 0.5) * CW, edge = (j) => L + j * CW, sy = (v) => T + (1 - v) * (H - T - B);
    for (const v of [0, 0.5, 1]) {
      el("line", Object.assign({ x1: L, x2: L + NCOL * CW, y1: sy(v), y2: sy(v), class: "grid" }, ns), svg);
      if (big) el("text", { x: L - 4, y: sy(v) + 4, "text-anchor": "end", class: "tk" }, svg).textContent = v;
    }
    for (const [j, lab] of [[0, "−100 aa"], [FL, "start"], [FL + 100, "end"], [NCOL - 1, "+100 aa"]]) {
      const x = j === FL || j === FL + 100 ? edge(j) : sx(j);
      if (!big && !(j === FL || j === FL + 100)) continue;
      el("line", Object.assign({ x1: x, x2: x, y1: T, y2: H - B, class: j === FL || j === FL + 100 ? "edge" : "grid" }, ns), svg);
      if (big) el("text", { x, y: H - 4, "text-anchor": j === 0 ? "start" : j === NCOL - 1 ? "end" : "middle", class: "tk" }, svg).textContent = lab;
    }
    let d = "";
    for (let j = 0; j < NCOL; j++) if (!isNaN(curve[j])) d += (d ? "L" : "M") + sx(j).toFixed(1) + "," + sy(curve[j]).toFixed(1);
    el("path", Object.assign({ d, fill: "none", stroke: color, "stroke-width": big ? 2 : 2 }, ns), svg);
    if (capped) {
      let d2 = "";
      for (let j = 0; j < NCOL; j++) if (!isNaN(capped[j])) d2 += (d2 ? "L" : "M") + sx(j).toFixed(1) + "," + sy(capped[j]).toFixed(1);
      el("path", Object.assign({ d: d2, fill: "none", stroke: tok("--ink") || "#1c211f", "stroke-width": big ? 1.2 : 1, class: "pu-cap" }, ns), svg);
    }
    if (big) el("text", { x: 10, y: (T + H - B) / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${(T + H - B) / 2})`, class: "lab" }, svg).textContent = "share covered";
  }

  // ---------- page ----------
  function panelTitle(p) {
    const ref = p.tool === "hmmscan" ? ' <span class="tag warn">reference: graded on its own models</span>' : "";
    const k = p.tool + "|" + p.setting;
    return `${META.names[p.tool]}, <span title="${p.setting}: ${META.words[k] || ""}">${META.labels[k] || p.setting}</span>${ref}`;
  }
  function render() {
    const ord = order();
    const grid = document.getElementById("pu-grid");
    grid.innerHTML = "";
    const shown = new Set(META.default_panels);
    if (!shown.has(P[state.focus].tool + "|" + P[state.focus].setting)) shown.add(P[state.focus].tool + "|" + P[state.focus].setting);
    P.forEach((p, k) => {
      if (!state.all && !shown.has(p.tool + "|" + p.setting)) return;
      const ro = readouts(COV[k]);
      const card = document.createElement("button");
      card.type = "button";
      card.className = "pu-panel" + (k === state.focus ? " on" : "");
      card.setAttribute("aria-pressed", k === state.focus);
      card.innerHTML = `<span class="pu-t"><span class="kind k-${p.kind}"></span>${panelTitle(p)}</span>`;
      const svg = el("svg", { class: "pu-curve", role: "img", "aria-label": "column average" });
      card.appendChild(svg);
      drawCurve(svg, ro.curve, kindColor(p.kind), false, readouts(COV[k], FAMILY_ROWS).curve);
      const cv = document.createElement("canvas");
      cv.width = NCOL; cv.height = 150; cv.className = "pu-ov";
      card.appendChild(cv);
      drawOverview(cv, COV[k], p, ord);
      const r = document.createElement("span");
      r.className = "pu-ro";
      r.dataset.missed = ro.missedFrac; r.dataset.body = ro.bodyMean; r.dataset.flank = ro.flankMean;
      r.innerHTML = `missed ${(100 * ro.missedFrac).toFixed(0)}% · inside ${ro.bodyMean.toFixed(2)} · flanks ${ro.flankMean.toFixed(2)}<br>` +
        `edge sharpness: start ${edgeShort(ro.startAbs)}, end ${edgeShort(ro.endAbs)}`;
      card.appendChild(r);
      card.addEventListener("click", () => { state.focus = k; render(); });
      grid.appendChild(card);
    });
    // focus panel
    const p = P[state.focus], m = COV[state.focus], s = buildScore(p), ro = readouts(m), rf = readouts(m, FAMILY_ROWS);
    const vals = []; for (let i = 0; i < s.length; i++) if (!isNaN(s[i])) vals.push(s[i]);
    vals.sort((a, b) => a - b);
    const sr = vals.length ? [vals[Math.floor(0.05 * (vals.length - 1))], vals[Math.floor(0.95 * (vals.length - 1))]] : [0, 1];
    document.getElementById("pu-ftitle").innerHTML = panelTitle(p);
    drawCurve(document.getElementById("pu-fcurve"), ro.curve, kindColor(p.kind), true, rf.curve);
    const fc = document.getElementById("pu-fheat");
    fc.width = STRIP + NCOL * 2; fc.height = N;
    fc.parentElement.style.marginLeft = (100 * FOCUS_L / FOCUS_W) + "%";
    drawFocus(fc, m, s, p, ord, sr);
    fc.dataset.order = ""; window.__pileup.order = ord;
    const famText = `With one domain per Pfam family (${FAMILY_ROWS.length.toLocaleString("en-US")} families): ${readoutText(rf)}.`;
    document.getElementById("pu-fro").textContent = "All domains: " + readoutText(ro) + ". " + famText;
    document.getElementById("pu-caption-numbers").textContent =
      `It misses ${(100 * ro.missedFrac).toFixed(0)}% of the domains and covers ${(100 * ro.bodyMean).toFixed(0)}% inside them on average, ` +
      `and ${(100 * ro.flankMean).toFixed(0)}% of the 100 residues either side.`;
    const allb = document.getElementById("pu-all");
    if (allb) { allb.setAttribute("aria-pressed", state.all); allb.textContent = state.all ? "One panel per tool" : `Show all ${P.length} settings`; }
    const sb = document.getElementById("pu-scorebar");
    sb.hidden = state.color !== "score";
    sb.querySelector(".name").textContent = `${META.names[p.tool]}: ${META.scores[p.tool]}, higher is better.`;
    sb.querySelector(".lo").textContent = sr[0].toPrecision(3);
    sb.querySelector(".hi").textContent = sr[1].toPrecision(3);
    sb.style.setProperty("--kind", kindColor(p.kind));
    document.querySelectorAll("button[data-sort]").forEach((b) => b.setAttribute("aria-pressed", b.dataset.sort === state.sort));
    document.querySelectorAll("button[data-color]").forEach((b) => b.setAttribute("aria-pressed", b.dataset.color === state.color));
    window.__pileup.readouts = P.map((_, k) => readouts(COV[k]));
    document.getElementById("pu-tip").hidden = true;
  }

  // ---------- hover ----------
  function hoverText(i) {
    const r = R[i], p = P[state.focus], row = p.rows[i];
    const calls = row.c.length ? row.c.slice(0, 6).map(([a, b]) => `${a}–${b}`).join(", ") + (row.c.length > 6 ? ", …" : "") : "no call";
    // Residues of the true region (up to 120), with a mark under each one a call covers.
    let res = "";
    if (r.seq) {
      const a = r.s, b = Math.min(r.e, r.plen, r.s + 119);
      const covered = (x) => row.c.some(([u, v]) => x >= u && x <= v);
      for (let x0 = a; x0 <= b; x0 += 60) {
        const x1 = Math.min(b, x0 + 59);
        const aa = r.seq.slice(x0 - r.off, x1 - r.off + 1);
        let mk = ""; for (let x = x0; x <= x1; x++) mk += covered(x) ? "|" : " ";
        res += `${String(x0).padStart(5)} ${aa}\n      ${mk}\n`;
      }
      if (r.e > b) res += `      … ${r.e - b} more residues\n`;
    }
    return `<b>${r.sym || r.acc}</b> <span class="mono">${r.acc}</span> · Pfam ${r.fam}<br>` +
      `Pfam domain at residues ${r.s}–${r.e} (${r.len} aa) of a ${r.plen} aa protein${r.past_end ? "; it runs past the end of this sequence" : ""}<br>` +
      `Identity to the closest ${r.fam} domain in zebrafish: ${r.pid === null ? "zebrafish has none" : r.pid + "%"} · mean pLDDT over the domain: ${r.plddt === null ? "no matching AlphaFold model" : r.plddt}<br>` +
      `${META.names[p.tool]}: ${row.n} call${row.n === 1 ? "" : "s"} of ${r.fam} within 100 residues of the domain. Residues they cover, merged: ${calls}<br>` +
      `<pre class="pu-res">${res}</pre>` +
      `<span class="muted">Top line: the human residues of the domain. Under each, | marks a residue inside a call; this is not an alignment. The zebrafish residues each call matched are not stored in the run's call tables.</span>`;
  }

  // ---------- histogram: one small panel per tool ----------
  // Each panel is drawn at its box's own width, so its 12 px text stays 12 px on a phone.
  function drawHist() {
    const box = document.getElementById("pu-hist"); box.innerHTML = "";
    const edges = [];
    for (let x = -100; x <= 100; x += 5) edges.push(x);
    const nb = edges.length + 1; // underflow, 40 bins of 5 aa, overflow
    const series = META.hist.map((h) => {
      const p = P.find((q) => q.tool === h.tool && q.setting === h.setting);
      const cnt = new Array(nb).fill(0); let n = 0;
      for (const row of p.rows) {
        if (row.d === null) continue; n++;
        const d = row.d; let b = d < -100 ? 0 : d >= 100 ? nb - 1 : 1 + Math.floor((d + 100) / 5);
        cnt[b]++;
      }
      return { h, p, frac: cnt.map((c) => c / Math.max(1, n)), n };
    });
    // One vertical scale for every tool that searched a relative. hmmscan's calls almost all
    // start on the domain (its peak is near 1), so it gets its own scale and is left out of
    // the grey lines in the other panels.
    const own = (s) => s.h.tool === "hmmscan";
    const scaleOf = (list) => {
      const top = Math.max(...list.flatMap((s) => s.frac));
      const ystep = top > 0.5 ? 0.2 : top > 0.2 ? 0.1 : 0.05;
      return [Math.ceil(top / ystep) * ystep, ystep];
    };
    const shared = scaleOf(series.filter((s) => !own(s)));
    for (const s of series) {
      const [ymax, ystep] = own(s) ? scaleOf([s]) : shared;
      const cell = document.createElement("div");
      cell.className = "pu-hcell";
      box.appendChild(cell);
      const W = Math.max(150, Math.round(cell.clientWidth || 240)), H = 150, L = 34, Rr = 6, T = 22, B = 22;
      const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
        "aria-label": `${META.names[s.h.tool]}: call start minus domain start` }, cell);
      const bw = (W - L - Rr) / nb, sx = (b) => L + b * bw, sy = (v) => T + (1 - v / ymax) * (H - T - B);
      el("text", { x: 0, y: 14, class: "lab" }, svg).textContent = `${META.names[s.h.tool]} (n = ${s.n.toLocaleString("en-US")})${own(s) ? ", own scale" : ""}`;
      for (let v = 0; v <= ymax + 1e-9; v += ystep) {
        el("line", { x1: L, x2: W - Rr, y1: sy(v), y2: sy(v), class: "grid" }, svg);
        el("text", { x: L - 4, y: sy(v) + 4, "text-anchor": "end", class: "tk" }, svg).textContent = v.toFixed(ystep < 0.1 ? 2 : 1);
      }
      for (const [b, lab] of [[0, "< −100"], [21, "0"], [nb - 1, "> +100"]]) {
        const x = b === 0 || b === nb - 1 ? sx(b) + bw / 2 : sx(b);
        el("text", { x, y: H - 6, "text-anchor": b === 0 ? "start" : b === nb - 1 ? "end" : "middle", class: "tk" }, svg).textContent = lab;
        if (b === 21) el("line", { x1: x, x2: x, y1: T, y2: H - B, class: "edge" }, svg);
      }
      const path = (frac) => frac.map((v, b) => (b ? "L" : "M") + sx(b).toFixed(1) + "," + sy(v).toFixed(1) + "H" + (sx(b) + bw).toFixed(1)).join("");
      for (const o of series) if (o !== s && !own(o) && !own(s)) el("path", { d: path(o.frac), fill: "none", stroke: tok("--hair") || "#d9ddda", "stroke-width": 1.2, class: "pu-ctx" }, svg);
      el("path", { d: path(s.frac), fill: "none", stroke: kindColor(s.p.kind), "stroke-width": 2, class: "pu-hs", "data-tool": s.h.tool }, svg);
    }
    window.__pileup.hist = series.map((s) => ({ tool: s.h.tool, n: s.n }));
  }

  // ---------- start ----------
  window.__pileup = {};
  fetch("data/pileup.json").then((r) => r.json()).then((d) => {
    D = d; R = d.regions; P = d.panels; N = R.length;
    ALL_ROWS = range(0, N, 1);
    // One row per family: the first by accession, then start.
    const seen = new Set();
    FAMILY_ROWS = ALL_ROWS.slice().sort((a, b) => (R[a].acc < R[b].acc ? -1 : R[a].acc > R[b].acc ? 1 : R[a].s - R[b].s))
      .filter((i) => (seen.has(R[i].fam) ? false : (seen.add(R[i].fam), true)));
    P.forEach((p) => COV.push(buildCov(p)));
    state.focus = Math.max(0, P.findIndex((p) => p.tool === META.focus_default));
    window.__pileup.cov = COV; window.__pileup.regions = R; window.__pileup.panels = P.map((p) => p.tool + "|" + p.setting);
    window.__pileup.setFocus = (k) => { state.focus = k; render(); };
    window.__pileup.setSort = (s) => { state.sort = s; render(); };
    window.__pileup.hoverText = (i) => hoverText(i);
    window.__pileup.readoutsOf = (k, fam) => { const ro = readouts(COV[k], fam ? FAMILY_ROWS : null); return { missedFrac: ro.missedFrac, bodyMean: ro.bodyMean, flankMean: ro.flankMean, rowsScored: ro.rowsScored,
      startAbs: ro.startAbs.text, endAbs: ro.endAbs.text, startRel: ro.startRel.text, endRel: ro.endRel.text }; };
    document.querySelectorAll("button[data-sort]").forEach((b) => b.addEventListener("click", () => { state.sort = b.dataset.sort; render(); }));
    document.querySelectorAll("button[data-color]").forEach((b) => b.addEventListener("click", () => { state.color = b.dataset.color; render(); }));
    const allb = document.getElementById("pu-all");
    if (allb) allb.addEventListener("click", () => { state.all = !state.all; render(); });
    let pending = null;
    window.addEventListener("resize", () => { clearTimeout(pending); pending = setTimeout(drawHist, 150); });
    const fc = document.getElementById("pu-fheat"), tip = document.getElementById("pu-tip");
    const show = (ev) => {
      const rect = fc.getBoundingClientRect();
      const y = Math.floor(((ev.clientY - rect.top) / rect.height) * N);
      if (y < 0 || y >= N) return;
      tip.innerHTML = hoverText(window.__pileup.order[y]);
      tip.hidden = false;
      const box = fc.parentElement.getBoundingClientRect();
      tip.style.top = Math.min(ev.clientY - box.top + fc.parentElement.scrollTop + 14, fc.parentElement.scrollHeight - 10) + "px";
    };
    fc.addEventListener("mousemove", show);
    fc.addEventListener("mouseleave", () => (tip.hidden = true));
    render();
    drawHist();
    // Redraw when the colour theme changes, so canvases and curves take the new kind colours.
    const redraw = () => { render(); drawHist(); };
    if (window.matchMedia) {
      const mq = window.matchMedia("(prefers-color-scheme: dark)");
      if (mq.addEventListener) mq.addEventListener("change", redraw);
    }
    new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
    window.__pileup.redraw = redraw;
    document.getElementById("pu-loading").hidden = true;
    window.__pileup.ready = true;
  });
})();
