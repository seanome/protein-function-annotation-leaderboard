// Tooltips, the "show all" redraw and the residue panel for one protein page.
// PAGE is set by the page: the protein sequence (q), the drawing scale (px, w), the hits
// (k: track, l: label, q: range on this protein, s: score values, t: target, r: range on the
// target, a: alignment if made at build time, m: kmerseek alphabet), score names per track,
// kmerseek's class tables (alpha) and the number of Swiss-Prot sequence files (shards).
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
    return `<div>${esc(PAGE.tracks[h.k].name)}</div><b>${esc(h.l)}</b><table>${rows}</table>` +
      (h.t ? "<div>Click to see the residues.</div>" : "");
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
      g.addEventListener("click", () => { tip.hidden = true; open(h); });
      g.addEventListener("keydown", (e) => { if (e.key === "Enter") open(h); });
    });
  }
  wire(document);

  // "Show all": redraw each tool's row with all of its hits, with the same lanes, boxes and
  // labels build.py draws for the best few, then move the rows below it down.
  const LANE_H = 18, BOX_H = 14, ROW_PAD = 8;
  function boxSvg(h, i, lane) {
    const tr = PAGE.tracks[h.k];
    const x0 = (h.q[0] - 1) * PAGE.px, w = Math.max(2, h.q[1] * PAGE.px - x0), y = ROW_PAD + lane * LANE_H;
    const rect = tr.open
      ? `<rect x="${(x0 + 1).toFixed(1)}" y="${y + 1}" width="${Math.max(1, w - 2).toFixed(1)}" height="${BOX_H - 2}" rx="2" fill="#fff" stroke="${tr.color}" stroke-width="2"/>`
      : `<rect x="${x0.toFixed(1)}" y="${y}" width="${w.toFixed(1)}" height="${BOX_H}" rx="2" fill="${tr.color}"/>`;
    let shown = h.l;
    if (h.l.length * 11 * 0.56 + 6 > w) { const n = Math.floor((w - 6) / (11 * 0.56)) - 1; shown = n >= 3 ? h.l.slice(0, n) + "…" : ""; }
    const lab = shown ? `<text x="${(x0 + 3).toFixed(1)}" y="${y + BOX_H - 3}" fill="${tr.open ? "#1a1a1a" : "#ffffff"}" class="bl">${esc(shown)}</text>` : "";
    return `<g class="hit" tabindex="0" data-i="${i}">${rect}${lab}</g>`;
  }
  const tracksDiv = document.getElementById("tracks");
  const original = tracksDiv.innerHTML;
  let allDrawn = null;
  function drawAll() {
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
      if (idx && idx.length > PAGE.top) {
        const lanes = [], parts = [];
        for (const i of idx) {
          const h = H[i];
          let li = lanes.findIndex((ends) => ends.every(([a, b]) => h.q[0] > b || h.q[1] < a));
          if (li < 0) { lanes.push([]); li = lanes.length - 1; }
          lanes[li].push(h.q);
          parts.push(boxSvg(h, i, li));
        }
        hgt = ROW_PAD + lanes.length * LANE_H + 4;
        const sep = rg.querySelector("line.sep");
        const line = `<line x1="0" y1="${ROW_PAD + BOX_H / 2}" x2="${PAGE.w}" y2="${ROW_PAD + BOX_H / 2}" stroke="#C8C8C8" stroke-width="1.5"/>`;
        rg.innerHTML = (sep ? sep.outerHTML : "") + line + parts.join("");
        const more = lg.querySelector(".more");
        if (more) more.textContent = `all ${idx.length}`;
      }
      lg.setAttribute("transform", `translate(0,${y})`);
      rg.setAttribute("transform", `translate(0,${y})`);
      y += hgt;
    });
    for (const s of [ls, rs]) { s.setAttribute("height", y); s.setAttribute("viewBox", `0 0 ${s === ls ? "100%" : PAGE.w} ${y}`); }
    ls.removeAttribute("viewBox");
    return tmp.innerHTML;
  }
  const tg = document.getElementById("toggle");
  if (tg) {
    const label = tg.textContent;
    let on = false;
    tg.addEventListener("click", () => {
      on = !on;
      if (on && allDrawn === null) allDrawn = drawAll();
      tracksDiv.innerHTML = on ? allDrawn : original;
      wire(tracksDiv);
      tg.textContent = on ? `Show the best ${PAGE.top} hits per tool` : label;
    });
  }
  document.getElementById("close").addEventListener("click", () => { panel.hidden = true; });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") panel.hidden = true; });

  // Swiss-Prot sequences live in ../targets/NNN.js, one file per group of accessions.
  function shardOf(acc) {
    let h = 0;
    for (const ch of acc) h = (Math.imul(h, 31) + ch.charCodeAt(0)) >>> 0;
    return h % PAGE.shards;
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

  function blocks(q, t, mid, qs, ts, qname, tname) {
    const width = 60, lw = Math.max(qname.length, tname.length);
    let out = [], qi = qs, ti = ts;
    for (let s = 0; s < q.length; s += width) {
      const qq = q.slice(s, s + width), tt = t.slice(s, s + width), mm = mid.slice(s, s + width);
      const qn = qq.replace(/-/g, "").length, tn = tt.replace(/-/g, "").length;
      const qe = qi + qn - 1, te = ti + tn - 1;
      const pad = (x) => x.padEnd(lw);
      out.push(`${pad(qname)} ${String(qi).padStart(6)} ${qq} ${qe}`);
      out.push(`${pad("")} ${"".padStart(6)} ${mm}`);
      out.push(`${pad(tname)} ${String(ti).padStart(6)} ${tt} ${te}`);
      out.push("");
      qi = qe + 1; ti = te + 1;
    }
    return out.join("\n");
  }

  function open(h) {
    if (!h.t) return;
    if (h.k === "hhsearch") { show(h, null); return; }
    pbody.innerHTML = "<p>Loading the residues...</p>"; panel.hidden = false;
    targetSeq(h.t).then((tseq) => show(h, tseq));
  }
  function show(h, tseq) {
    const q = PAGE.q;
    let html = `<h3>${esc(PAGE.tracks[h.k].name)}: ${esc(h.l)}</h3>`;
    if (h.k !== "hhsearch" && !tseq) {
      html += `<p>${esc(h.t)} is not in swissprot.fasta, so its residues cannot be shown.</p>`;
    } else if (h.k === "kmerseek") {
      const al = PAGE.alpha[h.m];
      const qs = q.slice(h.q[0] - 1, h.q[1]), ts = tseq.slice(h.r[0] - 1, h.r[1]);
      let mid = "", same = 0, ident = 0;
      const n = Math.max(qs.length, ts.length);
      for (let i = 0; i < n; i++) {
        const a = qs[i] || " ", b = ts[i] || " ";
        if (a === b) { mid += "|"; ident++; same++; }
        else if (al && al[a] && al[a] === al[b]) { mid += "+"; same++; }
        else mid += " ";
      }
      html += `<p>${ident} of ${n} identical` + (al ? `; ${same} of ${n} in the same class of ${esc(h.m)}.` : ".") + "</p>";
      html += `<p class="note">kmerseek places the two segments side by side without gaps. Match line: | the same residue; + a different residue in the same class of the alphabet; blank: different classes. Class table read from kmerseek's source file src/rust/alphabets.rs.</p>`;
      html += `<div class="aln">${esc(blocks(qs, ts, mid, h.q[0], h.r[0], "this protein", h.t))}</div>`;
      if (al) {
        const groups = {};
        for (const [r, c] of Object.entries(al)) (groups[c] = groups[c] || []).push(r);
        const enc = (s) => Array.from(s, (r) => al[r] || "x").join("");
        const qe = enc(qs), te = enc(ts);
        let cm = "";
        for (let i = 0; i < n; i++) cm += qe[i] === te[i] && qe[i] !== "x" ? "|" : " ";
        html += `<p class="note">The same segments written as classes of ${esc(h.m)}: ${Object.entries(groups).map(([c, g]) => "<code>" + esc(c) + "</code> = <code>" + g.sort().join("") + "</code>").join(", ")}; x: a residue outside the 20 standard ones. Match line: | the same class.</p>`;
        html += `<div class="aln">${esc(blocks(qe, te, cm, h.q[0], h.r[0], "this protein", h.t))}</div>`;
      } else {
        html += `<p class="note">No class table for ${esc(h.m || "this alphabet")} (kmerseek_alphabets.rs not given, or the alphabet is not in it), so only identical residues are marked.</p>`;
      }
    } else {
      let a = h.a, how;
      if (h.k === "hhsearch") {
        how = "Alignment as printed in the .hhr file by hhsearch.";
      } else if (a) {
        how = "Alignment recomputed when the site was built, not taken from the tool's output: the two coordinate ranges were aligned again with Biopython's PairwiseAligner (local, BLOSUM62, gap open -11, gap extend -1).";
      } else {
        const lo = Math.min(h.r[0], h.r[1]), hi = Math.max(h.r[0], h.r[1]);
        a = localAlign(q.slice(h.q[0] - 1, h.q[1]), tseq.slice(lo - 1, hi));
        a.qs += h.q[0] - 1; a.ts += lo - 1;
        how = "Alignment recomputed in your browser, not taken from the tool's output: the two coordinate ranges were aligned with the scores build.py gives Biopython for the best 5 hits per tool (local, BLOSUM62, gap open -11, gap extend -1).";
      }
      if (!a || !a.q) {
        html += "<p>No alignment could be made from these residues.</p>";
      } else {
        let mid = "", ident = 0;
        for (let i = 0; i < a.q.length; i++) {
          const same = a.q[i] !== "-" && a.q[i].toUpperCase() === (a.t[i] || "").toUpperCase();
          mid += same ? "|" : " "; if (same) ident++;
        }
        html += `<p>${ident} of ${a.q.length} aligned positions identical.</p><p class="note">${how} Match line: | the same residue. A - in a row is a gap.</p>`;
        html += `<div class="aln">${esc(blocks(a.q, a.t, mid, a.qs, a.ts, "this protein", h.t))}</div>`;
      }
    }
    pbody.innerHTML = html;
    panel.hidden = false;
  }
})();
