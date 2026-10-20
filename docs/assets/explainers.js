/* "How the scores work": five small widgets on a made-up example (#toy-data, written by
   build_site.py). Every number shown is computed here from that example; nothing comes
   from the benchmark tables except the default core-hours and prices in Figure 8.
   tests/check_explainers.py recomputes the same numbers in Python and compares.

   Marks, the same in every figure:
     true domain           box on the protein line, --ex-true
     call that counts      filled bar, --ink
     call that does not    dashed grey outline, --null (below the cutoff, or merged away)
     correct / wrong       the glyphs ✓ / × beside the call
     the cutoff you set    vertical line, --accent
     highest F (Fmax)      open ring, --ink */
(function () {
  "use strict";
  const TOY = JSON.parse(document.getElementById("toy-data").textContent);
  const SVGNS = "http://www.w3.org/2000/svg";
  const f3 = (x) => (x == null ? "—" : x.toFixed(3));

  // ---------------------------------------------------------------- the maths
  function iou(a, b) {
    // Residues are counted inclusively: 10–109 is 100 residues.
    const inter = Math.max(0, Math.min(a.end, b.end) - Math.max(a.start, b.start) + 1);
    const union = (a.end - a.start + 1) + (b.end - b.start + 1) - inter;
    return { inter, union, iou: inter / union };
  }
  // A call is correct when it names a true domain's family, overlaps it with IoU >= 0.5,
  // and no higher-scoring call has already found that domain.
  function judge(protein, calls) {
    const found = new Set();
    const out = calls.slice().sort((a, b) => b.score - a.score).map((c) => {
      let hit = -1;
      protein.true.forEach((t, i) => {
        if (hit < 0 && !found.has(i) && t.family === c.family && iou(c, t).iou >= 0.5) hit = i;
      });
      if (hit >= 0) found.add(hit);
      return Object.assign({}, c, { correct: hit >= 0, domain: hit });
    });
    return out;
  }
  // Dedup: keep the highest-scoring call; drop any later call of the same family that
  // overlaps a kept one with IoU >= 0.5.
  function dedup(calls) {
    const keep = [], merged = {};
    for (const c of calls.slice().sort((a, b) => b.score - a.score)) {
      const into = keep.find((k) => k.family === c.family && iou(k, c).iou >= 0.5);
      if (into) merged[c.id] = into.id; else keep.push(c);
    }
    return { keep, merged };
  }
  // Precision and recall at one cutoff, protein by protein, then averaged (CAFA style).
  function atCutoff(judged, proteins, t) {
    const precs = [], recs = [];
    let kept = 0, correct = 0, found = 0, nTrue = 0;
    judged.forEach((calls, i) => {
      const k = calls.filter((c) => c.score >= t);
      const ok = k.filter((c) => c.correct).length;
      kept += k.length; correct += ok; found += ok; nTrue += proteins[i].true.length;
      if (k.length) precs.push(ok / k.length);
      recs.push(ok / proteins[i].true.length);
    });
    const P = precs.length ? precs.reduce((a, b) => a + b, 0) / precs.length : null;
    const R = recs.reduce((a, b) => a + b, 0) / recs.length;
    const F = P == null || P + R === 0 ? 0 : (2 * P * R) / (P + R);
    return { t, P, R, F, kept, correct, found, nTrue };
  }
  function fmaxOf(judged, proteins) {
    const curve = [];
    for (let t = 0; t <= 100; t++) curve.push(atCutoff(judged, proteins, t));
    const best = Math.max(...curve.map((c) => c.F));
    const at = curve.filter((c) => Math.abs(c.F - best) < 1e-12).map((c) => c.t);
    return { curve, fmax: best, lo: Math.min(...at), hi: Math.max(...at) };
  }
  // AUPRC: every call of every protein in one list, best score first.
  function prCurve(judged, proteins) {
    const all = [];
    judged.forEach((calls, i) => calls.forEach((c) => all.push(Object.assign({ protein: i }, c))));
    all.sort((a, b) => b.score - a.score);
    const nTrue = proteins.reduce((a, p) => a + p.true.length, 0);
    let tp = 0, prevR = 0, area = 0;
    return {
      points: all.map((c, k) => {
        tp += c.correct ? 1 : 0;
        const P = tp / (k + 1), R = tp / nTrue;
        const piece = (R - prevR) * P;
        const p = { call: c, P, R, prevR, piece, areaSoFar: (area += piece) };
        prevR = R;
        return p;
      }),
      auprc: area,
    };
  }

  // ---------------------------------------------------------------- drawing helpers
  function svgEl(name, attrs, parent) {
    const e = document.createElementNS(SVGNS, name);
    for (const [k, v] of Object.entries(attrs || {})) e.setAttribute(k, v);
    if (parent) parent.appendChild(e);
    return e;
  }
  function text(parent, x, y, s, attrs) {
    const t = svgEl("text", Object.assign({ x, y }, attrs || {}), parent);
    t.textContent = s;
    return t;
  }
  function widthOf(box) {
    return Math.max(300, Math.min(760, box.clientWidth || 600));
  }
  // A protein: a thin line with a box for each true domain, then one row per call.
  // state(c) returns "kept" or "drop"; note(c) is text after the glyph.
  function drawTrack(box, protein, calls, state, note, opts) {
    opts = opts || {};
    box.innerHTML = "";
    const W = widthOf(box);
    const LW = Math.min(150, Math.round(W * 0.36)), R = 12, top = 30, rowH = 22;
    const H = top + 32 + calls.length * rowH + 6;
    const svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: "img",
      "aria-label": opts.label || "A protein with its true domains and the calls on it" }, box);
    const sx = (x) => LW + ((x - 1) / (protein.length - 1)) * (W - LW - R);
    text(svg, 0, top + 4, `${protein.name}, ${protein.length} aa`, { class: "ex-lab" });
    svgEl("line", { x1: sx(1), x2: sx(protein.length), y1: top, y2: top, class: "ex-backbone" }, svg);
    for (const t of protein.true) {
      svgEl("rect", { x: sx(t.start), y: top - 9, width: sx(t.end) - sx(t.start), height: 18, rx: 3, class: "ex-true" }, svg);
      text(svg, (sx(t.start) + sx(t.end)) / 2, top + 4, t.family, { class: "ex-domlab", "text-anchor": "middle" });
    }
    const tick = protein.length >= 400 ? 100 : 50;
    for (let x = 0; x <= protein.length; x += tick) {
      const xx = sx(Math.max(1, x));
      svgEl("line", { x1: xx, x2: xx, y1: top + 13, y2: top + 17, class: "ex-tick" }, svg);
      text(svg, xx, top + 28, x === 0 ? "1" : String(x), { class: "ex-tk", "text-anchor": "middle" });
    }
    calls.forEach((c, i) => {
      const y = top + 44 + i * rowH;
      const st = state(c);
      const g = svgEl("g", { "data-call": c.id, "data-state": st, "data-correct": c.correct ? "1" : "0" }, svg);
      text(g, 0, y + 4, `${c.family} ${c.score}`, { class: "ex-lab mono" + (st === "kept" ? "" : " ex-off") });
      text(g, LW - 30, y + 5, c.correct ? "✓" : "✗", { class: "ex-glyph" + (st === "kept" ? "" : " ex-off") });
      const n = note ? note(c) : "";
      // A note goes beside the bar, on whichever side has room, never over it.
      if (n) {
        const right = sx(c.end) + 6 + n.length * 6.2 <= W - R;
        text(g, right ? sx(c.end) + 6 : sx(c.start) - 6, y + 4, n, { class: "ex-note", "text-anchor": right ? "start" : "end" });
      }
      svgEl("rect", { x: sx(c.start), y: y - 5, width: Math.max(2, sx(c.end) - sx(c.start)), height: 10, rx: 2,
        class: st === "kept" ? "ex-kept" : "ex-drop" }, g);
    });
    return svg;
  }
  // A small x-y frame with grid lines at the ticks.
  function frame(box, o) {
    box.innerHTML = "";
    const W = widthOf(box), H = o.height || 230, L = 46, R = 14, T = 12, B = 40;
    const svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": o.label }, box);
    const sx = (x) => L + ((x - o.x0) / (o.x1 - o.x0)) * (W - L - R);
    const sy = (y) => T + (1 - (y - o.y0) / (o.y1 - o.y0)) * (H - T - B);
    for (const x of o.xticks) {
      svgEl("line", { x1: sx(x), x2: sx(x), y1: T, y2: H - B, class: "ex-grid" }, svg);
      text(svg, sx(x), H - B + 15, String(x), { class: "ex-tk", "text-anchor": "middle" });
    }
    for (const y of o.yticks) {
      svgEl("line", { x1: L, x2: W - R, y1: sy(y), y2: sy(y), class: "ex-grid" }, svg);
      text(svg, L - 6, sy(y) + 4, y.toFixed(o.ydp == null ? 1 : o.ydp), { class: "ex-tk", "text-anchor": "end" });
    }
    text(svg, (L + W - R) / 2, H - 6, o.xlab, { class: "ex-axlab", "text-anchor": "middle" });
    text(svg, 12, (T + H - B) / 2, o.ylab, { class: "ex-axlab", "text-anchor": "middle", transform: `rotate(-90 12 ${(T + H - B) / 2})` });
    return { svg, sx, sy, W, H, L, R, T, B };
  }
  function put(id, v) {
    const e = document.getElementById(id);
    if (!e) return;
    e.textContent = typeof v === "number" ? f3(v) : v;
    if (typeof v === "number") e.dataset.value = String(v);
  }
  function onResize(fn) {
    let w = window.innerWidth;
    window.addEventListener("resize", () => { if (window.innerWidth !== w) { w = window.innerWidth; fn(); } });
  }

  const A = TOY.proteins.A, B = TOY.proteins.B;
  const jA = judge(A, A.calls), jB = judge(B, B.calls);

  // ---------------------------------------------------------------- Figure 1: Fmax
  function fig1() {
    const slider = document.getElementById("f1-cut");
    const fm = fmaxOf([jA], [A]);
    function draw() {
      const t = +slider.value;
      const r = atCutoff([jA], [A], t);
      drawTrack(document.getElementById("f1-track"), A, jA, (c) => (c.score >= t ? "kept" : "drop"), null,
        { label: "Protein A with three true domains and seven calls; calls scoring below the cutoff are grey" });
      put("f1-t", String(t));
      put("f1-kept", String(r.kept));
      put("f1-correct", String(r.correct));
      put("f1-found", String(r.found));
      put("f1-ntrue", String(A.true.length));
      put("f1-p", r.P == null ? "not defined (no call kept)" : r.P);
      put("f1-r", r.R);
      put("f1-f", r.F);
      put("f1-fmax", fm.fmax);
      put("f1-best", fm.lo === fm.hi ? String(fm.hi) : `${fm.lo}–${fm.hi}`);
      document.getElementById("f1-f").dataset.value = String(r.F);
      document.getElementById("f1-p").dataset.value = r.P == null ? "" : String(r.P);
      // F against the cutoff.
      const fr = frame(document.getElementById("f1-curve"), {
        x0: 0, x1: 100, y0: 0, y1: 1, xticks: [0, 20, 40, 60, 80, 100], yticks: [0, 0.25, 0.5, 0.75, 1], ydp: 2,
        xlab: "score cutoff (keep calls scoring at least this)", ylab: "F", height: 210,
        label: "F at each score cutoff, with the highest F ringed and your cutoff as a vertical line" });
      // Step line: F holds its value until the cutoff passes the next call's score.
      let d = "";
      fm.curve.forEach((c, i) => {
        const x = fr.sx(c.t), y = fr.sy(c.F);
        d += i === 0 ? `M${x},${y}` : `H${x}V${y}`;
      });
      svgEl("path", { d, class: "ex-line" }, fr.svg);
      svgEl("circle", { cx: fr.sx(fm.hi), cy: fr.sy(fm.fmax), r: 7, class: "ex-ring" }, fr.svg);
      const lx = fr.sx(fm.hi) + 10 > fr.W - 120 ? fr.sx(fm.hi) - 10 : fr.sx(fm.hi) + 10;
      text(fr.svg, lx, fr.sy(fm.fmax) - 10, `Fmax ${f3(fm.fmax)}`, { class: "ex-note", "text-anchor": lx < fr.sx(fm.hi) ? "end" : "start" });
      svgEl("line", { x1: fr.sx(t), x2: fr.sx(t), y1: fr.T, y2: fr.H - fr.B, class: "ex-cut" }, fr.svg);
      svgEl("circle", { cx: fr.sx(t), cy: fr.sy(r.F), r: 3.5, class: "ex-dot-cut" }, fr.svg);
    }
    slider.addEventListener("input", draw);
    document.getElementById("f1-go").addEventListener("click", () => { slider.value = fm.hi; draw(); });
    onResize(draw);
    draw();
  }

  // ---------------------------------------------------------------- Figure 2: AUPRC
  function fig2() {
    const slider = document.getElementById("f2-cut");
    const chips = [...document.querySelectorAll("#f2-which button")];
    let withB = false;
    function draw() {
      const t = +slider.value;
      const judged = withB ? [jA, jB] : [jA];
      const prots = withB ? [A, B] : [A];
      const pr = prCurve(judged, prots);
      const fm = fmaxOf(judged, prots);
      const kept = pr.points.filter((p) => p.call.score >= t);
      const areaKept = kept.reduce((a, p) => a + p.piece, 0);
      put("f2-auprc", pr.auprc);
      put("f2-area", areaKept);
      put("f2-fmax", fm.fmax);
      put("f2-t", String(t));
      chips.forEach((b) => b.setAttribute("aria-pressed", (b.dataset.b === "1") === withB ? "true" : "false"));
      // Both proteins' scores, with and without B, for the side-by-side numbers.
      const a1 = prCurve([jA], [A]).auprc, a2 = prCurve([jA, jB], [A, B]).auprc;
      const f1 = fmaxOf([jA], [A]).fmax, f2 = fmaxOf([jA, jB], [A, B]).fmax;
      put("f2-a-only-auprc", a1); put("f2-ab-auprc", a2);
      put("f2-a-only-fmax", f1); put("f2-ab-fmax", f2);
      document.querySelectorAll("#f2-table [data-col]").forEach((c) =>
        c.classList.toggle("ex-now", (c.dataset.col === "ab") === withB));

      const fr = frame(document.getElementById("f2-curve"), {
        x0: 0, x1: 1, y0: 0, y1: 1.05, xticks: [0, 0.25, 0.5, 0.75, 1], yticks: [0, 0.25, 0.5, 0.75, 1], ydp: 2,
        xlab: "recall: share of true domains found", ylab: "precision", height: 260,
        label: "Precision against recall, one point per call taken in score order; the shaded area is AUPRC" });
      // Area: one rectangle per call that finds a new domain (width = recall gained).
      for (const p of pr.points) {
        if (p.piece <= 0) continue;
        const on = p.call.score >= t;
        svgEl("rect", { x: fr.sx(p.prevR), y: fr.sy(p.P), width: fr.sx(p.R) - fr.sx(p.prevR), height: fr.sy(0) - fr.sy(p.P),
          class: on ? "ex-area" : "ex-area-off" }, fr.svg);
      }
      // Outline the strips: a call that finds a domain adds a step to the right at its
      // precision; a wrong call drops straight down (no recall gained).
      let d = `M${fr.sx(0)},${fr.sy(pr.points[0].P)}`;
      pr.points.forEach((p) => { d += `V${fr.sy(p.P)}H${fr.sx(p.R)}`; });
      svgEl("path", { d, class: "ex-line-thin" }, fr.svg);
      pr.points.forEach((p) => {
        const on = p.call.score >= t;
        svgEl("circle", { cx: fr.sx(p.R), cy: fr.sy(p.P), r: 4, class: on ? "ex-dot" : "ex-dot-off" }, fr.svg);
      });
      // The ranked list beside it: every call of every protein in one list.
      const tb = document.querySelector("#f2-list tbody");
      tb.innerHTML = "";
      pr.points.forEach((p, k) => {
        const tr = document.createElement("tr");
        const on = p.call.score >= t;
        tr.className = on ? "" : "ex-off";
        tr.dataset.state = on ? "kept" : "drop";
        tr.innerHTML = `<td>${k + 1}</td><td>${prots[p.call.protein].short}</td><td class="mono">${p.call.family} ${p.call.score}</td>` +
          `<td>${p.call.correct ? "✓" : "✗"}</td><td class="mono">${f3(p.P)}</td><td class="mono">${f3(p.R)}</td>` +
          `<td class="mono">${p.piece > 0 ? f3(p.piece) : "0"}</td>`;
        tb.appendChild(tr);
      });
    }
    slider.addEventListener("input", draw);
    chips.forEach((b) => b.addEventListener("click", () => { withB = b.dataset.b === "1"; draw(); }));
    onResize(draw);
    draw();
  }

  // ---------------------------------------------------------------- Figure 3: IoU
  function fig3() {
    const P = TOY.iou;
    const s = document.getElementById("f3-start"), e = document.getElementById("f3-end");
    function draw(ev) {
      let a = +s.value, b = +e.value;
      // Keep start <= end: whichever slider moved pushes the other.
      if (a > b) { if (ev && ev.target === s) { b = a; e.value = b; } else { a = b; s.value = a; } }
      const call = { start: a, end: b };
      const r = iou(call, P.true);
      put("f3-start-v", String(a)); put("f3-end-v", String(b));
      put("f3-inter", String(r.inter)); put("f3-union", String(r.union)); put("f3-iou", r.iou);
      document.getElementById("f3-verdict").textContent = r.iou >= 0.5 ? "correct (IoU 0.5 or more)" : "not correct (IoU below 0.5)";
      const box = document.getElementById("f3-track");
      box.innerHTML = "";
      const W = widthOf(box), L = 14, R = 14, H = 150;
      const svg = svgEl("svg", { width: W, height: H, viewBox: `0 0 ${W} ${H}`, role: "img",
        "aria-label": `True domain ${P.true.start}–${P.true.end}, call ${a}–${b}, IoU ${f3(r.iou)}` }, box);
      const sx = (x) => L + ((x - 1) / (P.length - 1)) * (W - L - R);
      const y0 = 26;
      svgEl("line", { x1: sx(1), x2: sx(P.length), y1: y0, y2: y0, class: "ex-backbone" }, svg);
      svgEl("rect", { x: sx(P.true.start), y: y0 - 9, width: sx(P.true.end) - sx(P.true.start), height: 18, rx: 3, class: "ex-true" }, svg);
      text(svg, (sx(P.true.start) + sx(P.true.end)) / 2, y0 + 4, `true domain ${P.true.start}–${P.true.end}`, { class: "ex-domlab", "text-anchor": "middle" });
      svgEl("rect", { x: sx(a), y: y0 + 18, width: Math.max(2, sx(b) - sx(a)), height: 10, rx: 2, class: "ex-kept" }, svg);
      const cl = `call ${a}–${b}`;
      text(svg, Math.min(sx(a), W - R - cl.length * 6.6), y0 + 44, cl, { class: "ex-note" });
      // Brackets: shared residues, then residues covered by either.
      const lo = Math.min(a, P.true.start), hi = Math.max(b, P.true.end);
      const ia = Math.max(a, P.true.start), ib = Math.min(b, P.true.end);
      function bracket(x1, x2, y, label) {
        svgEl("path", { d: `M${sx(x1)},${y - 5}V${y}H${sx(x2)}V${y - 5}`, class: "ex-bracket" }, svg);
        const lab = label, w = lab.length * 6.6;
        const mid = (sx(x1) + sx(x2)) / 2;
        text(svg, Math.max(L, Math.min(W - R - w, mid - w / 2)), y + 15, lab, { class: "ex-note" });
      }
      if (r.inter > 0) bracket(ia, ib, y0 + 70, `shared: ${r.inter} residues`);
      else text(svg, L, y0 + 85, "shared: 0 residues", { class: "ex-note" });
      bracket(lo, hi, y0 + 102, `covered by either: ${r.union} residues`);

      // The median over correct calls: Figure 1's correct calls, plus yours if correct.
      const rows = jA.filter((c) => c.correct).map((c) => ({ name: `${c.family} ${c.score}`, v: iou(c, A.true[c.domain]).iou, mine: false }));
      if (r.iou >= 0.5) rows.push({ name: "your call", v: r.iou, mine: true });
      const vs = rows.map((x) => x.v).sort((x, y) => x - y);
      const n = vs.length, med = n % 2 ? vs[(n - 1) / 2] : (vs[n / 2 - 1] + vs[n / 2]) / 2;
      put("f3-median", med);
      put("f3-n", String(n));
      document.getElementById("f3-mine").textContent = r.iou >= 0.5
        ? "Your call is correct, so it is one of the values."
        : "Your call is not correct, so it is left out: the median is over correct calls only.";
      const fr = frame(document.getElementById("f3-strip"), {
        x0: 0.5, x1: 1, y0: 0, y1: 1, xticks: [0.5, 0.6, 0.7, 0.8, 0.9, 1], yticks: [], height: 150,
        xlab: "IoU of each correct call with its true domain", ylab: "",
        label: "IoU of each correct call, with the median as a vertical line" });
      svgEl("line", { x1: fr.sx(med), x2: fr.sx(med), y1: fr.T, y2: fr.H - fr.B, class: "ex-median" }, fr.svg);
      text(fr.svg, fr.sx(med) + (fr.sx(med) > fr.W - 110 ? -4 : 4), fr.T + 12, `median ${f3(med)}`,
        { class: "ex-note", "text-anchor": fr.sx(med) > fr.W - 110 ? "end" : "start" });
      rows.slice().sort((x, y) => x.v - y.v).forEach((x, i) => {
        // One row per call, lowest IoU at the top, each labelled with its call.
        const cx = fr.sx(x.v), cy = fr.T + 28 + i * 16;
        svgEl("circle", { cx, cy, r: 4.5, class: x.mine ? "ex-ring-small" : "ex-dot", "data-iou": x.v }, fr.svg);
        // Labels sit on the side away from the median line, so the line never crosses them.
        const lab = `${x.name}: ${f3(x.v)}`, w = lab.length * 6.2;
        const left = (x.v < med && cx - 8 - w > fr.L) || cx + 8 + w > fr.W - fr.R;
        text(fr.svg, left ? cx - 8 : cx + 8, cy + 4, lab, { class: "ex-note", "text-anchor": left ? "end" : "start" });
      });
    }
    s.addEventListener("input", draw);
    e.addEventListener("input", draw);
    onResize(draw);
    draw();
  }

  // ---------------------------------------------------------------- Figure 4: raw vs dedup
  function fig4() {
    const chips = [...document.querySelectorAll("#f4-which button")];
    let mode = "raw";
    const d = dedup(A.calls);
    const jD = judge(A, d.keep);
    const byId = Object.fromEntries(A.calls.map((c) => [c.id, c]));
    function scores(j) {
      const all = atCutoff([j], [A], 0);
      return { p0: all.P, fmax: fmaxOf([j], [A]).fmax, auprc: prCurve([j], [A]).auprc, n: j.length };
    }
    const raw = scores(jA), ded = scores(jD);
    function draw() {
      chips.forEach((b) => b.setAttribute("aria-pressed", b.dataset.mode === mode ? "true" : "false"));
      // In dedup, a merged call keeps its raw verdict on the drawing but no longer counts.
      drawTrack(document.getElementById("f4-track"), A, jA,
        (c) => (mode === "dedup" && d.merged[c.id] ? "drop" : "kept"),
        (c) => (mode === "dedup" && d.merged[c.id] ? `merged into ${byId[d.merged[c.id]].family} ${byId[d.merged[c.id]].score}` : ""),
        { label: "Protein A's seven calls; in dedup, calls merged into a higher-scoring one are grey" });
      for (const [k, v] of [["raw", raw], ["dedup", ded]]) {
        put(`f4-${k}-n`, String(v.n)); put(`f4-${k}-p0`, v.p0); put(`f4-${k}-fmax`, v.fmax); put(`f4-${k}-auprc`, v.auprc);
      }
      document.querySelectorAll("#f4-table [data-col]").forEach((c) => c.classList.toggle("ex-now", c.dataset.col === mode));
    }
    chips.forEach((b) => b.addEventListener("click", () => { mode = b.dataset.mode; draw(); }));
    onResize(draw);
    draw();
  }

  // ---------------------------------------------------------------- Figure 5: cost
  function fig5() {
    const ids = ["f5-cpu", "f5-gpu", "f5-rcpu", "f5-rgpu"];
    const sig = (x) => (x === 0 ? "0" : x >= 1 ? x.toFixed(3) : x.toPrecision(3));
    function draw() {
      const [cpu, gpu, rc, rg] = ids.map((id) => parseFloat(document.getElementById(id).value));
      if (![cpu, gpu, rc, rg].every((x) => isFinite(x) && x >= 0)) {
        put("f5-out", "enter four numbers of 0 or more");
        return;
      }
      const a = cpu * rc, b = gpu * rg, tot = a + b;
      put("f5-line", `${cpu} × ${rc} + ${gpu} × ${rg} = ${sig(a)} + ${sig(b)}`);
      const o = document.getElementById("f5-out");
      o.textContent = "$" + sig(tot);
      o.dataset.value = String(tot);
    }
    ids.forEach((id) => document.getElementById(id).addEventListener("input", draw));
    draw();
  }

  fig1(); fig2(); fig3(); fig4(); fig5();
  window.__explainers = { ready: true };
})();
