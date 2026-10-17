/* SCOPe40 page: Figure 19 (sensitivity curves), Figure 20 (by SCOPe class) and Table 3.
   Data: data/scope40.json, written by scripts/scope40_page.py from data/scope40_long.parquet.
   Hover text uses figures.js (any [data-tip] inside figure.fgfig). */
(function () {
  "use strict";
  const NS = "http://www.w3.org/2000/svg";
  const LEVELS = ["family", "superfamily", "fold"];
  const LETTER = { family: "a", superfamily: "b", fold: "c" };
  const f3 = (x) => x.toFixed(3);
  const pct = (x) => (100 * x).toFixed(1) + "%";
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const el = (tag, attrs, parent) => {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  };
  const text = (parent, x, y, s, attrs) => { const t = el("text", Object.assign({ x, y }, attrs || {}), parent); t.textContent = s; return t; };
  const mean = (a) => a.reduce((s, v) => s + v, 0) / a.length;
  const marker = (g, shape, x, y, color, tip, open) => {
    let m;
    if (shape === "square") m = el("rect", { x: x - 5, y: y - 5, width: 10, height: 10 }, g);
    else if (shape === "diamond") m = el("path", { d: `M${x},${y - 6.5} L${x + 6.5},${y} L${x},${y + 6.5} L${x - 6.5},${y} Z` }, g);
    else if (shape === "triangle") m = el("path", { d: `M${x},${y - 6} L${x + 6},${y + 5} L${x - 6},${y + 5} Z` }, g);
    else m = el("circle", { cx: x, cy: y, r: 5.5 }, g);
    m.setAttribute("style", open ? `fill:var(--surface);stroke:${color};stroke-width:2` : `fill:${color};stroke:var(--surface);stroke-width:1`);
    if (tip) { m.setAttribute("data-tip", tip); m.setAttribute("tabindex", "0"); }
    return m;
  };

  let D = null;
  const state = { on: {}, mode: "zero" };

  // Values of one tool at one level under the current mode (no-hit queries kept as 0, or dropped).
  const values = (key, lvl, mode) => {
    const v = D.sens[key][lvl], h = D.hits[key];
    return mode === "drop" ? v.filter((_, i) => h[i]) : v.slice();
  };
  const summary = (key, lvl) => {
    const t = D.tools.find((x) => x.key === key);
    return D.summary.overall.find((r) => r.tool === t.tool && r.setting === t.setting && r.level === lvl);
  };

  function drawLegend() {
    const box = document.getElementById("s40-curve-legend");
    const drop = state.mode === "drop";
    let h = '<div class="tablebox"><table class="lb s40-leg"><thead><tr><th class="txt">Tool</th>' +
      LEVELS.map((l) => `<th>${LETTER[l]}, ${l} AUC${drop ? ", no-hit queries dropped" : ""}</th>`).join("") +
      "<th>Coverage</th><th>n</th></tr></thead><tbody>";
    for (const t of D.tools) {
      if (!state.on[t.key]) continue;
      const s0 = summary(t.key, "family");
      h += `<tr data-key="${esc(t.key)}"><td class="txt"><svg width="26" height="12" aria-hidden="true"><line x1="1" x2="25" y1="6" y2="6" style="stroke:${t.color};stroke-width:3${t.dash ? `;stroke-dasharray:${t.dash}` : ""}"/></svg> ${esc(t.name)}${t.tool === "kmerseek" || t.chance ? "" : ", " + esc(t.setting)}</td>`;
      for (const l of LEVELS) {
        const a0 = mean(values(t.key, l, "zero"));
        if (!drop) h += `<td class="num" data-auc="${a0}">${f3(a0)}</td>`;
        else {
          const a1 = mean(values(t.key, l, "drop"));
          const d = a1 - a0;
          h += `<td class="num" data-auc="${a1}" data-delta="${d}">${f3(a1)} <span class="key">${d >= 0 ? "+" : "−"}${f3(Math.abs(d))} from ${f3(a0)}</span></td>`;
        }
      }
      const n = drop ? s0.n_with_hits : s0.n_queries;
      h += `<td class="num">${pct(s0.coverage)}</td><td class="num">${n.toLocaleString("en-US")}</td></tr>`;
    }
    h += "</tbody></table></div>";
    box.innerHTML = h;
  }

  function drawCurves() {
    const box = document.getElementById("s40-curves");
    box.innerHTML = "";
    const W = 360, H = 250, L = 54, R = 18, T = 12, B = 40;
    const X = (x) => L + x * (W - L - R), Y = (y) => T + (1 - y) * (H - T - B);
    for (const lvl of LEVELS) {
      const panel = document.createElement("div");
      panel.className = "fg-panel";
      panel.innerHTML = `<h4>${LETTER[lvl]}. ${lvl[0].toUpperCase() + lvl.slice(1)}</h4>`;
      box.appendChild(panel);
      const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, class: "fg", role: "img", "aria-label": `Sensitivity curves, ${lvl} level` }, panel);
      for (const v of [0, 0.25, 0.5, 0.75, 1]) {
        el("line", { x1: L, x2: W - R, y1: Y(v), y2: Y(v), class: "fg-grid" }, svg);
        text(svg, L - 6, Y(v) + 4, v.toFixed(2), { class: "fg-tick", "text-anchor": "end" });
        el("line", { x1: X(v), x2: X(v), y1: T, y2: H - B, class: "fg-grid" }, svg);
        text(svg, X(v), H - B + 15, v.toFixed(2), { class: "fg-tick", "text-anchor": "middle" });
      }
      el("line", { x1: L, x2: W - R, y1: H - B, y2: H - B, class: "fg-axis" }, svg);
      el("line", { x1: L, x2: L, y1: T, y2: H - B, class: "fg-axis" }, svg);
      text(svg, (L + W - R) / 2, H - 6, "fraction of queries, highest sensitivity first", { class: "fg-lab", "text-anchor": "middle" });
      text(svg, 15, (T + H - B) / 2, "sensitivity", { class: "fg-lab", "text-anchor": "middle", transform: `rotate(-90 15 ${(T + H - B) / 2})` });
      const order = D.tools.filter((t) => state.on[t.key]).sort((a, b) => (a.primary === b.primary ? 0 : a.primary ? 1 : -1));
      for (const t of order) {
        const v = values(t.key, lvl, state.mode).sort((a, b) => b - a);
        const n = v.length;
        let d = `M${X(0).toFixed(2)},${Y(v[0]).toFixed(2)}`;
        for (let i = 1; i < n; i++) {
          if (v[i] !== v[i - 1]) d += `H${X(i / n).toFixed(2)}V${Y(v[i]).toFixed(2)}`;
        }
        d += `H${X(1).toFixed(2)}`;
        const a = mean(v), s = summary(t.key, lvl);
        const tip = `${t.name}, ${t.setting}, ${lvl}: AUC ${f3(a)}` + (state.mode === "zero" ? ` (95% interval ${f3(s.lo)}–${f3(s.hi)})` : ", no-hit queries dropped") +
          `, coverage ${pct(s.coverage)}, n = ${n.toLocaleString("en-US")}`;
        const style = t.primary ? `fill:none;stroke:${t.color};stroke-width:2.2` : `fill:none;stroke:${t.color};stroke-width:1.2;opacity:0.35`;
        el("path", { d, style: style + (t.dash ? `;stroke-dasharray:${t.dash}` : ""), "data-tool": t.key }, svg);
        el("path", { d, style: "fill:none;stroke:transparent;stroke-width:10", "data-tip": tip, tabindex: "0" }, svg);
      }
    }
  }

  function drawClasses() {
    const leg = document.getElementById("s40-class-legend");
    const P = D.tools.filter((t) => t.primary);
    leg.innerHTML = '<span class="item">dot = mean, bar = 95% bootstrap interval:</span>' + P.map((t) => `<span class="item"><svg width="16" height="16" aria-hidden="true">` +
      (t.shape === "square" ? `<rect x="3" y="3" width="10" height="10" style="fill:${t.color}"/>`
        : t.chance ? `<circle cx="8" cy="8" r="5" style="fill:var(--surface);stroke:${t.color};stroke-width:2"/>` : `<circle cx="8" cy="8" r="5.5" style="fill:${t.color}"/>`) +
      `</svg>${esc(t.name)}</span>`).join("") +
      '<span class="item muted">Composition scan: not run</span>';
    const box = document.getElementById("s40-classes");
    box.innerHTML = "";
    const cls = Object.keys(D.classNames);
    const W = 820, H = 230, L = 54, R = 10, T = 12, B = 46;
    const slot = (W - L - R) / cls.length;
    const Y = (y) => T + (1 - y) * (H - T - B);
    for (const lvl of LEVELS) {
      const panel = document.createElement("div");
      panel.className = "fg-panel";
      panel.innerHTML = `<h4>${LETTER[lvl]}. ${lvl[0].toUpperCase() + lvl.slice(1)}</h4>`;
      box.appendChild(panel);
      const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, class: "fg", role: "img", "aria-label": `Sensitivity by SCOPe class, ${lvl} level` }, panel);
      for (const v of [0, 0.25, 0.5, 0.75, 1]) {
        el("line", { x1: L, x2: W - R, y1: Y(v), y2: Y(v), class: "fg-grid" }, svg);
        text(svg, L - 6, Y(v) + 4, v.toFixed(2), { class: "fg-tick", "text-anchor": "end" });
      }
      text(svg, 15, (T + H - B) / 2, "mean sensitivity", { class: "fg-lab", "text-anchor": "middle", transform: `rotate(-90 15 ${(T + H - B) / 2})` });
      cls.forEach((c, j) => {
        const x0 = L + j * slot, cx = x0 + slot / 2;
        const g = el("g", { class: "fg-col", "data-cls": c, "data-level": lvl, tabindex: "0", role: "button",
          "aria-label": `List the class ${c} queries` }, svg);
        el("rect", { x: x0 + 2, y: T, width: slot - 4, height: H - T - 4, style: "fill:transparent" }, g);
        const rows = D.summary.by_class.filter((r) => r.scope_class === c && r.level === lvl);
        P.forEach((t, k) => {
          const r = rows.find((x) => x.tool === t.tool && x.setting === t.setting);
          const x = cx + (k - (P.length - 1) / 2) * 16;
          el("line", { x1: x, x2: x, y1: Y(r.lo), y2: Y(r.hi), style: `stroke:${t.color};stroke-width:2` }, g);
          marker(g, t.shape, x, Y(r.mean), t.color,
            `${t.name}, class ${c} (${D.classNames[c]}), ${lvl}: mean ${f3(r.mean)}, 95% interval ${f3(r.lo)}–${f3(r.hi)}, n = ${r.n_queries.toLocaleString("en-US")}. Click to list the queries.`, t.chance);
        });
        text(g, cx, H - B + 16, `class ${c}`, { class: "fg-lab", "text-anchor": "middle" });
        text(g, cx, H - B + 31, `n = ${rows[0].n_queries.toLocaleString("en-US")}`, { class: "fg-tick", "text-anchor": "middle" });
      });
    }
  }

  function listClass(c, lvl) {
    document.querySelectorAll("#s40-classes g.fg-col.on").forEach((x) => x.classList.remove("on"));
    document.querySelectorAll(`#s40-classes g.fg-col[data-cls="${c}"][data-level="${lvl}"]`).forEach((x) => x.classList.add("on"));
    const idx = D.cls.map((x, i) => (x === c ? i : -1)).filter((i) => i >= 0);
    const P = D.tools.filter((t) => t.primary);
    const head = P.map((t) => `<th>${esc(t.name)}</th>`).join("");
    const rows = idx.map((i) => `<tr><td class="mono">${esc(D.queries[i])}</td><td class="mono">${esc(D.sccs[i])}</td>` +
      P.map((t) => `<td class="num">${f3(D.sens[t.key][lvl][i])}</td>`).join("") + "</tr>").join("");
    document.getElementById("s40-class-list").innerHTML =
      `<p><b>Class ${esc(c)} (${esc(D.classNames[c])}), ${esc(lvl)} level: ${idx.length.toLocaleString("en-US")} queries.</b> ` +
      `Sensitivity up to the first false positive per tool.</p><div class="tablebox" style="max-height:420px;overflow:auto">` +
      `<table class="fg-regions" id="s40-class-table"><thead><tr><th>query</th><th>SCOPe id</th>${head}</tr></thead><tbody>${rows}</tbody></table></div>`;
  }

  function sortTable(col) {
    const tab = document.getElementById("s40-tab");
    const th = tab.querySelectorAll("thead th")[col];
    const dir = th.getAttribute("aria-sort") === "descending" ? "ascending" : "descending";
    tab.querySelectorAll("thead th").forEach((x) => x.removeAttribute("aria-sort"));
    th.setAttribute("aria-sort", dir);
    const body = tab.tBodies[0];
    const rows = [...body.rows];
    const pinned = rows.filter((r) => r.classList.contains("pinned"));
    const toggles = rows.filter((r) => r.classList.contains("more-toggle") && !r.classList.contains("pinned"));
    const rest = rows.filter((r) => !r.classList.contains("pinned") && !r.classList.contains("more-toggle"));
    const key = (r) => {
      const td = r.cells[col];
      if (col === 0) return td.textContent.trim().toLowerCase();
      const v = td.getAttribute("data-v");
      return v === "" || v === null ? null : parseFloat(v);
    };
    rest.sort((a, b) => {
      const x = key(a), y = key(b);
      if (x === null && y === null) return 0;
      if (x === null) return 1;
      if (y === null) return -1;
      const c = typeof x === "string" ? x.localeCompare(y) : x - y;
      return dir === "ascending" ? c : -c;
    });
    for (const r of rest.concat(toggles, pinned)) body.appendChild(r);
  }

  function costs() {
    // Cost per 1,000 queries = core-hours x r_cpu + GPU-hours x r_gpu, where both are measured.
    const rc = parseFloat(document.getElementById("s40-rcpu").value), rg = parseFloat(document.getElementById("s40-rgpu").value);
    document.querySelectorAll("#s40-tab tbody tr").forEach((tr) => {
      const ch = tr.cells[6].getAttribute("data-v"), gh = tr.cells[7].getAttribute("data-v");
      const td = tr.cells[8];
      if (ch === "" || ch === null) { td.textContent = "—"; td.setAttribute("data-v", ""); return; }
      const v = parseFloat(ch) * rc + (gh ? parseFloat(gh) : 0) * rg;
      td.textContent = "$" + v.toFixed(4);
      td.setAttribute("data-v", String(v));
    });
  }

  function redraw() { drawLegend(); drawCurves(); }

  async function init() {
    try { D = await (await fetch("data/scope40.json")).json(); }
    catch (err) { document.getElementById("s40-curves").textContent = "Could not load data/scope40.json."; return; }
    const boxes = document.getElementById("s40-tool-boxes");
    const box = (t) => `<label><input type="checkbox"${t.primary ? " checked" : ""} data-key="${esc(t.key)}"> <svg width="20" height="10" aria-hidden="true"><line x1="1" x2="19" y1="5" y2="5" style="stroke:${t.color};stroke-width:3${t.dash ? `;stroke-dasharray:${t.dash}` : ""}"/></svg> ${esc(t.name)}</label>`;
    for (const t of D.tools) state.on[t.key] = t.primary;
    const extra = D.tools.filter((t) => !t.primary);
    boxes.innerHTML = D.tools.filter((t) => t.primary).map(box).join(" ") +
      (extra.length ? `<details class="s40-more-settings"><summary>+${extra.length} more settings (kmerseek alphabets and k-sizes, real and shuffled queries)</summary>${extra.map(box).join("<br>")}</details>` : "");
    boxes.addEventListener("change", (e) => { state.on[e.target.getAttribute("data-key")] = e.target.checked; redraw(); });
    document.querySelectorAll('input[name="s40-nohit"]').forEach((r) => r.addEventListener("change", (e) => { state.mode = e.target.value; redraw(); }));
    redraw();
    drawClasses();
    const cb = document.getElementById("s40-classes");
    const pick = (e) => { const g = e.target.closest && e.target.closest("g.fg-col"); if (g) listClass(g.getAttribute("data-cls"), g.getAttribute("data-level")); };
    cb.addEventListener("click", pick);
    cb.addEventListener("keydown", (e) => { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); pick(e); } });
    document.querySelectorAll("#s40-tab thead button").forEach((b) => b.addEventListener("click", () => sortTable(parseInt(b.getAttribute("data-col"), 10))));
    document.querySelectorAll("#s40-tab .s40-more").forEach((b) => b.addEventListener("click", () => {
      const open = b.getAttribute("aria-expanded") !== "true";
      b.setAttribute("aria-expanded", String(open));
      document.querySelectorAll(`#s40-tab tr.more-${b.getAttribute("data-tool")}`).forEach((r) => { r.hidden = !open; });
    }));
    document.querySelectorAll("#s40-tab tr.more").forEach((r) => { r.hidden = true; });
    ["s40-rcpu", "s40-rgpu"].forEach((id) => document.getElementById(id).addEventListener("input", costs));
    costs();
    document.body.setAttribute("data-s40-ready", "1");
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
