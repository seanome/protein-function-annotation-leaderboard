/* Leaderboard and score-against-cost pages. Reads the JSON that build_site.py inlines
   into the page (#lb-data); never computes a metric, only looks rows up. A (tool,
   setting) with no row for the chosen filters is drawn as "—" and the words "not run";
   a missing compute row is "—" and "not measured". */
(function () {
  "use strict";
  const DATA = JSON.parse(document.getElementById("lb-data").textContent);
  const M = DATA.meta;
  const RESULTS = DATA.results; // key -> {tool|setting: [fmax, auprc, iou, n]}
  const COMPUTE = DATA.compute; // tool|setting|species -> {...}
  const page = document.body.dataset.page;

  // ---------- state, kept in the URL hash so both pages share the filters ----------
  const state = Object.assign({}, M.defaults);
  function readHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    for (const k of ["truth", "calls", "species", "axis", "stratum", "metric", "sort", "dir"]) {
      if (h.has(k)) state[k] = h.get(k);
    }
    for (const k of ["r_cpu", "r_gpu"]) {
      if (h.has(k) && isFinite(parseFloat(h.get(k)))) state[k] = parseFloat(h.get(k));
    }
  }
  function writeHash() {
    const h = new URLSearchParams();
    for (const k of ["truth", "calls", "species", "axis", "stratum", "metric", "sort", "dir", "r_cpu", "r_gpu"]) {
      if (state[k] !== undefined && state[k] !== null) h.set(k, state[k]);
    }
    history.replaceState(null, "", "#" + h.toString());
    // Keep the filters when moving between the leaderboard and the plot.
    document.querySelectorAll("a[data-keep-filters]").forEach((a) => {
      a.href = a.getAttribute("href").split("#")[0] + "#" + h.toString();
    });
  }

  function resultKey(s) {
    return [s.truth, s.calls, s.species, s.axis, s.stratum].join("|");
  }
  function lookup(s, row) {
    const block = RESULTS[resultKey(s)];
    return block ? block[row.tool + "|" + row.setting] || null : null;
  }
  function computeRow(s, row) {
    return COMPUTE[row.tool + "|" + row.setting + "|" + s.species] || null;
  }
  function cost(c) {
    // cost_usd_est = cpu_h * r_cpu + gpu_h * r_gpu, per 1,000 queries.
    if (!c) return null;
    const cpu = c.cpu_h_per_1000 != null ? c.cpu_h_per_1000 : c.cpu_h_reported_per_1000;
    if (cpu == null) return null;
    return cpu * state.r_cpu + (c.gpu_h_per_1000 || 0) * state.r_gpu;
  }

  // ---------- number formats ----------
  const fmt3 = (x) => x.toFixed(3);
  const fmtInt = (x) => x.toLocaleString("en-US");
  function sig(x, n) {
    if (x === 0) return "0";
    const d = Math.max(0, n - 1 - Math.floor(Math.log10(Math.abs(x))));
    return x.toFixed(Math.min(d, 6));
  }

  // ---------- chips ----------
  function chip(label, pressed, onClick, title) {
    const b = document.createElement("button");
    b.type = "button";
    b.className = "chip";
    b.textContent = label;
    b.setAttribute("aria-pressed", pressed ? "true" : "false");
    if (title) b.title = title;
    b.addEventListener("click", onClick);
    return b;
  }
  function renderFilters() {
    const box = document.getElementById("filters");
    box.innerHTML = "";
    for (const row of M.chiprows) {
      const r = document.createElement("div");
      r.className = "chiprow";
      r.dataset.row = row.id;
      const a = document.createElement("div");
      a.className = "axis";
      a.innerHTML = row.label + (row.help ? "<small>" + row.help + "</small>" : "");
      r.appendChild(a);
      const chips = document.createElement("div");
      chips.className = "chips";
      for (const g of row.groups) {
        const grp = document.createElement("span");
        grp.className = "chipgroup";
        if (g.label) {
          const gl = document.createElement("span");
          gl.className = "glabel";
          gl.textContent = g.label;
          grp.appendChild(gl);
        }
        for (const c of g.chips) {
          let pressed;
          if (row.kind === "stratum") {
            pressed = c.value === "all" ? state.axis !== row.axis : state.axis === row.axis && state.stratum === c.value;
          } else {
            pressed = state[row.id] === c.value;
          }
          grp.appendChild(
            chip(c.label, pressed, () => {
              if (row.kind === "stratum") {
                // The run cut the queries along one stratum axis at a time, so picking a
                // value here clears the other three stratum rows.
                if (c.value === "all") {
                  if (state.axis === row.axis) { state.axis = "all"; state.stratum = "all"; }
                } else { state.axis = row.axis; state.stratum = c.value; }
              } else {
                state[row.id] = c.value;
              }
              update();
            }, c.title)
          );
        }
        chips.appendChild(grp);
      }
      r.appendChild(chips);
      box.appendChild(r);
    }
  }

  // ---------- leaderboard table ----------
  const COLS = [
    { id: "tool", label: "Tool", txt: true },
    { id: "setting", label: "Setting", txt: true },
    { id: "fmax", label: "Fmax", key: "fmax" },
    { id: "auprc", label: "AUPRC", key: "auprc" },
    { id: "iou", label: "Median IoU", key: "median_iou_tp" },
    { id: "n", label: "True instances", key: "n_truth_instances" },
    { id: "cpu", label: 'CPU-h per 1,000 queries <span class="prov measured">measured</span>', key: "cpu_h" },
    { id: "cpu_rep", label: 'CPU-h per 1,000 queries <span class="prov reported">reported</span>', key: "cpu_h_reported" },
    { id: "ram", label: 'Peak RAM, GB <span class="prov measured">measured</span>', key: "peak_rss_gb" },
    { id: "gpu", label: "GPU in search", key: "gpu" },
    { id: "cost", label: "Estimated cost, USD per 1,000 queries", key: "cost_usd_est" },
  ];
  function rowValues(s, row) {
    const r = lookup(s, row);
    const c = computeRow(s, row);
    return {
      row, r, c,
      fmax: r ? r[0] : null, auprc: r ? r[1] : null, iou: r ? r[2] : null, n: r ? r[3] : null,
      cpu: c ? c.cpu_h_per_1000 : null,
      cpu_rep: c ? c.cpu_h_reported_per_1000 : null,
      ram: c ? c.peak_rss_gb : null,
      gpu: c ? c.gpu : null,
      cost: cost(c),
    };
  }
  function sortRows(vals) {
    const k = state.sort, dir = state.dir === "asc" ? 1 : -1;
    return vals.sort((a, b) => {
      let x = k === "tool" ? M.tools[a.row.tool].name : k === "setting" ? a.row.setting : a[k];
      let y = k === "tool" ? M.tools[b.row.tool].name : k === "setting" ? b.row.setting : b[k];
      // Missing values always sort last, whichever way the column is sorted.
      if (x == null && y == null) return a.row.order - b.row.order;
      if (x == null) return 1;
      if (y == null) return -1;
      if (x < y) return -dir;
      if (x > y) return dir;
      return a.row.order - b.row.order;
    });
  }
  function missingCell(kind, words, span) {
    const td = document.createElement("td");
    td.className = "missing";
    td.dataset.missing = kind;
    if (span) td.colSpan = span;
    td.innerHTML = "— " + words;
    return td;
  }
  function numCell(text) {
    const td = document.createElement("td");
    td.className = "num";
    td.textContent = text;
    return td;
  }
  function renderTable() {
    const thead = document.querySelector("#lb thead tr");
    thead.innerHTML = "";
    for (const col of COLS) {
      const th = document.createElement("th");
      th.scope = "col";
      if (col.txt) th.className = "txt";
      if (state.sort === col.id) th.setAttribute("aria-sort", state.dir === "asc" ? "ascending" : "descending");
      const b = document.createElement("button");
      b.innerHTML = col.label + '<span class="arrow"></span>' + (col.key ? '<span class="key">' + col.key + "</span>" : "");
      b.addEventListener("click", () => {
        if (state.sort === col.id) state.dir = state.dir === "asc" ? "desc" : "asc";
        else { state.sort = col.id; state.dir = col.txt || col.id === "cpu" || col.id === "cpu_rep" || col.id === "ram" || col.id === "cost" ? "asc" : "desc"; }
        update();
      });
      th.appendChild(b);
      thead.appendChild(th);
    }
    const tbody = document.querySelector("#lb tbody");
    tbody.innerHTML = "";
    const vals = sortRows(M.rows.map((row) => rowValues(state, row)));
    for (const v of vals) {
      const tr = document.createElement("tr");
      tr.dataset.tool = v.row.tool;
      tr.dataset.setting = v.row.setting;
      const t = M.tools[v.row.tool];
      const td1 = document.createElement("td");
      td1.className = "txt";
      td1.innerHTML = '<span class="kind k-' + t.kind + '" aria-hidden="true"></span><a href="tools/' + v.row.tool + '.html">' + t.name + "</a>";
      tr.appendChild(td1);
      const td2 = document.createElement("td");
      td2.className = "txt mono";
      td2.textContent = v.row.setting;
      tr.appendChild(td2);
      if (v.r) {
        tr.appendChild(numCell(fmt3(v.fmax)));
        tr.appendChild(numCell(fmt3(v.auprc)));
        tr.appendChild(numCell(v.iou == null ? "" : fmt3(v.iou)));
        tr.appendChild(numCell(fmtInt(v.n)));
      } else {
        tr.appendChild(missingCell("results", "not run", 4));
      }
      if (v.c) {
        tr.appendChild(numCell(sig(v.cpu, 3)));
        tr.appendChild(v.cpu_rep == null ? missingCell("reported", "none reported") : numCell(sig(v.cpu_rep, 3)));
        tr.appendChild(numCell(v.ram.toFixed(1)));
        const g = numCell(v.gpu);
        if (v.row.tool === "prostt5") g.innerHTML = v.gpu + ' <a href="tools/prostt5.html#failure-modes" title="Predicting the 3Di letters used one GPU; not counted">*</a>';
        tr.appendChild(g);
        tr.appendChild(numCell("$" + sig(v.cost, 3)));
      } else {
        // "not measured" only when this tool searched this relative (it has result rows
        // under some filter) but no trace row holds the search task.
        const why = M.searched.includes(v.row.tool + "|" + v.row.setting + "|" + state.species) ? "not measured" : "not run";
        tr.appendChild(missingCell("compute", why));
        tr.appendChild(missingCell("reported", "none reported"));
        tr.appendChild(missingCell("compute", why));
        tr.appendChild(missingCell("compute", why));
        tr.appendChild(missingCell("compute", why));
      }
      tbody.appendChild(tr);
    }
    const nRun = vals.filter((v) => v.r).length;
    document.getElementById("count").textContent =
      nRun + " of " + vals.length + " tool and setting rows have results for these filters.";
  }

  // ---------- prices ----------
  function renderPrices() {
    for (const k of ["r_cpu", "r_gpu"]) {
      const el = document.getElementById(k);
      if (!el) continue;
      el.value = state[k];
      el.oninput = () => {
        const v = parseFloat(el.value);
        if (isFinite(v) && v >= 0) { state[k] = v; update(false); }
      };
    }
  }

  // ---------- score against cost ----------
  const SVGNS = "http://www.w3.org/2000/svg";
  function el(name, attrs, parent) {
    const e = document.createElementNS(SVGNS, name);
    for (const [k, v] of Object.entries(attrs || {})) e.setAttribute(k, v);
    if (parent) parent.appendChild(e);
    return e;
  }
  function marker(shape, x, y, r, color, filled, parent) {
    const common = { stroke: color, "stroke-width": 1.6, fill: filled ? color : "var(--surface)" };
    if (shape === "square") return el("rect", Object.assign({ x: x - r * 0.9, y: y - r * 0.9, width: r * 1.8, height: r * 1.8 }, common), parent);
    if (shape === "triangle") return el("path", Object.assign({ d: `M${x},${y - r * 1.1} L${x + r},${y + r * 0.75} L${x - r},${y + r * 0.75} Z` }, common), parent);
    return el("circle", Object.assign({ cx: x, cy: y, r: r }, common), parent);
  }
  function frontier(points) {
    // Points no other point beats on both axes: no other point has lower or equal cost
    // AND a higher or equal score with at least one strictly better.
    return points.filter((p) => !points.some((q) => q !== p && q.x <= p.x && q.y >= p.y && (q.x < p.x || q.y > p.y)))
      .sort((a, b) => a.x - b.x);
  }
  function renderPlot() {
    const metric = M.metrics.find((m) => m.id === state.metric) || M.metrics[0];
    const box = document.getElementById("plot");
    box.innerHTML = "";
    const pts = [], notPlotted = [];
    for (const row of M.rows) {
      const v = rowValues(state, row);
      const y = v[metric.field];
      if (!v.r) { notPlotted.push([row, "not run"]); continue; }
      if (!v.c || v.cpu == null) { notPlotted.push([row, "CPU-hours not measured"]); continue; }
      if (y == null) { notPlotted.push([row, "no value"]); continue; }
      pts.push({ x: v.cpu, y, v, row });
    }
    const W = 760, H = 440, L = 64, R = 18, T = 14, B = 52;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, role: "img", "aria-label": "Score against CPU-hours, one marker per tool and setting" }, box);
    const tip = document.createElement("div");
    tip.className = "tip";
    box.appendChild(tip);
    if (!pts.length) {
      el("text", { x: W / 2, y: H / 2, "text-anchor": "middle" }, svg).textContent = "No tool has both a score and a measured cost for these filters.";
    } else {
      const xs = pts.map((p) => p.x), ys = pts.map((p) => p.y);
      let lo = Math.floor(Math.log10(Math.min(...xs))), hi = Math.ceil(Math.log10(Math.max(...xs)));
      if (hi === lo) hi = lo + 1;
      const ymax = Math.max(...ys), ymin = Math.min(...ys);
      const pad = (ymax - ymin) * 0.08 || 0.02;
      const y0 = Math.max(0, ymin - pad), y1 = ymax + pad;
      const sx = (x) => L + ((Math.log10(x) - lo) / (hi - lo)) * (W - L - R);
      const sy = (y) => T + (1 - (y - y0) / (y1 - y0)) * (H - T - B);
      // x ticks at every power of ten, with 2 and 5 in between.
      const g = el("g", { class: "tick" }, svg);
      for (let e = lo; e <= hi; e++) {
        for (const m of [1, 2, 5]) {
          const x = m * Math.pow(10, e);
          if (Math.log10(x) > hi + 1e-9) continue;
          el("line", { x1: sx(x), x2: sx(x), y1: T, y2: H - B, class: "grid", "stroke-width": m === 1 ? 1 : 0.5 }, g);
          if (m === 1 || hi - lo <= 2) el("text", { x: sx(x), y: H - B + 16, "text-anchor": "middle" }, g).textContent = sig(x, 1);
        }
      }
      const step = niceStep((y1 - y0) / 5);
      for (let y = Math.ceil(y0 / step) * step; y <= y1 + 1e-12; y += step) {
        el("line", { x1: L, x2: W - R, y1: sy(y), y2: sy(y), class: "grid" }, g);
        el("text", { x: L - 6, y: sy(y) + 4, "text-anchor": "end" }, g).textContent = y.toFixed(step < 0.01 ? 3 : 2);
      }
      el("text", { x: (L + W - R) / 2, y: H - 10, "text-anchor": "middle", "font-size": 13 }, svg).textContent =
        "CPU-hours per 1,000 queries, measured on Sherlock (less is cheaper)";
      el("text", { x: 14, y: (T + H - B) / 2, "text-anchor": "middle", "font-size": 13, transform: `rotate(-90 14 ${(T + H - B) / 2})` }, svg).textContent =
        metric.label + " (higher is better)";
      const fr = frontier(pts);
      if (fr.length > 1) el("polyline", { points: fr.map((p) => sx(p.x) + "," + sy(p.y)).join(" "), class: "frontier" }, svg);
      for (const p of pts) {
        const t = M.tools[p.row.tool];
        const m = marker(t.shape, sx(p.x), sy(p.y), 6, M.kinds[t.kind].color, true, svg);
        m.setAttribute("tabindex", "0");
        m.setAttribute("role", "link");
        m.dataset.tool = p.row.tool;
        m.dataset.setting = p.row.setting;
        m.style.cursor = "pointer";
        m.setAttribute("aria-label", `${t.name} ${p.row.setting}: ${metric.label} ${fmt3(p.y)}, ${sig(p.x, 3)} CPU-hours per 1,000 queries`);
        const show = (ev) => {
          tip.style.display = "block";
          tip.innerHTML = `<b>${t.name}</b> <span class="mono">${p.row.setting}</span><br>${metric.label} ${fmt3(p.y)}<br>` +
            `Fmax ${fmt3(p.v.fmax)} · AUPRC ${fmt3(p.v.auprc)} · median IoU ${p.v.iou == null ? "—" : fmt3(p.v.iou)}<br>` +
            `${fmtInt(p.v.n)} true instances<br>${sig(p.x, 3)} CPU-h per 1,000 queries · peak RAM ${p.v.ram.toFixed(1)} GB` +
            (fr.includes(p) ? "<br>On the line: no other tool is both cheaper and better here." : "");
          const bb = box.getBoundingClientRect();
          const mx = (ev.clientX || bb.left + (sx(p.x) / W) * bb.width) - bb.left;
          const my = (ev.clientY || bb.top + (sy(p.y) / H) * bb.height) - bb.top;
          tip.style.left = Math.min(mx + 12, bb.width - 290) + "px";
          tip.style.top = Math.max(0, my - 10) + "px";
        };
        m.addEventListener("mousemove", show);
        m.addEventListener("focus", show);
        m.addEventListener("mouseleave", () => (tip.style.display = "none"));
        m.addEventListener("blur", () => (tip.style.display = "none"));
        const go = () => (location.href = "tools/" + p.row.tool + ".html");
        m.addEventListener("click", go);
        m.addEventListener("keydown", (e) => { if (e.key === "Enter") go(); });
      }
      box.dataset.nPlotted = pts.length;
      box.dataset.nFrontier = fr.length;
    }
    const np = document.getElementById("notplotted");
    np.innerHTML = notPlotted.length
      ? "<b>Not plotted:</b> " + notPlotted.map(([row, why]) => `${M.tools[row.tool].name} <span class="mono">${row.setting}</span> (— ${why})`).join("; ") + "."
      : "";
    document.querySelectorAll("button[data-metric]").forEach((b) => b.setAttribute("aria-pressed", b.dataset.metric === metric.id ? "true" : "false"));
  }
  function niceStep(raw) {
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= raw) return m * p;
    return 10 * p;
  }

  function describe() {
    const d = document.getElementById("filterline");
    if (!d) return;
    const sp = M.species[state.species];
    let strat = "all proteins";
    if (state.axis !== "all") {
      const row = M.chiprows.find((r) => r.axis === state.axis);
      const c = row.groups.flatMap((g) => g.chips).find((c) => c.value === state.stratum);
      strat = row.label.toLowerCase() + ": " + (c ? c.label : state.stratum);
    }
    d.textContent = `${M.truthNames[state.truth]}, ${state.calls === "dedup" ? "overlapping calls merged" : "raw calls"}, ` +
      `relative searched: ${sp.name}${sp.mya != null ? " (" + sp.mya + " Mya)" : ""}, ${strat}.`;
  }

  function update(rebuildFilters = true) {
    writeHash();
    if (rebuildFilters) renderFilters();
    describe();
    if (page === "leaderboard") renderTable();
    if (page === "cost") renderPlot();
  }

  // Exposed for the headless check: render one filter combination and return what the
  // page shows, so every combination can be compared against results.parquet.
  window.__lb = {
    meta: M,
    set(s) { Object.assign(state, s); update(); },
    state: () => Object.assign({}, state),
  };

  readHash();
  if (page === "cost") {
    document.querySelectorAll("button[data-metric]").forEach((b) =>
      b.addEventListener("click", () => { state.metric = b.dataset.metric; update(false); }));
  }
  renderPrices();
  update();
})();
