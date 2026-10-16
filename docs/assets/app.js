/* Leaderboard and compare-two-columns pages. Reads the JSON that build_site.py inlines
   into the page (#lb-data); never computes a metric, only looks rows up. A (tool,
   setting) with no row for the chosen filters is drawn as "—" and the words "not run";
   a missing compute row is "—" and "not measured". */
(function () {
  "use strict";
  const DATA = JSON.parse(document.getElementById("lb-data").textContent);
  const M = DATA.meta;
  const RESULTS = DATA.results; // key -> {tool|setting: [fmax, auprc, median IoU, true instances]}
  const COMPUTE = DATA.compute; // tool|setting|species -> {...}
  const page = (document.querySelector("main[data-page]") || document.body).dataset.page;

  // ---------- state, kept in the URL hash so both pages share the filters ----------
  const state = Object.assign({}, M.defaults);
  // The values each hash key may take, from the page's own data. Anything else in the hash
  // is ignored and the default kept, so a hand-edited or stale link cannot break the page.
  function chipValues(id) {
    const row = M.chiprows.find((r) => r.id === id);
    return row ? row.groups.flatMap((g) => g.chips.map((c) => c.value)) : [];
  }
  function readHash() {
    const h = new URLSearchParams(location.hash.slice(1));
    const ok = {
      truth: (v) => chipValues("truth").includes(v) && v in M.truthNames,
      calls: (v) => chipValues("calls").includes(v),
      species: (v) => chipValues("species").includes(v) && v in M.species,
      x: (v) => M.axes.some((a) => a.id === v),
      y: (v) => M.axes.some((a) => a.id === v && a.kind === "score"),
      sort: (v) => COLS.some((c) => c.id === v),
      dir: (v) => v === "asc" || v === "desc",
    };
    for (const k of Object.keys(ok)) {
      if (h.has(k) && ok[k](h.get(k))) state[k] = h.get(k);
    }
    // Axis and stratum go together: the stratum must be a chip of that axis's row.
    const axis = h.get("axis"), stratum = h.get("stratum");
    const row = M.chiprows.find((r) => r.kind === "stratum" && r.axis === axis);
    if (axis === "all" || (h.has("axis") && !row)) { state.axis = "all"; state.stratum = "all"; }
    else if (row && stratum !== "all" && row.groups.some((g) => g.chips.some((c) => c.value === stratum))) {
      state.axis = axis; state.stratum = stratum;
    }
    // Links made before both axes could be chosen name the score as "metric".
    if (h.has("metric") && !h.has("y") && ok.y(h.get("metric"))) state.y = h.get("metric");
    for (const k of ["r_cpu", "r_gpu"]) {
      const v = Number(h.get(k));
      if (h.has(k) && h.get(k) !== "" && isFinite(v) && v >= 0) state[k] = v;
    }
  }
  function writeHash() {
    const h = new URLSearchParams();
    for (const k of ["truth", "calls", "species", "axis", "stratum", "x", "y", "sort", "dir", "r_cpu", "r_gpu"]) {
      if (state[k] !== undefined && state[k] !== null) h.set(k, state[k]);
    }
    try { history.replaceState(null, "", "#" + h.toString()); } catch (e) { /* framed viewers may refuse */ }
    // Keep the filters when moving between the leaderboard and the plot.
    document.querySelectorAll("a[data-keep-filters]").forEach((a) => {
      a.href = a.getAttribute("href").split("#")[0] + "#" + h.toString();
    });
  }

  // truth|axis|stratum for every stratum that has a row under some relative and dedup choice.
  const STRATA = new Set(Object.keys(RESULTS).map((k) => { const p = k.split("|"); return p[0] + "|" + p[3] + "|" + p[4]; }));
  function stratumExists(truth, axis, stratum) {
    return axis === "all" || STRATA.has(truth + "|" + axis + "|" + stratum);
  }
  // Colour of a kind of method, read from the page's CSS tokens at draw time so that it
  // matches the legends in light and dark mode.
  const rootStyle = getComputedStyle(document.documentElement);
  function kindColor(kind) {
    return rootStyle.getPropertyValue("--k-" + kind).trim() || M.kinds[kind].color;
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
    return x.toFixed(d);
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
      a.innerHTML = row.label + (M.terms && M.terms.filters[row.id] ? " " + M.terms.filters[row.id] : "") +
        (row.help ? "<small>" + row.help + "</small>" : "");
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
            chip(c.empty ? c.label + " (not run)" : c.label, pressed, () => {
              if (row.kind === "stratum") {
                // The run cut the queries along one stratum axis at a time, so picking a
                // value here clears the other three stratum rows.
                if (c.value === "all") {
                  if (state.axis === row.axis) { state.axis = "all"; state.stratum = "all"; }
                } else { state.axis = row.axis; state.stratum = c.value; }
              } else {
                state[row.id] = c.value;
                // A stratum the new truth set does not have (a Swiss-Prot feature type under
                // Pfam) would leave every row "not run"; go back to all proteins instead.
                if (row.id === "truth" && !stratumExists(state.truth, state.axis, state.stratum)) {
                  state.axis = "all"; state.stratum = "all";
                }
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
    { id: "fmax", label: "Fmax", key: "fmax", term: "fmax" },
    { id: "auprc", label: "AUPRC", key: "auprc", term: "auprc" },
    { id: "iou", label: "Median IoU", key: "median_iou_tp", term: "iou" },
    { id: "n", label: "True instances", key: "n_truth_instances", term: "n" },
    { id: "cpu", label: 'Core-hours per 1,000 queries, per relative <span class="prov measured">measured</span>', key: "cpu_h", term: "cpu" },
    { id: "cpu_rep", label: 'CPU-h per 1,000 queries <span class="prov reported">reported</span>', key: "cpu_h_reported", term: "cpu_rep" },
    { id: "ram", label: 'Peak RAM, GB <span class="prov measured">measured</span>', key: "peak_rss_gb", term: "ram" },
    { id: "gpu", label: "GPU in search", key: "gpu", term: "gpu" },
    { id: "cost", label: "Estimated cost, USD per 1,000 queries", key: "cost_usd_est", term: "cost" },
  ];
  function rowValues(s, row) {
    const r = lookup(s, row);
    const c = computeRow(s, row);
    return {
      row, r, c,
      // A stored median IoU of 0 means no correct call on a region longer than one residue
      // (any such call has IoU >= 0.5), so it is shown as words and sorts with the missing.
      fmax: r ? r[0] : null, auprc: r ? r[1] : null, iou: r && r[2] ? r[2] : null, n: r ? r[3] : null,
      cpu: c ? c.cpu_h_per_1000 : null,
      cpu_rep: c ? c.cpu_h_reported_per_1000 : null,
      ram: c ? c.peak_rss_gb : null,
      gpu: c ? c.gpu : null,
      cost: cost(c),
    };
  }
  // Words for a median IoU of 0. Fmax above 0 means at least one call was correct; with no
  // IoU those correct calls were all on single-residue annotations.
  function iouWords(v) {
    return v.fmax > 0 ? "correct calls only on single residues" : "no correct calls";
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
  // The caption's first words, from the current sort: "Tools ranked by Fmax, highest first".
  const SORT_NAMES = {
    tool: "tool name", setting: "setting name", fmax: "Fmax", auprc: "AUPRC", iou: "median IoU",
    n: "true instances", cpu: "measured core-hours per 1,000 queries", cpu_rep: "reported CPU-hours per 1,000 queries",
    ram: "peak RAM", gpu: "GPU in search", cost: "estimated cost",
  };
  function sortWords() {
    const asc = state.dir === "asc";
    if (state.sort === "tool" || state.sort === "setting" || state.sort === "gpu") {
      return `Tools sorted by ${SORT_NAMES[state.sort]}, ${asc ? "A to Z" : "Z to A"}`;
    }
    return `Tools ranked by ${SORT_NAMES[state.sort] || state.sort}, ${asc ? "lowest" : "highest"} first`;
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
      // The "?" sits beside the sort button, not inside it, so it never sorts.
      if (col.term && M.terms) th.insertAdjacentHTML("beforeend", M.terms.th[col.term]);
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
      const tags = (t.tags || []).map((x) => '<span class="tag">' + x + "</span>");
      if (v.row.tool === "hmmscan" && state.truth === "pfam") tags.push('<span class="tag warn">graded on its own models</span>');
      td1.innerHTML = '<span class="kind k-' + t.kind + '" aria-hidden="true"></span><a href="tools/' + v.row.tool + '.html">' + t.name + "</a>" +
        (tags.length ? "<br>" + tags.join(" ") : "");
      tr.appendChild(td1);
      const td2 = document.createElement("td");
      td2.className = "txt mono";
      td2.textContent = v.row.setting;
      td2.title = v.row.words;
      tr.appendChild(td2);
      if (v.r) {
        tr.appendChild(numCell(fmt3(v.fmax)));
        tr.appendChild(numCell(fmt3(v.auprc)));
        if (v.iou == null) {
          const td = document.createElement("td");
          td.className = "none";
          td.dataset.none = "iou";
          td.innerHTML = "<i>" + iouWords(v) + "</i>";
          td.title = v.fmax > 0 ? "an IoU against one residue says nothing about where a call's ends are" : "no call was correct";
          tr.appendChild(td);
        } else tr.appendChild(numCell(fmt3(v.iou)));
        tr.appendChild(numCell(fmtInt(v.n)));
      } else {
        tr.appendChild(missingCell("results", "not run", 4));
      }
      if (v.c) {
        const cc = numCell(sig(v.cpu, 3));
        if (v.c.n_runs > 1) cc.innerHTML += ` <span class="muted runs">(${sig(v.c.cpu_min, 3)}–${sig(v.c.cpu_max, 3)}, ${v.c.n_runs} runs)</span>`;
        tr.appendChild(cc);
        tr.appendChild(v.cpu_rep == null ? missingCell("reported", "none reported") : numCell(sig(v.cpu_rep, 3)));
        tr.appendChild(numCell(v.ram.toFixed(1)));
        const g = numCell(v.c.gpu_h_per_1000 ? sig(v.c.gpu_h_per_1000, 3) + " GPU-h" : "none");
        if (v.c.gpu_h_per_1000) g.title = "one GPU, while predicting 3Di letters for the queries";
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
    const rank = document.getElementById("rankline");
    if (rank) rank.textContent = sortWords();
    const nRun = vals.filter((v) => v.r).length;
    document.getElementById("count").textContent = nRun
      ? nRun + " of " + vals.length + " tool and setting rows have results for these filters."
      : "No tool has results for these filters: either this stratum holds no regions for this relative, or it was not scored.";
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

  // ---------- Figure 1: any two columns against each other ----------
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
    if (shape === "diamond") return el("path", Object.assign({ d: `M${x},${y - r * 1.2} L${x + r * 1.2},${y} L${x},${y + r * 1.2} L${x - r * 1.2},${y} Z` }, common), parent);
    if (shape === "triangle") return el("path", Object.assign({ d: `M${x},${y - r * 1.1} L${x + r},${y + r * 0.75} L${x - r},${y + r * 0.75} Z` }, common), parent);
    return el("circle", Object.assign({ cx: x, cy: y, r: r }, common), parent);
  }
  function frontier(points) {
    // Only for a lower-is-better x against a higher-is-better y. A point is on the line when
    // no other point has lower or equal x AND higher or equal y, with one of them strictly.
    return points.filter((p) => !points.some((q) => q !== p && q.x <= p.x && q.y >= p.y && (q.x < p.x || q.y > p.y)))
      .sort((a, b) => a.x - b.x);
  }
  function axisDef(id, fallback) {
    return M.axes.find((a) => a.id === id) || M.axes.find((a) => a.id === fallback);
  }
  // The value a row has on one axis, or the reason it has none ("— not run" and so on).
  function axisValue(ax, v) {
    if (ax.kind === "score") {
      if (!v.r) return { why: "not run" };
      const y = v[ax.id];
      if (y == null) return { why: ax.id === "iou" ? "no correct calls" : "no " + ax.short };
      return { val: y };
    }
    if (!v.c) {
      const searched = M.searched.includes(v.row.tool + "|" + v.row.setting + "|" + state.species);
      return { why: searched ? ax.short + " not measured" : "not run" };
    }
    const val = { cpu: v.cpu, cost: v.cost, wall: v.c.wall_h_per_1000, ram: v.ram, db: v.c.db_cpu_h_per_db }[ax.id];
    if (val == null) return { why: ax.id === "db" ? "no database build in the run's records" : ax.short + " not measured" };
    // An axis where each step is ten times the last has no place for zero.
    if (val <= 0) return { why: ax.short + " is zero, which this axis cannot show" };
    return { val };
  }
  function fmtAxis(ax, x) {
    if (ax.kind === "score") return fmt3(x);
    if (ax.id === "cost") return "$" + sig(x, 3);
    return sig(x, 3);
  }
  function unitWords(ax) {
    return { cpu: "core-hours", cost: "cost", wall: "time", ram: "memory", db: "build time" }[ax.id] || "";
  }
  // Split a label into lines of at most maxChars characters, at spaces.
  function wrap(text, maxChars) {
    const lines = [];
    let cur = "";
    for (const w of text.split(" ")) {
      if (cur && (cur + " " + w).length > maxChars) { lines.push(cur); cur = w; } else cur = cur ? cur + " " + w : w;
    }
    if (cur) lines.push(cur);
    return lines;
  }
  function scale(ax, vals, lo0, hi0) {
    // Cost, time and memory: each step is ten times the last. Scores: linear, padded.
    if (ax.kind !== "score") {
      let lo = Math.floor(Math.log10(Math.min(...vals))), hi = Math.ceil(Math.log10(Math.max(...vals)));
      if (hi === lo) hi = lo + 1;
      const f = (x) => lo0 + ((Math.log10(x) - lo) / (hi - lo)) * (hi0 - lo0);
      const ticks = [];
      for (let e = lo; e <= hi; e++) {
        for (const m of [1, 2, 5]) {
          const x = m * Math.pow(10, e);
          if (Math.log10(x) > hi + 1e-9) continue;
          ticks.push({ x, major: m === 1, label: m === 1 || hi - lo <= 2 ? (ax.id === "cost" ? "$" : "") + sig(x, 1) : null });
        }
      }
      return { f, ticks };
    }
    const vmax = Math.max(...vals), vmin = Math.min(...vals);
    const pad = (vmax - vmin) * 0.08 || 0.02;
    const a = Math.max(0, vmin - pad), b = vmax + pad;
    const f = (x) => lo0 + ((x - a) / (b - a)) * (hi0 - lo0);
    const step = niceStep((b - a) / 5);
    const ticks = [];
    for (let x = Math.ceil(a / step - 1e-9) * step; x <= b + 1e-12; x += step) {
      ticks.push({ x, major: true, label: x.toFixed(step < 0.01 ? 3 : 2) });
    }
    return { f, ticks };
  }
  let plotted = [];
  function renderAxisPick() {
    const box = document.getElementById("axispick");
    if (!box) return;
    box.innerHTML = "";
    for (const [key, label, opts] of [
      ["y", "Vertical axis, a score:", M.axes.filter((a) => a.kind === "score")],
      ["x", "Horizontal axis:", M.axes],
    ]) {
      const r = document.createElement("div");
      r.className = "pickrow";
      r.setAttribute("role", "group");
      r.setAttribute("aria-label", label);
      r.dataset.axis = key;
      const s = document.createElement("span");
      s.className = "picklabel";
      s.textContent = label;
      r.appendChild(s);
      const chips = document.createElement("span");
      chips.className = "chips";
      for (const ax of opts) {
        const other = key === "x" ? state.y : state.x;
        const b = chip(ax.short.charAt(0).toUpperCase() + ax.short.slice(1), state[key] === ax.id, () => {
          state[key] = ax.id;
          update(false);
        }, ax.label);
        b.dataset.pick = ax.id;
        // The same column on both axes would put every marker on the diagonal.
        if (ax.id === other) { b.disabled = true; b.title = "already on the other axis"; }
        chips.appendChild(b);
      }
      r.appendChild(chips);
      box.appendChild(r);
    }
  }
  function renderPlot() {
    if (!M.axes.some((a) => a.id === state.y && a.kind === "score")) state.y = "fmax";
    const ya = axisDef(state.y, "fmax");
    let xa = axisDef(state.x, "cpu");
    if (xa.id === ya.id) xa = axisDef("cpu");
    state.x = xa.id;
    renderAxisPick();
    const box = document.getElementById("plot");
    box.innerHTML = "";
    const pts = [], notPlotted = [];
    for (const row of M.rows) {
      const v = rowValues(state, row);
      const py = axisValue(ya, v), px = axisValue(xa, v);
      // A row with no result has no score and no cost to speak of: one reason, "not run".
      const why = !v.r ? "not run" : py.why || px.why;
      if (why) { notPlotted.push([row, why]); continue; }
      pts.push({ x: px.val, y: py.val, v, row });
    }
    const pareto = xa.better === "lower" && ya.better === "higher";
    // Draw at the box's own width so the 12-13 px text stays that size on a phone.
    const W = Math.max(300, Math.min(900, (box.clientWidth || 760) - 14));
    const narrow = W < 560;
    const H = narrow ? Math.round(W * 1.0) : 440;
    const maxChars = Math.floor((W - 40) / 7);
    const xlines = wrap(xa.label + ` (${xa.better} is better)`, maxChars);
    const L = narrow ? 52 : 64, R = 14, T = 12, B = 34 + 16 * xlines.length;
    const svg = el("svg", { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: "img",
      "aria-label": `${ya.short} against ${xa.short}, one marker per tool and setting` }, box);
    const tip = document.createElement("div");
    tip.className = "tip";
    box.appendChild(tip);
    let fr = [];
    if (!pts.length) {
      for (const [i, line] of wrap(`No tool has both ${ya.short} and ${xa.short} for these filters.`, maxChars).entries())
        el("text", { x: W / 2, y: H / 2 + i * 18, "text-anchor": "middle" }, svg).textContent = line;
    } else {
      const sx = scale(xa, pts.map((p) => p.x), L, W - R);
      const sy = scale(ya, pts.map((p) => p.y), H - B, T);
      const g = el("g", { class: "tick" }, svg);
      let lastRight = -1e9;
      for (const t of sx.ticks) {
        const X = sx.f(t.x);
        el("line", { x1: X, x2: X, y1: T, y2: H - B, class: "grid", "stroke-width": t.major ? 1 : 0.5 }, g);
        // Skip a tick label that would run into the one before it.
        if (t.label && X - t.label.length * 3.6 > lastRight + 4) {
          el("text", { x: X, y: H - B + 16, "text-anchor": "middle" }, g).textContent = t.label;
          lastRight = X + t.label.length * 3.6;
        }
      }
      for (const t of sy.ticks) {
        const Y = sy.f(t.x);
        el("line", { x1: L, x2: W - R, y1: Y, y2: Y, class: "grid", "stroke-width": t.major ? 1 : 0.5 }, g);
        if (t.label) el("text", { x: L - 6, y: Y + 4, "text-anchor": "end" }, g).textContent = t.label;
      }
      xlines.forEach((line, i) => {
        el("text", { x: (L + W - R) / 2, y: H - B + 36 + i * 16, "text-anchor": "middle", "font-size": 13 }, svg).textContent = line;
      });
      const ymid = (T + H - B) / 2;
      el("text", { x: 14, y: ymid, "text-anchor": "middle", "font-size": 13, transform: `rotate(-90 14 ${ymid})` }, svg).textContent =
        ya.label + " (higher is better)";
      fr = pareto ? frontier(pts) : [];
      if (fr.length > 1) el("polyline", { points: fr.map((p) => sx.f(p.x) + "," + sy.f(p.y)).join(" "), class: "frontier" }, svg);
      for (const p of pts) {
        const t = M.tools[p.row.tool];
        const X = sx.f(p.x), Y = sy.f(p.y);
        const m = marker(p.row.shape, X, Y, narrow ? 5 : 6, kindColor(t.kind), true, svg);
        m.setAttribute("tabindex", "0");
        m.setAttribute("role", "link");
        m.dataset.tool = p.row.tool;
        m.dataset.setting = p.row.setting;
        m.dataset.x = String(p.x);
        m.dataset.y = String(p.y);
        m.style.cursor = "pointer";
        m.setAttribute("aria-label", `${t.name} ${p.row.setting}: ${ya.short} ${fmtAxis(ya, p.y)}, ${xa.short} ${fmtAxis(xa, p.x)}`);
        const show = (ev) => {
          const c = p.v.c;
          tip.style.display = "block";
          tip.innerHTML = `<b>${t.name}</b> <span class="mono">${p.row.setting}</span><br>` +
            `Fmax ${fmt3(p.v.fmax)} · AUPRC ${fmt3(p.v.auprc)} · median IoU ${p.v.iou == null ? "— no correct calls" : fmt3(p.v.iou)}<br>` +
            `${fmtInt(p.v.n)} true instances<br>` +
            (c ? `${sig(c.cpu_h_per_1000, 3)} core-hours · $${sig(p.v.cost, 3)} · ${sig(c.wall_h_per_1000, 3)} wall-clock hours on ${c.cores} cores, per 1,000 queries per relative<br>` +
              `peak RAM ${p.v.ram.toFixed(1)} GB · database build ` +
              (c.db_cpu_h_per_db == null ? "— no record" : `${sig(c.db_cpu_h_per_db, 3)} core-hours per database (${c.db_n} recorded)`) +
              (c.gpu_h_per_1000 ? `<br>plus ${sig(c.gpu_h_per_1000, 3)} GPU-h` : "")
              : "cost — not measured") +
            (fr.includes(p) ? `<br>On the line: no other tool is both ${xa.less} and better.` : "");
          const bb = box.getBoundingClientRect();
          const mx = (ev.clientX || bb.left + (X / W) * bb.width) - bb.left;
          const my = (ev.clientY || bb.top + (Y / H) * bb.height) - bb.top;
          tip.style.left = Math.max(0, Math.min(mx + 12, bb.width - 290)) + "px";
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
    }
    box.dataset.nPlotted = pts.length;
    box.dataset.nFrontier = fr.length;
    plotted = pts.map((p) => ({ tool: p.row.tool, setting: p.row.setting, x: p.x, y: p.y, frontier: fr.includes(p) }));

    // Caption and key follow the chosen axes.
    const fk = document.getElementById("frontier-key");
    if (fk) fk.hidden = !(pareto && fr.length > 1);
    // Filled or open says where a cost came from; with two scores no cost is drawn.
    document.querySelectorAll("[data-cost-key]").forEach((e) => (e.hidden = xa.kind === "score"));
    const fw = document.getElementById("frontier-words");
    if (fw && pareto) fw.textContent = `the line: no other tool is both ${xa.less} and better on ${ya.short}`;
    const cap = document.getElementById("plotcaption");
    if (cap) {
      let words = `${ya.short.charAt(0).toUpperCase() + ya.short.slice(1)} against ${xa.short}. Up is better`;
      words += xa.better === "lower" ? `; left is ${xa.less}.` : "; right is better.";
      if (xa.kind !== "score") words += ` Each step along the horizontal axis is ten times the ${unitWords(xa)} of the last.`;
      if (!pareto) words += " No line is drawn: both axes are scores, where higher is better.";
      else if (fr.length === 1) {
        const p = fr[0];
        words += ` No line is drawn: one marker, ${M.tools[p.row.tool].name} ${p.row.setting}, beats or ties every other on both axes.`;
      }
      if (xa.kind === "score") words += " The two scores can rank tools in opposite orders.";
      cap.textContent = words;
    }
    const np = document.getElementById("notplotted");
    np.innerHTML = notPlotted.length
      ? "<b>Not plotted:</b> " + notPlotted.map(([row, why]) =>
        `<span class="np" data-tool="${row.tool}" data-setting="${row.setting}">${M.tools[row.tool].name} <span class="mono">${row.setting}</span> (— ${why})</span>`).join("; ") + "."
      : "";
  }
  function niceStep(raw) {
    const p = Math.pow(10, Math.floor(Math.log10(raw)));
    for (const m of [1, 2, 2.5, 5, 10]) if (m * p >= raw) return m * p;
    return 10 * p;
  }

  function describe() {
    const d = document.getElementById("filterline");
    const short = document.getElementById("filtersummary");
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
    if (short) short.textContent = "now: " + d.textContent;
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
    // What Figure 1 drew: one entry per marker, with the raw x and y values.
    plotted: () => plotted.map((p) => Object.assign({}, p)),
  };

  readHash();
  // A link with a new hash, or back and forward, changes the filters without a reload.
  window.addEventListener("hashchange", () => {
    Object.assign(state, M.defaults);
    readHash();
    renderPrices();
    update();
  });
  // Redraw when the colour theme changes, so the plot takes the new kind colours.
  const redraw = () => update(false);
  if (window.matchMedia) {
    const mq = window.matchMedia("(prefers-color-scheme: dark)");
    if (mq.addEventListener) mq.addEventListener("change", redraw);
  }
  new MutationObserver(redraw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
  if (page === "cost") {
    // Redraw at the new width, so the text keeps its size on a phone.
    let pending = null;
    window.addEventListener("resize", () => {
      clearTimeout(pending);
      pending = setTimeout(() => renderPlot(), 150);
    });
  }
  renderPrices();
  update();
})();
