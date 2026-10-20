/* Figure 1 of the "Why this project" page: one tile per protein, lit tiles first, then dark,
   COLS tiles per row in both maps so a tile is the same size in both. Hover or tap a tile to
   see the protein. Data: data/why_map.json, written by scripts/why_page.py. */
(function () {
  "use strict";
  const canvases = [...document.querySelectorAll("figure.maps canvas[data-species]")];
  if (!canvases.length) return;
  const tip = document.getElementById("maptip");
  let data = null;

  const colour = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  function draw() {
    if (!data) return;
    const cols = data.cols;
    // One tile size for both maps: the narrower canvas sets it.
    const w = Math.min(...canvases.map((c) => c.parentElement.clientWidth));
    const tile = w / cols;
    const dpr = window.devicePixelRatio || 1;
    const lit = colour("--lit"), unlit = colour("--unlit");
    for (const c of canvases) {
      const sp = data.species[c.dataset.species];
      const rows = Math.ceil(sp.proteins.length / cols);
      c.style.width = w + "px";
      c.style.height = rows * tile + "px";
      c.width = Math.round(w * dpr);
      c.height = Math.round(rows * tile * dpr);
      const g = c.getContext("2d");
      g.setTransform(dpr, 0, 0, dpr, 0, 0);
      g.clearRect(0, 0, w, rows * tile);
      sp.proteins.forEach((p, i) => {
        g.fillStyle = p[3] || p[4] ? lit : unlit;
        const x = (i % cols) * tile, y = Math.floor(i / cols) * tile;
        // Snap to device pixels so neighbouring tiles leave no seams.
        g.fillRect(x, y, tile + 0.5 / dpr, tile + 0.5 / dpr);
      });
      c._tile = tile;
    }
  }

  function at(c, ev) {
    const r = c.getBoundingClientRect();
    const col = Math.floor((ev.clientX - r.left) / c._tile), row = Math.floor((ev.clientY - r.top) / c._tile);
    if (col < 0 || col >= data.cols || row < 0) return null;
    return data.species[c.dataset.species].proteins[row * data.cols + col] || null;
  }

  function show(c, ev) {
    if (!data) return;
    const p = at(c, ev);
    if (!p) { tip.hidden = true; return; }
    const labels = [p[3] ? "an InterPro family" : null, p[4] ? "another Swiss-Prot protein at least 30% identical" : null].filter(Boolean);
    tip.textContent = "";
    const b = document.createElement("b");
    b.textContent = p[1] ? `${p[1]} (${p[0]})` : p[0];
    const name = document.createElement("div");
    name.textContent = p[2] || "no name in UniProt";
    const lab = document.createElement("div");
    lab.className = "muted";
    lab.textContent = labels.length ? "lit: " + labels.join(", ") : "dark: no InterPro family, no other Swiss-Prot protein at least 30% identical";
    tip.append(b, name, lab);
    tip.hidden = false;
    const fig = c.closest("figure").getBoundingClientRect();
    const x = Math.min(ev.clientX - fig.left + 12, fig.width - tip.offsetWidth - 4);
    tip.style.left = Math.max(0, x) + "px";
    tip.style.top = ev.clientY - fig.top + 14 + "px";
  }

  for (const c of canvases) {
    c.addEventListener("mousemove", (ev) => show(c, ev));
    c.addEventListener("click", (ev) => show(c, ev));
    c.addEventListener("mouseleave", () => { tip.hidden = true; });
  }
  window.addEventListener("resize", draw);
  window.matchMedia("(prefers-color-scheme: dark)").addEventListener("change", draw);
  new MutationObserver(draw).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });

  fetch("data/why_map.json").then((r) => r.json()).then((d) => {
    data = d;
    draw();
    window.__whyMap = { ready: true, cols: d.cols, counts: Object.fromEntries(Object.entries(d.species).map(([k, v]) => [k, { n: v.proteins.length, lit: v.proteins.filter((p) => p[3] || p[4]).length }])) };
  });
})();
