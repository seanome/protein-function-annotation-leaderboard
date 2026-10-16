/* Hover and tap text for Figures 11 to 17: any mark with a data-tip attribute inside a
   figure shows that text in a box beside it. Tap shows it on a phone; Escape or a tap
   elsewhere hides it. */
(function () {
  "use strict";
  let tip = null;
  const box = () => {
    if (!tip) {
      tip = document.createElement("div");
      tip.className = "tip fg-tip";
      tip.setAttribute("role", "status");
      document.body.appendChild(tip);
    }
    return tip;
  };
  const show = (el, x, y) => {
    const t = box();
    t.textContent = el.getAttribute("data-tip");
    t.style.display = "block";
    const w = t.offsetWidth, h = t.offsetHeight;
    const left = Math.min(window.scrollX + x + 12, window.scrollX + document.documentElement.clientWidth - w - 8);
    t.style.left = Math.max(window.scrollX + 8, left) + "px";
    t.style.top = (window.scrollY + y - h - 10 < window.scrollY ? window.scrollY + y + 14 : window.scrollY + y - h - 10) + "px";
  };
  const hide = () => { if (tip) tip.style.display = "none"; };
  const target = (e) => e.target.closest && e.target.closest("figure.fgfig [data-tip]");
  document.addEventListener("mouseover", (e) => { const el = target(e); if (el) show(el, e.clientX, e.clientY); });
  document.addEventListener("mouseout", (e) => { if (target(e)) hide(); });
  document.addEventListener("click", (e) => { const el = target(e); if (el) show(el, e.clientX, e.clientY); else hide(); });
  document.addEventListener("focusin", (e) => {
    const el = target(e);
    if (!el) return;
    const r = el.getBoundingClientRect();
    show(el, r.right, r.top);
  });
  document.addEventListener("focusout", (e) => { if (target(e)) hide(); });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") hide(); });

  /* Figure 15: a click (or Enter) on a column lists its regions under the figure, from
     data/upset.json: combination -> [[symbol, accession, Pfam, start, end, residues], ...]. */
  let lists = null;
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const listCol = async (g) => {
    const box = document.getElementById("upset-list");
    if (!box) return;
    document.querySelectorAll("g.fg-col.on").forEach((x) => x.classList.remove("on"));
    g.classList.add("on");
    if (!lists) {
      box.innerHTML = '<p class="muted">Loading…</p>';
      try { lists = await (await fetch("data/upset.json")).json(); }
      catch (err) { box.innerHTML = '<p class="muted">Could not load data/upset.json.</p>'; return; }
    }
    const rows = lists[g.getAttribute("data-combo")] || [];
    const head = g.getAttribute("data-tip").replace(/ Click to list them\.$/, "");
    box.innerHTML = `<p><b>${esc(head)}</b></p><div class="tablebox"><table class="fg-regions"><thead><tr>` +
      "<th>gene</th><th>UniProt</th><th>Pfam</th><th>residues</th><th>residues of the true region</th></tr></thead><tbody>" +
      rows.map((r) => `<tr><td>${esc(r[0])}</td><td class="mono">${esc(r[1])}</td><td class="mono">${esc(r[2])}</td>` +
        `<td class="num">${r[3].toLocaleString("en-US")}–${r[4].toLocaleString("en-US")}</td>` +
        `<td class="mono seq">${esc(r[5])}</td></tr>`).join("") + "</tbody></table></div>";
  };
  document.addEventListener("click", (e) => {
    const g = e.target.closest && e.target.closest("g.fg-col");
    if (g) listCol(g);
  });
  document.addEventListener("keydown", (e) => {
    const g = e.target.closest && e.target.closest("g.fg-col");
    if (g && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); listCol(g); }
  });
})();
