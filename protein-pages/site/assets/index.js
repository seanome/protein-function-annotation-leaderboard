// Sorting and filters for the index table. The filter lives in the page address (#none,
// #one, #only-kmerseek ...), so a link to a filtered table can be shared.
(function () {
  const tbody = document.querySelector("#t tbody");
  const rows = Array.from(tbody.rows);
  const count = document.getElementById("count");
  const buttons = document.querySelectorAll(".filters button");
  function keep(r, f) {
    if (f === "all") return true;
    if (f.startsWith("only-")) return r.dataset.s === "one" && r.dataset.f === f.slice(5);
    return r.dataset.s === f;
  }
  function apply() {
    const f = decodeURIComponent(location.hash.slice(1)) || "all";
    let n = 0;
    for (const r of rows) { const k = keep(r, f); r.hidden = !k; if (k) n++; }
    buttons.forEach((b) => b.classList.toggle("on", b.dataset.f === f));
    const b = Array.from(buttons).find((b) => b.dataset.f === f);
    count.textContent = `Showing ${n.toLocaleString()} of ${rows.length.toLocaleString()} proteins` + (b && f !== "all" ? `: ${b.textContent.toLowerCase()}.` : ".");
  }
  buttons.forEach((b) => b.addEventListener("click", () => { location.hash = b.dataset.f; }));
  addEventListener("hashchange", apply);
  apply();
  let last = null, dir = 1;
  document.querySelectorAll("#t th[data-k]").forEach((th) => th.addEventListener("click", () => {
    const k = +th.dataset.k, str = th.dataset.t === "s";
    dir = last === k ? -dir : 1; last = k;
    const val = (r) => { const c = r.cells[k]; return str ? c.textContent : parseFloat(c.dataset.v); };
    rows.sort((a, b) => { const x = val(a), y = val(b); return (x < y ? -1 : x > y ? 1 : 0) * dir; });
    for (const r of rows) tbody.appendChild(r);
  }));
})();
