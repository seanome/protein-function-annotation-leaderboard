// Sorting, filters and search for the index table. The filter and the search words live in the
// page address (#none, #only-kmerseek, #agree&q=MHC ...), so a link to a filtered table can be shared.
(function () {
  const tbody = document.querySelector("#t tbody");
  const rows = Array.from(tbody.rows);
  const count = document.getElementById("count");
  const buttons = document.querySelectorAll(".filters button, .side a[data-f]");
  function keep(r, f) {
    if (f === "all") return true;
    if (f.startsWith("only-")) return r.dataset.s === "one" && r.dataset.f === f.slice(5);
    return r.dataset.s === f;
  }
  const S = window.SEARCH || { programs: [], text: [], rows: [] };
  const lower = S.text.map((t) => t.toLowerCase());
  rows.forEach((r, i) => { r._pairs = S.rows[i] || []; r._id = r.cells[0].querySelector("a").textContent.toLowerCase(); });
  const box = document.getElementById("q");
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  // Address: #<filter>&q=<search words>&e=<log10 of the E-value limit>, each part optional.
  function state() {
    const h = location.hash.slice(1), i = h.indexOf("&");
    const f = decodeURIComponent(i < 0 ? h : h.slice(0, i)) || "all";
    const ps = new URLSearchParams(i < 0 ? "" : h.slice(i + 1));
    return { f, q: ps.get("q") || "", e: ((v) => (Number.isNaN(v) ? 1 : Math.min(1, v)))(parseInt(ps.get("e") ?? "1", 10)), u: ps.get("u") === "1" };
  }
  function setState(f, q, e, u) {
    const ps = new URLSearchParams();
    if (q) ps.set("q", q);
    if (e !== 1) ps.set("e", e);
    if (u) ps.set("u", "1");
    const rest = ps.toString();
    history.replaceState(null, "", "#" + encodeURIComponent(f) + (rest ? "&" + rest : "")); apply();
  }
  const slider = document.getElementById("emax"), ev = document.getElementById("ev"), unnamed = document.getElementById("unnamed");
  const bestE = (r) => parseFloat(r.cells[4].dataset.v);
  const progs = (mask) => S.programs.filter((_, b) => mask & (1 << b)).join(", ");
  // Every word must start a word in the protein id or in some text of its hits ("ring" finds
  // "RING finger", not "during"). Returns the texts
  // that matched (with the programs whose hits carry them), [] for an id-only match, or null.
  function matches(r, words) {
    const found = [];
    for (const w of words) {
      const re = new RegExp("(^|[^a-z0-9])" + w.replace(/[.*+?^${}()|[\]\\]/g, "\\$&"));
      let any = re.test(r._id);
      for (let j = 0; j < r._pairs.length; j += 2) {
        const t = lower[r._pairs[j]];
        if (t.includes(w) && re.test(t)) {
          any = true;
          if (found.length < 40 && !found.some((x) => x[0] === r._pairs[j])) found.push([r._pairs[j], r._pairs[j + 1]]);
        }
      }
      if (!any) return null;
    }
    return found;
  }
  function apply() {
    const { f, q, e, u } = state();
    if (unnamed) unnamed.checked = u;
    if (box && box.value !== q && document.activeElement !== box) box.value = q;
    if (slider && +slider.value !== e) slider.value = e;
    if (ev) ev.textContent = e === 1 ? "10" : e === 0 ? "1" : "1e" + e;
    const maxE = Math.pow(10, e);
    const words = q.toLowerCase().split(/\s+/).filter(Boolean);
    document.getElementById("t").classList.toggle("searching", words.length > 0);
    let n = 0;
    for (const r of rows) {
      let k = keep(r, f) && (e === 1 || bestE(r) <= maxE) && (!u || r.dataset.n === "0");
      const cell = r.cells[r.cells.length - 1];
      if (k && words.length) {
        const m = matches(r, words);
        if (m === null) k = false;
        else cell.innerHTML = m.length ? m.slice(0, 4).map(([t, mask]) => `<b>${esc(S.text[t])}</b> (${mask ? esc(progs(mask)) : "the protein's name"})`).join("<br>") +
          (m.length > 4 ? `<br>and ${m.length - 4}${m.length >= 40 ? "+" : ""} more` : "") : "the protein id";
      }
      r.hidden = !k; if (k) n++;
    }
    buttons.forEach((b) => b.classList.toggle("on", b.dataset.f === f));
    const b = Array.from(buttons).find((b) => b.dataset.f === f);
    count.textContent = `Showing ${n.toLocaleString()} of ${rows.length.toLocaleString()} proteins` +
      (b && f !== "all" ? `: ${f.startsWith("only-") ? "only " + b.querySelector(".pn").textContent + " has a hit" : b.textContent.toLowerCase()}` : "") + (e !== 1 ? `${b && f !== "all" ? "," : ":"} best E-value at most ${e === 0 ? "1" : "1e" + e}` : "") +
      (u ? `${(b && f !== "all") || e !== 1 ? "," : ":"} no gene name in the ID` : "") +
      (q ? `${(b && f !== "all") || e !== 1 || u ? "," : ":"} matching “${q}”.` : ".");
  }
  buttons.forEach((b) => b.addEventListener("click", (e) => { e.preventDefault(); setState(b.dataset.f, state().q, state().e, state().u); }));
  if (box) box.addEventListener("input", () => setState(state().f, box.value.trim(), state().e, state().u));
  if (slider) slider.addEventListener("input", () => setState(state().f, state().q, +slider.value, state().u));
  if (unnamed) unnamed.addEventListener("change", () => setState(state().f, state().q, state().e, unnamed.checked));
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
