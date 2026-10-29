// The index of every protein (build_full.py): search, the count buttons and bars as filters, sorting,
// and the table shown 100 rows at a time.
(function () {
  var F = window.FULL, R = F.rows, CATS = F.cats;
  // row: 0 id, 1 length, 2 existing, 3 proposed, 4 name check, 5 Swiss-Prot name, 6 accession,
  // 7 identity, 8 coverage, 9 E-value, 10 Pfam names, 11 category bits, 12 Fester/BHF match,
  // 13 cluster size, 14 drawn page
  function catNames(m) { var o = []; for (var i = 0; i < CATS.length; i++) if (m >> i & 1) o.push(CATS[i]); return o; }
  var text = R.map(function (r) {
    return [r[0], r[2], r[3], r[5], r[6], r[10], catNames(r[11]).join(" "), r[12]].join(" ").toLowerCase();
  });
  var TILES = {
    all: function () { return true; }, sp: function (r) { return !!r[6]; }, pfam: function (r) { return !!r[10]; },
    none: function (r) { return !r[6] && !r[10]; }, immune: function (r) { return r[11] > 0; },
    allo: function (r) { return !!r[12]; }, drawn: function (r) { return !!r[14]; }
  };
  var TILE_WORDS = { all: "", sp: "with a Swiss-Prot hit", pfam: "with a Pfam domain", none: "with neither a Swiss-Prot hit nor a Pfam domain",
    immune: "that are immune-gene candidates", allo: "with a Fester or BHF hit", drawn: "with every program's hits drawn" };
  var st = { q: "", tile: "all", cat: -1, chk: -1, key: 0, dir: 1, shown: 100 };
  var tbody = document.querySelector("#t tbody"), count = document.getElementById("count"), more = document.getElementById("more");
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function fmt(n) { return n.toLocaleString("en-US"); }
  var list = [];
  function filter() {
    var words = st.q.toLowerCase().split(/\s+/).filter(Boolean), tf = TILES[st.tile];
    list = [];
    for (var i = 0; i < R.length; i++) {
      var r = R[i];
      if (!tf(r)) continue;
      if (st.cat >= 0 && !(r[11] >> st.cat & 1)) continue;
      if (st.chk >= 0 && r[4] !== st.chk) continue;
      var t = text[i], ok = true;
      for (var w = 0; w < words.length; w++) if (t.indexOf(words[w]) < 0) { ok = false; break; }
      if (ok) list.push(i);
    }
    var k = st.key, d = st.dir, s = k === 0 || k === 2 || k === 3 || k === 5 || k === 10 || k === 12;
    list.sort(function (a, b) {
      var x = R[a][k], y = R[b][k];
      var ex = x === null || x === "", ey = y === null || y === "";
      if (ex !== ey) return ex ? 1 : -1;  // empty cells last, whichever the direction
      if (s) return d * String(x).localeCompare(String(y));
      return d * (x - y) || R[a][0].localeCompare(R[b][0]);
    });
    st.shown = 100;
    draw();
  }
  function draw() {
    var h = [];
    for (var j = 0; j < Math.min(st.shown, list.length); j++) {
      var r = R[list[j]];
      h.push("<tr><td><a class=\"mono\" href=\"protein.html#" + encodeURIComponent(r[0]) + "\">" + esc(r[0]) + "</a></td>" +
        "<td class=\"r\">" + fmt(r[1]) + "</td><td>" + esc(r[2]) + "</td><td>" + esc(r[3]) + "</td>" +
        "<td>" + esc(r[5]) + "</td><td class=\"r\">" + (r[7] === null ? "" : r[7].toFixed(1)) + "</td>" +
        "<td class=\"r\">" + (r[8] === null ? "" : Math.round(r[8])) + "</td><td>" + esc(r[10]) + "</td>" +
        "<td>" + esc(catNames(r[11]).join(", ")) + "</td><td>" + esc(r[12]) + "</td><td class=\"r\">" + fmt(r[13]) + "</td>" +
        "<td>" + (r[14] ? "<a href=\"" + esc(r[14]) + "\">drawing</a>" : "") + "</td></tr>");
    }
    tbody.innerHTML = h.join("");
    var what = [TILE_WORDS[st.tile]];
    if (st.cat >= 0) what.push("in the immune category " + CATS[st.cat]);
    if (st.chk >= 0) what.push("where " + F.check[st.chk]);
    if (st.q) what.push("matching “" + st.q + "”");
    count.textContent = fmt(list.length) + " of " + fmt(R.length) + " proteins" +
      (what.filter(Boolean).length ? " " + what.filter(Boolean).join(", ") : "") +
      "; showing " + fmt(Math.min(st.shown, list.length)) + ".";
    more.hidden = st.shown >= list.length;
  }
  function mark() {
    document.querySelectorAll(".num").forEach(function (b) { b.classList.toggle("on", b.dataset.t === st.tile && st.tile !== "all"); });
    document.querySelectorAll(".bar[data-cat]").forEach(function (b) { b.classList.toggle("on", +b.dataset.cat === st.cat); });
    document.querySelectorAll(".bar[data-chk]").forEach(function (b) { b.classList.toggle("on", +b.dataset.chk === st.chk); });
  }
  function go() { mark(); filter(); document.getElementById("t").scrollIntoView({ behavior: "smooth", block: "start" }); }
  document.querySelectorAll(".num").forEach(function (b) {
    b.addEventListener("click", function () { st.tile = b.dataset.t; st.cat = -1; st.chk = -1; go(); });
  });
  document.querySelectorAll(".bar[data-cat]").forEach(function (b) {
    b.addEventListener("click", function () { st.cat = st.cat === +b.dataset.cat ? -1 : +b.dataset.cat; go(); });
  });
  document.querySelectorAll(".bar[data-chk]").forEach(function (b) {
    b.addEventListener("click", function () { st.chk = st.chk === +b.dataset.chk ? -1 : +b.dataset.chk; go(); });
  });
  var timer;
  document.getElementById("q").addEventListener("input", function (e) {
    clearTimeout(timer);
    timer = setTimeout(function () { st.q = e.target.value.trim(); filter(); }, 150);
  });
  document.querySelectorAll("#t th[data-k]").forEach(function (th) {
    th.addEventListener("click", function () {
      var k = +th.dataset.k;
      st.dir = st.key === k ? -st.dir : (th.dataset.t === "s" ? 1 : -1);
      st.key = k;
      document.querySelectorAll("#t th").forEach(function (x) { x.removeAttribute("aria-sort"); });
      th.setAttribute("aria-sort", st.dir > 0 ? "ascending" : "descending");
      filter();
    });
  });
  more.addEventListener("click", function () { st.shown += 100; draw(); });
  filter();
})();
