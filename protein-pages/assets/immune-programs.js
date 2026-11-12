// The index page's "Possible immune genes, by program" panel (build_full.py immune_programs_panel).
// The rows are drawn by build_full.py; this adds, once full/immune_programs.js has loaded:
//  - the ladder: count only proteins that at least k programs call immune, and redraw every row;
//  - pointing at a program: every other row shows, in vermillion, how many of its proteins the
//    pointed-at program also calls immune;
//  - clicking a program: its proteins open under the row, most agreed first, 50 at a time, each
//    with one square per program (filled = that program calls it immune) and this program's evidence;
//  - a search over the open list, and a TSV of what it keeps.
(function () {
  var box = document.getElementById("ip");
  if (!box) return;
  var V = (document.currentScript.src.split("?")[1] || ""), D = null, loading = null;
  var st = { k: 1, open: -1, hover: -1, shown: 50, q: [] };
  var rowsEl = box.querySelectorAll(".ip-row"), f = function (v) { return v.toLocaleString("en-US"); };
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function bits(m) { var n = 0; while (m) { n += m & 1; m >>= 1; } return n; }

  function load() {
    if (D) return Promise.resolve(D);
    return loading || (loading = loadScript("full/immune_programs.js" + (V ? "?" + V : "")).then(function () {
      return gunzipB64(window.IMMUNE_PROG_Z);
    }).then(function (d) {
      window.IMMUNE_PROG_Z = null;
      d.rows.forEach(function (r) { r.n = bits(r[2]); r.text = (r[0] + " " + r[4] + " " + r[3].join(" ")).toLowerCase(); });
      D = d; return d;
    }));
  }

  function kept() { return D.rows.filter(function (r) { return r.n >= st.k; }); }

  function draw() {
    var rs = kept(), P = D.programs.length, per = [], top = 1, ov = [];
    for (var i = 0; i < P; i++) { per.push([0, 0]); ov.push(0); }
    rs.forEach(function (r) {
      for (var i = 0; i < P; i++) if (r[2] >> i & 1) {
        per[i][r[1] ? 0 : 1]++;
        if (st.hover >= 0 && r[2] >> st.hover & 1) ov[i]++;
      }
    });
    per.forEach(function (p) { top = Math.max(top, p[0] + p[1]); });
    rowsEl.forEach(function (el, i) {
      var a = per[i][0], b = per[i][1], t = a + b;
      el.querySelector("i.new").style.width = (100 * a / top) + "%";
      el.querySelector("i.old").style.width = (100 * b / top) + "%";
      var o = el.querySelector("b.ov") || el.querySelector(".ip-bar").appendChild(document.createElement("b"));
      o.className = "ov";
      var show = st.hover >= 0 && i !== st.hover;
      o.style.width = show ? (100 * ov[i] / top) + "%" : "0";
      el.querySelector(".ip-n").innerHTML = show
        ? f(t) + ' <span class="ovn" title="of these ' + f(t) + ", " + f(ov[i]) + " are also called immune by " + esc(D.programs[st.hover][1]) + '">' + f(ov[i]) + " shared</span>"
        : f(t);
      el.querySelector(".ip-split").textContent = f(a) + " new · " + f(b) + " old";
      el.classList.toggle("dim", show && ov[i] === 0);
      el.classList.toggle("pt", i === st.hover);
    });
    var note = box.querySelector(".ip-hovernote");
    note.textContent = st.hover < 0 ? "" : "Vermillion: of each program's proteins, those " + D.programs[st.hover][1] +
      " also calls immune (" + f(per[st.hover][0] + per[st.hover][1]) + " proteins" + (st.k > 1 ? ", " + st.k + " or more programs" : "") + ").";
    drawList();
  }

  function squares(m) {
    var s = "";
    for (var i = 0; i < D.programs.length; i++)
      s += '<i class="' + (m >> i & 1 ? "on" : "") + '" title="' + esc(D.programs[i][1]) + (m >> i & 1 ? ": calls it immune" : ": does not") + '"></i>';
    return s;
  }

  function listRows() {
    var p = st.open;
    return kept().filter(function (r) {
      if (!(r[2] >> p & 1)) return false;
      for (var j = 0; j < st.q.length; j++) if (r.text.indexOf(st.q[j]) < 0) return false;
      return true;
    });
  }

  function evidence(r, p) {  // the open program's evidence: r[3] holds one entry per set bit, in program order
    var j = 0;
    for (var i = 0; i < p; i++) if (r[2] >> i & 1) j++;
    return r[3][j];
  }

  function drawList() {
    var old = box.querySelector(".ip-list");
    if (old) old.remove();
    rowsEl.forEach(function (el) { el.classList.remove("open"); el.setAttribute("aria-expanded", "false"); });
    document.getElementById("ip-dl").hidden = st.open < 0;
    if (st.open < 0) return;
    var p = st.open, el = rowsEl[p], rs = listRows(), n = Math.min(st.shown, rs.length);
    el.classList.add("open"); el.setAttribute("aria-expanded", "true");
    var html = '<div class="ip-list" role="list"><div class="ip-lhead"><span>protein · gene name</span><span title="one square per program, in the order of the rows above; filled = it calls the protein immune">programs that agree, in row order</span><span>' +
      esc(D.programs[p][1]) + ": its evidence</span></div>";
    for (var i = 0; i < n; i++) {
      var r = rs[i];
      var id = r[1] ? '<a class="mono" href="protein.html#' + encodeURIComponent(r[0]) + '">' + esc(r[0]) + "</a>"
        : '<span class="mono" title="old proteome; no identical protein in this one, so no page here">' + esc(r[0]) + "</span>";
      id = '<span class="pid">' + id + (r[1] ? (r[4] ? ' <small>' + esc(r[4]) + "</small>" : "") : " <small>old proteome</small>") + "</span>";
      html += '<div class="ip-li" role="listitem">' + id + '<span class="sq">' + squares(r[2]) + '<b class="mono">' + r.n + "</b></span><span class=\"ev\">" + esc(evidence(r, p)) + "</span></div>";
    }
    html += '<p class="ip-more">' + f(n) + " of " + f(rs.length) + " proteins" + (st.q.length ? " that match the search" : "") +
      (n < rs.length ? ' <button class="btn" id="ip-more">show 50 more</button>' : "") + "</p></div>";
    el.insertAdjacentHTML("afterend", html);
    var more = document.getElementById("ip-more");
    if (more) more.onclick = function () { st.shown += 50; drawList(); };
  }

  function busy(fn) { return function (e) { var a = arguments, t = this; load().then(function () { fn.apply(t, a); }); }; }

  rowsEl.forEach(function (el, i) {
    el.addEventListener("mouseenter", busy(function () { st.hover = i; draw(); }));
    el.addEventListener("mouseleave", busy(function () { if (st.hover === i) { st.hover = -1; draw(); } }));
    var toggle = busy(function () { st.open = st.open === i ? -1 : i; st.shown = 50; draw(); });
    el.addEventListener("click", toggle);
    el.addEventListener("keydown", function (e) {
      if (e.key === "Enter" || e.key === " ") { e.preventDefault(); e.stopPropagation(); toggle(); }
      if (e.key === "Escape" && st.open >= 0) { e.stopPropagation(); st.open = -1; draw(); }
    });
    el.addEventListener("focus", busy(function () { st.hover = i; draw(); }));
    el.addEventListener("blur", busy(function () { if (st.hover === i) { st.hover = -1; draw(); } }));
  });
  box.querySelectorAll(".ip-ladder .chip").forEach(function (b) {
    b.addEventListener("click", busy(function () {
      st.k = +b.dataset.k; st.shown = 50;
      box.querySelectorAll(".ip-ladder .chip").forEach(function (c) { c.classList.toggle("on", c === b); });
      draw();
    }));
  });
  var q = document.getElementById("ip-q");
  q.addEventListener("input", busy(function () {
    st.q = q.value.toLowerCase().split(/\s+/).filter(Boolean); st.shown = 50;
    if (st.open < 0 && st.q.length) st.open = 0;
    drawList();
  }));
  q.addEventListener("keydown", function (e) { e.stopPropagation(); });
  document.getElementById("ip-dl").addEventListener("click", function () {
    var p = st.open, lines = ["protein_id\tproteome\tn_programs\t" + D.programs.map(function (x) { return x[0]; }).join("\t")];
    listRows().forEach(function (r) {
      var cols = [], j = 0;
      for (var i = 0; i < D.programs.length; i++) cols.push(r[2] >> i & 1 ? r[3][j++] : "");
      lines.push([r[0], r[1] ? "new" : "old", r.n].concat(cols).join("\t"));
    });
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }));
    a.download = "immune_" + D.programs[p][0] + (st.k > 1 ? "_at_least_" + st.k : "") + ".tsv";
    document.body.appendChild(a); a.click(); a.remove();
  });
  // start loading as soon as the fold opens (or now, if it is open)
  var fold = box.closest("details");
  if (!fold || fold.open) load().then(draw);
  else fold.addEventListener("toggle", function () { if (fold.open) load().then(draw); });
})();
