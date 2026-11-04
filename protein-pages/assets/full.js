// The first page (build_full.py): one grid row per protein, 100 at a time, under a sticky header.
// Filters (the chips, the two lists and the overview bars), search, sorting by column, a TSV of
// every protein the filters keep, Help, and keys (j/k rows, Enter opens, / search, Esc).
(function () {
  var F = window.FULL, R = F.rows, CATS = F.cats, PFAM = F.pfam || [];
  // row: 0 id, 1 length, 2 existing, 3 proposed, 4 name check, 5 Swiss-Prot name, 6 accession,
  // 7 identity, 8 coverage, 9 E-value, 10 Pfam names, 11 category bits, 12 empty, 13 cluster size,
  // 14 drawn page, 15 InterProScan's Pfam matches [start, end, index into PFAM], 16 MMseqs2's best
  // hit [start, end, identity %, accession] or null, 17 the other programs with a hit, as bits of F.other, 19 kmerseek's best hit
  // [E-value, setting, accession, name] when no other program found anything and E is at most 0.001, else 0,
  // 18 1 when the two names conflict (not only their spelling; build_full.names_conflict)
  var grid = document.getElementById("grid"), AX = +grid.dataset.ax, SEQ = grid.dataset.seq, DOM = grid.dataset.dom;
  function catNames(m) { var o = []; for (var i = 0; i < CATS.length; i++) if (m >> i & 1) o.push(CATS[i]); return o; }
  var text = R.map(function (r) {
    return [r[0], r[2], r[3], r[5], r[6], r[10], catNames(r[11]).join(" ")].join(" ").toLowerCase();
  });
  var OTHER = F.other || [], CUTS = [1e-3, 1e-10, 1e-50, 1e-100];
  // a Swiss-Prot match counts when DIAMOND's E-value is at or below the cutoff on the ladder
  function isMatch(r) { var e = ev(r[9]); return !!r[6] && e !== null && e <= st.e; }
  var TILES = {
    all: function () { return true; }, sp: isMatch, pfam: function (r) { return !!r[10]; },
    none: function (r) { return !isMatch(r) && !r[10]; }, immune: function (r) { return r[11] > 0; },
    drawn: function (r) { return !!r[14]; }, names: function (r) { return r[18] === 1; },
    km3: function (r) { return !!r[19]; }, km5: function (r) { return r[19] && r[19][0] <= 1e-5; }
  };
  function others(b) { var o = []; for (var i = 0; i < OTHER.length; i++) if (b >> i & 1) o.push(OTHER[i][0]); return o; }
  function nbits(b) { var n = 0; while (b) { n += b & 1; b >>= 1; } return n; }
  var TILE_WORDS = { all: "", sp: "with a Swiss-Prot match", pfam: "with a Pfam domain", none: "with neither a Swiss-Prot match nor a Pfam domain",
    immune: "that are immune-gene candidates", drawn: "with every program's hits drawn", names: "whose FASTA name and proposed name conflict",
    km3: "that only kmerseek found, E-value 0.001 or below", km5: "that only kmerseek found, E-value 0.00001 or below" };
  var st = { q: "", tile: "all", cat: -1, chk: -1, key: -1, dir: 1, shown: 100, cur: -1, e: 1e-3 };
  var count = document.getElementById("count"), more = document.getElementById("more"), lbl = document.querySelector(".ix .ruler .lbl");
  var catSel = document.getElementById("cat"), chkSel = document.getElementById("chk"), q = document.getElementById("q");
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function fmt(n) { return n.toLocaleString("en-US"); }
  function href(r) { return r[14] ? r[14] : "protein.html#" + encodeURIComponent(r[0]); }
  function ev(x) { var v = parseFloat(x); return isNaN(v) ? null : v; }
  var list = [];
  // the value a column sorts by; null sorts last whichever the direction
  function val(r, k) {
    if (k === 9) return r[9] === "" ? null : ev(r[9]);
    if (k === 17) return nbits(r[17]);
    var v = r[k];
    return v === null || v === "" ? null : v;
  }
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
    if (st.key >= 0) {
      var k = st.key, d = st.dir;
      list.sort(function (a, b) {
        var x = val(R[a], k), y = val(R[b], k);
        if ((x === null) !== (y === null)) return x === null ? 1 : -1;
        if (x === null) return R[a][0].localeCompare(R[b][0]);
        var c = typeof x === "string" ? x.localeCompare(y) : x - y;
        return d * c || R[a][0].localeCompare(R[b][0]);
      });
    }
    else if (st.tile === "km3" || st.tile === "km5") list.sort(function (a, b) { return R[a][19][0] - R[b][19][0]; });  // best first
    st.shown = 100; st.cur = -1;
    draw();
  }
  function bin(p) { return p >= 60 ? 1 : p >= 40 ? 0.75 : p >= 20 ? 0.5 : 0.3; }
  var TW = +grid.dataset.tw;  // the protein scale takes this % of the track (build_full.py TRACK_PCT); the rest holds "› 16,045 aa"
  function px(res) { return (TW * (Math.min(res, AX) - 1) / (AX - 1)).toFixed(2) + "%"; }
  function track(r) {
    var L = r[1], g = ['<line x1="0" x2="' + px(L) + '" y1="7" y2="7" stroke="#c8c8c8" stroke-width="2"/>'];
    (r[15] || []).forEach(function (d) {
      if (d[0] > AX) return;
      var p = PFAM[d[2]] || ["", ""], w = (TW * (Math.min(d[1], AX) - d[0]) / (AX - 1));
      g.push('<rect x="' + px(d[0]) + '" y="2" width="' + Math.max(w, 0.25).toFixed(2) + '%" height="10" rx="1.5" fill="#fff" stroke="' + DOM +
        '" stroke-width="1.2"><title>Pfam ' + esc(p[0]) + " " + esc(p[1]) + ", residues " + d[0] + "–" + d[1] + " (InterProScan)</title></rect>");
    });
    var m = r[16];
    if (m && m[0] <= AX) {
      var w2 = (TW * (Math.min(m[1], AX) - m[0]) / (AX - 1));
      g.push('<rect x="' + px(m[0]) + '" y="14" width="' + Math.max(w2, 0.25).toFixed(2) + '%" height="4" fill="' + SEQ + '" fill-opacity="' + bin(m[2]) +
        '"><title>MMseqs2\'s best Swiss-Prot hit ' + esc(m[3]) + (m[3] !== r[6] ? (r[6] ? " (DIAMOND\'s best is " + esc(r[6]) + ")" : " (DIAMOND found no match)") : " (DIAMOND\'s best too)") +
        ", residues " + m[0] + "–" + m[1] + ", " + m[2] + "% identical</title></rect>");
    }
    if (L > AX) g.push('<text x="' + TW + '%" y="11" dx="3" class="ax">&#8250; ' + fmt(L) + " aa</text>");
    return '<svg role="img" aria-label="' + esc(r[0]) + ' along its length">' + g.join("") + "</svg>";
  }
  function rowHtml(r, j) {
    var same = r[2] && r[3] && r[2].toLowerCase() === r[3].toLowerCase();
    var names = same ? esc(r[3]) + " (both)" : r[2] || r[3] ? esc(r[2] || "—") + " &#8594; " + esc(r[3] || "—") : '<span class="dim">no name</span>';
    var ips = (r[15] || []).length, mmo = r[16] && r[16][3] !== r[6];
    return '<div class="row g ixg' + (j === st.cur ? " cur" : "") + '" role="listitem" data-j="' + j + '" tabindex="-1">' +
      '<div><a class="id mono" href="' + esc(href(r)) + '">' + esc(r[0]) + "</a>" +
      '<span class="s">' + names + (r[14] ? '<span class="dr">· every program drawn</span>' : "") + "</span></div>" +
      '<div class="num">' + fmt(r[1]) + "</div>" +
      '<div class="txt" title="' + esc(r[5]) + '">' + (r[6] ? (isMatch(r) ? esc(r[5]) : '<span class="dim">above the cutoff: ' + esc(r[5]) + "</span>") : '<span class="dim">no match</span>') +
        '<span class="s">' + (r[6] ? '<span class="mono">' + esc(r[6]) + "</span>" : "") +
        (mmo ? (r[6] ? " · " : "") + "MMseqs2's best: " + '<span class="mono">' + esc(r[16][3]) + "</span>" : "") + "</span></div>" +
      '<div class="num">' + (r[7] === null ? "" : r[7].toFixed(1) + "%") + (r[8] === null ? "" : '<span class="s">covers ' + Math.round(r[8]) + "%</span>") + "</div>" +
      '<div class="num"' + (r[9] && ev(r[9]) === 0 ? ' title="DIAMOND writes 0 when the E-value is smaller than the smallest number it stores">0' : ">" + esc(r[9] || "")) + "</div>" +
      '<div class="txt" title="' + esc(r[10]) + '">' + (r[10] ? esc(r[10]) : '<span class="dim">none</span>') +
        '<span class="s">InterProScan: ' + (ips ? ips + " Pfam match" + (ips > 1 ? "es" : "") : "none") + "</span></div>" +
      '<div class="txt" title="' + esc(catNames(r[11]).join(", ")) + '">' + esc(catNames(r[11]).join(", ")) + "</div>" +
      (r[19] && (st.tile === "km3" || st.tile === "km5") ? '<div class="txt" title="' + esc(r[19][3] + " (" + r[19][2] + "), " + r[19][1]) + '">kmerseek ' + r[19][0].toExponential(1) + " · " + esc(r[19][3]) + "</div>" :
      '<div class="txt" title="' + esc(others(r[17]).join(", ")) + '">' + (r[17] ? esc(others(r[17]).join(" · ")) : '<span class="dim">none</span>') + "</div>") +
      '<div class="num">' + fmt(r[13]) + (r[13] === 1 ? '<span class="s">alone</span>' : "") + "</div>" +
      '<div class="track">' + track(r) + "</div></div>";
  }
  function draw() {
    var h = [], n = Math.min(st.shown, list.length);
    for (var j = 0; j < n; j++) h.push(rowHtml(R[list[j]], j));
    grid.innerHTML = h.join("");
    var what = [TILE_WORDS[st.tile]];
    if (st.cat >= 0) what.push("in the immune category " + CATS[st.cat]);
    if (st.chk >= 0) what.push("where " + F.check[st.chk]);
    if (st.q) what.push("matching “" + st.q + "”");
    what = what.filter(Boolean);
    count.innerHTML = "<b>" + fmt(list.length) + "</b> of " + fmt(R.length) + " proteins" + (what.length ? " " + esc(what.join(", ")) : "");
    lbl.textContent = list.length ? "rows 1–" + fmt(n) + " of " + fmt(list.length) + " shown" + (n < list.length ? "; more at the bottom" : "") : "no protein matches";
    more.hidden = n >= list.length;
    more.textContent = "Show " + fmt(Math.min(100, list.length - n)) + " more (of " + fmt(list.length - n) + " not shown yet)";
    mark();
  }
  // the ladder: what counts as a Swiss-Prot match at each DIAMOND cutoff; the chips' counts follow it
  var ladder = document.getElementById("eladder");
  function fmtE(e) { return e >= 1e-3 ? String(e) : "1e" + Math.round(Math.log10(e)); }
  function drawLadder() {
    ladder.innerHTML = CUTS.map(function (c) {
      var n = 0; for (var i = 0; i < R.length; i++) { var e = ev(R[i][9]); if (R[i][6] && e !== null && e <= c) n++; }
      return '<button class="chip' + (c === st.e ? " on" : "") + '" data-e="' + c + '" title="' + fmt(n) + " proteins have a DIAMOND match at E-value " + fmtE(c) +
        ' or below">≤ ' + fmtE(c) + ': <b class="mono">' + fmt(n) + "</b></button>";
    }).join("");
    ladder.querySelectorAll(".chip").forEach(function (b) {
      b.addEventListener("click", function () { st.e = +b.dataset.e; drawLadder(); filter(); });
    });
    var c = { sp: 0, none: 0 };
    for (var i = 0; i < R.length; i++) { if (TILES.sp(R[i])) c.sp++; if (TILES.none(R[i])) c.none++; }
    ["sp", "none"].forEach(function (k) { document.querySelector('#tiles .chip[data-t="' + k + '"] b').textContent = fmt(c[k]); });
  }
  function mark() {
    document.querySelectorAll("#tiles .chip").forEach(function (b) { b.classList.toggle("on", b.dataset.t === st.tile); });
    document.querySelectorAll(".bar[data-cat]").forEach(function (b) { b.classList.toggle("on", +b.dataset.cat === st.cat); });
    document.querySelectorAll(".bar[data-chk]").forEach(function (b) { b.classList.toggle("on", +b.dataset.chk === st.chk); });
    catSel.value = String(st.cat); chkSel.value = String(st.chk);
    document.querySelectorAll(".colhead button[data-k]").forEach(function (b) {
      var on = +b.dataset.k === st.key;
      b.classList.toggle("on", on);
      b.querySelector("i").textContent = on ? (st.dir > 0 ? "▲" : "▼") : "↕";
    });
  }
  function top() { var g = grid.getBoundingClientRect().top + window.scrollY - hh(); if (window.scrollY > g) window.scrollTo(0, g); }
  document.querySelectorAll("#tiles .chip").forEach(function (b) {
    b.addEventListener("click", function () { st.tile = b.dataset.t; filter(); top(); });
  });
  catSel.addEventListener("change", function () { st.cat = +catSel.value; filter(); top(); });
  chkSel.addEventListener("change", function () { st.chk = +chkSel.value; filter(); top(); });
  document.querySelectorAll(".bar[data-cat]").forEach(function (b) {
    b.addEventListener("click", function () { st.cat = st.cat === +b.dataset.cat ? -1 : +b.dataset.cat; filter(); grid.scrollIntoView(); });
  });
  document.querySelectorAll(".bar[data-chk]").forEach(function (b) {
    b.addEventListener("click", function () { st.chk = st.chk === +b.dataset.chk ? -1 : +b.dataset.chk; filter(); grid.scrollIntoView(); });
  });
  var timer;
  q.addEventListener("input", function () {
    clearTimeout(timer);
    timer = setTimeout(function () { st.q = q.value.trim(); filter(); }, 150);
  });
  // first click on a column = best first: numbers high to low (E-value low to high), words A to Z
  document.querySelectorAll(".colhead button[data-k]").forEach(function (b) {
    b.addEventListener("click", function () {
      var k = +b.dataset.k, first = b.dataset.t === "s" || b.dataset.t === "e" ? 1 : -1;
      st.dir = st.key === k ? -st.dir : first;
      st.key = k;
      filter();
    });
  });
  more.addEventListener("click", function () { st.shown += 100; draw(); });
  document.getElementById("clear").addEventListener("click", function () {
    st = { q: "", tile: "all", cat: -1, chk: -1, key: -1, dir: 1, shown: 100, cur: -1, e: 1e-3 }; q.value = ""; drawLadder(); filter();
  });
  document.getElementById("dl").addEventListener("click", function () {
    var head = ["protein", "length_aa", "existing_name", "proposed_name", "swissprot_name", "swissprot_accession", "identity_pct",
      "coverage_pct", "evalue", "pfam_domains", "immune_categories", "other_programs_with_a_hit", "cluster_size", "only_kmerseek_evalue", "only_kmerseek_setting", "only_kmerseek_hit"];
    var lines = [head.join("\t")].concat(list.map(function (i) {
      var r = R[i];
      return [r[0], r[1], r[2], r[3], r[5], r[6], r[7] === null ? "" : r[7], r[8] === null ? "" : r[8], r[9], r[10],
        catNames(r[11]).join(", "), others(r[17]).join(", "), r[13]].concat(r[19] ? [r[19][0], r[19][1], r[19][2] + " " + r[19][3]] : ["", "", ""]).join("\t");
    }));
    var a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }));
    a.download = "botryllus_proteins.tsv"; a.click();
    setTimeout(function () { URL.revokeObjectURL(a.href); }, 1000);
  });
  var help = document.getElementById("help"), hb = document.getElementById("help-btn");
  function showHelp(on) { help.classList.toggle("open", on); hb.setAttribute("aria-expanded", on ? "true" : "false"); }
  hb.addEventListener("click", function () { showHelp(!help.classList.contains("open")); });
  grid.addEventListener("click", function (e) {
    if (e.target.closest("a")) return;
    var row = e.target.closest(".row");
    if (row) location.href = href(R[list[+row.dataset.j]]);
  });
  function hh() { return document.querySelector("header.top").getBoundingClientRect().height; }
  function setHH() { document.documentElement.style.setProperty("--hh", hh() + "px"); }
  window.addEventListener("resize", setHH);
  function go(j) {
    var n = grid.children.length;
    if (!n) return;
    st.cur = Math.max(0, Math.min(n - 1, j));
    grid.querySelectorAll(".row.cur").forEach(function (x) { x.classList.remove("cur"); });
    var row = grid.children[st.cur];
    row.classList.add("cur"); row.focus({ preventScroll: true }); row.scrollIntoView({ block: "nearest" });
  }
  document.addEventListener("keydown", function (e) {
    var t = e.target.tagName;
    if (e.key === "Escape") {
      if (help.classList.contains("open")) showHelp(false);
      else if (st.q) { q.value = ""; st.q = ""; filter(); }
      else if (st.cur >= 0) { grid.querySelectorAll(".row.cur").forEach(function (x) { x.classList.remove("cur"); }); st.cur = -1; }
      return;
    }
    if (t === "INPUT" || t === "SELECT" || t === "TEXTAREA" || e.metaKey || e.ctrlKey || e.altKey) return;
    if (e.key === "/") { e.preventDefault(); q.focus(); }
    else if (e.key === "j") go(st.cur + 1);
    else if (e.key === "k") go(st.cur - 1);
    else if ((e.key === "Enter" || e.key === " ") && st.cur >= 0) { e.preventDefault(); location.href = href(R[list[st.cur]]); }
  });
  // A phone: the grid scrolls sideways inside main; the column names and the ruler move with it
  var main = document.querySelector("main"), heads = document.querySelectorAll(".ix .colhead, .ix .ruler");
  main.addEventListener("scroll", function () {
    heads.forEach(function (h) { h.style.transform = main.scrollLeft ? "translateX(" + -main.scrollLeft + "px)" : ""; });
  }, { passive: true });
  drawLadder();
  filter();
  setHH();
})();
