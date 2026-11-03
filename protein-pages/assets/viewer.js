// protein.html: one page for any protein. The ID after "#" picks the protein; its data are in
// full/NNN.js, where NNN is the same hash of the ID that build_full.py's shard_of computes.
// A protein with a drawn page (build.py) goes straight to it: build_full.py adds an element
// #biomni to that page, and this script fills it with the same tables.
(function () {
  var F = window.FULL, CATS = F.cats, N_SHARDS = 96;
  var box = document.getElementById("biomni"), BASE = box ? "../" : "";
  var byId = {};
  F.rows.forEach(function (r) { byId[r[0]] = r; });
  function shardOf(id) {
    var h = 0;
    for (var i = 0; i < id.length; i++) h = (Math.imul(h, 31) + id.charCodeAt(i)) >>> 0;
    return h % N_SHARDS;
  }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function fmt(n) { return n.toLocaleString("en-US"); }
  var main = box || document.getElementById("main"), loaded = {};
  function load(s, then) {
    if (loaded[s]) return then();
    var el = document.createElement("script");
    el.src = BASE + "full/" + ("00" + s).slice(-3) + ".js";
    el.onload = function () { loaded[s] = true; then(); };
    el.onerror = function () { main.innerHTML = "<p>Could not load " + esc(el.src) + ".</p>"; };
    document.head.appendChild(el);
  }
  function table(head, rows) {
    return "<div class=\"tablewrap\"><table class=\"small\"><thead><tr>" + head.map(function (h) {
      return "<th" + (h[1] ? " class=\"r\"" : "") + ">" + esc(h[0]) + "</th>";
    }).join("") + "</tr></thead><tbody>" + rows.map(function (r) {
      return "<tr>" + r.map(function (c, i) { return "<td" + (head[i][1] ? " class=\"r\"" : "") + ">" + c + "</td>"; }).join("") + "</tr>";
    }).join("") + "</tbody></table></div>";
  }
  function link(id) { return "<a class=\"mono\" href=\"" + BASE + "protein.html#" + encodeURIComponent(id) + "\">" + esc(id) + "</a>"; }
  function show() {
    var id = location.hash.slice(1), r;
    try { id = decodeURIComponent(id); } catch (e) { /* malformed escape: keep the raw text */ }
    r = byId[id];
    document.getElementById("pid").textContent = id || "No protein named";
    if (!r) { main.innerHTML = "<p>No protein " + esc(id) + " in this proteome. <a href=\"index.html\">Search all proteins</a>.</p>"; return; }
    if (r[14]) { location.replace(r[14]); return; }  // the drawn page is this protein's page
    document.title = id + (r[3] || r[2] ? " " + (r[3] || r[2]) : "") + " · Botryllus protein";
    load(shardOf(id), function () { render(id, r, window.FULL_SHARD[id]); window.scrollTo(0, 0); });
  }
  // Links to other sites open in a new tab. uni: a UniProt accession (P12345, P12345-2). pfam: a
  // Pfam accession, with or without its version (PF00179.32 links to PF00179). accHtml: s escaped,
  // with every Pfam (PF00069) and InterPro (IPR000719) accession in it made a link; accessions
  // are letters and digits, so escaping leaves them as they were.
  function ext(url, text) { return "<a class=\"mono\" href=\"" + esc(url) + "\" target=\"_blank\" rel=\"noopener\">" + esc(text) + "</a>"; }
  function uni(acc) { return ext("https://www.uniprot.org/uniprotkb/" + encodeURIComponent(acc) + "/entry", acc); }
  function pfam(acc) { return ext("https://www.ebi.ac.uk/interpro/entry/pfam/" + encodeURIComponent(acc.split(".")[0]) + "/", acc); }
  function ipr(acc) { return ext("https://www.ebi.ac.uk/interpro/entry/InterPro/" + encodeURIComponent(acc) + "/", acc); }
  function accHtml(s) {
    return esc(s).replace(/\b(PF\d{5}|IPR\d{6})(\.\d+)?\b/g, function (m, a) { return a.indexOf("PF") === 0 ? pfam(m) : ipr(m); });
  }
  // The top of every protein page, in the same order as a drawn page: name line, one sentence on
  // which programs report a hit, then one row per program with its best hit.
  function top(id, r, d) {
    var cats = []; for (var i = 0; i < CATS.length; i++) if (r[11] >> i & 1) cats.push(CATS[i]);
    var progs = [], rows = [];
    // how: one short grey line under the program's name saying how it was run
    function row(name, found, hit, score, where, how) {
      if (found) progs.push(name);
      rows.push(["<b>" + esc(name) + "</b>" + (how ? "<br><span class=\"how\">" + esc(how) + "</span>" : ""),
        found ? hit : "<i>no hit</i>", found ? score : "", found ? where : ""]);
    }
    var nowhere = "not in the table";
    row("DIAMOND against Swiss-Prot 2026_03", !!r[6], r[6] ? uni(r[6]) + " " + esc(r[5]) + (d.sp_gene ? " (" + esc(d.sp_gene) + ")" : "") + (d.sp_org ? ", <i>" + esc(d.sp_org) + "</i>" : "") : "",
      r[6] ? (r[7] === null ? "" : r[7].toFixed(1) + "% identity over " + Math.round(r[8]) + "% of your protein, ") + "E-value " + esc(d.sp_e) : "", nowhere,
      "run by Biomni; E-value 0.001 or below");
    var pf = {}; d.pfam.forEach(function (p) { pf[p[1]] = (pf[p[1]] || 0) + 1; });
    row("hmmsearch with Pfam", d.pfam.length > 0, Object.keys(pf).map(function (n) { return esc(n) + (pf[n] > 1 ? " ×" + pf[n] : ""); }).join(", "),
      "gathering threshold", nowhere, "run by Biomni; Pfam's own cutoffs");
    var j = d.jh[0];
    row("jackhmmer from immune seed proteins", d.jh.length > 0, j ? esc(j[0]) + " (seed " + esc(j[1]) + ")" + (d.jh.length > 1 ? " and " + (d.jh.length - 1) + " more" : "") : "",
      j ? "E-value " + esc(j[3]) + ", round " + j[2] : "", nowhere, "run by Biomni from 24 immune seed proteins");
    var f = d.fa;
    // FANTASIA always names its most similar protein; it counts as a hit below FANTASIA_MAX_DIST
    // (build.py; figures/fantasia_vs_diamond.png)
    var FANTASIA_MAX_DIST = 0.15, fName = "FANTASIA (most similar protein by a protein language model)",
      fHow = "run by Biomni; picks the most similar of 124,397 well-studied UniProt proteins, so it always names one";
    var fHit = f ? uni(f[0]) + " " + esc(f[1] || "") + (f[2] ? " (" + esc(f[2]) + ")" : "") + (f[3] ? ", <i>" + esc(f[3]) + "</i>" : "") : "";
    if (f && f[4] >= FANTASIA_MAX_DIST) rows.push(["<b>" + esc(fName) + "</b><br><span class=\"how\">" + esc(fHow) + "</span>",
      "<i>not close enough to count</i>; most similar: " + fHit,
      "distance " + f[4].toFixed(2) + " (it counts below " + FANTASIA_MAX_DIST.toFixed(2) + "; 0 is identical)", "whole protein, no position"]);
    else row(fName, !!f, fHit, f ? "distance " + f[4].toFixed(2) + ", close enough to count (below " + FANTASIA_MAX_DIST.toFixed(2) + ")" : "",
      f ? "whole protein, no position" : "", fHow);
    if (d.mm !== undefined && d.mm !== null) {
      var m = d.mm;
      row("MMseqs2 against Swiss-Prot 2026_03", !!m, m ? uni(m[0]) + " " + esc(m[6] || "") + (m[7] ? ", <i>" + esc(m[7]) + "</i>" : "") : "",
        m ? (100 * m[1]).toFixed(1) + "% identity, E-value " + esc(m[4]) : "", m ? "residues " + m[2] + "-" + m[3] : "",
        "E-value 0.001 or below");
    }
    if (d.fs !== undefined && d.fs !== null) {
      var g = d.fs;
      row("ProstT5 + Foldseek against AlphaFold models of Swiss-Prot (590,183)", !!g, g ? uni(g[0]) + " " + esc(g[6] || "") + (g[7] ? ", <i>" + esc(g[7]) + "</i>" : "") : "",
        g ? (100 * g[1]).toFixed(1) + "% identity, E-value " + esc(g[4]) : "", g ? "residues " + g[2] + "-" + g[3] : "",
        "run by Biomni; shape letters predicted from the sequence; E-value 0.001 or below");
    }
    if (d.fst === "not run") {
      rows.push(["<b>Foldseek on the predicted 3-D structure</b>", "<i>not run</i>: structures were predicted only for proteins of 200 aa or less", "", ""]);
    } else if (d.fst !== undefined && d.fst !== null) {
      var s = d.fst;
      row("Foldseek on the predicted 3-D structure, against AlphaFold models of Swiss-Prot (590,183)", !!s,
        s ? uni(s[0]) + " " + esc(s[6] || "") + (s[7] ? ", <i>" + esc(s[7]) + "</i>" : "") : "",
        s ? (100 * s[1]).toFixed(1) + "% identity, E-value " + esc(s[4]) : "", s ? "residues " + s[2] + "-" + s[3] : "",
        "run by Biomni on structures it predicted for proteins of 200 aa or less; E-value 0.001 or below");
    }
    if (d.ips !== undefined && d.ips !== null) {
      var IPS_FEAT = { Coils: 1, MobiDBLite: 1, Phobius: 1, TMHMM: 1, AntiFam: 1 };
      if (d.ips === "not run") rows.push(["<b>InterProScan</b>", "<i>not run (over 5,000 aa)</i>", "", ""]);
      else {
        var dom = d.ips.filter(function (x) { return !IPS_FEAT[x[0]]; }), feat = d.ips.filter(function (x) { return IPS_FEAT[x[0]]; });
        var ent = {}; dom.forEach(function (x) { var k = x[7] || x[2] || x[1]; if (!ent[k]) ent[k] = [x[3], x[4]]; });
        var dbs = {}; d.ips.forEach(function (x) { dbs[x[0]] = 1; });
        var ek = Object.keys(ent);
        row("InterProScan", d.ips.length > 0, ek.length ? ek.slice(0, 4).map(accHtml).join(", ") + (ek.length > 4 ? " and " + (ek.length - 4) + " more" : "") :
          "only " + Object.keys(dbs).map(esc).join(", ") + " (disorder, coils or membrane)",
          d.ips.length + " matches from " + Object.keys(dbs).length + " member databases",
          ek.length ? ek.slice(0, 4).map(function (k) { return ent[k][0] + "-" + ent[k][1]; }).join(", ") :
          feat.map(function (x) { return x[3] + "-" + x[4]; }).filter(function (v, i, a) { return a.indexOf(v) === i; }).join(", "),
          "version 5.59, all member databases");
      }
    }
    if (d.km !== undefined && d.km !== null) {
      var k = d.km;
      row("kmerseek against Swiss-Prot 2026_03", !!k, k ? uni(k[2]) + " " + esc(k[3]) + (k[4] ? ", <i>" + esc(k[4]) + "</i>" : "") : "",
        k ? "E-value " + k[1].toPrecision(3) + ", best of " + k[5] + " settings (" + esc(k[0]) + ")" : "", k ? (box ? "drawn above" : "not in this table") : "",
        "best of its alphabet and k-mer size settings; E-value 10 or below");
    }
    if (d.uf === "not run") rows.push(["<b>UniFIRE (UniProt's annotation rules)</b>", "<i>not run (over 5,000 aa)</i>", "", ""]);
    else if (d.uf) {
      var u = d.uf, pick = function (t) { return u.filter(function (x) { return x[2] === t; }).map(function (x) { return x[3]; }); };
      var nm = pick("protein.recommendedName.fullName")[0] || pick("comment.similarity")[0] ||
        pick("keyword").filter(function (k, i, a) { return a.indexOf(k) === i; }).slice(0, 5).join(", ");
      var src = {}; u.forEach(function (x) { src[x[0]] = src[x[0]] || {}; src[x[0]][x[1]] = 1; });
      var feats = u.filter(function (x) { return x[4] !== null; });
      row("UniFIRE (UniProt's annotation rules)", u.length > 0, esc(nm),
        u.length + " predictions from " + Object.keys(src).map(function (k) { var n = Object.keys(src[k]).length; return n + " " + k + " rule" + (n > 1 ? "s" : ""); }).join(", "),
        feats.length ? feats.length + " sites or regions; see below" : "whole protein",
        "UniProt's rules (ARBA, UniRule, PIRSR) applied to InterProScan's matches");
    }
    var N = rows.length;
    var name = (d.card && d.card.card_name) || r[3] || r[2] || "";
    var line = d.card ? d.card.one_line : (cats.length ? "immune-gene candidate: " + cats.join(", ") : "");
    if (!box) document.getElementById("sub").textContent = (name ? name + " · " : "") + (line ? line + " · " : "") + fmt(r[1]) + " amino acids";
    var n = progs.length, sentence = n === 0 ? "None of the " + N + " programs reports a hit." :
      n + " of the " + N + " programs report" + (n === 1 ? "s" : "") + " a hit: " + progs.join(", ") + ".";
    var h = box ? ["<h2>What the whole-proteome annotation found</h2><p class=\"lead\">" + esc(sentence) + "</p>",
      "<p class=\"note\">Each program's best hit for this protein, from runs on the whole proteome. " +
      "Weaker hits are left out. The slider above does not change this table.</p>"] :
      ["<p class=\"lead\">" + esc(sentence) + "</p>", "<h2>What each program found</h2>",
      "<p class=\"note\">The best hit of each program. Nothing is drawn along the protein on this page; the last column gives the residues where a program reports them" +
      (d.uf ? "; UniFIRE gives positions only for the sites it predicts, listed in its table below." : ".") + "</p>"];
    h.push(table([["Program"], ["Best hit"], ["Score"], ["Where on your protein"]], rows));
    return box ? h : h.concat(card(r, d));
  }
  function card(r, d) {
    var h = [];
    if (d.card) {
      var c = d.card, v = c.name_verdict;
      h.push("<h2>Gene card</h2><p class=\"note\">Written by a language model (" + esc(c.model) + ") from the evidence on this page only; " +
        "check it against the tables below.</p><p>" + accHtml(c.summary) + "</p><p><b>Confidence: " + esc(c.confidence) + ".</b> " + accHtml(c.confidence_reason) + "</p>");
      if (v && v.choice !== "not_applicable") h.push("<p><b>Header name " + esc(r[2]) + " or proposed name " + esc(r[3]) + "?</b> " +
        esc((F.verdicts[v.choice] || v.choice).replace(/^./, function (x) { return x.toUpperCase(); })) + "; suggested name <b>" + esc(v.recommended_name) + "</b>. " + accHtml(v.reason) + "</p>");
    }
    return h;
  }
  function render(id, r, d) {
    var h = top(id, r, d);
    if (box) {  // on a drawn page each program is a row of the drawing; under it, only what the rows do not show
      box.innerHTML = "";
      main = document.getElementById("biomni-more");
      h = card(r, d).concat(["<h2>More about this protein</h2>"]);
    }
    // Names
    h.push("<h2>Name</h2>" + table([["Existing name (FASTA header)"], ["Proposed name"], ["Comparison"]],
      [[esc(r[2]) || "none", esc(r[3]) || "none", esc(F.check[r[4]])]]));
    // Best Swiss-Prot match (on a drawn page its row in the drawing has the hit; only its GO terms here)
    h.push(box ? (r[6] ? "<h2>GO terms and keywords of DIAMOND's best hit, " + uni(r[6]) + "</h2>" : "") : "<h2>Best Swiss-Prot match (DIAMOND, Swiss-Prot 2026_03)</h2>");
    if (r[6]) {
      if (!box) h.push(table([["Accession"], ["Protein"], ["Identity (%)", 1], ["Coverage of your protein (%)", 1], ["E-value", 1]],
        [[uni(r[6]),
          esc(r[5]), r[7] === null ? "" : r[7].toFixed(1), r[8] === null ? "" : Math.round(r[8]), esc(d.sp_e)]]));
      if (!box) h.push("<p class=\"note\">Where on your protein the match lies is not in the table this page was built from.</p>");
      if (d.go.length) h.push("<h3>GO terms of the matched Swiss-Prot entry (not measured on your protein)</h3><p class=\"chips\">" + d.go.map(function (g) {
        return "<a href=\"https://www.ebi.ac.uk/QuickGO/term/" + encodeURIComponent(g[0]) + "\">" + esc(g[1] || g[0]) + "</a>"; }).join(" · ") + "</p>");
      if (d.kw.length) h.push("<h3>UniProt keywords of that entry</h3><p>" + d.kw.map(esc).join(", ") + "</p>");
    } else if (!box) h.push("<p>No Swiss-Prot match at E-value 0.001 or below.</p>");
    // Pfam (a drawn page has hmmscan's Pfam domains along the protein)
    if (!box) h.push("<h2>Pfam domains (hmmsearch, gathering thresholds)</h2>");
    if (box) { /* drawn: see the hmmscan row */ } else if (d.pfam.length) {
      var seen = {}, rows = [];
      d.pfam.forEach(function (p) { seen[p[0]] = seen[p[0]] || { n: p[1], c: 0 }; seen[p[0]].c++; });
      Object.keys(seen).forEach(function (a) {
        rows.push([pfam(a),
          esc(seen[a].n), seen[a].c]);
      });
      h.push(table([["Accession"], ["Name"], ["Times listed", 1]], rows));
      h.push("<p class=\"note\">Listed in the table's order; where each domain lies is not in the table.</p>");
    } else h.push("<p>No Pfam domain.</p>");
    // Immune
    var cats = []; for (var i = 0; i < CATS.length; i++) if (r[11] >> i & 1) cats.push(CATS[i]);
    h.push("<h2>Immune-gene candidate</h2>");
    if (cats.length) h.push(table([["Category"], ["Evidence"]], d.ev.map(function (e) { return [esc(e[0]), accHtml(e[1])]; })));
    else h.push("<p>Not in the immune-gene table.</p>");
    // jackhmmer (a drawn page lists its best hit in a row of the drawing; every hit here)
    h.push("<h2>jackhmmer hits from immune seed proteins</h2>");
    if (box && d.jh.length <= 1) h.pop();
    else if (d.jh.length) h.push(table([["Seed family"], ["Seed protein"], ["Round", 1], ["E-value", 1], ["Score (bits)", 1]],
      d.jh.map(function (j) { return [esc(j[0]), esc(j[1]), j[2], esc(j[3]), esc(j[4])]; })));
    else h.push("<p>None.</p>");
    // InterProScan (a drawn page has every match along the protein)
    if (d.ips && !box) {
      h.push("<h2>InterProScan matches</h2>");
      if (d.ips === "not run") h.push("<p>Not run: InterProScan was not run on proteins over 5,000 amino acids.</p>");
      else if (!d.ips.length) h.push("<p>No match.</p>");
      else h.push(table([["Member database"], ["Signature"], ["Description"], ["Residues", 1], ["Score or E-value", 1], ["InterPro entry"]],
        d.ips.map(function (x) {
          return [esc(x[0]), x[0] === "Pfam" ? pfam(x[1]) : "<span class=\"mono\">" + esc(x[1]) + "</span>", accHtml(x[2]), x[3] + "-" + x[4], esc(x[5]),
            x[6] ? ipr(x[6]) + " " + esc(x[7]) : ""];
        })) + "<p class=\"note\">InterProScan 5.59, all member databases, run on the whole proteome; in order along your protein.</p>");
    }
    // UniFIRE
    if (d.uf) {
      h.push("<h2>UniFIRE predictions (UniProt's annotation rules)</h2>");
      if (d.uf === "not run") h.push("<p>Not run: UniFIRE was not run on proteins over 5,000 amino acids.</p>");
      else if (!d.uf.length) h.push("<p>No prediction.</p>");
      else {
        var LAB = { "keyword": "UniProt keyword", "xref.GO": "GO term", "comment.subcellular_location": "where in the cell",
          "comment.similarity": "family", "comment.catalytic_activity": "reaction it catalyses", "comment.function": "function",
          "comment.cofactor": "cofactor", "comment.pathway": "pathway", "comment.subunit": "subunit", "comment.domain": "domain",
          "comment.PTM": "modification", "protein.recommendedName.fullName": "name", "protein.recommendedName.shortName": "short name",
          "protein.recommendedName.ecNumber": "EC number", "protein.alternativeName.fullName": "other name",
          "protein.alternativeName.shortName": "other short name", "protein.alternativeName.ecNumber": "other EC number",
          "gene.name.primary": "gene name" };
        var ORDER = Object.keys(LAB);
        function rank(t) { var i = ORDER.indexOf(t); return i < 0 ? 100 + (t.indexOf("feature.") === 0 ? 0 : 1) : i; }
        function label(t) { return LAB[t] || (t.indexOf("feature.") === 0 ? "site: " + t.slice(8).toLowerCase().replace(/_/g, " ") : t); }
        function rule(x) {
          var url = x[1].indexOf("ARBA") === 0 ? "https://www.uniprot.org/arba/" : "https://www.uniprot.org/unirule/";
          return esc(x[0]) + " <a class=\"mono\" href=\"" + url + encodeURIComponent(x[1]) + "\">" + esc(x[1]) + "</a>";
        }
        var seenV = {}, urows = [];
        d.uf.slice().sort(function (a, b) { return rank(a[2]) - rank(b[2]) || (a[4] || 0) - (b[4] || 0); }).forEach(function (x) {
          var k = x[2] + "|" + x[3] + "|" + x[4] + "|" + x[5];
          if (seenV[k]) { seenV[k].push(rule(x)); return; }
          seenV[k] = [rule(x)];
          urows.push([x, seenV[k]]);
        });
        h.push(table([["What"], ["Prediction"], ["Where on your protein"], ["Rule"]], urows.map(function (z) {
          var x = z[0];
          return [esc(label(x[2])), esc(x[3]), x[4] === null ? "" : x[4] + (x[5] !== x[4] ? "-" + x[5] : ""), z[1].join("<br>")];
        })));
        h.push("<p class=\"note\">A prediction several rules make is listed once, with each rule. Positions are residue numbers on your protein.</p>");
      }
    }
    // Cluster
    h.push("<h2>MMseqs2 cluster</h2>");
    if (r[13] > 1) {
      h.push("<p>" + fmt(r[13]) + " proteins, including this one; its representative, the protein MMseqs2 keeps to stand for the cluster, is " + link(d.rep) + ".</p>");
      h.push(table([["Other member"], ["Name"]], d.mem.map(function (m) { return [link(m[0]), esc(m[1])]; })));
      if (d.more) h.push("<p class=\"note\">and " + fmt(d.more) + " more.</p>");
    } else h.push("<p>Alone in its cluster.</p>");
    main.innerHTML = h.join("");
  }
  if (box) {
    var id = box.dataset.id;
    load(shardOf(id), function () { render(id, byId[id], window.FULL_SHARD[id]); });
  } else {
    window.addEventListener("hashchange", show);
    show();
  }
})();
