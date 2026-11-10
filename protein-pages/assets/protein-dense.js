// protein.html: the dense page (assets/dense-bhf.js) for any of the 45,339 proteins, from the
// whole-proteome tables (build_full.py's data files, full/NNN.js). A protein with a drawn page goes
// straight to it. The others have each program's best hit only: MMseqs2, ProstT5 + Foldseek and
// Foldseek on the predicted structure with their positions, InterProScan's matches, and the
// programs that name a hit with no position (kmerseek's best hit, DIAMOND, jackhmmer from immune
// seeds, FANTASIA, UniFIRE). With the protein's sequence in its data file (build_full.py
// --proteome-fasta) a hit row gives its residue composition and opens the two proteins aligned in the
// browser; without it, neither (PAGE.noseq). This file builds PAGE and the globals
// dense-bhf.js reads (WHOLE, NOTRUN, HOW, ALN), then loads dense-bhf.js and viewer.js.
(function () {
  var C = window.DENSE_CONST, N_SHARDS = 96, r;
  function shardOf(id) {
    var h = 0;
    for (var i = 0; i < id.length; i++) h = (Math.imul(h, 31) + id.charCodeAt(i)) >>> 0;
    return h % N_SHARDS;
  }
  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; }); }
  function script(src, then) {
    var el = document.createElement("script");
    el.src = src; el.onload = then;
    el.onerror = function () { document.getElementById("grid").innerHTML = "<p class=\"note\">Could not load " + esc(src) + ".</p>"; };
    document.body.appendChild(el);
  }
  function uni(acc) {
    return "<a class=\"mono\" href=\"https://www.uniprot.org/uniprotkb/" + encodeURIComponent(acc) + "/entry\" target=\"_blank\" rel=\"noopener\">" + esc(acc) + "</a>";
  }
  var id = location.hash.slice(1);
  try { id = decodeURIComponent(id); } catch (e) { /* malformed escape: keep the raw text */ }
  window.addEventListener("hashchange", function () { location.reload(); });
  document.querySelector(".titlebar h1").textContent = id || "No protein named";
  function nothing() { document.getElementById("grid").innerHTML = "<p class=\"note\">No protein " + esc(id) + " in this proteome. <a href=\"index.html\">Search all proteins</a>.</p>"; }
  if (!id) return nothing();
  loadShard(shardOf(id), "").then(function () {
    var d = window.FULL_SHARD[id];
    if (!d) return nothing();
    r = d.r;
    if (r[14]) location.replace(r[14]);  // the drawn page is this protein's page
    else build(d);
  }).catch(function (e) { document.getElementById("grid").innerHTML = "<p class=\"note\">" + esc(e.message) + "</p>"; });

  function build(d) {
    var L = r[1], name = r[3] || r[2] || "";
    document.title = id + (name ? " " + name : "") + ": what every annotation program found";
    if (name) document.querySelector(".titlebar .pn").textContent = name;
    document.querySelector(".titlebar .len").textContent = L.toLocaleString("en-US") + " amino acids" +
      (d.iso ? " · isoform " + (d.iso.map(function (m) { return m[0]; }).indexOf(id) + 1) + " of " + d.iso.length + " (listed below)" : "");
    var hits = [], tracksRun = {};
    function best(k, x, extra) {  // a best-hit table row: [accession, identity 0-1, start, end, E-value, bits, name, species, ..., target start, target end]
      tracksRun[k] = true;
      if (!x) return;
      var h = { k: k, l: x[0], t: x[0], q: [x[2], x[3]], s: [x[4], String(x[1]), x[5]].concat(extra || []), e: parseFloat(x[4]),
                n: x[6] || undefined, o: x[7] || undefined, f: null };
      if (x.length >= 10 && isFinite(x[x.length - 2])) h.r = [x[x.length - 2], x[x.length - 1]];  // the range on the other protein
      if (d.tl && d.tl[x[0]]) h.tl = d.tl[x[0]];
      hits.push(h);
    }
    best("mmseqs2", d.mm);
    best("prostt5", d.fs, ["not in the file"]);
    var shortEnough = L <= C.struct_max_aa;
    if (shortEnough) best("foldseek", d.fst, d.fst && d.fst.length > 11 ? [String(d.fst[8]), String(d.fst[9])] : ["", ""]);
    var ipsRun = d.ips !== "not run" && d.ips !== null && d.ips !== undefined;
    if (ipsRun) {
      tracksRun.interproscan = true;
      d.ips.forEach(function (m) {  // [analysis, signature, description, start, end, score, InterPro entry, its name]
        var word = (C.ips_drawn[m[0]] || {})[m[1]] || (C.ips_drawn[m[0]] || {})["*"];
        var h = { k: "interproscan", l: word || m[2] || m[1], q: [m[3], m[4]], s: [m[0], m[1], m[5] || "-", m[6] || "none"],
                  e: C.ips_evalue.indexOf(m[0]) >= 0 && isFinite(parseFloat(m[5])) ? parseFloat(m[5]) : null,
                  f: m[0] === "Pfam" ? [(m[2] || m[1]) + " (" + m[1] + ")"] : [] };
        if (word) h.x = 1;  // a prediction (disorder, coil, membrane span, signal peptide): drawn dashed, never counted
        hits.push(h);
      });
    }
    if (d.km) {  // [setting, E-value, accession, name, species, settings, moltype, k, start, end, target start, target end, target residues]
      tracksRun.kmerseek = true;
      if (d.km.length > 6) {
        var kq = d.q ? d.q.slice(d.km[8] - 1, d.km[9]) : "", same = 0, kt = d.km[12];
        for (var ki = 0; ki < kq.length; ki++) if (kq[ki] === kt[ki]) same++;
        var kh = { k: "kmerseek", l: d.km[2], t: d.km[2], q: [d.km[8], d.km[9]], r: [d.km[10], d.km[11]], tsub: kt, e: d.km[1], s: [String(d.km[1]), d.km[0], d.km[5] + " settings"],
                   n: d.km[3] || undefined, o: d.km[4] || undefined, f: null, m: d.km[6], kk: d.km[7] };
        if (kq.length === kt.length && kq.length) kh.id = [same, kq.length];
        if (d.tl && d.tl[d.km[2]]) kh.tl = d.tl[d.km[2]];
        hits.push(kh);
      }
    }
    if (d.dm) {  // [start, end, target start, target end, identical, columns]: DIAMOND's best hit, realigned
      tracksRun.diamond = true;
      hits.push({ k: "diamond", l: r[6], t: r[6], q: [d.dm[0], d.dm[1]], r: [d.dm[2], d.dm[3]], id: [d.dm[4], d.dm[5]], e: parseFloat(d.sp_e),
                  s: [d.sp_e, r[7] === null ? "" : String(r[7]), r[8] === null ? "" : String(r[8])], n: r[5] || undefined, o: d.sp_org || undefined, f: null, tl: d.tl && d.tl[r[6]] });
    } else if (r[6] === "" || r[6] === undefined) tracksRun.diamond = true;
    var ufNo = 0;
    if (d.uf && d.uf !== "not run") {  // one hit for each rule and place: a site's own position, else where the rule's signatures lie
      tracksRun.unifire = true;
      var byRule = {};
      d.uf.forEach(function (u) {
        var at = u[4] ? [[u[4], u[5] || u[4], "the rule's own position"]] : (u[6] || []).map(function (i) { var m = d.ips[i]; return [m[3], m[4], m[1] + (m[6] ? " / " + m[6] : "")]; });
        if (!at.length) { ufNo++; at = [[1, L, "whole protein: UniProt's current version of this rule names no signature that InterProScan found here"]]; }
        at.forEach(function (p) {
          var key = u[1] + ":" + p[0] + "-" + p[1], o = byRule[key] || (byRule[key] = { k: "unifire", l: u[3], q: [p[0], p[1]], s: [u[0], u[1], u[2], p[2]], e: null, f: [], v: [] });
          if (o.v.indexOf(u[3]) < 0) o.v.push(u[3]);
        });
      });
      Object.keys(byRule).forEach(function (key) { var o = byRule[key]; o.l = o.v[0] + (o.v.length > 1 ? " +" + (o.v.length - 1) : ""); o.n = o.v.join("; "); hits.push(o); });
    }
    // the order dense-bhf.js expects: by track, best E-value first
    var order = Object.keys(C.tracks);
    hits.sort(function (a, b) { return order.indexOf(a.k) - order.indexOf(b.k) || ((a.e === null ? 1e300 : a.e) - (b.e === null ? 1e300 : b.e)) || a.q[0] - b.q[0]; });

    // the programs that name a hit with no position: [key, name, how, text, flag, second line]
    var W = [];
    if (d.km !== null && d.km !== undefined && !(d.km && d.km.length > 6)) W.push(["w_kmerseek", "kmerseek, best hit", "run on Sherlock, not by Biomni: best hit against Swiss-Prot 2026_03 over 35 alphabet+k settings; this site keeps only that best hit, so no position",
      d.km ? uni(d.km[2]) + " " + esc(d.km[3]) + (d.km[4] ? ", <i>" + esc(d.km[4]) + "</i>" : "") : "no hit (E-value " + C.kmerseek_max_e + " or below)", d.km ? 1 : 0,
      d.km ? "E-value " + esc(String(+d.km[1].toPrecision(2))) + " with " + esc(d.km[0]) + ", the best of the " + d.km[5] + " settings whose E-value fit worked · whole protein, no position" : undefined]);
    if (!d.dm && r[6]) W.push(["w_diamond", "DIAMOND against Swiss-Prot", "run by Biomni; E-value 0.001 or below",
      r[6] ? uni(r[6]) + " " + esc(r[5]) + (d.sp_gene ? " (" + esc(d.sp_gene) + ")" : "") + (d.sp_org ? ", <i>" + esc(d.sp_org) + "</i>" : "") : "no hit (E-value 0.001 or below)",
      r[6] ? 1 : 0, r[6] ? "E-value " + esc(d.sp_e) + (r[7] === null ? "" : ", " + Math.round(r[7]) + "% identity over " + Math.round(r[8] || 0) + "% of your protein") + " · whole protein, no position" : undefined]);
    var jh = d.jh || [];
    W.push(["w_jackhmmer", "jackhmmer with immune-gene seeds", "run by Biomni in the reverse direction: immune-gene seed proteins (24 in Biomni's run; 20 have a hit in its table) as queries against the Botryllus proteome",
      jh.length ? "found from the seed " + esc(jh[0][1]) + " (" + esc(jh[0][0]) + ")" + (jh.length > 1 ? " and " + (jh.length - 1) + " more" : "") : "no hit", jh.length ? 1 : 0,
      jh.length ? "E-value " + esc(jh[0][3]) + ", round " + jh[0][2] + " · whole protein, no position" : undefined]);
    if (d.fa !== null && d.fa !== undefined) {  // FANTASIA ran only on proteins with no Swiss-Prot match or in the immune-gene table
      var fa = d.fa, counts = fa[4] < C.fantasia_max;
      W.push(["w_fantasia", "FANTASIA (language model)", "run by Biomni; names its most similar protein at any distance",
        uni(fa[0]) + " " + esc(fa[1]) + (fa[2] ? " (" + esc(fa[2]) + ")" : "") + (fa[3] ? ", <i>" + esc(fa[3]) + "</i>" : "") + " · its most similar protein", counts ? 1 : 2,
        "embedding distance " + fa[4].toFixed(2) + (counts ? "" : ", weak: above " + C.fantasia_max.toFixed(2)) + " · an embedding match, so it has no position or alignment"]);
    }
    // UniFIRE skipped proteins over 5,000 aa

    // why each other row has nothing
    var NR = {};
    var tracks = C.tracks;
    Object.keys(tracks).forEach(function (k) {
      if (tracksRun[k]) return;
      NR[k] = k === "foldseek" && !shortEnough ? "not run: structures were predicted only for proteins of " + C.struct_max_aa + " aa or less" :
        k === "interproscan" || k === "unifire" ? "not run: proteins over 5,000 amino acids were skipped" :
        k === "kmerseek" ? "no position kept for this hit" :
        k === "hmmscan" ? "run only on the 51 test proteins; Pfam domains here: InterProScan's Pfam matches on its row, and hmmsearch's names under the grid" :
        "run only on the 51 test proteins";
    });

    // the sentence at each E-value cutoff: programs with a hit at that cutoff, of the programs run here
    var ran = Object.keys(tracksRun).map(function (k) { return C.tracks[k].p; }).filter(function (p, i, a) { return a.indexOf(p) === i; });
    var wholeRan = W;
    var M = ran.length + wholeRan.length;
    var note = "<p class=\"note\">This protein is not one of the 51 test proteins, so only each program's best hit is known (tables of the whole proteome, from Biomni's runs and our kmerseek run on Sherlock); " +
      "whether two programs agree on a Pfam family is not judged here. " + (d.q ? "Click a hit for the residues: the two proteins aligned in your browser " +
      "(BLOSUM62, gap open -11, extend -1), since the tables keep no alignment. " : "This protein's sequence is not on this site, so a hit opens no alignment. ") +
      "phmmer, jackhmmer against Swiss-Prot, hmmscan, hhsearch and MMseqs2 three rounds were run only on those 51.</p>";
    var sent = [1, 0, -1, -2, -3].reverse().map(function (v) {
      var progs = {};
      hits.forEach(function (h) { if (!h.x && (h.e === null || h.e <= Math.pow(10, v))) progs[C.tracks[h.k].p] = 1; });
      // kmerseek's best hit counts at the cutoff like the drawn rows; the other whole-protein programs by their own cutoffs
      var n = Object.keys(progs).length + W.filter(function (w) { return w[4] >= 1 && (w[0] !== "w_kmerseek" || d.km[1] <= Math.pow(10, v)); }).length;
      return [v, n === 0 ? "None of the " + M + " programs run reports a hit (best hits only)." :
        n + " of the " + M + " programs run report" + (n === 1 ? "s" : "") + " a hit (best hits only).", note];
    });
    var names = { kmerseek: ["region_evalue", "alphabet+k", "best of"], diamond: ["evalue", "% identity (DIAMOND's own)", "% of your protein (DIAMOND's own)"],
                  unifire: ["source", "rule", "what it predicts", "where it points"], mmseqs2: ["evalue", "fident", "bits"], prostt5: ["evalue", "fident", "bits", "alntmscore"],
                  foldseek: ["evalue", "fident", "bits", "alntmscore", "lddt"], interproscan: ["member database", "signature", "its score or E-value", "InterPro entry"] };
    window.PAGE = {
      id: id, pn: name, q: d.q || new Array(L + 1).join("X"), noseq: d.q ? 0 : 1, px: 640 / L, w: 640, top: 10, m: ran.length, sent: sent, hits: hits, td: {}, ds: {},
      names: names, alpha: C.alpha || {}, nohit: { mmseqs2: "E ≤ 0.001: its table keeps nothing weaker", diamond: "E ≤ 0.001: its table keeps nothing weaker", kmerseek: "E ≤ " + C.kmerseek_max_e, prostt5: "E ≤ 10", foldseek: "E ≤ 10", unifire: "no rule fired" }, tracks: tracks, pfam_color: C.pfam_color, tseq_files: C.tseq_files || 1, tbase: "", agree: [],
      twelve: "the " + ran.length + " programs with positions run here (" + ran.map(function (p) { return C.program_names[p]; }).join(", ") +
        ") plus " + wholeRan.length + " that name a hit with no position: " + wholeRan.map(function (w) { return w[1]; }).join(", ")
    };
    window.WHOLE = W; window.NOTRUN = NR; window.ALN = {};
    window.HOW = { kmerseek: "best hit only: run on Sherlock, not by Biomni, against Swiss-Prot 2026_03 over 35 alphabet+k settings", diamond: "best hit only: run by Biomni; its table has no positions, so this pair was aligned again (BLOSUM62, local, gap open -11, extend -1)",
                   unifire: "run by Biomni on InterProScan's matches. Rules work on signatures, so each prediction is placed where the rule's signatures lie on this protein; a site uses the rule's own position",
                   mmseqs2: "best hit only: Biomni's run against all of Swiss-Prot 2026_03 keeps one hit per protein here, at E-value 0.001 or below",
                   prostt5: "best hit only: ProstT5 3Di letters, Foldseek version 10; one hit per protein, E-value 10 or below",
                   foldseek: "best hit only: Foldseek on the structure ESMFold2 predicted (proteins of " + C.struct_max_aa + " aa or less); E-value 10 or below",
                   interproscan: "InterProScan 5.59, all member databases; every match" };
    var b = document.getElementById("biomni");
    b.dataset.id = id;
    script("assets/dense-bhf.js" + C.ver, function () { script("assets/viewer.js" + C.ver, function () {}); });
  }
})();
