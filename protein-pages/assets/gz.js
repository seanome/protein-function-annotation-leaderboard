// Data files are gzip compressed and base64 encoded, so a page downloads only what it shows.
// loadShard fetches one protein file (full/NNN.js); bootPage starts a drawn page from the data
// build_full.py packed into it, then fetches Mol* once the rows are on the screen.
(function () {
  var v = (document.currentScript.src.split("?")[1] || ""), V = v ? "?" + v : "", shards = {};
  window.FULL_SHARD = window.FULL_SHARD || {};
  window.FULL_GZ = window.FULL_GZ || {};
  window.gunzipB64 = function (b64) {
    var s = atob(b64), a = new Uint8Array(s.length);
    for (var i = 0; i < s.length; i++) a[i] = s.charCodeAt(i);
    return new Response(new Blob([a]).stream().pipeThrough(new DecompressionStream("gzip"))).json();
  };
  window.loadScript = function (src) {
    return new Promise(function (ok, no) {
      var el = document.createElement("script");
      el.src = src; el.async = false; el.onload = ok;
      el.onerror = function () { no(new Error("Could not load " + src)); };
      document.body.appendChild(el);
    });
  };
  window.loadShard = function (n, base) {
    return shards[n] || (shards[n] = loadScript(base + "full/" + ("00" + n).slice(-3) + ".js" + V).then(function () {
      return gunzipB64(window.FULL_GZ[n]);
    }).then(function (d) { Object.assign(window.FULL_SHARD, d); delete window.FULL_GZ[n]; }));
  };
  function idle(f) { (window.requestIdleCallback || setTimeout)(f, { timeout: 800 }); }
  window.bootPage = function (c) {
    gunzipB64(c.z).then(function (d) {
      window.PAGE = d.page; window.WHOLE = d.whole; window.NOTRUN = d.notrun; window.HOW = d.how;
      var a = d.add;
      if (a) {
        Object.assign(PAGE.tracks, a.tracks); PAGE.m += Object.keys(a.tracks).length; PAGE.hits.push.apply(PAGE.hits, a.hits);
        PAGE.names = Object.assign(PAGE.names || {}, a.names); PAGE.nohit = a.nohit; PAGE.twelve += a.twelve; PAGE.ss = a.ss || null;
        WHOLE = window.WHOLE = WHOLE.filter(function (w) { return a.move.indexOf(w[0]) < 0; }); Object.assign(HOW, a.how);
      }
      return loadScript(c.dense + V);
    }).then(function () {
      return loadScript(c.base + "full/meta.js" + V);
    }).then(function () {
      return loadScript(c.base + "assets/viewer.js" + V);
    }).then(function () {
      if (!c.molstar) return;
      idle(function () {
        var l = document.createElement("link");
        l.rel = "stylesheet"; l.href = c.molstar[0]; document.head.appendChild(l);
        loadScript(c.molstar[1]).then(function () { return loadScript(c.molstar[2]); }).then(function () { return loadScript(c.molstar[3]); });
      });
    }).catch(function (e) { document.getElementById("grid").innerHTML = "<p class=\"note\">" + e.message + "</p>"; });
  };
})();
