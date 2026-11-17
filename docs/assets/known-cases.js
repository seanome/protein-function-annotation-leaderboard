// The Known cases page: the cost columns follow the two price boxes, and a mark in Figure 18
// opens its case. Everything the page shows comes from the HTML the builder wrote; this file
// only recomputes cost and moves the view.
(function () {
  "use strict";

  function money(x) {
    if (!isFinite(x)) return "—";
    if (x === 0) return "$0";
    if (x < 0.01) return "$" + x.toFixed(4);
    if (x < 1) return "$" + x.toFixed(3);
    return "$" + x.toFixed(2);
  }

  function costs() {
    var rcpu = parseFloat(document.getElementById("kc-rcpu").value);
    var rgpu = parseFloat(document.getElementById("kc-rgpu").value);
    document.querySelectorAll(".cost1k").forEach(function (el) {
      var cpu = parseFloat(el.dataset.cpu) || 0;
      var gpu = parseFloat(el.dataset.gpu) || 0;
      el.textContent = money(cpu * rcpu + gpu * rgpu);
      el.title = cpu.toFixed(4) + " core-hours and " + gpu.toFixed(4) +
        " GPU-hours per 1,000 queries, measured in the Table 1 run of this tool against one " +
        "relative's proteome, not in a search of Swiss-Prot";
    });
  }

  function open(id) {
    var card = document.getElementById("case-" + id);
    if (card) card.scrollIntoView({ behavior: "smooth", block: "start" });
  }

  document.addEventListener("DOMContentLoaded", function () {
    var a = document.getElementById("kc-rcpu");
    var b = document.getElementById("kc-rgpu");
    if (a && b) {
      a.addEventListener("input", costs);
      b.addEventListener("input", costs);
      costs();
    }
    document.querySelectorAll(".kcmatrix .mk").forEach(function (g) {
      g.addEventListener("click", function () { open(g.dataset.case); });
      g.addEventListener("keydown", function (e) {
        if (e.key === "Enter" || e.key === " ") { e.preventDefault(); open(g.dataset.case); }
      });
    });
  });
})();
