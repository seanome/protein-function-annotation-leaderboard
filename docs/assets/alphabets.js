/* Alphabet page: the chips above the table pick which identity bin the usable-k column shows.
   All three bins are in the page; this only hides two of them. */
(function () {
  "use strict";
  const chips = [...document.querySelectorAll(".kbins .chip")];
  const nOf = {};
  for (const kb of document.querySelectorAll("td.uk .kb")) nOf[kb.dataset.bin] = kb.dataset.n;
  function show(bin) {
    for (const c of chips) c.setAttribute("aria-pressed", String(c.dataset.kbin === bin));
    for (const kb of document.querySelectorAll("td.uk .kb")) kb.hidden = kb.dataset.bin !== bin;
    const n = document.getElementById("kbin-n");
    if (n && nOf[bin]) n.textContent = "n = " + Number(nOf[bin]).toLocaleString("en-US") + " pairs";
  }
  for (const c of chips) c.addEventListener("click", () => show(c.dataset.kbin));
  window.__al = {show};
})();
