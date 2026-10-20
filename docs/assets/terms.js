/* Definitions beside column names and labels. Each is a small "?" button followed by its
   definition (role="tooltip", linked by aria-describedby), which ends with a link to the
   explainer. The definition shows on mouse hover, on keyboard focus and on tap, and hides
   on Escape, on a tap or click elsewhere, or when focus leaves both. It is placed with
   position:fixed so a table that scrolls sideways does not cut it off.
   Works for buttons drawn later (Table 1 is redrawn on every filter change): all events
   are caught on the document. */
(function () {
  "use strict";
  let open = null, leaveTimer = null, quiet = false;
  const tipOf = (b) => document.getElementById(b.getAttribute("aria-controls"));

  function place(b, t) {
    const r = b.getBoundingClientRect();
    const w = Math.min(320, window.innerWidth - 16);
    t.style.width = w + "px";
    t.style.left = Math.min(Math.max(8, r.left + r.width / 2 - w / 2), window.innerWidth - w - 8) + "px";
    const h = t.offsetHeight;
    let top = r.bottom + 6;
    if (top + h > window.innerHeight - 8 && r.top - h - 6 > 8) top = r.top - h - 6;
    t.style.top = top + "px";
  }
  function show(b) {
    clearTimeout(leaveTimer);
    if (open && open !== b) hide(open);
    const t = tipOf(b);
    if (!t) return;
    t.hidden = false;
    place(b, t);
    b.setAttribute("aria-expanded", "true");
    open = b;
  }
  function hide(b) {
    const t = tipOf(b);
    if (t) t.hidden = true;
    b.setAttribute("aria-expanded", "false");
    if (open === b) open = null;
  }
  const inOpen = (node) => open && (open.contains(node) || (tipOf(open) && tipOf(open).contains(node)));

  document.addEventListener("focusin", (e) => {
    const b = e.target.closest && e.target.closest("button.term");
    if (b && !quiet) show(b);
    else if (open && !inOpen(e.target)) hide(open);
  });
  document.addEventListener("click", (e) => {
    const b = e.target.closest && e.target.closest("button.term");
    if (b) { show(b); return; }
    if (open && !inOpen(e.target)) hide(open);
  });
  document.addEventListener("mouseover", (e) => {
    const b = e.target.closest && e.target.closest("button.term");
    if (b) show(b);
    else if (inOpen(e.target)) clearTimeout(leaveTimer);
  });
  document.addEventListener("mouseout", (e) => {
    if (!open || !inOpen(e.target) || inOpen(e.relatedTarget)) return;
    const b = open;
    // Leave time to move the pointer from the button onto the definition's link.
    leaveTimer = setTimeout(() => { if (open === b && document.activeElement !== b) hide(b); }, 250);
  });
  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape" && open) { const b = open; hide(b); quiet = true; b.focus({ preventScroll: true }); quiet = false; }
  });
  const replace = () => { if (open) place(open, tipOf(open)); };
  window.addEventListener("scroll", replace, true);
  window.addEventListener("resize", replace);
})();
