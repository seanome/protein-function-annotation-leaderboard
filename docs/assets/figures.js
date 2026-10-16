/* Hover and tap text for Figures 11 to 17: any mark with a data-tip attribute inside a
   figure shows that text in a box beside it. Tap shows it on a phone; Escape or a tap
   elsewhere hides it. */
(function () {
  "use strict";
  let tip = null;
  const box = () => {
    if (!tip) {
      tip = document.createElement("div");
      tip.className = "tip fg-tip";
      tip.setAttribute("role", "status");
      document.body.appendChild(tip);
    }
    return tip;
  };
  const show = (el, x, y) => {
    const t = box();
    t.textContent = el.getAttribute("data-tip");
    t.style.display = "block";
    const w = t.offsetWidth, h = t.offsetHeight;
    const left = Math.min(window.scrollX + x + 12, window.scrollX + document.documentElement.clientWidth - w - 8);
    t.style.left = Math.max(window.scrollX + 8, left) + "px";
    t.style.top = (window.scrollY + y - h - 10 < window.scrollY ? window.scrollY + y + 14 : window.scrollY + y - h - 10) + "px";
  };
  const hide = () => { if (tip) tip.style.display = "none"; };
  const target = (e) => e.target.closest && e.target.closest("figure.fgfig [data-tip]");
  document.addEventListener("mouseover", (e) => { const el = target(e); if (el) show(el, e.clientX, e.clientY); });
  document.addEventListener("mouseout", (e) => { if (target(e)) hide(); });
  document.addEventListener("click", (e) => { const el = target(e); if (el) show(el, e.clientX, e.clientY); else hide(); });
  document.addEventListener("focusin", (e) => {
    const el = target(e);
    if (!el) return;
    const r = el.getBoundingClientRect();
    show(el, r.right, r.top);
  });
  document.addEventListener("focusout", (e) => { if (target(e)) hide(); });
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") hide(); });
})();
