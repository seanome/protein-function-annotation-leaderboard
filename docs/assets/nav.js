/* Header menus ("More", and the one-line page menu on a phone): close on a click
   elsewhere, on Escape, or when the other one opens. */
(function () {
  "use strict";
  const menus = () => [...document.querySelectorAll("header.site details")];
  document.addEventListener("click", (e) => {
    for (const d of menus()) if (d.open && !d.contains(e.target)) d.open = false;
  });
  document.addEventListener("keydown", (e) => {
    if (e.key !== "Escape") return;
    for (const d of menus()) if (d.open) { d.open = false; d.querySelector("summary").focus(); }
  });
})();
