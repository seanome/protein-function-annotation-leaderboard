// The 3-D structure on a drawn page: ribbons from window.PDB (proteins/<id>.pdb.js), drawn with
// 3Dmol.js and coloured by the predictor's confidence per residue (the B-factor column; data-scale
// turns a 0 to 1 file into 0 to 100), in the pLDDT strip's colours (build.py PLDDT).
// Linked highlighting: window.STRUCT.show(a, b) colours residues a to b in the highlight colour
// (data-hl) on the structure and outlines them on the pLDDT strip, with the rest grey and faded;
// protein.js calls it when a hit box or a table row is pointed at. Dragging across the strip
// (#pstrip), or clicking a hit box or table row, selects residues the same way and keeps them
// until the strip is clicked again.
// Pointing at the structure shows the residue's number and pLDDT.
(function () {
  const el = document.getElementById("mol");
  if (!el || !window.PDB || !window.$3Dmol) return;
  const scale = +el.dataset.scale || 1, n = +el.dataset.n || 0, HL = el.dataset.hl || "#E6007E";
  const colour = (b) => { const p = b * scale; return p >= 90 ? "#0053D6" : p >= 70 ? "#65CBF3" : p >= 50 ? "#FFDB13" : "#FF7D45"; };
  const v = $3Dmol.createViewer(el, { backgroundColor: "white" });
  v.addModel(window.PDB, "pdb");
  const base = { cartoon: { colorfunc: (atom) => colour(atom.b) } };
  v.setStyle({}, base);
  v.zoomTo();
  v.render();

  const strip = document.getElementById("pstrip");
  const box = document.getElementById("pbrush"), m1 = document.getElementById("pmask1"), m2 = document.getElementById("pmask2");
  const sel = document.getElementById("molsel");
  let pinned = null;  // [a, b] selected on the strip

  function mark(a, b, label) {
    if (!a) {
      v.setStyle({}, base);
      [box, m1, m2].forEach((r) => r && r.setAttribute("width", 0));
      if (sel) sel.hidden = true;
    } else {
      a = Math.max(1, Math.min(a, b)); b = Math.min(n || b, Math.max(a, b));
      v.setStyle({}, { cartoon: { color: "#D9D9D9" } });
      v.setStyle({ resi: `${a}-${b}` }, { cartoon: { color: HL, thickness: 0.6 } });
      if (box) { box.setAttribute("x", a - 1); box.setAttribute("width", b - a + 1); }
      if (m1) { m1.setAttribute("x", 0); m1.setAttribute("width", a - 1); }
      if (m2) { m2.setAttribute("x", b); m2.setAttribute("width", Math.max(0, n - b)); }
      if (sel) { sel.textContent = label || `residues ${a} to ${b}`; sel.hidden = false; }
    }
    v.render();
  }
  window.STRUCT = {
    show(a, b) {
      if (a) mark(a, b);
      else if (pinned) mark(pinned[0], pinned[1], `residues ${pinned[0]} to ${pinned[1]}, selected (click the strip to clear)`);
      else mark(null);
    },
    pin(a, b) { pinned = a ? [a, b] : null; this.show(null); },  // a click on a hit box or table row
  };

  // Drag across the strip to select residues; a click without a drag clears the selection.
  if (strip && n) {
    const resAt = (e) => { const r = strip.getBoundingClientRect(); return Math.max(1, Math.min(n, Math.floor((e.clientX - r.left) / r.width * n) + 1)); };
    let start = null, moved = false;
    strip.addEventListener("pointerdown", (e) => { start = resAt(e); moved = false; strip.setPointerCapture(e.pointerId); e.preventDefault(); });
    strip.addEventListener("pointermove", (e) => {
      if (start === null) return;
      const r = resAt(e);
      if (r !== start) moved = true;
      if (moved) mark(Math.min(start, r), Math.max(start, r));
    });
    strip.addEventListener("pointerup", (e) => {
      if (start === null) return;
      const r = resAt(e);
      pinned = moved ? [Math.min(start, r), Math.max(start, r)] : null;
      start = null;
      window.STRUCT.show(null);
    });
  }

  // Point at the structure: the residue's number and pLDDT, and its place on the strip.
  v.setHoverable({}, true, (atom) => {
    if (!atom || !sel) return;
    sel.textContent = `residue ${atom.resi} (${atom.resn}), pLDDT ${Math.round(atom.b * scale)}`;
    sel.hidden = false;
    if (box && !pinned) { box.setAttribute("x", atom.resi - 1); box.setAttribute("width", 1); }
  }, () => { window.STRUCT.show(null); });
  v.render();
})();
