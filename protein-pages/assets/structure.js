// The 3-D structure on a drawn page: ribbons from window.PDB (proteins/<id>.pdb.js), drawn with
// 3Dmol.js and coloured by the predictor's confidence per residue (the B-factor column; data-scale
// turns a 0 to 1 file into 0 to 100), in the pLDDT strip's colours (build.py PLDDT).
// window.STRUCT.show(a, b) draws residues a to b in black; protein.js calls it on hover.
(function () {
  const el = document.getElementById("mol");
  if (!el || !window.PDB || !window.$3Dmol) return;
  const scale = +el.dataset.scale || 1;
  const colour = (b) => { const p = b * scale; return p >= 90 ? "#0053D6" : p >= 70 ? "#65CBF3" : p >= 50 ? "#FFDB13" : "#FF7D45"; };
  const v = $3Dmol.createViewer(el, { backgroundColor: "white" });
  v.addModel(window.PDB, "pdb");
  const base = { cartoon: { colorfunc: (atom) => colour(atom.b) } };
  v.setStyle({}, base);
  v.zoomTo();
  v.render();
  window.STRUCT = {
    show(a, b) {
      v.setStyle({}, base);
      if (a) v.setStyle({ resi: `${a}-${b}` }, { cartoon: { color: "#1a1a1a" } });
      v.render();
    },
  };
})();
