// The 3-D structure on a drawn page, in Mol* (PDBe's build, pdbe-molstar), with Mol*'s own control
// panels and sequence panel. The model comes from window.PDB (proteins/<id>.pdb.js) and is coloured
// by the predictor's confidence per residue (pLDDT, from the B-factor column; data-scale turns a 0
// to 1 file into 0 to 100) in the pLDDT strip's colours (build.py PLDDT).
// Linked highlighting: window.STRUCT.show(a, b) colours residues a to b in the highlight colour
// (data-hl) on the structure, the rest grey, and outlines them on the pLDDT strip with the rest
// faded; protein.js calls it when a hit box or a table row is pointed at. Dragging across the strip
// (#pstrip), or clicking a hit box or table row, keeps a selection until the strip is clicked again.
// Pointing at a residue in Mol* marks it on the strip.
(function () {
  const el = document.getElementById("mol");
  if (!el || !window.PDB || !window.PDBeMolstarPlugin) return;
  const scale = +el.dataset.scale || 1, n = +el.dataset.n || 0, HL = el.dataset.hl || "#E6007E";
  const bin = (p) => (p >= 90 ? "#0053D6" : p >= 70 ? "#65CBF3" : p >= 50 ? "#FFDB13" : "#FF7D45");

  // pLDDT of each residue, by residue number, from the CA atoms of the PDB text
  const plddt = {};
  for (const line of window.PDB.split("\n")) {
    if (line.startsWith("ATOM") && line.slice(12, 16).trim() === "CA") plddt[+line.slice(22, 26)] = +line.slice(60, 66) * scale;
  }
  // runs of neighbouring residues in one pLDDT colour: the base colouring, as Mol* selections
  const base = [];
  Object.keys(plddt).map(Number).sort((x, y) => x - y).forEach((r) => {
    const c = bin(plddt[r]), last = base[base.length - 1];
    if (last && last.color === c && last.end_residue_number === r - 1) last.end_residue_number = r;
    else base.push({ start_residue_number: r, end_residue_number: r, color: c });
  });

  const viewer = new PDBeMolstarPlugin();
  const url = URL.createObjectURL(new Blob([window.PDB], { type: "text/plain" }));
  let ready = false, want = null;
  viewer.render(el, {
    customData: { url, format: "pdb", binary: false },
    bgColor: "white", sequencePanel: true, hideControls: true, pdbeLink: false, reactive: true, hideStructure: ["water"],
  });
  viewer.events.loadComplete.subscribe(() => { ready = true; paint(want); });

  function paint(range) {
    if (!ready) { want = range; return; }
    if (range) viewer.visual.select({ data: [{ start_residue_number: range[0], end_residue_number: range[1], color: HL }], nonSelectedColor: "#D9D9D9" });
    else viewer.visual.select({ data: base });
  }

  const strip = document.getElementById("pstrip");
  const box = document.getElementById("pbrush"), m1 = document.getElementById("pmask1"), m2 = document.getElementById("pmask2");
  const said = document.getElementById("molsel");
  let pinned = null, drawn = "";  // pinned: [a, b] kept until the strip is clicked; drawn: what Mol* shows now

  function mark(a, b, kept) {
    const key = a ? `${a}-${b}` : "";
    if (!a) {
      [box, m1, m2].forEach((r) => r && r.setAttribute("width", 0));
      if (said) said.textContent = "";
    } else {
      a = Math.max(1, Math.min(a, b)); b = Math.min(n || b, Math.max(a, b));
      if (box) { box.setAttribute("x", a - 1); box.setAttribute("width", b - a + 1); }
      if (m1) m1.setAttribute("width", a - 1);
      if (m2) { m2.setAttribute("x", b); m2.setAttribute("width", Math.max(0, n - b)); }
      if (said) said.textContent = `Shown in magenta: residues ${a} to ${b}` + (kept ? " (click the strip to clear)." : ".");
    }
    if (key !== drawn) { drawn = key; paint(a ? [a, b] : null); }
  }
  window.STRUCT = {
    show(a, b) {
      if (a) mark(a, b, false);
      else if (pinned) mark(pinned[0], pinned[1], true);
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
      if (moved) mark(Math.min(start, r), Math.max(start, r), false);
    });
    strip.addEventListener("pointerup", (e) => {
      if (start === null) return;
      const r = resAt(e);
      pinned = moved ? [Math.min(start, r), Math.max(start, r)] : null;
      start = null;
      window.STRUCT.show(null);
    });
  }

  // Pointing at a residue in Mol*: its place on the strip (Mol* shows its own label)
  document.addEventListener("PDB.molstar.mouseover", (e) => {
    const d = (e && e.eventData) || {};
    const r = d.auth_seq_id || d.residueNumber || d.seq_id;
    if (r && box && !pinned && !drawn) { box.setAttribute("x", r - 1); box.setAttribute("width", 1); }
  });
  document.addEventListener("PDB.molstar.mouseout", () => { if (!pinned && !drawn && box) box.setAttribute("width", 0); });
})();
