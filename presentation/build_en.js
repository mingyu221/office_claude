const PptxGenJS = require("pptxgenjs");

const INK="16183A", DEEP="2F3596", MID="6A72CE", ACCENT="C0692B",
      TINT="EEF0FA", TINT2="F7F8FC", WHITE="FFFFFF",
      GRAY="5C6180", TEXT="1C1E33", LINE="D6DAEE", INKSUB="A8AECF";
const F = "Calibri";
const W=13.3, H=7.5, M=0.65, CW=W-2*M;

const pres = new PptxGenJS();
pres.layout = "LAYOUT_WIDE";
pres.author = "CODH project";
pres.title  = "Docking-distilled pLM — progress report 2026-09-16";

let pageNo = 0;
function newSlide(dark){ const s = pres.addSlide(); if (dark) s.background = { color: INK }; return s; }
function head(s, kicker, title, dark){
  s.addText(kicker, { x:M, y:0.40, w:CW, h:0.26, fontFace:F, fontSize:11.5, bold:true,
    color:ACCENT, charSpacing:1.5, isTextBox:true, margin:0 });
  s.addText(title, { x:M, y:0.68, w:CW, h:1.02, fontFace:F, fontSize:26, bold:true,
    color: dark?WHITE:TEXT, isTextBox:true, margin:0, valign:"top", lineSpacingMultiple:1.02 });
}
function foot(s, dark){
  pageNo += 1;
  s.addText(String(pageNo), { x:W-M-0.6, y:H-0.52, w:0.6, h:0.25, fontFace:F, fontSize:9.5,
    color: dark?INKSUB:GRAY, align:"right", isTextBox:true, margin:0 });
}
function bullets(s, items, box){
  const paras = [];
  items.forEach((runs, i) => {
    const rs = Array.isArray(runs) ? runs : [{ text: runs }];
    rs.forEach((r, j) => {
      paras.push({ text: r.text, options: {
        bullet: false, bold: !!r.bold, color: r.color || box.color || TEXT,
        breakLine: (j === rs.length-1) && (i < items.length-1),
        paraSpaceAfter: (j === rs.length-1) ? (box.gap===undefined?9:box.gap) : 0
      }});
    });
  });
  s.addText(paras, { x:box.x, y:box.y, w:box.w, h:box.h, fontFace:F,
    fontSize: box.size||14, color: box.color||TEXT, isTextBox:true, margin:0,
    valign:"top", lineSpacingMultiple:1.16 });
}
function plain(s, txt, box){
  s.addText(txt, { x:box.x, y:box.y, w:box.w, h:box.h, fontFace:F, fontSize:box.size||14,
    bold:!!box.bold, color:box.color||TEXT, align:box.align||"left", isTextBox:true,
    margin:box.margin===undefined?0:box.margin, valign:box.valign||"top",
    lineSpacingMultiple:box.ls||1.13, italic:!!box.italic });
}
function card(s, x,y,w,h, fill, line){
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h, fill:{ color: fill },
    line: line ? { color: line, width:1 } : { type:"none" }, rectRadius:0.06 });
}
function badge(s, x, y, d, txt, fill, col){
  s.addShape(pres.ShapeType.ellipse, { x, y, w:d, h:d, fill:{ color: fill }, line:{ type:"none" } });
  s.addText(txt, { x, y, w:d, h:d, fontFace:F, fontSize: d>=0.44?14:12, bold:true,
    color: col||WHITE, align:"center", valign:"middle", isTextBox:true, margin:0 });
}
function chip(s, x,y,w,h, txt, fill, col, size){
  s.addShape(pres.ShapeType.roundRect, { x, y, w, h, fill:{ color: fill }, line:{ type:"none" }, rectRadius:0.05 });
  s.addText(txt, { x:x+0.06, y, w:w-0.12, h, fontFace:F, fontSize:size||12, bold:true,
    color:col, align:"center", valign:"middle", isTextBox:true, margin:0, lineSpacingMultiple:1.03 });
}
function arrow(s, x, y, txt){
  s.addText(txt||"→", { x, y, w:0.3, h:0.3, fontFace:F, fontSize:15, bold:true,
    color:MID, align:"center", valign:"middle", isTextBox:true, margin:0 });
}
function table(s, rows, box){
  s.addTable(rows, { x:box.x, y:box.y, w:box.w, colW:box.colW, fontFace:F,
    fontSize:box.size||11.5, color:TEXT, border:{ type:"solid", pt:0.75, color:LINE },
    align:"center", valign:"middle", rowH:box.rowH||0.3, autoPage:false });
}
function hrow(cells, fill){
  return cells.map(c => ({ text:c, options:{ bold:true, color:WHITE, fill:{ color: fill||DEEP }, fontSize:11 } }));
}
function drow(cells, opts){
  return cells.map((c,i) => ({ text:c, options: Object.assign({ align: i===0?"left":"center" }, (opts&&opts[i])||{}) }));
}
function stat(s, x,y,w, big, label, col, size){
  const fs = size||34, off = fs*0.0195;
  s.addText(big, { x, y, w, h:off+0.06, fontFace:F, fontSize:fs, bold:true, color:col||DEEP,
    align:"center", isTextBox:true, margin:0 });
  s.addText(label, { x, y:y+off, w, h:0.4, fontFace:F, fontSize: fs>=30?11:10, color:GRAY,
    align:"center", isTextBox:true, margin:0, lineSpacingMultiple:1.06 });
}

/* ── 1. Title ── */
{
  const s = newSlide(true);
  s.addText("PROGRESS REPORT  ·  2026-09-16", { x:M, y:1.55, w:CW, h:0.3, fontFace:F, fontSize:12.5,
    bold:true, color:ACCENT, charSpacing:1.7, isTextBox:true, margin:0 });
  s.addText("Engineering CO dehydrogenases with\nenhanced electron mediator affinity through\ndocking-distilled protein language models",
    { x:M, y:2.05, w:CW-0.6, h:2.1, fontFace:F, fontSize:32, bold:true, color:WHITE,
      isTextBox:true, margin:0, lineSpacingMultiple:1.16 });
  s.addText("Ranking the IPR010047 family by mediator affinity from sequence alone",
    { x:M, y:4.32, w:CW-0.6, h:0.4, fontFace:F, fontSize:14, color:INKSUB, isTextBox:true, margin:0 });
  const tags = ["①  Receptor = dimer", "②  Center = residue-41 CA", "③  157 train / 55 eval", "④  321 dimers assembled"];
  tags.forEach((t,i) => chip(s, M+i*3.0, 5.35, 2.82, 0.5, t, "242764", WHITE, 11.5));
  s.addNotes("Four design decisions were frozen today; the rest of the deck gives the evidence for each.");
  foot(s, true);
}

/* ── 2. Summary ── */
{
  const s = newSlide();
  head(s, "SUMMARY", "Four design decisions were frozen today");
  const D = [
    ["①","Receptor","Dimer","Monomer rejected — correlation with measured Km collapses from 0.88 to 0.55."],
    ["②","Grid center","F41 CE2 → residue-41 CA","Equivalent (ρ 0.968), and applicable sequences rise from 1,613 to 1,750."],
    ["③","Dataset","157 train / 55 eval","212 representatives split by cluster. Zero leakage, frozen before any label exists."],
    ["④","Structures","321 dimers assembled","Assembly bias 0.080 kcal/mol, below the measurement IQR of 0.09."]
  ];
  const cw=5.85, ch=1.95;
  D.forEach((d,i) => {
    const x = M + (i%2)*(cw+0.30), y = 1.80 + Math.floor(i/2)*(ch+0.28);
    card(s, x, y, cw, ch, TINT2, LINE);
    badge(s, x+0.30, y+0.28, 0.44, d[0], DEEP);
    plain(s, d[1], { x:x+0.86, y:y+0.36, w:cw-1.2, size:12.5, bold:true, color:GRAY });
    plain(s, d[2], { x:x+0.30, y:y+0.86, w:cw-0.6, size:21, bold:true, color:DEEP });
    plain(s, d[3], { x:x+0.30, y:y+1.36, w:cw-0.6, size:12, color:TEXT, ls:1.13 });
  });
  plain(s, "IN PROGRESS   Dock 321 structures  →  measure the noise floor (SNR)  →  attach labels  →  train A / B / C / D1 / D2",
    { x:M, y:6.35, w:CW, size:12.5, bold:true, color:ACCENT });
  s.addNotes("This slide alone should carry the meeting; everything after it is supporting evidence.");
  foot(s);
}

/* ── 3. Objective & scope strategy (NEW) ── */
{
  const s = newSlide();
  head(s, "OBJECTIVE  &  SCOPE", "Rank the CODH family by mediator affinity from sequence alone");
  card(s, M, 1.78, CW, 0.90, TINT, LINE);
  plain(s, "Identify CO dehydrogenase sequences with enhanced electron-mediator (EV) affinity by ranking candidates with a protein language model distilled from docking scores.",
    { x:M+0.32, y:1.98, w:CW-0.64, size:14.5, bold:true, color:DEEP, ls:1.18 });

  const C = [
    ["Why it is hard", ACCENT, "Only 21 Km values have been measured — too few to supervise a sequence model directly. The family also has almost no close relatives above 90% identity."],
    ["Approach", DEEP, "Docking generates labels at scale: 288 variants and 321 homologs. A regression head distills them into ESM embeddings, so inference needs no structure and extends to all 1,827 sequences."],
    ["Success criterion", DEEP, "Spearman ρ on the frozen 55-sequence evaluation set, required to exceed the docking noise floor measured from the 109 within-cluster pairs."]
  ];
  C.forEach((c,i) => {
    const x = M + i*4.12;
    card(s, x, 2.82, 3.85, 1.68, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    plain(s, c[0], { x:x+0.26, y:3.00, w:3.33, size:13, bold:true, color:c[1] });
    plain(s, c[2], { x:x+0.26, y:3.34, w:3.33, size:11, color:TEXT, ls:1.18 });
  });

  plain(s, "Scope strategy — two regimes, each answering what the other cannot",
    { x:M, y:4.64, w:CW, size:14, bold:true, color:TEXT });
  const R = [
    ["Variant regime", "local  ·  dense", "288 variants on a single 1SU6 backbone, 15 positions × 19 substitutions. Side chains only, no backbone relaxation. This regime anchors the docking score to measured Km (ρ = 0.90)."],
    ["Homolog regime", "global  ·  sparse", "212 cluster representatives spanning 28.5–74.8% identity, each with its own predicted and assembled structure. This regime tests whether the score survives a change of backbone."]
  ];
  R.forEach((r,i) => {
    const x = M + i*6.15;
    card(s, x, 4.98, 5.85, 1.32, TINT2, LINE);
    plain(s, r[0], { x:x+0.28, y:5.14, w:3.2, size:13.5, bold:true, color:DEEP });
    plain(s, r[1], { x:x+3.40, y:5.18, w:2.2, size:11, bold:true, color:ACCENT, align:"right" });
    plain(s, r[2], { x:x+0.28, y:5.50, w:5.3, size:11, color:TEXT, ls:1.18 });
  });
  plain(s, "OUT OF SCOPE   absolute binding affinity  ·  Fe-S cofactor electrostatics  ·  experimental validation beyond the 21 existing Km values",
    { x:M, y:6.48, w:CW, size:10.5, bold:true, color:GRAY });
  s.addNotes("The two regimes are complementary: variants tie the score to experiment, homologs test generalization. Neither alone is sufficient.");
  foot(s);
}

/* ── 4. Label bottleneck ── */
{
  const s = newSlide();
  head(s, "STRATEGY", "The bottleneck is labels, not models");
  bullets(s, [
    [{text:"Only 21 Km values exist. ", bold:true}, {text:"That is too few to supervise a protein language model directly."}],
    [{text:"Docking scores serve as surrogate labels. ", bold:true}, {text:"Across 288 variants they correlate with measured Km at ρ = 0.89."}],
    [{text:"Strategy — a docking-distilled pLM. ", bold:true}, {text:"Dock into the EV binding site to generate binding-affinity predictions at scale, distill them into a sequence model, and at inference time predict without any structure."}]
  ], { x:M, y:1.95, w:6.3, h:2.6, size:14, gap:14 });
  card(s, M, 4.78, 6.3, 1.5, TINT2, LINE);
  plain(s, "Why distill at all", { x:M+0.3, y:4.98, w:5.7, size:12.5, bold:true, color:ACCENT });
  plain(s, "Docking needs a structure and costs minutes per sequence. A distilled model ranks from sequence alone in milliseconds, which is what makes the full 1,827-sequence family reachable.",
    { x:M+0.3, y:5.32, w:5.7, size:12, color:TEXT, ls:1.18 });

  const bx=7.55, bw=5.1;
  const steps = [
    ["Measured Km", "21 values", ACCENT],
    ["Docking score generation", "288 variants + 321 homologs", DEEP],
    ["ESM embedding + regression head", "SeekRank", DEEP],
    ["Sequence-only ranking", "extends to all 1,827", MID]
  ];
  steps.forEach((st,i) => {
    const y = 1.95 + i*1.13;
    card(s, bx, y, bw, 0.82, i===0?"FBEEE2":TINT, i===0?"E8C9A8":LINE);
    plain(s, st[0], { x:bx+0.25, y:y+0.12, w:bw-0.5, size:13, bold:true, color:st[2] });
    plain(s, st[1], { x:bx+0.25, y:y+0.45, w:bw-0.5, size:11, color:GRAY });
    if (i<3) s.addShape(pres.ShapeType.downArrow, { x:bx+bw/2-0.11, y:y+0.87, w:0.22, h:0.22,
      fill:{ color: MID }, line:{ type:"none" } });
  });
  s.addNotes("Docking is the label generator; the pLM is the predictor that ships.");
  foot(s);
}

/* ── 5. Pipeline ── */
{
  const s = newSlide();
  head(s, "OVERVIEW", "Pipeline");
  const rows = [
    ["Sequences", [["IPR010047\n1,827"],["70% clustering\n219"],["dockability gate\n1,722"],["medoid representatives\n212"],["157 train / 55 eval\n+ 109 D2 expansion"]]],
    ["Structures", [["AFDB\n298"],["ColabFold\n23"],["architecture-matched\nSF4→1SU6 · SF2→6OND"],["dimer assembly\n321 / 321"]]],
  ];
  let y = 1.85;
  rows.forEach(([lab, chips]) => {
    chip(s, M, y, 1.42, 0.95, lab, INK, WHITE, 12.5);
    const n = chips.length, ax = 0.32;
    const cwid = (CW - 1.42 - 0.28 - (n-1)*ax) / n;
    chips.forEach((c,i) => {
      const x = M + 1.42 + 0.28 + i*(cwid+ax);
      chip(s, x, y, cwid, 0.95, c[0], TINT, DEEP, 11);
      if (i<n-1) arrow(s, x+cwid+0.01, y+0.33);
    });
    y += 1.18;
  });
  const wide = [
    ["Docking", "residue-41 CA of each structure  ·  15.000 Å half-width (npts 80 × 0.375)  ·  AutoDock-GPU, 50 runs  ·  primary metric med"],
    ["Training", "SeekRank regression head  ·  experiments A sequential-by-identity / B all-at-once / C random order / D1 data volume / D2 within-cluster"]
  ];
  wide.forEach(([lab, txt]) => {
    chip(s, M, y, 1.42, 0.95, lab, INK, WHITE, 12.5);
    chip(s, M+1.70, y, CW-1.70, 0.95, txt, TINT, DEEP, 11.5);
    y += 1.18;
  });
  plain(s, "Sequence selection, structure retrieval and dimer assembly are settled. Docking is the next step.",
    { x:M, y:6.55, w:CW, size:11.5, italic:true, color:GRAY });
  s.addNotes("Everything above the docking row was frozen today.");
  foot(s);
}

/* ── 6. Dimer ── */
{
  const s = newSlide();
  head(s, "DECISION ①  RECEPTOR", "The receptor must be a dimer — the monomer loses predictive power");
  table(s, [
    hrow(["Metric (n = 21)","Dimer","Monomer","Δ"]),
    drow(["ρ with measured Km, med_lig","0.903","0.531","−0.372"], [null,{bold:true,color:DEEP},{bold:true,color:ACCENT},null]),
    drow(["ρ with measured Km, med_ds","0.856","0.512","−0.344"]),
    drow(["ρ with measured Km, best_ds","0.876","0.558","−0.319"]),
    drow(["ρ with dimer ranking, n = 288","—","0.702","—"])
  ], { x:M, y:1.95, w:6.6, colW:[2.9,1.2,1.25,1.25], rowH:0.42, size:11.5 });
  plain(s, "288 truncated-monomer controls · identical center, half-width, nrun and seed",
    { x:M, y:4.15, w:6.6, size:11, italic:true, color:GRAY });
  bullets(s, [
    [{text:"This is not a convergence failure. ", bold:true}, {text:"Ligand-swap correlation in the monomer runs is 0.995, higher than the dimer's. The pocket simply became simpler; the search converged normally."}],
    [{text:"It is the consequence of changing the structure. ", bold:true}, {text:"The SU2 contribution that was removed was not random noise but signal that matched experiment."}],
    [{text:"Consequence — multimer prediction is a prerequisite. ", bold:true, color:ACCENT}, {text:"The plan to build homologs with AF2 monomer was abandoned.", color:ACCENT}]
  ], { x:M, y:4.58, w:6.6, h:2.3, size:12, gap:10 });

  const bx=7.65;
  card(s, bx, 1.95, 5.0, 2.55, TINT2, LINE);
  plain(s, "Correlation with measured Km  (med_lig)", { x:bx+0.3, y:2.15, w:4.4, size:12, bold:true, color:GRAY });
  stat(s, bx+0.25, 2.60, 2.0, "0.903", "dimer", DEEP);
  s.addText("→", { x:bx+2.30, y:2.68, w:0.45, h:0.5, fontFace:F, fontSize:22, bold:true,
    color:GRAY, align:"center", valign:"middle", isTextBox:true, margin:0 });
  stat(s, bx+2.75, 2.60, 2.0, "0.531", "monomer", ACCENT);
  plain(s, "41% of the predictive power lost", { x:bx+0.3, y:4.02, w:4.4, size:12.5, bold:true, color:ACCENT, align:"center" });

  card(s, bx, 4.72, 5.0, 1.85, TINT, LINE);
  plain(s, "Why the dimer is required", { x:bx+0.3, y:4.92, w:4.4, size:12.5, bold:true, color:DEEP });
  plain(s, "The EV binding site straddles the C2 interface of the two protomers. Truncating to a monomer deletes half of that site.",
    { x:bx+0.3, y:5.28, w:4.4, size:12, color:TEXT, ls:1.18 });
  s.addNotes("The loss in Km predictive power is far larger than the loss in rank correlation. That is the argument.");
  foot(s);
}

/* ── 7. Interface ── */
{
  const s = newSlide();
  head(s, "EVIDENCE", "The EV binding site straddles the C2 interface");
  bullets(s, [
    [{text:"Atoms inside the box split across both protomers. ", bold:true}, {text:"SU1 contributes 436 (54.6%), SU2 contributes 362 (45.4%)."}],
    [{text:"Both protomers face each other with the same 34–67 loop. ", bold:true}, {text:"Within 12 Å of the center, the pocket wall is 28 residues from SU1 and 18 from SU2."}],
    [{text:"The box encloses the symmetry axis. ", bold:true}, {text:"The grid center sits 4.76 Å from the C2 axis, and the two F41 CE2 atoms form a symmetric pair 9.53 Å apart."}],
    [{text:"The C-cluster lies outside the box. ", bold:true, color:ACCENT}, {text:"The CO oxidation site is 26.4 Å away. The box sits near the D-cluster, on the electron transfer path.", color:ACCENT}]
  ], { x:M, y:1.95, w:6.2, h:3.4, size:12.5, gap:13 });
  card(s, M, 5.45, 6.2, 1.15, TINT2, LINE);
  plain(s, "Both protomers share chain A with identical residue numbering, so residue-based checks cannot tell them apart. The dimer was confirmed from atom count (11,370) and a median duplicate-pair distance of 51.61 Å.",
    { x:M+0.28, y:5.63, w:5.64, size:11, color:TEXT, ls:1.18 });

  const cx=9.95, cy=3.35;
  card(s, 7.15, 1.95, 5.5, 4.65, TINT2, LINE);
  plain(s, "15 Å half-width box", { x:cx-1.25, y:2.05, w:2.5, size:11, bold:true, color:ACCENT, align:"center" });
  s.addShape(pres.ShapeType.line, { x:cx, y:2.32, w:0, h:2.78, line:{ color:GRAY, width:1.25, dashType:"dash" } });
  s.addShape(pres.ShapeType.ellipse, { x:cx-2.55, y:cy-1.25, w:2.5, h:2.5,
    fill:{ color:"D5DAF2" }, line:{ color:MID, width:1.25 } });
  s.addShape(pres.ShapeType.ellipse, { x:cx+0.05, y:cy-1.25, w:2.5, h:2.5,
    fill:{ color:"E4E7F7" }, line:{ color:MID, width:1.25 } });
  plain(s, "SU1", { x:cx-2.30, y:cy-0.18, w:1.2, size:14, bold:true, color:DEEP, align:"center" });
  plain(s, "SU2", { x:cx+1.10, y:cy-0.18, w:1.2, size:14, bold:true, color:DEEP, align:"center" });
  s.addShape(pres.ShapeType.roundRect, { x:cx-0.95, y:cy-0.95, w:1.9, h:1.9,
    fill:{ type:"none" }, line:{ color:ACCENT, width:2.25 }, rectRadius:0.04 });
  plain(s, "C2 axis", { x:cx+0.14, y:4.72, w:1.2, size:10.5, bold:true, color:GRAY });
  plain(s, "Atoms in box    SU1 54.6%  /  SU2 45.4%", { x:7.45, y:5.42, w:5.0, size:12, bold:true, color:DEEP, align:"center" });
  plain(s, "grid center → C2 axis  4.76 Å      ·      C-cluster  26.4 Å (outside)",
    { x:7.45, y:5.80, w:5.0, size:10.5, color:GRAY, align:"center" });
  plain(s, "schematic — to be replaced by a PyMOL render for external use",
    { x:7.45, y:6.16, w:5.0, size:9.5, italic:true, color:GRAY, align:"center" });
  s.addNotes("The numbers alone do not convey that this is an interface site, hence the schematic.");
  foot(s);
}

/* ── 8. CA switch ── */
{
  const s = newSlide();
  head(s, "DECISION ②  GRID CENTER", "Move the grid center to the residue-41 CA — it is defined in every lineage");
  table(s, [
    hrow(["Residue at position 41","Sequences","CE2"]),
    drow(["F / Y / W  (aromatic)","1,613","present"], [null,null,{color:DEEP,bold:true}]),
    drow(["I  (Ile lineage)","118","absent"], [null,{bold:true,color:ACCENT},{color:ACCENT,bold:true}]),
    drow(["M · T · H · L · S","19","absent"], [null,null,{color:ACCENT}]),
    drow(["not annotated","77","absent"], [null,null,{color:ACCENT}])
  ], { x:M, y:1.95, w:5.8, colW:[3.0,1.4,1.4], rowH:0.42, size:11.5 });
  card(s, M, 4.20, 5.8, 1.05, "FBEEE2", "E8C9A8");
  plain(s, "214 sequences (11.7%) have no CE2", { x:M+0.3, y:4.38, w:5.2, size:17, bold:true, color:ACCENT });
  plain(s, "His is aromatic but carries only CE1 / NE2, so it does not help", { x:M+0.3, y:4.78, w:5.2, size:11, color:GRAY });
  plain(s, "A CE2-based definition cannot cover the whole family.", { x:M, y:5.45, w:5.8, size:12.5, bold:true, color:TEXT });
  plain(s, "A non-aromatic residue at position 41 is not an annotation error — it is a real feature of these homologs.",
    { x:M, y:5.78, w:5.8, size:12, color:GRAY, ls:1.18 });

  const bx=7.05, bw=5.6;
  bullets(s, [
    [{text:"Position 41 itself is kept. ", bold:true}, {text:"Of the 15 mutated positions (40, 42, 43, 44, 46, 57–66), only 41 is excluded. It was left out because it is the key residue exchanging electrons with EV, so centering on it is mechanistic rather than convenient."}],
    [{text:"The definition follows from sequence alone. ", bold:true}, {text:"f41_pos = anchor_pos + 2 holds for 1,750 of 1,750 sequences. Position 41 is the third residue of the FeS-binding Cys motif."}]
  ], { x:bx, y:1.95, w:bw, h:2.5, size:12.5, gap:13 });
  card(s, bx, 4.45, bw, 1.62, TINT, LINE);
  const motif = [["C39","",INK,WHITE],["G40","mut",TINT2,DEEP],["F41","center",ACCENT,WHITE],
                 ["G42","mut",TINT2,DEEP],["E43","mut",TINT2,DEEP],["T44","mut",TINT2,DEEP],
                 ["x45","",TINT2,GRAY],["L46","mut",TINT2,DEEP],["C47","",INK,WHITE]];
  motif.forEach((m,i) => {
    const x = bx + 0.03 + i*0.62;
    chip(s, x, 4.68, 0.58, 0.44, m[0], m[2], m[3], 11);
    if (m[1]) plain(s, m[1], { x:x-0.02, y:5.16, w:0.62, size:9, color:GRAY, align:"center" });
  });
  plain(s, "Mutations 40 · 42 · 43 · 44 · 46 are this motif's non-Cys positions",
    { x:bx+0.25, y:5.58, w:bw-0.5, size:11, italic:true, color:GRAY, align:"center" });
  plain(s, "The Cys-pair midpoint is unusable: the two lineages space their Cys differently (8 vs 3 positions), so the geometry does not match.",
    { x:bx, y:6.22, w:bw, size:11, color:GRAY, ls:1.15 });
  s.addNotes("Why CE2 fails is the crux. Position 41 is retained on mechanistic grounds.");
  foot(s);
}

/* ── 9. Equivalence ── */
{
  const s = newSlide();
  head(s, "DECISION ②  VALIDATION", "The switch changes nothing — 4.661 Å of movement, ranking preserved");
  table(s, [
    hrow(["Metric","ρ vs CE2","ρ within position","ρ with measured Km"]),
    drow(["med_ds","0.968","0.928","0.856 → 0.889"], [null,null,null,{bold:true,color:DEEP}]),
    drow(["med_lig","0.971","0.929","0.903 → 0.896"]),
    drow(["best_ds","0.978","0.956","0.876 → 0.884"]),
    drow(["best_lig","0.970","0.954","0.882 → 0.892"])
  ], { x:M, y:1.95, w:6.5, colW:[1.55,1.4,1.75,1.8], rowH:0.40, size:11.5 });
  plain(s, "288 structures re-docked · half-width fixed at 15, center only replaced · baseline = CE2 + 15 Å",
    { x:M, y:4.05, w:6.5, size:11, italic:true, color:GRAY });
  plain(s, "Why moving 4.661 Å changes nothing", { x:M, y:4.50, w:6.5, size:13.5, bold:true, color:DEEP });
  bullets(s, [
    [{text:"The two boxes overlap by 76.8%. ", bold:true}, {text:"Per-axis overlap is 0.920 × 0.875 × 0.954."}],
    [{text:"Docking lands on the same site. ", bold:true}, {text:"The top two clusters differ by 2.86 and 2.83 Å — within grid resolution, so it is the same physical site."}],
    [{text:"What remains is grid discretization. ", bold:true}, {text:"The mean difference is +0.01 to +0.02, so there is no bias, only per-structure scatter."}]
  ], { x:M, y:4.88, w:6.5, h:1.9, size:12, gap:9 });

  const bx=7.55, bw=5.1;
  card(s, bx, 1.95, bw, 2.15, TINT2, LINE);
  plain(s, "Applicable sequences", { x:bx+0.3, y:2.15, w:bw-0.6, size:12, bold:true, color:GRAY });
  stat(s, bx+0.25, 2.55, 1.95, "1,613", "CE2-based", GRAY);
  s.addText("→", { x:bx+2.25, y:2.63, w:0.45, h:0.5, fontFace:F, fontSize:22, bold:true,
    color:ACCENT, align:"center", valign:"middle", isTextBox:true, margin:0 });
  stat(s, bx+2.70, 2.55, 1.95, "1,750", "CA-based", DEEP);
  card(s, bx, 4.30, bw, 1.25, "FBEEE2", "E8C9A8");
  plain(s, "med_ds actually improved", { x:bx+0.3, y:4.52, w:bw-0.6, size:13.5, bold:true, color:ACCENT });
  plain(s, "0.856 → 0.889, while med_lig held at 0.903 → 0.896.",
    { x:bx+0.3, y:4.88, w:bw-0.6, size:12, color:TEXT });
  card(s, bx, 5.72, bw, 0.9, TINT, LINE);
  plain(s, "Coverage increases and nothing is paid for it.",
    { x:bx+0.3, y:5.94, w:bw-0.6, size:12.5, bold:true, color:DEEP });
  s.addNotes("The point is that the switch is free. med_ds even rose.");
  foot(s);
}

/* ── 10. Box size ── */
{
  const s = newSlide();
  head(s, "DECISION ②  BOX SIZE", "15 Å is not negotiable — a larger box lets the ligand slide down a groove");
  table(s, [
    hrow(["Half-width (Å)","Poses within 15 Å","Per-variant minimum","ρ within position"]),
    drow(["15.000","0.999","0.82","0.928"], [{bold:true},{bold:true,color:DEEP},{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["16.125","0.301","0.00","0.833"], [null,null,{bold:true,color:ACCENT},null]),
    drow(["17.250","0.176","0.00","0.781"]),
    drow(["18.000","0.168","0.00","0.779"]),
    drow(["20.250","0.136","0.00","0.774"])
  ], { x:M, y:1.95, w:6.2, colW:[1.35,1.8,1.7,1.35], rowH:0.365, size:11 });
  card(s, M, 4.50, 6.2, 1.0, "FBEEE2", "E8C9A8");
  plain(s, "At 16.125 Å the per-variant minimum is already 0.00", { x:M+0.28, y:4.66, w:5.6, size:13, bold:true, color:ACCENT });
  plain(s, "Some variants send all 50 poses outside 15 Å — a non-uniform bias.",
    { x:M+0.28, y:5.03, w:5.6, size:11.5, color:TEXT });
  bullets(s, [
    [{text:"Not a cliff but a descent from 15.0 to 17.25, then a plateau. ", bold:true}, {text:"60% of the drop occurs in the first 1.1 Å."}],
    [{text:"15 Å cuts the groove off at its mouth, ", bold:true, color:DEEP}, {text:"and that stretch overlaps the mutated residues, which is where the signal comes from.", color:DEEP}]
  ], { x:M, y:5.70, w:6.2, h:1.3, size:11.5, gap:8 });

  const bx=7.35, bw=5.3;
  card(s, bx, 1.95, bw, 4.4, TINT2, LINE);
  plain(s, "Cluster centers migrate along one groove", { x:bx+0.3, y:2.15, w:bw-0.6, size:13, bold:true, color:DEEP });
  plain(s, "Direction (1, 0, 1)/√2 — y stays pinned at 3 while x and z grow",
    { x:bx+0.3, y:2.50, w:bw-0.6, size:11, color:GRAY });
  const ax0=bx+0.55, axw=4.2, ay=4.30;
  s.addShape(pres.ShapeType.rect, { x:ax0, y:ay-0.85, w:1.35, h:1.7, fill:{ color:"E2E6F6" }, line:{ type:"none" } });
  plain(s, "mutated residues", { x:ax0-0.05, y:ay-1.14, w:1.45, size:9, bold:true, color:DEEP, align:"center" });
  s.addShape(pres.ShapeType.line, { x:ax0, y:ay, w:axw, h:0, line:{ color:GRAY, width:1.25 } });
  s.addShape(pres.ShapeType.line, { x:ax0+1.35, y:ay-0.95, w:0, h:1.9, line:{ color:ACCENT, width:1.5, dashType:"dash" } });
  plain(s, "15 Å boundary", { x:ax0+0.95, y:ay+1.00, w:1.6, size:10, bold:true, color:ACCENT, align:"center" });
  const pts = [[0.28,"[3,3,6]","r15"],[1.02,"[6,3,6]","r15"],[2.55,"[12,3,12]","r16–18"],[3.75,"[15,3,15]","r20–22"]];
  pts.forEach(([dx,lab,rad],i) => {
    const x = ax0 + dx;
    s.addShape(pres.ShapeType.ellipse, { x:x-0.10, y:ay-0.10, w:0.20, h:0.20,
      fill:{ color: i<2?DEEP:ACCENT }, line:{ type:"none" } });
    plain(s, lab, { x:x-0.45, y:ay-0.52, w:0.90, size:10, bold:true, color:i<2?DEEP:ACCENT, align:"center" });
    plain(s, rad, { x:x-0.45, y:ay+0.16, w:0.90, size:9.5, color:GRAY, align:"center" });
  });
  plain(s, "Once the ligand slips past the mutated residues and into the groove, the score is measuring something other than what we intend.",
    { x:bx+0.3, y:5.58, w:bw-0.6, size:11.5, bold:true, color:TEXT, ls:1.18 });
  s.addNotes("15 A is a physically motivated boundary, not a compromise. Do not change it.");
  foot(s);
}

/* ── 11. Principles ── */
{
  const s = newSlide();
  head(s, "DECISION ③  DATASET", "The evaluation set was frozen before any label existed — four principles");
  const P = [
    ["1","Cluster before gating","Sequence partitioning is a property of the protein family, not of our pipeline. Applying technical constraints first would make cluster boundaries depend on those constraints and stop being reproducible."],
    ["2","One representative per cluster (CAP = 1)","Sequences in a cluster are at least 70% identical, so their labels are nearly the same. Adding more only lets that lineage cast extra votes in the loss without contributing new information."],
    ["3","Split by cluster, not by sequence","If members of one cluster land in both train and evaluation, the evaluation measures recall rather than generalization. With one representative per cluster, this leak is structurally impossible."],
    ["4","Fix the split before seeing labels","Partitioning after inspecting docking scores would justify any outcome. The 212 representatives and their train/eval assignment were fixed before docking began."]
  ];
  P.forEach((p,i) => {
    const y = 1.88 + i*1.13;
    badge(s, M+0.05, y+0.10, 0.46, p[0], i<2?DEEP:MID);
    plain(s, p[1], { x:M+0.70, y:y+0.00, w:11.2, size:14, bold:true, color:TEXT });
    plain(s, p[2], { x:M+0.70, y:y+0.36, w:11.2, size:11.5, color:GRAY, ls:1.16 });
  });
  card(s, M, 6.25, CW, 0.62, "FBEEE2", "E8C9A8");
  plain(s, "Evaluating on the training set would overstate performance — these four principles are the answer to that.",
    { x:M+0.3, y:6.41, w:CW-0.6, size:12.5, bold:true, color:ACCENT });
  s.addNotes("These principles govern the whole selection protocol.");
  foot(s);
}

/* ── 12. Funnel ── */
{
  const s = newSlide();
  head(s, "DECISION ③  SELECTION", "From 1,827 to 212 — gating applies only when picking representatives");
  const F2 = [
    ["IPR010047, all sequences", "1,827 sequences", 6.9, DEEP],
    ["mmseqs clustering  --min-seq-id 0.70", "219 clusters", 5.9, DEEP],
    ["dockability gate", "1,722 sequences", 4.9, MID],
    ["one medoid per cluster", "212 representatives", 3.9, MID],
    ["bin-stratified 25% split", "157 train / 55 eval", 3.0, ACCENT]
  ];
  const ccx = 4.35;
  F2.forEach((f,i) => {
    const y = 1.92 + i*0.92, w = f[2];
    card(s, ccx-w/2, y, w, 0.74, i===4?"FBEEE2":TINT, i===4?"E8C9A8":LINE);
    plain(s, f[0], { x:ccx-w/2+0.2, y:y+0.08, w:w-0.4, size:11, color:GRAY, align:"center" });
    plain(s, f[1], { x:ccx-w/2+0.2, y:y+0.36, w:w-0.4, size:13, bold:true, color:f[3], align:"center" });
    if (i<4) s.addShape(pres.ShapeType.downArrow, { x:ccx-0.10, y:y+0.76, w:0.20, h:0.15,
      fill:{ color:MID }, line:{ type:"none" } });
  });
  const bx=8.35, bw=4.3;
  plain(s, "Excluded at the gate", { x:bx, y:1.92, w:bw, size:13.5, bold:true, color:TEXT });
  table(s, [
    hrow(["Reason","Seqs","Clusters"]),
    drow(["no equivalent of residue 41","80","4"]),
    drow(["quality filter","25","3"]),
    drow(["total","105","7"], [{bold:true},{bold:true},{bold:true}])
  ], { x:bx, y:2.30, w:bw, colW:[2.3,1.0,1.0], rowH:0.36, size:10.5 });
  bullets(s, [
    [{text:"A missing residue 41 is biology, not failure. ", bold:true}, {text:"In 4.4% of the family the position does not exist, so no binding site can be defined."}],
    [{text:"AFDB membership is not a gate. ", bold:true}, {text:"It is a property of the database, not the protein. The 23 sequences needing prediction were all A0AB* / A0AC* — recent deposits that gating would have discarded."}]
  ], { x:bx, y:4.05, w:bw, h:1.7, size:11, gap:10 });
  card(s, bx, 5.82, bw, 0.98, TINT, LINE);
  plain(s, "One cluster has all 13 members missing residue 41 — possibly a clade that lost the F41 site, flagged for follow-up.",
    { x:bx+0.26, y:5.97, w:bw-0.52, size:10.5, bold:true, color:DEEP, ls:1.12 });
  s.addNotes("This is principle 1 in action: cluster first, gate only at representative selection.");
  foot(s);
}

/* ── 13. Bins ── */
{
  const s = newSlide();
  head(s, "DECISION ③  BINS", "Bin boundaries follow the data — the intent is unchanged");
  card(s, M, 1.95, 6.1, 1.85, TINT2, LINE);
  plain(s, "Identity to ChCODH2 across 212 clusters", { x:M+0.3, y:2.13, w:5.5, size:12, bold:true, color:GRAY });
  const dd = [["0.285","minimum"],["0.419","median"],["0.550","95th pct"],["0.748","maximum"]];
  dd.forEach((d,i) => stat(s, M+0.20+i*1.44, 2.55, 1.38, d[0], d[1], DEEP, 21));
  card(s, M, 4.00, 6.1, 1.35, "FBEEE2", "E8C9A8");
  plain(s, "Only 3–4 sequences exceed 90%, and 62.5–72.5% is completely empty.",
    { x:M+0.3, y:4.22, w:5.5, size:12.5, bold:true, color:ACCENT, ls:1.18 });
  plain(s, "This family has almost no close relatives of ChCODH2.",
    { x:M+0.3, y:4.88, w:5.5, size:11.5, color:TEXT });
  plain(s, "Under the original proposal the top four bins would each hold single digits and 1,780 sequences would pile into the last one.",
    { x:M, y:5.50, w:6.1, size:12, bold:true, color:TEXT, ls:1.18 });
  plain(s, "The intent — training in order of identity — is kept; only the boundaries moved to fit the data.",
    { x:M, y:6.12, w:6.1, size:12, bold:true, color:DEEP, ls:1.18 });

  const bx=7.35, bw=5.3;
  plain(s, "Original", { x:bx, y:1.95, w:2.4, size:12.5, bold:true, color:GRAY });
  plain(s, "Adopted", { x:bx+2.9, y:1.95, w:2.4, size:12.5, bold:true, color:DEEP });
  const pairs = [["90 – 100","≥ 60","6"],["80 – 90","50 – 60","14"],["70 – 80","40 – 50","107"],["60 – 70","30 – 40","75"],["< 60","< 30","10"]];
  pairs.forEach((p,i) => {
    const y = 2.35 + i*0.82;
    chip(s, bx, y, 2.4, 0.62, p[0], "E9EAEF", GRAY, 12);
    arrow(s, bx+2.48, y+0.16);
    chip(s, bx+2.9, y, 2.4, 0.62, p[1] + "      " + p[2] + " reps", TINT, DEEP, 12);
  });
  plain(s, "10% width · cluster median fident", { x:bx, y:6.50, w:bw, size:11, italic:true, color:GRAY });
  s.addNotes("Explains why we deviate from the proposed boundaries while keeping the intent.");
  foot(s);
}

/* ── 14. Selection result ── */
{
  const s = newSlide();
  head(s, "DECISION ③  RESULT", "212 representatives stand in for 98.8% of the family");
  table(s, [
    hrow(["Bin","Reps","Seqs covered","Train","Eval","SF2","Non-aromatic 41"]),
    drow(["≥ 60","6","31","4","2","0","0"]),
    drow(["50 – 60","14","96","10","4","0","0"]),
    drow(["40 – 50","107","619","80","27","13","23"]),
    drow(["30 – 40","75","897","56","19","17","16"]),
    drow(["< 30","10","163","7","3","3","4"]),
    drow(["Total","212","1,806","157","55","33","43"],
      [{bold:true,fill:{color:TINT}},{bold:true,color:DEEP,fill:{color:TINT}},{bold:true,fill:{color:TINT}},
       {bold:true,color:DEEP,fill:{color:TINT}},{bold:true,color:DEEP,fill:{color:TINT}},
       {bold:true,fill:{color:TINT}},{bold:true,fill:{color:TINT}}])
  ], { x:M, y:1.95, w:7.55, colW:[1.25,0.85,1.35,0.9,0.85,0.85,1.5], rowH:0.40, size:11 });
  bullets(s, [
    [{text:"Both SF2 and non-aromatic position 41 appear only below 50% identity. ", bold:true}, {text:"The high-identity bins are not only small, they carry little information — and that region is already covered densely by the 288 variants."}],
    [{text:"Zero shared clusters, zero duplicate sequences. ", bold:true}, {text:"With one representative per cluster, leakage cannot occur."}]
  ], { x:M, y:5.15, w:7.55, h:1.6, size:11.5, gap:9 });

  const bx=8.55, bw=4.1;
  card(s, bx, 1.95, bw, 1.65, TINT2, LINE);
  plain(s, "Family coverage", { x:bx+0.28, y:2.12, w:bw-0.56, size:11.5, bold:true, color:GRAY, align:"center" });
  plain(s, "98.8%", { x:bx+0.28, y:2.45, w:bw-0.56, size:38, bold:true, color:DEEP, align:"center" });
  plain(s, "1,806 of 1,827", { x:bx+0.28, y:3.15, w:bw-0.56, size:11, color:GRAY, align:"center" });
  card(s, bx, 3.78, bw, 1.15, TINT, LINE);
  plain(s, "Primary metric", { x:bx+0.28, y:3.94, w:bw-0.56, size:11, bold:true, color:GRAY });
  plain(s, "Spearman over all 55", { x:bx+0.28, y:4.24, w:bw-0.56, size:13.5, bold:true, color:DEEP });
  card(s, bx, 5.08, bw, 1.67, "FBEEE2", "E8C9A8");
  plain(s, "Per-bin cannot be primary", { x:bx+0.28, y:5.26, w:bw-0.56, size:12, bold:true, color:ACCENT });
  plain(s, "The ≥ 60 bin has 2 evaluation sequences and < 30 has 3. Rank correlations at that size carry intervals far too wide. Only 40–50 (27) and 30–40 (19) are reported separately; the tails are pooled.",
    { x:bx+0.28, y:5.58, w:bw-0.56, size:10, color:TEXT, ls:1.16 });
  s.addNotes("Coverage and zero leakage are the headline. The per-bin caveat is stated up front.");
  foot(s);
}

/* ── 15. Experiments ── */
{
  const s = newSlide();
  head(s, "DECISION ③  EXPERIMENTS", "Same training data, same evaluation set — only order, mode and volume differ");
  table(s, [
    hrow(["","Mode","Order","Volume","Question answered","Reps"]),
    drow(["A","sequential (warm-start)","by identity","4 → 14 → 94 → 150 → 157","does identity-ordered training help","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["B","all at once","—","157","upper baseline, endpoint for A","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["C","sequential","random","same sizes as A","isolates the contribution of order","5"],
      [{bold:true,color:ACCENT},null,null,null,{align:"left"},null]),
    drow(["D1","all at once","random (stratified)","20 · 40 · … · 157","where data volume saturates","5"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["D2","all at once","—","157 + 109 = 266","was one representative enough","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null])
  ], { x:M, y:1.95, w:CW, colW:[0.7,2.35,1.85,2.55,3.45,1.1], rowH:0.44, size:11 });
  const nb = [
    ["Why C is indispensable", "Without C, the gap between A and B cannot be attributed to sequencing versus identity ordering.", ACCENT],
    ["A requires warm-start training", "Cumulative retraining would make A's final point identical to B, erasing the contrast. Ridge, SVR and RF have no warm-start notion.", DEEP],
    ["Embeddings computed once", "72 training runs but only 321 unique sequences. ESM embeddings are cached per sequence.", DEEP]
  ];
  nb.forEach((n,i) => {
    const x = M + i*4.12;
    card(s, x, 4.85, 3.85, 1.55, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    plain(s, n[0], { x:x+0.25, y:5.02, w:3.35, size:12, bold:true, color:n[2] });
    plain(s, n[1], { x:x+0.25, y:5.36, w:3.35, size:11, color:TEXT, ls:1.16 });
  });
  plain(s, "OPEN — the structure of the SeekRank regression head must be confirmed; whether experiment A is feasible depends on it.",
    { x:M, y:6.58, w:CW, size:11.5, bold:true, color:ACCENT });
  s.addNotes("C's purpose and A's warm-start requirement are the likely questions.");
  foot(s);
}

/* ── 16. Assembly ── */
{
  const s = newSlide();
  head(s, "DECISION ④  STRUCTURES", "All 321 dimers assembled — references split by architecture");
  plain(s, "Structure sources", { x:M, y:1.95, w:5.9, size:13.5, bold:true, color:TEXT });
  table(s, [
    hrow(["Source","Count","Note"]),
    drow(["AFDB","298","downloaded"]),
    drow(["ColabFold","23","A0AB* / A0AC* — not in AFDB"]),
    drow(["Total","321","core 212 + D2 expansion 109"], [{bold:true},{bold:true,color:DEEP},{bold:true}])
  ], { x:M, y:2.33, w:5.9, colW:[1.5,0.95,3.45], rowH:0.38, size:11 });
  card(s, M, 4.05, 5.9, 1.2, "FBEEE2", "E8C9A8");
  plain(s, "Two prediction methods are mixed", { x:M+0.28, y:4.22, w:5.34, size:12, bold:true, color:ACCENT });
  plain(s, "Both are AlphaFold2, but the MSA pipelines differ. Recorded in a struct_source column and checked against score after docking.",
    { x:M+0.28, y:4.55, w:5.34, size:11, color:TEXT, ls:1.16 });
  plain(s, "Assembly quality", { x:M, y:5.42, w:5.9, size:13.5, bold:true, color:TEXT });
  const q = [["321 / 321","succeeded"],["0.80 – 1.65","RMSD by bin"],["0 – 1.22","steric clashes"],["94.6 – 96.6","pLDDT"]];
  q.forEach((x,i) => {
    const cx2 = M + i*1.50;
    card(s, cx2, 5.78, 1.42, 0.82, TINT2, LINE);
    plain(s, x[0], { x:cx2+0.05, y:5.90, w:1.32, size:12.5, bold:true, color:i===0?DEEP:TEXT, align:"center" });
    plain(s, x[1], { x:cx2+0.05, y:6.24, w:1.32, size:9, color:GRAY, align:"center" });
  });

  const bx=7.05, bw=5.6;
  plain(s, "Reference chosen by D-cluster architecture", { x:bx, y:1.95, w:bw, size:13.5, bold:true, color:TEXT });
  const refs = [["SF4 · SF4v · SFx","ChCODH2  (1SU6)"],["SF2","DvCODH  (6OND)"]];
  refs.forEach((r,i) => {
    const y = 2.38 + i*0.95;
    chip(s, bx, y, 2.5, 0.75, r[0], TINT, DEEP, 12);
    arrow(s, bx+2.58, y+0.22);
    chip(s, bx+3.0, y, 2.6, 0.75, r[1], INK, WHITE, 12);
  });
  card(s, bx, 4.32, bw, 0.95, TINT2, LINE);
  plain(s, "Both global and local superposition modes were tried; the one minimising (clashes, RMSD) was kept.",
    { x:bx+0.28, y:4.55, w:bw-0.56, size:11.5, color:TEXT, ls:1.18 });
  card(s, bx, 5.42, bw, 1.3, TINT, LINE);
  plain(s, "Steric clashes stay flat across bins", { x:bx+0.28, y:5.62, w:bw-0.56, size:12.5, bold:true, color:DEEP });
  plain(s, "RMSD rises with decreasing identity but clashes do not — the physical plausibility of the assembly does not track identity.",
    { x:bx+0.28, y:5.96, w:bw-0.56, size:11, color:TEXT, ls:1.16 });
  s.addNotes("The point is that SF2 homologs must not be assembled against 1SU6; 6OND is used instead.");
  foot(s);
}

/* ── 17. Assembly validation ── */
{
  const s = newSlide();
  head(s, "DECISION ④  VALIDATION", "The assembly procedure does not distort docking scores");
  card(s, M, 2.05, 5.75, 2.1, TINT2, LINE);
  plain(s, "Wildtype re-dock — pipeline consistency check", { x:M+0.3, y:2.28, w:5.15, size:12, bold:true, color:GRAY });
  plain(s, "−4.380", { x:M+0.3, y:2.62, w:2.6, size:36, bold:true, color:DEEP });
  plain(s, "Δ vs reference  0.000", { x:M+2.95, y:2.95, w:2.5, size:13, bold:true, color:DEEP });
  plain(s, "kcal/mol", { x:M+0.3, y:3.42, w:2.6, size:10.5, color:GRAY });

  card(s, 6.85, 2.05, 5.8, 2.1, "FBEEE2", "E8C9A8");
  plain(s, "Assembled ChCODH2 dimer", { x:7.15, y:2.28, w:5.2, size:12, bold:true, color:GRAY });
  plain(s, "−4.30", { x:7.15, y:2.62, w:2.6, size:36, bold:true, color:ACCENT });
  plain(s, "Δ  +0.080", { x:9.95, y:2.95, w:2.4, size:13, bold:true, color:ACCENT });
  plain(s, "kcal/mol", { x:7.15, y:3.42, w:2.6, size:10.5, color:GRAY });

  card(s, M, 4.40, CW, 1.15, TINT, LINE);
  plain(s, "The 0.080 assembly effect is smaller than the crystal's own 50-run range of 0.11.",
    { x:M+0.35, y:4.62, w:CW-0.7, size:16.5, bold:true, color:DEEP });
  plain(s, "In other words, the difference introduced by assembly is buried in the run-to-run variation of docking itself.",
    { x:M+0.35, y:5.02, w:CW-0.7, size:12, color:TEXT });
  bullets(s, [
    [{text:"The reproduction check also passed. ", bold:true}, {text:"D63A_A559W was re-docked with the CE2 center at npts 80 and every sorted energy matched the stored result to atol = 1e-3 before the batch was started."}],
    [{text:"The reference score is Wildtype ad_med_ds = −4.38 with IQR 0.09. ", bold:true}, {text:"All homolog scores are read against it."}]
  ], { x:M, y:5.70, w:CW, h:1.3, size:11.5, gap:8 });
  s.addNotes("Without this check the homolog scores would not be trustworthy.");
  foot(s);
}

/* ── 18. Next steps ── */
{
  const s = newSlide();
  head(s, "IN PROGRESS ⑤", "Next: dock 321 structures and measure the noise floor");
  plain(s, "Why a noise floor is needed — the noise sources changed", { x:M, y:1.92, w:6.2, size:13.5, bold:true, color:TEXT });
  table(s, [
    hrow(["","288 variants","321 homologs"]),
    drow(["Backbone","one shared 1SU6","predicted separately"]),
    drow(["Assembly","not required","assembled separately"]),
    drow(["Noise","run-to-run only (IQR 0.09)","run + prediction + assembly"],
      [{bold:true},null,{bold:true,color:ACCENT}])
  ], { x:M, y:2.32, w:6.2, colW:[1.35,2.35,2.5], rowH:0.40, size:11 });
  card(s, M, 4.10, 6.2, 1.0, "FBEEE2", "E8C9A8");
  plain(s, "The SNR of the 288 set cannot be carried over.", { x:M+0.28, y:4.28, w:5.64, size:12.5, bold:true, color:ACCENT });
  plain(s, "No difference smaller than the instrument's own noise can be claimed.", { x:M+0.28, y:4.62, w:5.64, size:11.5, color:TEXT });
  plain(s, "How it is measured — each of the 109 expansion sequences is the second member of a training cluster, paired with its representative.",
    { x:M, y:5.22, w:6.2, size:11.5, bold:true, color:TEXT, ls:1.18 });
  card(s, M, 5.78, 6.2, 0.95, TINT, LINE);
  plain(s, "within(c) = | med(representative) − med(second member) |", { x:M+0.28, y:5.92, w:5.64, size:12, bold:true, color:DEEP });
  plain(s, "SNR = Var(between clusters) / mean(within²)", { x:M+0.28, y:6.28, w:5.64, size:12, bold:true, color:DEEP });

  const bx=7.35, bw=5.3;
  plain(s, "Decision thresholds fixed in advance", { x:bx, y:1.92, w:bw, size:13.5, bold:true, color:TEXT });
  table(s, [
    hrow(["Median within","SNR","Verdict"]),
    drow(["~ 0.05","40","use labels as they are"]),
    drow(["~ 0.30","7","comparable to 288. proceed"], [null,null,{bold:true,color:DEEP}]),
    drow(["~ 0.80","2.5","borderline. average repeat assemblies"]),
    drow(["~ 1.50","1.3","labels buried in noise. redesign"], [null,null,{color:ACCENT}])
  ], { x:bx, y:2.32, w:bw, colW:[1.5,0.8,3.0], rowH:0.365, size:10.5 });
  plain(s, "The D2 expansion is docked in the same batch, so the added cost is zero.",
    { x:bx, y:4.20, w:bw, size:11.5, bold:true, color:DEEP });
  plain(s, "Schedule", { x:bx, y:4.70, w:bw, size:13.5, bold:true, color:TEXT });
  const plan = ["Dock 321  (core 212 + expansion 109 together)","Measure the noise floor → SNR verdict",
    "Covariate regression  struct_source · mode · clashes · pLDDT","Attach labels  FASTA header fident → med score",
    "Train  A / B / C / D1 / D2"];
  plan.forEach((p,i) => {
    const y = 5.08 + i*0.36;
    badge(s, bx, y, 0.26, String(i+1), i===0?ACCENT:MID);
    plain(s, p, { x:bx+0.38, y:y+0.01, w:bw-0.38, size:11, color:TEXT });
  });
  s.addNotes("Without this number the claim of predicting homolog ranking cannot be defended to a reviewer.");
  foot(s);
}

/* ── Backup divider ── */
{
  const s = newSlide(true);
  s.addText("BACKUP", { x:M, y:3.05, w:CW, h:0.35, fontFace:F, fontSize:13, bold:true,
    color:ACCENT, charSpacing:2.2, isTextBox:true, margin:0 });
  s.addText("Backup slides", { x:M, y:3.45, w:CW, h:0.7, fontFace:F, fontSize:34, bold:true,
    color:WHITE, isTextBox:true, margin:0 });
  s.addText("B1 orientation dependence  ·  B2 choice of primary metric  ·  B3 simple descriptors  ·  B4 fixed constants  ·  B5 limitations  ·  B6 interpretive scope",
    { x:M, y:4.30, w:CW-0.6, h:0.4, fontFace:F, fontSize:12.5, color:INKSUB, isTextBox:true, margin:0 });
  foot(s, true);
}

/* ── B1 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B1", "At 15 Å the ranking is independent of orientation");
  table(s, [
    hrow(["Half-width","Rotation control","Seed control","Gap","Verdict"]),
    drow(["15.000","0.905","0.927","−0.022","orientation-free"],
      [{bold:true},null,null,{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["18.000","0.897","0.979","−0.082","orientation-dependent"], [null,null,null,{color:ACCENT},{color:ACCENT}])
  ], { x:M, y:1.95, w:7.3, colW:[1.3,1.5,1.4,1.1,2.0], rowH:0.42, size:11 });
  plain(s, "med, within position · receptor rotated about the grid center (zyx 37/23/51°) and re-docked",
    { x:M, y:3.42, w:7.3, size:11, italic:true, color:GRAY });
  bullets(s, [
    [{text:"A seed control is mandatory. ", bold:true}, {text:"Rotation changes the search trajectory itself, so without a yardstick for “a different search under equivalent conditions” the result cannot be interpreted."}],
    [{text:"The 0.520 maximum shift is not negligible. ", bold:true, color:ACCENT}, {text:"That is 1.7σ against the 0.298 standard deviation across 288. The mean difference is 0.005, so there is no bias — a few structures simply move a lot. Ranking survives; absolute values do not.", color:ACCENT}]
  ], { x:M, y:3.88, w:7.3, h:2.2, size:12, gap:10 });
  const bx=8.25, bw=4.4;
  card(s, bx, 1.95, bw, 4.4, TINT2, LINE);
  plain(s, "Paradox", { x:bx+0.3, y:2.15, w:bw-0.6, size:12, bold:true, color:ACCENT });
  plain(s, "High reproducibility is not correctness", { x:bx+0.3, y:2.48, w:bw-0.6, size:14, bold:true, color:TEXT, ls:1.13 });
  plain(s, "In the seed control, 18 Å (0.979) reproduces better than 15 Å (0.927). At 18 Å, 79% of poses are stuck in a single cluster, so changing the seed barely moves the median. At 15 Å several clusters compete at nearly equal energy, so it wobbles.",
    { x:bx+0.3, y:3.05, w:bw-0.6, size:11, color:TEXT, ls:1.18 });
  card(s, bx+0.3, 4.90, bw-0.6, 1.25, "FBEEE2", "E8C9A8");
  plain(s, "Measuring the wrong thing stably is worse than measuring the right thing with some wobble.",
    { x:bx+0.52, y:5.12, w:bw-1.04, size:11.5, bold:true, color:ACCENT, ls:1.18 });
  foot(s);
}

/* ── B2 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B2", "med is the primary metric because it survives box perturbation");
  table(s, [
    hrow(["Metric","ρ within position  (15 → 18 Å, CE2 center)"]),
    drow(["med  (median of 50 runs)","0.918  /  0.943"], [{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["best  (single extreme value)","0.674  /  0.691"], [null,{color:ACCENT}])
  ], { x:M, y:2.05, w:7.0, colW:[3.0,4.0], rowH:0.45, size:12 });
  bullets(s, [
    [{text:"best jumps the moment a new minimum appears. ", bold:true}, {text:"med is a robust statistic over the 50-run ensemble, so it holds its ranking when the landscape shifts slightly."}],
    [{text:"Choosing a metric on n = 27 would be overfitting. ", bold:true}, {text:"Instead of selecting by correlation with measured Km, the question was answered with an independent criterion: stability under box perturbation."}],
    [{text:"med = −median(50 runs), reported with IQR. ", bold:true, color:DEEP}, {text:"Spread is quoted alongside the value everywhere.", color:DEEP}]
  ], { x:M, y:3.85, w:7.0, h:2.4, size:12.5, gap:12 });
  const bx=8.05, bw=4.6;
  card(s, bx, 2.05, bw, 2.35, TINT2, LINE);
  plain(s, "50-run energy range", { x:bx+0.3, y:2.28, w:bw-0.6, size:12, bold:true, color:GRAY, align:"center" });
  plain(s, "0.22", { x:bx+0.3, y:2.62, w:bw-0.6, size:40, bold:true, color:DEEP, align:"center" });
  plain(s, "kcal/mol  at 15 Å", { x:bx+0.3, y:3.32, w:bw-0.6, size:11.5, color:GRAY, align:"center" });
  plain(s, "The landscape is flat — no basis for trusting a single extreme.",
    { x:bx+0.3, y:3.70, w:bw-0.6, size:11, italic:true, color:GRAY, align:"center", ls:1.13 });
  card(s, bx, 4.60, bw, 1.65, TINT, LINE);
  plain(s, "Run flags (required for reproducibility)", { x:bx+0.3, y:4.80, w:bw-0.6, size:12, bold:true, color:DEEP });
  plain(s, "--nrun 50 --seed 1\n--heuristics 0 --autostop 0", { x:bx+0.3, y:5.15, w:bw-0.6, size:12, color:TEXT, ls:1.22 });
  plain(s, "Drop the last two and the run count becomes variable, collapsing reproducibility.", { x:bx+0.3, y:5.85, w:bw-0.6, size:10, color:GRAY });
  foot(s);
}

/* ── B3 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B3", "Docking leads the simple descriptors, but at n = 21 the margin is not decisive");
  table(s, [
    hrow(["Descriptor","ρ with measured Km","Applies to homologs"]),
    drow(["ad_med_lig  (docking)","0.903","yes"], [{bold:true,color:DEEP},{bold:true,color:DEEP},null]),
    drow(["ad_best_ds  (docking)","0.876","yes"]),
    drow(["ad_med_ds  (docking)","0.856","yes"]),
    drow(["dchg  (mutation charge change)","0.815","no"], [{bold:true,color:ACCENT},{bold:true,color:ACCENT},{color:ACCENT}]),
    drow(["cdock  (C-Docker)","0.719","yes"]),
    drow(["chg","0.594","no"]),
    drow(["ev_dist","0.465","yes"]),
    drow(["dvol / vol","0.458 / 0.439","no"])
  ], { x:M, y:1.95, w:6.8, colW:[3.1,2.0,1.7], rowH:0.365, size:11 });
  const bx=7.75, bw=4.9;
  bullets(s, [
    [{text:"dchg cannot be applied to homologs. ", bold:true}, {text:"It is defined only for variants, so extending to the whole family requires docking."}],
    [{text:"Whether docking carries information beyond charge is unverified. ", bold:true, color:ACCENT}, {text:"Its incremental explanatory power after controlling for dchg should be measured across all 288 — the computation is straightforward but has not been run.", color:ACCENT}]
  ], { x:bx, y:1.95, w:bw, h:2.6, size:12, gap:12 });
  card(s, bx, 4.55, bw, 1.85, TINT2, LINE);
  plain(s, "8X9F geometry (for reference)", { x:bx+0.3, y:4.75, w:bw-0.6, size:12, bold:true, color:DEEP });
  plain(s, "The EV centroid of chain A sits 12.65 Å from the grid center and chain A-2 sits at 16.89 Å. The two EV copies are 28.03 Å apart, placing the grid center between them. CA RMSD on superposition to 1SU6 protomer 1 is 0.28 Å.",
    { x:bx+0.3, y:5.10, w:bw-0.6, size:11, color:TEXT, ls:1.18 });
  foot(s);
}

/* ── B4 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B4", "Fixed constants — the interface between sessions");
  table(s, [
    hrow(["Item","Value","Basis"]),
    drow(["Cluster threshold","0.70","chosen from a 0.30–0.90 sweep; reproduces the original 219 clusters"]),
    drow(["Box center","residue-41 CA of each structure","CE2 is absent in 11.7% of homologs (Ile and others)"], [null,{bold:true,color:DEEP},null]),
    drow(["Box size","npts 80 × 0.375 = 15.000 Å half-width","above 18 Å the within-position correlation falls 0.93 → 0.78"], [null,{bold:true,color:DEEP},null]),
    drow(["Ligand","ev_ds.pdbqt","16 atoms / charge sum 1.998 / TORSDOF 3"]),
    drow(["Run flags","--nrun 50 --seed 1 --heuristics 0 --autostop 0","without the last two the run count varies, breaking reproducibility"]),
    drow(["Primary metric","med (median of 50 runs)","more robust than best under box perturbation"]),
    drow(["Reference score","Wildtype ad_med_ds = −4.38, IQR 0.09","engines_288_CA_r15.csv"]),
    drow(["Seed","20260916","reproduces selection, split and experiment sets"])
  ], { x:M, y:1.95, w:CW, colW:[1.85,4.15,6.0], rowH:0.42, size:10.5 });
  card(s, M, 5.88, CW, 1.12, "FBEEE2", "E8C9A8");
  plain(s, "Regenerating the ligand with obabel or prepare_ligand4 drops its charge sum from +2 to 0.000.",
    { x:M+0.32, y:6.04, w:CW-0.64, size:12, bold:true, color:ACCENT });
  plain(s, "The signal disappears entirely while the output numbers still look normal, so it is hard to notice. GPF generators are likewise forbidden — a parameter_file None line makes autogrid4 segfault.",
    { x:M+0.32, y:6.38, w:CW-0.64, size:11, color:TEXT, ls:1.15 });
  foot(s);
}

/* ── B5 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B5", "Four limitations to state explicitly");
  const L = [
    ["Fe-S clusters are absent","FES (2Fe-2S) lies 5.13 Å from the residue-41 CA, right beside the binding site, but receptor_types carries no metal so it is missing from both sets. The comparison is internally consistent while the electrostatic environment is incomplete. AlphaFold cannot predict cofactors, so there is no alternative — and the fact that ρ = 0.88 against measured Km is reached even so suggests the cluster is not essential for ranking."],
    ["No high-identity homologs","Only 3–4 sequences exceed 90% identity and 62.5–72.5% is empty. The model has to work in the 30–50% regime from the outset."],
    ["Structural confidence in the low-identity bin","The < 30 bin has the least accurate assembly RMSD and only 10 clusters, so the final increment of the learning curve rests on 7 sequences."],
    ["Mixed structure prediction methods","23 of 321 (7%) come from ColabFold. Recorded as struct_source and checked after docking together with mode, clashes and pLDDT."]
  ];
  let y = 1.85;
  L.forEach((l,i) => {
    const h = i===0 ? 1.40 : 0.95;
    card(s, M, y, CW, h, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    badge(s, M+0.28, y+0.20, 0.40, String(i+1), i===0?ACCENT:MID);
    plain(s, l[0], { x:M+0.82, y:y+0.18, w:CW-1.2, size:13, bold:true, color:i===0?ACCENT:TEXT });
    plain(s, l[1], { x:M+0.82, y:y+0.52, w:CW-1.2, size:11, color:TEXT, ls:1.16 });
    y += h + 0.15;
  });
  foot(s);
}

/* ── B6 ── */
{
  const s = newSlide();
  head(s, "BACKUP  B6  ·  FOR EXTERNAL USE", "Outside this group, state the interpretive scope of the docking score");
  card(s, M, 1.95, 6.1, 0.95, TINT, LINE);
  plain(s, "Internal — direct phrasing", { x:M+0.28, y:2.10, w:5.54, size:11.5, bold:true, color:GRAY });
  plain(s, "“dock into the EV binding site to predict binding affinity”",
    { x:M+0.28, y:2.42, w:5.54, size:12, bold:true, color:DEEP });
  card(s, 7.05, 1.95, 5.6, 0.95, "FBEEE2", "E8C9A8");
  plain(s, "External manuscript or talk — needs adjusting", { x:7.33, y:2.10, w:5.04, size:11.5, bold:true, color:GRAY });
  plain(s, "“scores from a box sweeping the mutated region around residue 41 correlate with measured Km”",
    { x:7.33, y:2.42, w:5.04, size:11, bold:true, color:ACCENT, ls:1.13 });
  plain(s, "Why — this is not a direct measurement of site-specific binding affinity", { x:M, y:3.20, w:CW, size:13.5, bold:true, color:TEXT });
  bullets(s, [
    [{text:"The 50-run energy range at 15 Å is 0.22 kcal/mol. ", bold:true}, {text:"The energy landscape is flat."}],
    [{text:"Enlarging the box always reveals a better site farther away. ", bold:true}, {text:"17.23 → 21.42 → 23.62 Å — there is no converging binding site."}],
    [{text:"The 8X9F crystallographic EV centroid sits 12.65 Å from the grid center.", bold:true}],
    [{text:"EV occupies markedly different positions in the MDW and GLW variant crystal structures. ", bold:true}, {text:"The crystallographic EV site is not conserved, so it cannot anchor the box design."}]
  ], { x:M, y:3.60, w:CW, h:2.2, size:12, gap:11 });
  card(s, M, 5.95, CW, 0.95, TINT2, LINE);
  plain(s, "In the manuscript title, enhanced electron mediator affinity states the objective; the Results should describe ranking throughout.",
    { x:M+0.32, y:6.20, w:CW-0.64, size:12, bold:true, color:DEEP });
  foot(s);
}

const out = __dirname + "/20260916_progress_en.pptx";
pres.writeFile({ fileName: out }).then(() => console.log("WROTE " + out + "  slides=" + pageNo));
