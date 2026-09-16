const path = require("path");
const PptxGenJS = require("pptxgenjs");

const INK="16183A", DEEP="2F3596", MID="6A72CE", ACCENT="C0692B",
      TINT="EEF0FA", TINT2="F7F8FC", WHITE="FFFFFF",
      GRAY="5C6180", TEXT="1C1E33", LINE="D6DAEE", INKSUB="A8AECF";
const F = "Malgun Gothic";
const W=13.3, H=7.5, M=0.65, CW=W-2*M;

const pres = new PptxGenJS();
pres.layout = "LAYOUT_WIDE";
pres.author = "CODH project";
pres.title  = "도킹 기반 pLM 증류 — 진척 보고 2026-09-16";

let pageNo = 0;
function newSlide(dark){
  const s = pres.addSlide();
  if (dark) s.background = { color: INK };
  return s;
}
function head(s, kicker, title, dark){
  s.addText(kicker, { x:M, y:0.40, w:CW, h:0.26, fontFace:F, fontSize:11, bold:true,
    color: dark?ACCENT:ACCENT, charSpacing:1.4, isTextBox:true, margin:0 });
  s.addText(title, { x:M, y:0.68, w:CW, h:0.92, fontFace:F, fontSize:27, bold:true,
    color: dark?WHITE:TEXT, isTextBox:true, margin:0, valign:"top", lineSpacingMultiple:1.05 });
}
function foot(s, dark){
  pageNo += 1;
  s.addText(String(pageNo), { x:W-M-0.6, y:H-0.52, w:0.6, h:0.25, fontFace:F, fontSize:9,
    color: dark?INKSUB:GRAY, align:"right", isTextBox:true, margin:0 });
}
// bullets: items = [ [ {text,bold?,color?}, ... ], ... ]
function bullets(s, items, box){
  const paras = [];
  items.forEach((runs, i) => {
    const rs = Array.isArray(runs) ? runs : [{ text: runs }];
    rs.forEach((r, j) => {
      paras.push({ text: r.text, options: {
        bullet: false,
        bold: !!r.bold,
        color: r.color || box.color || TEXT,
        breakLine: (j === rs.length-1) && (i < items.length-1),
        paraSpaceAfter: (j === rs.length-1) ? (box.gap===undefined?9:box.gap) : 0
      }});
    });
  });
  s.addText(paras, { x:box.x, y:box.y, w:box.w, h:box.h, fontFace:F,
    fontSize: box.size||13, color: box.color||TEXT, isTextBox:true, margin:0,
    valign:"top", lineSpacingMultiple:1.18 });
}
function plain(s, txt, box){
  s.addText(txt, { x:box.x, y:box.y, w:box.w, h:box.h, fontFace:F, fontSize:box.size||13,
    bold:!!box.bold, color:box.color||TEXT, align:box.align||"left", isTextBox:true,
    margin:box.margin===undefined?0:box.margin, valign:box.valign||"top",
    lineSpacingMultiple:box.ls||1.15, italic:!!box.italic });
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
  s.addText(txt, { x:x+0.06, y, w:w-0.12, h, fontFace:F, fontSize:size||11.5, bold:true,
    color:col, align:"center", valign:"middle", isTextBox:true, margin:0, lineSpacingMultiple:1.05 });
}
function arrow(s, x, y, txt){
  s.addText(txt||"→", { x, y, w:0.3, h:0.3, fontFace:F, fontSize:14, bold:true,
    color:MID, align:"center", valign:"middle", isTextBox:true, margin:0 });
}
function table(s, rows, box){
  s.addTable(rows, { x:box.x, y:box.y, w:box.w, colW:box.colW, fontFace:F,
    fontSize:box.size||11, color:TEXT, border:{ type:"solid", pt:0.75, color:LINE },
    align:"center", valign:"middle", rowH:box.rowH||0.3, autoPage:false });
}
function hrow(cells, fill){
  return cells.map(c => ({ text:c, options:{ bold:true, color:WHITE, fill:{ color: fill||DEEP }, fontSize:10.5 } }));
}
function drow(cells, opts){
  return cells.map((c,i) => ({ text:c, options: Object.assign({ align: i===0?"left":"center" }, (opts&&opts[i])||{}) }));
}
function stat(s, x,y,w, big, label, col, size){
  const fs = size||34, off = fs*0.0195;
  s.addText(big, { x, y, w, h:off+0.06, fontFace:F, fontSize:fs, bold:true, color:col||DEEP,
    align:"center", isTextBox:true, margin:0 });
  s.addText(label, { x, y:y+off, w, h:0.4, fontFace:F, fontSize: fs>=30?10.5:9.5, color:GRAY,
    align:"center", isTextBox:true, margin:0, lineSpacingMultiple:1.08 });
}
function note(s, x, y, w, txt, col){
  s.addText(txt, { x, y, w, h:0.5, fontFace:F, fontSize:11.5, bold:true, color:col||DEEP,
    isTextBox:true, margin:0, lineSpacingMultiple:1.15 });
}

/* ───────────────────────── S1 표지 ───────────────────────── */
{
  const s = newSlide(true);
  s.addText("진척 보고  ·  2026-09-16", { x:M, y:1.55, w:CW, h:0.3, fontFace:F, fontSize:12,
    bold:true, color:ACCENT, charSpacing:1.6, isTextBox:true, margin:0 });
  s.addText("Engineering CO dehydrogenases with\nenhanced electron mediator affinity through\ndocking-distilled protein language models",
    { x:M, y:2.05, w:CW-0.6, h:2.1, fontFace:F, fontSize:30, bold:true, color:WHITE,
      isTextBox:true, margin:0, lineSpacingMultiple:1.18 });
  s.addText("도킹으로 증류한 단백질 언어모델을 이용한 CO 탈수소효소의 전자매개체 결합친화도 개량",
    { x:M, y:4.30, w:CW-0.6, h:0.4, fontFace:F, fontSize:13, color:INKSUB, isTextBox:true, margin:0 });
  const tags = ["① 수용체 = 이량체", "② 박스 중심 = 41번 CA", "③ 학습 157 / 평가 55", "④ 조립 321개"];
  tags.forEach((t,i) => chip(s, M+i*3.0, 5.35, 2.82, 0.5, t, "242764", WHITE, 11));
  s.addNotes("오늘 확정한 네 가지 설계 결정을 먼저 제시하고, 각각의 근거를 순서대로 전개한다.");
  foot(s, true);
}

/* ───────────────────────── S2 요약 ───────────────────────── */
{
  const s = newSlide();
  head(s, "SUMMARY", "오늘 확정한 것 — 도킹 파이프라인의 네 가지 설계 결정");
  const D = [
    ["①","수용체","이량체 확정","단량체 기각. 실측 Km 상관이 0.88 → 0.55 로 붕괴한다."],
    ["②","박스 중심","F41 CE2 → 41번 CA","등가 ρ 0.968. 적용 가능 서열 1,613 → 1,750 으로 증가."],
    ["③","데이터셋","학습 157 / 평가 55","212 대표를 클러스터 단위로 분할. 누수 0, 라벨 생성 전 고정."],
    ["④","구조","이량체 321개 조립","조립 편향 0.080 kcal/mol < 측정 IQR 0.09. 왜곡 없음."]
  ];
  const cw=5.85, ch=1.95;
  D.forEach((d,i) => {
    const x = M + (i%2)*(cw+0.30), y = 1.80 + Math.floor(i/2)*(ch+0.28);
    card(s, x, y, cw, ch, TINT2, LINE);
    badge(s, x+0.30, y+0.28, 0.44, d[0], DEEP);
    plain(s, d[1], { x:x+0.86, y:y+0.36, w:cw-1.2, size:12, bold:true, color:GRAY });
    plain(s, d[2], { x:x+0.30, y:y+0.86, w:cw-0.6, size:21, bold:true, color:DEEP });
    plain(s, d[3], { x:x+0.30, y:y+1.38, w:cw-0.6, size:11.5, color:TEXT, ls:1.15 });
  });
  plain(s, "진행 중   도킹 321개 실행  →  노이즈 바닥(SNR) 측정  →  라벨 부착  →  학습 A / B / C / D1 / D2",
    { x:M, y:6.35, w:CW, size:12, bold:true, color:ACCENT });
  s.addNotes("이 한 장만 보고도 오늘 회의가 끝나야 한다. 이후 슬라이드는 각 결정의 근거 전개.");
  foot(s);
}

/* ───────────────────────── S3 Objective & scope ───────────────────────── */
{
  const s = newSlide();
  head(s, "OBJECTIVE  &  SCOPE", "서열만으로 CODH 족의 매개체 결합친화도 순위를 매긴다");
  card(s, M, 1.82, CW, 0.88, TINT, LINE);
  plain(s, "도킹 점수로 증류한 단백질 언어모델로 후보 서열의 순위를 매겨, 전자매개체(EV) 결합친화도가 향상된 CO 탈수소효소 서열을 찾아낸다.",
    { x:M+0.32, y:2.02, w:CW-0.64, size:14, bold:true, color:DEEP, ls:1.18 });

  const OC = [
    ["무엇이 어려운가", ACCENT, "실측 Km이 21건뿐이라 서열 모델을 직접 지도할 수 없다. 이 족에는 90%를 넘는 근연 상동체도 거의 없다."],
    ["접근", DEEP, "도킹으로 라벨을 대량 생성한다 — 변이체 288개와 상동체 321개. 회귀 헤드가 이를 ESM 임베딩으로 증류하므로 추론에는 구조가 필요 없고 족 전체 1,827개로 확장된다."],
    ["성공 기준", DEEP, "고정된 55개 평가셋의 Spearman ρ. 109개 클러스터 내 쌍에서 측정한 도킹 노이즈 바닥을 넘어야 한다."]
  ];
  OC.forEach((c,i) => {
    const x = M + i*4.12;
    card(s, x, 2.88, 3.85, 1.62, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    plain(s, c[0], { x:x+0.26, y:3.06, w:3.33, size:12.5, bold:true, color:c[1] });
    plain(s, c[2], { x:x+0.26, y:3.40, w:3.33, size:10.5, color:TEXT, ls:1.2 });
  });

  plain(s, "Scope strategy — 두 레짐이 서로 답하지 못하는 것을 답한다",
    { x:M, y:4.70, w:CW, size:13.5, bold:true, color:TEXT });
  const RG = [
    ["변이체 레짐", "local  ·  dense", "1SU6 골격 하나 위의 변이체 288개, 15위치 × 19치환. 곁사슬만 바뀌고 골격 이완이 없다. 이 레짐이 도킹 점수를 실측 Km에 고정시킨다 (ρ = 0.90)."],
    ["상동체 레짐", "global  ·  sparse", "identity 28.5–74.8%에 걸친 클러스터 대표 212개, 구조마다 별개로 예측하고 조립했다. 이 레짐이 골격이 바뀌어도 점수가 살아남는지를 검정한다."]
  ];
  RG.forEach((r,i) => {
    const x = M + i*6.15;
    card(s, x, 5.04, 5.85, 1.30, TINT2, LINE);
    plain(s, r[0], { x:x+0.28, y:5.20, w:3.2, size:13, bold:true, color:DEEP });
    plain(s, r[1], { x:x+3.40, y:5.24, w:2.2, size:10.5, bold:true, color:ACCENT, align:"right" });
    plain(s, r[2], { x:x+0.28, y:5.56, w:5.3, size:10.5, color:TEXT, ls:1.2 });
  });
  plain(s, "SCOPE 밖   절대 결합친화도  ·  Fe-S 보조인자 정전기 환경  ·  기존 21건을 넘는 실험 검증",
    { x:M, y:6.48, w:CW, size:10.5, bold:true, color:GRAY });
  s.addNotes("두 레짐은 상호 보완이다. 변이체는 점수를 실험에 묶고, 상동체는 일반화를 검정한다. 어느 한쪽만으로는 주장이 성립하지 않는다.");
  foot(s);
}

/* ───────────────────────── S3 문제 정의 ───────────────────────── */
{
  const s = newSlide();
  head(s, "전략", "문제는 라벨 병목이다");
  bullets(s, [
    [{text:"실측 Km 은 21건뿐이다. ", bold:true}, {text:"단백질 언어모델 회귀를 직접 지도하기에 부족하다."}],
    [{text:"도킹 점수가 대리 라벨이 된다. ", bold:true}, {text:"변이체 288개에서 실측 Km 과 ρ = 0.89 를 보인다."}],
    [{text:"전략 — docking-distilled pLM. ", bold:true}, {text:"EV 결합부위 도킹으로 결합친화도 예측값을 대량 생성해 서열 모델로 증류하고, 추론 단계에서는 구조 없이 서열만으로 예측한다."}]
  ], { x:M, y:1.95, w:6.3, h:2.6, size:13.5, gap:14 });
  card(s, M, 4.75, 6.3, 1.55, TINT2, LINE);
  plain(s, "왜 증류인가", { x:M+0.3, y:4.95, w:5.7, size:12, bold:true, color:ACCENT });
  plain(s, "도킹은 구조가 있어야 하고 한 서열당 수 분이 든다. 증류된 모델은 구조 없이 서열만으로 즉시 순위를 낸다 — 족 전체 1,827개로 확장 가능해진다.",
    { x:M+0.3, y:5.32, w:5.7, size:11.5, color:TEXT, ls:1.2 });

  const bx=7.55, bw=5.1;
  const steps = [
    ["실측 Km", "21 건", ACCENT],
    ["도킹 점수 생성", "변이체 288 + 상동체 321", DEEP],
    ["ESM 임베딩 + 회귀 헤드", "SeekRank", DEEP],
    ["서열만으로 순위 예측", "족 전체 1,827 로 확장", MID]
  ];
  steps.forEach((st,i) => {
    const y = 1.95 + i*1.13;
    card(s, bx, y, bw, 0.82, i===0?"FBEEE2":TINT, i===0?"E8C9A8":LINE);
    plain(s, st[0], { x:bx+0.25, y:y+0.13, w:bw-0.5, size:13, bold:true, color:st[2] });
    plain(s, st[1], { x:bx+0.25, y:y+0.45, w:bw-0.5, size:10.5, color:GRAY });
    if (i<3) s.addShape(pres.ShapeType.downArrow, { x:bx+bw/2-0.11, y:y+0.87, w:0.22, h:0.22,
      fill:{ color: MID }, line:{ type:"none" } });
  });
  s.addNotes("라벨 병목이 이 프로젝트의 출발점. 도킹은 라벨 생성기이고 pLM이 최종 예측기다.");
  foot(s);
}

/* ───────────────────────── S4 파이프라인 ───────────────────────── */
{
  const s = newSlide();
  head(s, "OVERVIEW", "파이프라인 전체도");
  const rows = [
    ["서열", [["IPR010047\n1,827"],["70% 클러스터링\n219"],["도킹 가능성 게이트\n1,722"],["medoid 대표\n212"],["학습 157 / 평가 55\n+ D2 확장 109"]]],
    ["구조", [["AFDB\n298"],["ColabFold\n23"],["아키텍처별 참조 중첩\nSF4→1SU6 · SF2→6OND"],["이량체 조립\n321 / 321"]]],
  ];
  let y = 1.85;
  rows.forEach(([lab, chips]) => {
    chip(s, M, y, 1.25, 0.95, lab, INK, WHITE, 13);
    const n = chips.length, ax = 0.34;
    const cwid = (CW - 1.25 - 0.28 - (n-1)*ax) / n;
    chips.forEach((c,i) => {
      const x = M + 1.25 + 0.28 + i*(cwid+ax);
      chip(s, x, y, cwid, 0.95, c[0], TINT, DEEP, 11);
      if (i<n-1) arrow(s, x+cwid+0.02, y+0.33);
    });
    y += 1.18;
  });
  const wide = [
    ["도킹", "각 구조의 41번 CA 중심  ·  반변 15.000 Å (npts 80 × 0.375)  ·  AutoDock-GPU 50런  ·  주 지표 med"],
    ["학습", "SeekRank 회귀 헤드  ·  실험 A 유사도 순차 / B 한번에 / C 무작위 순차 / D1 데이터 양 / D2 클러스터 내 확장"]
  ];
  wide.forEach(([lab, txt]) => {
    chip(s, M, y, 1.25, 0.95, lab, INK, WHITE, 13);
    chip(s, M+1.53, y, CW-1.53, 0.95, txt, TINT, DEEP, 11.5);
    y += 1.18;
  });
  plain(s, "서열 선정 · 구조 확보 · 이량체 조립까지 확정되었다. 도킹부터가 다음 단계다.", { x:M, y:6.55, w:CW, size:11, italic:true, color:GRAY });
  s.addNotes("서열 선정 → 구조 확보/조립 → 도킹 → 학습. 오늘은 앞의 세 단계가 확정됐다.");
  foot(s);
}

/* ───────────────────────── S5 이량체 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ①  수용체", "수용체는 이량체여야 한다 — 단량체는 예측력을 잃는다");
  table(s, [
    hrow(["지표 (n=21)","이량체","단량체","차이"]),
    drow(["ρ(실측 Km), med_lig","0.903","0.531","−0.372"], [null,{bold:true,color:DEEP},{bold:true,color:ACCENT},null]),
    drow(["ρ(실측 Km), med_ds","0.856","0.512","−0.344"]),
    drow(["ρ(실측 Km), best_ds","0.876","0.558","−0.319"]),
    drow(["ρ(순위, n=288), med_ds","—","0.702","—"])
  ], { x:M, y:1.95, w:6.6, colW:[2.7,1.3,1.3,1.3], rowH:0.42, size:11 });
  plain(s, "단량체 절단 대조군 288개 · 중심 · 반변 · nrun · seed 모두 동일",
    { x:M, y:4.15, w:6.6, size:10.5, italic:true, color:GRAY });
  bullets(s, [
    [{text:"수렴 실패가 아니다. ", bold:true}, {text:"단량체 런의 리간드 교체 상관은 0.995 로 오히려 더 높다. 포켓이 단순해져 경쟁 최소점이 줄었을 뿐 계산은 정상 수렴했다."}],
    [{text:"구조를 바꾼 결과다. ", bold:true}, {text:"잘려나간 SU2 의 기여가 무작위 잡음이 아니라 실측과 맞아떨어지던 신호였다."}],
    [{text:"귀결 — 상동체도 multimer 예측이 전제. ", bold:true, color:ACCENT}, {text:"AF2 단량체로 상동체를 만드는 계획은 폐기했다.", color:ACCENT}]
  ], { x:M, y:4.60, w:6.6, h:2.3, size:12, gap:10 });

  const bx=7.65;
  card(s, bx, 1.95, 5.0, 2.55, TINT2, LINE);
  plain(s, "실측 Km 상관  (med_lig)", { x:bx+0.3, y:2.15, w:4.4, size:11.5, bold:true, color:GRAY });
  stat(s, bx+0.25, 2.60, 2.0, "0.903", "이량체", DEEP);
  s.addText("→", { x:bx+2.30, y:2.68, w:0.45, h:0.5, fontFace:F, fontSize:22, bold:true,
    color:GRAY, align:"center", valign:"middle", isTextBox:true, margin:0 });
  stat(s, bx+2.75, 2.60, 2.0, "0.531", "단량체", ACCENT);
  plain(s, "예측력의 41% 손실", { x:bx+0.3, y:4.02, w:4.4, size:12, bold:true, color:ACCENT, align:"center" });

  card(s, bx, 4.72, 5.0, 1.85, TINT, LINE);
  plain(s, "왜 이량체가 필요한가", { x:bx+0.3, y:4.92, w:4.4, size:11.5, bold:true, color:DEEP });
  plain(s, "EV 결합부위가 두 protomer 의 C2 계면에 걸쳐 있기 때문이다. 단량체로 자르면 결합부위의 절반이 사라진다.",
    { x:bx+0.3, y:5.28, w:4.4, size:11.5, color:TEXT, ls:1.2 });
  s.addNotes("가장 중요한 결정. 순위 상관 손실(0.70)보다 Km 예측력 손실이 훨씬 크다는 점이 핵심.");
  foot(s);
}

/* ───────────────────────── S6 계면 근거 ───────────────────────── */
{
  const s = newSlide();
  head(s, "근거", "EV 결합부위가 C2 계면에 걸쳐 있다");
  bullets(s, [
    [{text:"박스 내 원자가 두 protomer 에 나뉜다. ", bold:true}, {text:"SU1 436개(54.6%) / SU2 362개(45.4%)."}],
    [{text:"양쪽이 같은 34–67 루프를 마주 댄다. ", bold:true}, {text:"중심 12 Å 이내 포켓 벽이 SU1 28잔기, SU2 18잔기로 구성된다."}],
    [{text:"박스가 대칭축을 품는다. ", bold:true}, {text:"그리드 중심에서 C2 축까지 4.76 Å. 두 protomer 의 F41 CE2 는 9.53 Å 떨어진 대칭 쌍이다."}],
    [{text:"C-클러스터는 박스 밖이다. ", bold:true, color:ACCENT}, {text:"CO 산화 활성자리까지 26.4 Å. 박스는 D-클러스터 주변, 즉 전자 전달 경로에 앉아 있다.", color:ACCENT}]
  ], { x:M, y:1.95, w:6.2, h:3.4, size:12.5, gap:13 });
  card(s, M, 5.45, 6.2, 1.15, TINT2, LINE);
  plain(s, "이량체가 같은 체인 A 에 같은 잔기번호로 들어 있어 잔기번호 기준 점검으로는 구분되지 않는다. 원자 수(11,370)와 중복쌍 거리 중앙값(51.61 Å)으로 확인했다.",
    { x:M+0.28, y:5.63, w:5.64, size:10.5, color:TEXT, ls:1.2 });

  const cx=9.95, cy=3.35;
  card(s, 7.15, 1.95, 5.5, 4.65, TINT2, LINE);
  plain(s, "반변 15 Å 박스", { x:cx-1.10, y:2.05, w:2.2, size:10.5, bold:true, color:ACCENT, align:"center" });
  s.addShape(pres.ShapeType.line, { x:cx, y:2.32, w:0, h:2.78, line:{ color:GRAY, width:1.25, dashType:"dash" } });
  s.addShape(pres.ShapeType.ellipse, { x:cx-2.55, y:cy-1.25, w:2.5, h:2.5,
    fill:{ color:"D5DAF2" }, line:{ color:MID, width:1.25 } });
  s.addShape(pres.ShapeType.ellipse, { x:cx+0.05, y:cy-1.25, w:2.5, h:2.5,
    fill:{ color:"E4E7F7" }, line:{ color:MID, width:1.25 } });
  plain(s, "SU1", { x:cx-2.30, y:cy-0.20, w:1.2, size:14, bold:true, color:DEEP, align:"center" });
  plain(s, "SU2", { x:cx+1.10, y:cy-0.20, w:1.2, size:14, bold:true, color:DEEP, align:"center" });
  s.addShape(pres.ShapeType.roundRect, { x:cx-0.95, y:cy-0.95, w:1.9, h:1.9,
    fill:{ type:"none" }, line:{ color:ACCENT, width:2.25 }, rectRadius:0.04 });
  plain(s, "C2 축", { x:cx+0.14, y:4.72, w:1.0, size:10, bold:true, color:GRAY });
  plain(s, "박스 내 원자   SU1 54.6%  /  SU2 45.4%", { x:7.45, y:5.42, w:5.0, size:11.5, bold:true, color:DEEP, align:"center" });
  plain(s, "그리드 중심 → C2 축  4.76 Å      ·      C-클러스터  26.4 Å (박스 밖)",
    { x:7.45, y:5.82, w:5.0, size:10, color:GRAY, align:"center" });
  plain(s, "모식도 — 외부 발표 시 PyMOL 렌더로 교체", { x:7.45, y:6.18, w:5.0, size:9, italic:true, color:GRAY, align:"center" });
  s.addNotes("숫자만으로는 계면이라는 게 안 와닿아서 모식도를 넣었다. 실제 렌더는 외부용.");
  foot(s);
}

/* ───────────────────────── S7 CA 전환 사유 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ②  박스 중심", "박스 중심을 41번 CA 로 전환한다 — 상동체 전 계통에 정의된다");
  table(s, [
    hrow(["41번 잔기","서열 수","CE2"]),
    drow(["F / Y / W  (방향족)","1,613","있음"], [null,null,{color:DEEP,bold:true}]),
    drow(["I  (Ile 계통)","118","없음"], [null,{bold:true,color:ACCENT},{color:ACCENT,bold:true}]),
    drow(["M · T · H · L · S","19","없음"], [null,null,{color:ACCENT}]),
    drow(["주석 없음","77","없음"], [null,null,{color:ACCENT}])
  ], { x:M, y:1.95, w:5.8, colW:[2.9,1.45,1.45], rowH:0.42, size:11 });
  card(s, M, 4.20, 5.8, 1.05, "FBEEE2", "E8C9A8");
  plain(s, "214 개 (11.7%) 에 CE2 가 없다", { x:M+0.3, y:4.38, w:5.2, size:17, bold:true, color:ACCENT });
  plain(s, "His 는 방향족이지만 CE1 / NE2 만 가져 역시 불가", { x:M+0.3, y:4.78, w:5.2, size:10.5, color:GRAY });
  plain(s, "전 계통을 포함하려면 CE2 정의는 사용할 수 없다.", { x:M, y:5.45, w:5.8, size:12, bold:true, color:TEXT });
  plain(s, "비방향족 41번은 오류가 아니라 상동체에 실재하는 생물학적 사실이다.",
    { x:M, y:5.78, w:5.8, size:11.5, color:GRAY, ls:1.2 });

  const bx=7.05, bw=5.6;
  bullets(s, [
    [{text:"41번 자리 자체는 유지한다. ", bold:true}, {text:"변이 15위치(40, 42, 43, 44, 46, 57–66)에서 41번만 제외되어 있다. EV 와 전자를 주고받는 키 잔기라 변이 대상에서 뺀 것이므로, 중심을 41번에 두는 것은 편의가 아니라 기전적 선택이다."}],
    [{text:"정의가 서열에서 자동 결정된다. ", bold:true}, {text:"f41_pos = anchor_pos + 2 가 1,750 / 1,750 예외 없이 성립한다. 41번은 FeS 결합 Cys 모티프의 3번째 잔기다."}]
  ], { x:bx, y:1.95, w:bw, h:2.5, size:12.5, gap:13 });
  card(s, bx, 4.45, bw, 1.62, TINT, LINE);
  const motif = [["C39","",INK,WHITE],["G40","\uBCC0\uC774",TINT2,DEEP],["F41","\uC911\uC2EC",ACCENT,WHITE],
                 ["G42","\uBCC0\uC774",TINT2,DEEP],["E43","\uBCC0\uC774",TINT2,DEEP],["T44","\uBCC0\uC774",TINT2,DEEP],
                 ["x45","",TINT2,GRAY],["L46","\uBCC0\uC774",TINT2,DEEP],["C47","",INK,WHITE]];
  motif.forEach((m,i) => {
    const x = bx + 0.03 + i*0.62;
    chip(s, x, 4.68, 0.58, 0.44, m[0], m[2], m[3], 10.5);
    if (m[1]) plain(s, m[1], { x:x-0.02, y:5.16, w:0.62, size:8.5, color:GRAY, align:"center" });
  });
  plain(s, "\uBCC0\uC774 40 \u00B7 42 \u00B7 43 \u00B7 44 \u00B7 46 \uC740 \uC774 \uBAA8\uD2F0\uD504\uC758 \uBE44-Cys \uC790\uB9AC\uB2E4",
    { x:bx+0.25, y:5.58, w:bw-0.5, size:10.5, italic:true, color:GRAY, align:"center" });
  plain(s, "Cys 쌍 중점은 두 계통의 Cys 간격이 달라(8칸 vs 3칸) 기하가 어긋나므로 사용 금지.",
    { x:bx, y:6.20, w:bw, size:10.5, color:GRAY });
  s.addNotes("CE2를 못 쓰는 이유가 핵심. 41번 자리 자체는 기전적 근거로 유지한다.");
  foot(s);
}

/* ───────────────────────── S8 등가성 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ②  검증", "전환해도 결과는 같다 — 중심 4.661 Å 이동, 순위 보존");
  table(s, [
    hrow(["지표","ρ (CE2 대비)","ρ 위치내","ρ (실측 Km)"]),
    drow(["med_ds","0.968","0.928","0.856 → 0.889"], [null,null,null,{bold:true,color:DEEP}]),
    drow(["med_lig","0.971","0.929","0.903 → 0.896"]),
    drow(["best_ds","0.978","0.956","0.876 → 0.884"]),
    drow(["best_lig","0.970","0.954","0.882 → 0.892"])
  ], { x:M, y:1.95, w:6.5, colW:[1.7,1.6,1.4,1.8], rowH:0.40, size:11 });
  plain(s, "288개 재도킹 · 반변 15 고정, 중심만 교체 · 기준선 = CE2 + 반변 15",
    { x:M, y:4.08, w:6.5, size:10.5, italic:true, color:GRAY });
  plain(s, "왜 4.661 Å 를 옮겨도 안 바뀌는가", { x:M, y:4.55, w:6.5, size:13, bold:true, color:DEEP });
  bullets(s, [
    [{text:"두 박스가 76.8% 겹친다. ", bold:true}, {text:"축별 겹침 0.920 × 0.875 × 0.954."}],
    [{text:"같은 자리에 도킹한다. ", bold:true}, {text:"1·2위 군집 좌표 차이가 2.86 / 2.83 Å — 격자 해상도 이내, 동일한 물리적 자리다."}],
    [{text:"남는 차이는 격자 이산화뿐. ", bold:true}, {text:"평균차 +0.01~0.02 로 편향이 없고 개별 값만 흔들린다."}]
  ], { x:M, y:4.92, w:6.5, h:1.9, size:12, gap:9 });

  const bx=7.55, bw=5.1;
  card(s, bx, 1.95, bw, 2.15, TINT2, LINE);
  plain(s, "적용 가능 서열", { x:bx+0.3, y:2.15, w:bw-0.6, size:11.5, bold:true, color:GRAY });
  stat(s, bx+0.25, 2.55, 1.95, "1,613", "CE2 기준", GRAY);
  s.addText("→", { x:bx+2.25, y:2.63, w:0.45, h:0.5, fontFace:F, fontSize:22, bold:true,
    color:ACCENT, align:"center", valign:"middle", isTextBox:true, margin:0 });
  stat(s, bx+2.70, 2.55, 1.95, "1,750", "41번 CA 기준", DEEP);
  card(s, bx, 4.30, bw, 1.25, "FBEEE2", "E8C9A8");
  plain(s, "med_ds 는 오히려 상승했다", { x:bx+0.3, y:4.52, w:bw-0.6, size:13, bold:true, color:ACCENT });
  plain(s, "0.856 → 0.889.  med_lig 는 0.903 → 0.896 으로 동등.",
    { x:bx+0.3, y:4.88, w:bw-0.6, size:11.5, color:TEXT });
  card(s, bx, 5.72, bw, 0.9, TINT, LINE);
  plain(s, "결론 — 커버리지는 늘고 성능 손실은 없다.",
    { x:bx+0.3, y:5.94, w:bw-0.6, size:12.5, bold:true, color:DEEP });
  s.addNotes("전환의 대가가 없다는 점이 요지. 오히려 med_ds는 올랐다.");
  foot(s);
}

/* ───────────────────────── S9 박스 크기 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ②  박스 크기", "반변 15 Å 은 변경 불가 — 키우면 리간드가 홈을 따라 이탈한다");
  table(s, [
    hrow(["반변 (Å)","15 Å 내 포즈","변이체별 최소","ρ 위치내"]),
    drow(["15.000","0.999","0.82","0.928"], [{bold:true},{bold:true,color:DEEP},{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["16.125","0.301","0.00","0.833"], [null,null,{bold:true,color:ACCENT},null]),
    drow(["17.250","0.176","0.00","0.781"]),
    drow(["18.000","0.168","0.00","0.779"]),
    drow(["20.250","0.136","0.00","0.774"])
  ], { x:M, y:1.95, w:6.2, colW:[1.55,1.6,1.6,1.45], rowH:0.365, size:11 });
  card(s, M, 4.50, 6.2, 1.0, "FBEEE2", "E8C9A8");
  plain(s, "r16 에서 이미 변이체별 최소가 0.00 이다", { x:M+0.28, y:4.66, w:5.6, size:13, bold:true, color:ACCENT });
  plain(s, "50개 포즈가 전부 15 Å 밖으로 나가는 변이체가 존재한다 — 비균질 편향.",
    { x:M+0.28, y:5.03, w:5.6, size:11, color:TEXT });
  bullets(s, [
    [{text:"절벽이 아니라 15.0 – 17.25 에 걸친 하강 후 평탄. ", bold:true}, {text:"하강의 60% 가 첫 1.1 Å 에서 일어난다."}],
    [{text:"반변 15 는 홈을 입구에서 끊는 크기다. ", bold:true, color:DEEP}, {text:"그 구간이 변이 잔기 영역과 겹쳐서 신호가 나온다.", color:DEEP}]
  ], { x:M, y:5.70, w:6.2, h:1.3, size:11.5, gap:8 });

  const bx=7.35, bw=5.3;
  card(s, bx, 1.95, bw, 4.4, TINT2, LINE);
  plain(s, "군집 중심이 같은 홈을 따라 밀려난다", { x:bx+0.3, y:2.15, w:bw-0.6, size:12.5, bold:true, color:DEEP });
  plain(s, "방향 (1, 0, 1)/√2 — y 는 3 에 고정된 채 x · z 만 늘어난다",
    { x:bx+0.3, y:2.50, w:bw-0.6, size:10.5, color:GRAY });
  const ax0=bx+0.55, axw=4.2, ay=4.30;
  s.addShape(pres.ShapeType.rect, { x:ax0, y:ay-0.85, w:1.35, h:1.7, fill:{ color:"E2E6F6" }, line:{ type:"none" } });
  plain(s, "변이 잔기 영역", { x:ax0, y:ay-1.12, w:1.35, size:9.5, bold:true, color:DEEP, align:"center" });
  s.addShape(pres.ShapeType.line, { x:ax0, y:ay, w:axw, h:0, line:{ color:GRAY, width:1.25 } });
  s.addShape(pres.ShapeType.line, { x:ax0+1.35, y:ay-0.95, w:0, h:1.9, line:{ color:ACCENT, width:1.5, dashType:"dash" } });
  plain(s, "반변 15 경계", { x:ax0+0.95, y:ay+1.00, w:1.6, size:9.5, bold:true, color:ACCENT, align:"center" });
  const pts = [[0.28,"[3,3,6]","r15"],[1.02,"[6,3,6]","r15"],[2.55,"[12,3,12]","r16–18"],[3.75,"[15,3,15]","r20–22"]];
  pts.forEach(([dx,lab,rad],i) => {
    const x = ax0 + dx;
    s.addShape(pres.ShapeType.ellipse, { x:x-0.10, y:ay-0.10, w:0.20, h:0.20,
      fill:{ color: i<2?DEEP:ACCENT }, line:{ type:"none" } });
    plain(s, lab, { x:x-0.45, y:ay-0.52, w:0.90, size:9.5, bold:true, color:i<2?DEEP:ACCENT, align:"center" });
    plain(s, rad, { x:x-0.45, y:ay+0.16, w:0.90, size:9, color:GRAY, align:"center" });
  });
  plain(s, "리간드가 변이 잔기를 지나쳐 홈 안쪽으로 들어가면, 우리가 재려는 것과 다른 것을 재게 된다.",
    { x:bx+0.3, y:5.62, w:bw-0.6, size:11, bold:true, color:TEXT, ls:1.2 });
  s.addNotes("반변 15는 타협의 결과가 아니라 물리적 근거가 있는 경계. 절대 바꾸지 말 것.");
  foot(s);
}

/* ───────────────────────── S10 원칙 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ③  데이터셋", "평가셋은 라벨을 보기 전에 고정했다 — 네 가지 원칙");
  const P = [
    ["1","클러스터링을 게이트보다 먼저 한다","서열 분할은 단백질 족의 성질이지 우리 파이프라인의 성질이 아니다. 기술적 제약을 먼저 걸면 클러스터 경계가 그 제약에 따라 달라져 재현되지 않는다."],
    ["2","클러스터당 대표 1개 (CAP = 1)","같은 클러스터의 서열은 70% 이상 닮아 라벨이 거의 같다. 여럿 넣으면 그 계통이 손실함수에서 표를 더 행사할 뿐 새 정보가 없다."],
    ["3","분할 단위는 클러스터, 서열이 아니다","같은 클러스터 멤버가 학습과 평가로 갈리면 외운 것을 다시 맞히는 평가가 된다. 대표가 1개이므로 이 누수는 구조적으로 발생할 수 없다."],
    ["4","라벨을 보기 전에 분할을 확정한다","도킹 점수를 본 뒤 나누면 어떤 결과가 나와도 정당화할 수 있다. 212개와 학습·평가 배정은 도킹 시작 전에 고정되었다."]
  ];
  P.forEach((p,i) => {
    const y = 1.88 + i*1.13;
    badge(s, M+0.05, y+0.10, 0.46, p[0], i<2?DEEP:MID);
    plain(s, p[1], { x:M+0.70, y:y+0.02, w:11.2, size:14, bold:true, color:TEXT });
    plain(s, p[2], { x:M+0.70, y:y+0.38, w:11.2, size:11.5, color:GRAY, ls:1.18 });
  });
  card(s, M, 6.25, CW, 0.62, "FBEEE2", "E8C9A8");
  plain(s, "학습셋으로 평가하면 성능이 과대평가된다 — 이 네 가지가 그에 대한 답이다.",
    { x:M+0.3, y:6.41, w:CW-0.6, size:12, bold:true, color:ACCENT });
  s.addNotes("선정 프로토콜 전체를 지배하는 원칙. 이후 슬라이드는 이 원칙의 실행 결과다.");
  foot(s);
}

/* ───────────────────────── S11 퍼널 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ③  선정", "1,827 에서 212 까지 — 게이트는 대표를 고를 때만 적용한다");
  const F2 = [
    ["IPR010047 전체", "1,827 서열", 6.9, DEEP],
    ["mmseqs 클러스터링  --min-seq-id 0.70", "219 클러스터", 5.9, DEEP],
    ["도킹 가능성 게이트", "1,722 서열", 4.9, MID],
    ["클러스터당 medoid 1개", "212 대표", 3.9, MID],
    ["구간 층화 25% 분할", "학습 157 / 평가 55", 3.0, ACCENT]
  ];
  const ccx = 4.35;
  F2.forEach((f,i) => {
    const y = 1.92 + i*0.92, w = f[2];
    card(s, ccx-w/2, y, w, 0.74, i===4?"FBEEE2":TINT, i===4?"E8C9A8":LINE);
    plain(s, f[0], { x:ccx-w/2+0.2, y:y+0.09, w:w-0.4, size:10.5, color:GRAY, align:"center" });
    plain(s, f[1], { x:ccx-w/2+0.2, y:y+0.36, w:w-0.4, size:13, bold:true, color:f[3], align:"center" });
    if (i<4) s.addShape(pres.ShapeType.downArrow, { x:ccx-0.10, y:y+0.76, w:0.20, h:0.15,
      fill:{ color:MID }, line:{ type:"none" } });
  });
  const bx=8.35, bw=4.3;
  plain(s, "게이트에서 탈락한 것", { x:bx, y:1.92, w:bw, size:13, bold:true, color:TEXT });
  table(s, [
    hrow(["사유","서열","클러스터"]),
    drow(["41번 등가 잔기 부재","80","4"]),
    drow(["품질 필터","25","3"]),
    drow(["합계","105","7"], [{bold:true},{bold:true},{bold:true}])
  ], { x:bx, y:2.30, w:bw, colW:[2.1,1.1,1.1], rowH:0.36, size:10.5 });
  bullets(s, [
    [{text:"41번 부재는 생물학적 사실이다. ", bold:true}, {text:"족의 4.4% 는 해당 자리가 아예 없어 결합 부위를 정의할 수 없다."}],
    [{text:"AFDB 등재 여부는 게이트가 아니다. ", bold:true}, {text:"단백질의 성질이 아니라 데이터베이스의 수록 상태다. 예측 필요분 23개가 전부 A0AB* / A0AC* 였다 — 게이트로 썼다면 최근 등록된 환경·메타게놈 계통을 통째로 잃었다.", }]
  ], { x:bx, y:4.05, w:bw, h:2.4, size:11, gap:10 });
  card(s, bx, 5.72, bw, 1.08, TINT, LINE);
  plain(s, "멤버 13개가 전원 41번 결실인 클러스터가 하나 있다 — F41 자리를 잃은 분류군일 가능성이 있어 별도 확인 대상.",
    { x:bx+0.26, y:5.90, w:bw-0.52, size:10.5, bold:true, color:DEEP, ls:1.18 });
  s.addNotes("클러스터링을 먼저 하고 게이트는 대표 선정 때만 적용한다는 원칙 1의 실행.");
  foot(s);
}

/* ───────────────────────── S12 구간 경계 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ③  구간", "구간 경계는 데이터에 맞췄다 — 의도는 그대로 유지했다");
  card(s, M, 1.95, 6.1, 1.85, TINT2, LINE);
  plain(s, "ChCODH2 대비 identity 분포  (212 클러스터)", { x:M+0.3, y:2.13, w:5.5, size:11.5, bold:true, color:GRAY });
  const dd = [["0.285","최소"],["0.419","중앙"],["0.550","95백분위"],["0.748","최대"]];
  dd.forEach((d,i) => stat(s, M+0.20+i*1.44, 2.55, 1.38, d[0], d[1], DEEP, 21));
  card(s, M, 4.00, 6.1, 1.35, "FBEEE2", "E8C9A8");
  plain(s, "90% 를 넘는 서열은 3~4개뿐이고, 62.5 – 72.5% 구간은 완전히 비어 있다.",
    { x:M+0.3, y:4.22, w:5.5, size:12, bold:true, color:ACCENT, ls:1.2 });
  plain(s, "이 족에는 ChCODH2 와 가까운 친척이 거의 없다.",
    { x:M+0.3, y:4.88, w:5.5, size:11, color:TEXT });
  plain(s, "원안대로면 위 네 구간이 모두 한 자릿수가 되고 마지막 하나에 1,780개가 몰린다.",
    { x:M, y:5.55, w:6.1, size:12, bold:true, color:TEXT, ls:1.2 });
  plain(s, "유사도 순차 학습이라는 의도는 그대로 두고 경계만 데이터에 맞췄다.",
    { x:M, y:6.05, w:6.1, size:12, bold:true, color:DEEP, ls:1.2 });

  const bx=7.35, bw=5.3;
  plain(s, "원안", { x:bx, y:1.95, w:2.4, size:12, bold:true, color:GRAY });
  plain(s, "확정", { x:bx+2.9, y:1.95, w:2.4, size:12, bold:true, color:DEEP });
  const pairs = [["90 – 100","≥ 60","6"],["80 – 90","50 – 60","14"],["70 – 80","40 – 50","107"],["60 – 70","30 – 40","75"],["< 60","< 30","10"]];
  pairs.forEach((p,i) => {
    const y = 2.35 + i*0.82;
    chip(s, bx, y, 2.4, 0.62, p[0], "E9EAEF", GRAY, 12);
    arrow(s, bx+2.48, y+0.16);
    chip(s, bx+2.9, y, 2.4, 0.62, p[1] + "      " + p[2] + " 대표", TINT, DEEP, 12);
  });
  plain(s, "10% 폭 · 클러스터 median fident 기준", { x:bx, y:6.50, w:bw, size:10.5, italic:true, color:GRAY });
  s.addNotes("교수님 원안과 다른 경계를 쓰는 이유. 의도는 유지했다는 점을 강조.");
  foot(s);
}

/* ───────────────────────── S13 결과표 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ③  결과", "212 대표가 족의 98.8% 를 대리한다");
  table(s, [
    hrow(["구간","대표","커버 서열","학습","평가","SF2","비방향족 41"]),
    drow(["≥ 60","6","31","4","2","0","0"]),
    drow(["50 – 60","14","96","10","4","0","0"]),
    drow(["40 – 50","107","619","80","27","13","23"]),
    drow(["30 – 40","75","897","56","19","17","16"]),
    drow(["< 30","10","163","7","3","3","4"]),
    drow(["합계","212","1,806","157","55","33","43"],
      [{bold:true,fill:{color:TINT}},{bold:true,color:DEEP,fill:{color:TINT}},{bold:true,fill:{color:TINT}},
       {bold:true,color:DEEP,fill:{color:TINT}},{bold:true,color:DEEP,fill:{color:TINT}},
       {bold:true,fill:{color:TINT}},{bold:true,fill:{color:TINT}}])
  ], { x:M, y:1.95, w:7.55, colW:[1.35,0.95,1.35,0.95,0.95,0.95,1.05], rowH:0.40, size:11 });
  bullets(s, [
    [{text:"SF2 도 비방향족 41번도 50% 아래에서만 나타난다. ", bold:true}, {text:"고유사도 구간은 개수뿐 아니라 정보량 자체가 낮으며, 그 영역은 이미 변이체 288개가 촘촘히 덮고 있다."}],
    [{text:"공유 클러스터 0 · 중복 서열 0. ", bold:true}, {text:"대표가 클러스터당 1개이므로 누수가 구조적으로 불가능하다."}]
  ], { x:M, y:5.15, w:7.55, h:1.6, size:11.5, gap:9 });

  const bx=8.55, bw=4.1;
  card(s, bx, 1.95, bw, 1.65, TINT2, LINE);
  plain(s, "족 서열 포괄률", { x:bx+0.28, y:2.12, w:bw-0.56, size:11, bold:true, color:GRAY, align:"center" });
  plain(s, "98.8%", { x:bx+0.28, y:2.45, w:bw-0.56, size:38, bold:true, color:DEEP, align:"center" });
  plain(s, "1,827 중 1,806", { x:bx+0.28, y:3.15, w:bw-0.56, size:10.5, color:GRAY, align:"center" });
  card(s, bx, 3.78, bw, 1.15, TINT, LINE);
  plain(s, "주 지표", { x:bx+0.28, y:3.94, w:bw-0.56, size:10.5, bold:true, color:GRAY });
  plain(s, "전체 55개의 Spearman", { x:bx+0.28, y:4.24, w:bw-0.56, size:13, bold:true, color:DEEP });
  card(s, bx, 5.08, bw, 1.67, "FBEEE2", "E8C9A8");
  plain(s, "구간별은 주 지표가 될 수 없다", { x:bx+0.28, y:5.26, w:bw-0.56, size:11.5, bold:true, color:ACCENT });
  plain(s, "≥ 60 의 평가가 2개, < 30 이 3개다. 이 크기의 순위상관은 신뢰구간이 지나치게 넓다. 40–50(27) · 30–40(19) 만 개별 보고하고 양 끝은 합쳐서 기술한다.",
    { x:bx+0.28, y:5.58, w:bw-0.56, size:10, color:TEXT, ls:1.18 });
  s.addNotes("포괄률과 누수 0이 핵심. 구간별 평가의 한계도 미리 못박아 둔다.");
  foot(s);
}

/* ───────────────────────── S14 실험 설계 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ③  실험", "같은 학습셋 · 같은 평가셋 — 다른 것은 순서 · 방식 · 양뿐이다");
  table(s, [
    hrow(["","학습 방식","순서","양","답하는 것","반복"]),
    drow(["A","순차 (이어 학습)","유사도","4 → 14 → 94 → 150 → 157","유사도 순차 학습의 효과","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["B","한번에","—","157","상한 기준선 · A 의 종착점 비교 대상","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["C","순차","무작위","A 와 동일 크기","A 와의 차이 = 순서의 기여","5"],
      [{bold:true,color:ACCENT},null,null,null,{align:"left"},null]),
    drow(["D1","한번에","무작위 (층화)","20 · 40 · … · 157","데이터 양의 포화점","5"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null]),
    drow(["D2","한번에","—","157 + 109 = 266","대표 1개로 충분했는가","1"],
      [{bold:true,color:DEEP},null,null,null,{align:"left"},null])
  ], { x:M, y:1.95, w:CW, colW:[0.75,2.1,1.5,2.7,3.75,1.2], rowH:0.44, size:11 });
  const nb = [
    ["C 가 반드시 필요한 이유", "C 가 빠지면 A 와 B 의 차이가 순차성 때문인지 유사도 순서 때문인지 가릴 수 없다.", ACCENT],
    ["A 는 '이어 학습' 이어야 성립한다", "누적 재학습이면 마지막 점이 B 와 같아져 대비가 사라진다. Ridge · SVR · RF 계열에는 이어 학습 개념이 없다.", DEEP],
    ["임베딩은 한 번만 계산한다", "학습 실행은 72회지만 고유 서열은 321개. ESM 임베딩은 서열당 1회 캐시한다.", DEEP]
  ];
  nb.forEach((n,i) => {
    const x = M + i*4.12;
    card(s, x, 4.85, 3.85, 1.55, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    plain(s, n[0], { x:x+0.25, y:5.03, w:3.35, size:11.5, bold:true, color:n[2] });
    plain(s, n[1], { x:x+0.25, y:5.38, w:3.35, size:10.5, color:TEXT, ls:1.18 });
  });
  plain(s, "미결 — SeekRank 회귀 헤드의 구조 확인 필요. 실험 A 의 가능 여부가 여기에 달려 있다.",
    { x:M, y:6.60, w:CW, size:11, bold:true, color:ACCENT });
  s.addNotes("C의 존재 이유와 A의 이어 학습 조건이 질문 나올 지점.");
  foot(s);
}

/* ───────────────────────── S15 조립 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ④  구조", "321개 전량 이량체 조립 완료 — 아키텍처별 참조를 나눴다");
  plain(s, "구조 출처", { x:M, y:1.95, w:5.6, size:13, bold:true, color:TEXT });
  table(s, [
    hrow(["출처","개수","비고"]),
    drow(["AFDB","298","내려받기 완료"]),
    drow(["ColabFold","23","A0AB* / A0AC* — AFDB 미등재"]),
    drow(["합계","321","core 212 + D2 확장 109"], [{bold:true},{bold:true,color:DEEP},{bold:true}])
  ], { x:M, y:2.33, w:5.9, colW:[1.5,0.95,3.45], rowH:0.38, size:10.5 });
  card(s, M, 4.05, 5.9, 1.2, "FBEEE2", "E8C9A8");
  plain(s, "예측 방법이 섞인다", { x:M+0.28, y:4.22, w:5.34, size:11.5, bold:true, color:ACCENT });
  plain(s, "같은 AlphaFold2 지만 MSA 파이프라인이 다르다. struct_source 컬럼으로 기록하고 도킹 후 점수와의 상관을 확인한다.",
    { x:M+0.28, y:4.55, w:5.34, size:10.5, color:TEXT, ls:1.18 });
  plain(s, "조립 품질", { x:M, y:5.42, w:5.9, size:13, bold:true, color:TEXT });
  const q = [["321 / 321","성공 (실패 0)"],["0.80 – 1.65","RMSD (구간별)"],["0 – 1.22","임계충돌"],["94.6 – 96.6","pLDDT"]];
  q.forEach((x,i) => {
    const cx2 = M + i*1.50;
    card(s, cx2, 5.78, 1.42, 0.82, TINT2, LINE);
    plain(s, x[0], { x:cx2+0.05, y:5.90, w:1.32, size:13, bold:true, color:i===0?DEEP:TEXT, align:"center" });
    plain(s, x[1], { x:cx2+0.05, y:6.24, w:1.32, size:8.5, color:GRAY, align:"center" });
  });

  const bx=7.05, bw=5.6;
  plain(s, "참조 선택 — D-cluster 아키텍처 기준, 예외 없음", { x:bx, y:1.95, w:bw, size:13, bold:true, color:TEXT });
  const refs = [["SF4 · SF4v · SFx","ChCODH2  (1SU6)"],["SF2","DvCODH  (6OND)"]];
  refs.forEach((r,i) => {
    const y = 2.38 + i*0.95;
    chip(s, bx, y, 2.5, 0.75, r[0], TINT, DEEP, 12);
    arrow(s, bx+2.58, y+0.22);
    chip(s, bx+3.0, y, 2.6, 0.75, r[1], INK, WHITE, 12);
  });
  card(s, bx, 4.32, bw, 0.95, TINT2, LINE);
  plain(s, "모드 — global · local 둘 다 시도해 (임계충돌, RMSD) 최소를 채택했다.",
    { x:bx+0.28, y:4.55, w:bw-0.56, size:11.5, color:TEXT, ls:1.2 });
  card(s, bx, 5.42, bw, 1.3, TINT, LINE);
  plain(s, "임계충돌은 구간과 무관하게 평평하다", { x:bx+0.28, y:5.62, w:bw-0.56, size:12.5, bold:true, color:DEEP });
  plain(s, "RMSD 는 유사도를 따라 오르지만 충돌은 그렇지 않다 — 조립의 물리적 타당성이 유사도에 끌려가지 않는다.",
    { x:bx+0.28, y:5.98, w:bw-0.56, size:11, color:TEXT, ls:1.18 });
  s.addNotes("SF2 계통은 1SU6로 맞추면 안 된다는 점이 요지. 6OND를 별도 참조로 썼다.");
  foot(s);
}

/* ───────────────────────── S16 조립 검증 ───────────────────────── */
{
  const s = newSlide();
  head(s, "결정 ④  검증", "조립 절차는 도킹 점수를 왜곡하지 않는다");
  card(s, M, 2.05, 5.75, 2.1, TINT2, LINE);
  plain(s, "Wildtype 재도킹 — 파이프라인 일치 확인", { x:M+0.3, y:2.28, w:5.15, size:11.5, bold:true, color:GRAY });
  plain(s, "−4.380", { x:M+0.3, y:2.62, w:2.6, size:36, bold:true, color:DEEP });
  plain(s, "정본 대비  0.000", { x:M+3.0, y:2.95, w:2.45, size:13, bold:true, color:DEEP });
  plain(s, "kcal/mol", { x:M+0.3, y:3.42, w:2.6, size:10, color:GRAY });

  card(s, 6.85, 2.05, 5.8, 2.1, "FBEEE2", "E8C9A8");
  plain(s, "조립한 ChCODH2 이량체", { x:7.15, y:2.28, w:5.2, size:11.5, bold:true, color:GRAY });
  plain(s, "−4.30", { x:7.15, y:2.62, w:2.6, size:36, bold:true, color:ACCENT });
  plain(s, "차이  +0.080", { x:9.85, y:2.95, w:2.5, size:13, bold:true, color:ACCENT });
  plain(s, "kcal/mol", { x:7.15, y:3.42, w:2.6, size:10, color:GRAY });

  card(s, M, 4.40, CW, 1.15, TINT, LINE);
  plain(s, "조립 효과 0.080 은 결정구조 자신의 50런 범위 0.11 보다 작다.",
    { x:M+0.35, y:4.62, w:CW-0.7, size:17, bold:true, color:DEEP });
  plain(s, "즉, 조립에서 생긴 차이가 도킹 자체의 런 간 변동에 묻힌다.",
    { x:M+0.35, y:5.02, w:CW-0.7, size:11.5, color:TEXT });
  bullets(s, [
    [{text:"재현 검증도 통과했다. ", bold:true}, {text:"D63A_A559W 를 CE2 중심 npts 80 으로 재도킹해 기존 결과와 정렬 에너지 전량 일치(atol = 1e-3)를 확인한 뒤 배치를 시작했다."}],
    [{text:"기준 점수는 Wildtype ad_med_ds = −4.38, IQR 0.09. ", bold:true}, {text:"이후 모든 상동체 점수는 이 값을 기준으로 읽는다."}]
  ], { x:M, y:5.70, w:CW, h:1.3, size:11.5, gap:8 });
  s.addNotes("조립이 신호를 만들어내지 않았다는 확인. 이게 없으면 상동체 점수를 못 믿는다.");
  foot(s);
}

/* ───────────────────────── S17 다음 단계 ───────────────────────── */
{
  const s = newSlide();
  head(s, "진행 ⑤  다음 단계", "다음은 도킹 321개와 노이즈 바닥 측정이다");
  plain(s, "왜 노이즈 바닥이 필요한가 — 잡음원이 바뀌었다", { x:M, y:1.92, w:6.2, size:13, bold:true, color:TEXT });
  table(s, [
    hrow(["","변이체 288","상동체 321"]),
    drow(["골격","1SU6 하나 공유","구조마다 별개 예측"]),
    drow(["조립","불필요","구조마다 별개 조립"]),
    drow(["잡음원","런 간만 (IQR 0.09)","런 간 + 구조 예측 + 조립"],
      [{bold:true},null,{bold:true,color:ACCENT}])
  ], { x:M, y:2.32, w:6.2, colW:[1.2,2.4,2.6], rowH:0.40, size:10.5 });
  card(s, M, 4.10, 6.2, 1.0, "FBEEE2", "E8C9A8");
  plain(s, "288 세트의 SNR 을 그대로 옮겨 쓸 수 없다.", { x:M+0.28, y:4.28, w:5.64, size:12.5, bold:true, color:ACCENT });
  plain(s, "자기 잡음보다 작은 차이는 주장할 수 없다.", { x:M+0.28, y:4.62, w:5.64, size:11, color:TEXT });
  plain(s, "어떻게 재는가 — 확장 109개는 학습 클러스터의 2번째 멤버로, 대표와 짝을 이룬다.",
    { x:M, y:5.30, w:6.2, size:11.5, bold:true, color:TEXT, ls:1.2 });
  card(s, M, 5.78, 6.2, 0.95, TINT, LINE);
  plain(s, "within(c) = | med(대표) − med(2번째 멤버) |", { x:M+0.28, y:5.92, w:5.64, size:12, bold:true, color:DEEP });
  plain(s, "SNR = Var(클러스터 간) / mean(within²)", { x:M+0.28, y:6.28, w:5.64, size:12, bold:true, color:DEEP });

  const bx=7.35, bw=5.3;
  plain(s, "판정 기준은 미리 고정한다", { x:bx, y:1.92, w:bw, size:13, bold:true, color:TEXT });
  table(s, [
    hrow(["within 중앙값","SNR","판정"]),
    drow(["~ 0.05","40","라벨 그대로 사용"]),
    drow(["~ 0.30","7","288 과 비슷. 진행 가능"], [null,null,{bold:true,color:DEEP}]),
    drow(["~ 0.80","2.5","경계. 반복 조립 평균 필요"]),
    drow(["~ 1.50","1.3","라벨이 잡음에 묻힘. 재설계"], [null,null,{color:ACCENT}])
  ], { x:bx, y:2.32, w:bw, colW:[1.6,0.9,2.8], rowH:0.365, size:10.5 });
  plain(s, "D2 확장분을 함께 도킹하므로 추가 비용은 0 이다.",
    { x:bx, y:4.20, w:bw, size:11.5, bold:true, color:DEEP });
  plain(s, "일정", { x:bx, y:4.70, w:bw, size:13, bold:true, color:TEXT });
  const plan = ["도킹 321개  (core 212 + expansion 109 동시)","노이즈 바닥 측정 → SNR 판정",
    "공변량 회귀  struct_source · mode · 임계충돌 · pLDDT","라벨 부착  FASTA 헤더 fident → med 점수",
    "학습  A / B / C / D1 / D2"];
  plan.forEach((p,i) => {
    const y = 5.08 + i*0.36;
    badge(s, bx, y, 0.26, String(i+1), i===0?ACCENT:MID);
    plain(s, p, { x:bx+0.38, y:y+0.02, w:bw-0.38, size:10.5, color:TEXT });
  });
  s.addNotes("이 수치 없이는 상동체 순위를 예측했다는 주장을 방어할 수 없다. 심사자 질문 대비.");
  foot(s);
}

/* ───────────────────────── 백업 표지 ───────────────────────── */
{
  const s = newSlide(true);
  s.addText("BACKUP", { x:M, y:3.05, w:CW, h:0.35, fontFace:F, fontSize:13, bold:true,
    color:ACCENT, charSpacing:2.2, isTextBox:true, margin:0 });
  s.addText("백업 슬라이드", { x:M, y:3.45, w:CW, h:0.7, fontFace:F, fontSize:34, bold:true,
    color:WHITE, isTextBox:true, margin:0 });
  s.addText("B1 방향 의존성  ·  B2 주 지표 선정  ·  B3 단순 기술자 대비  ·  B4 고정 상수  ·  B5 한계  ·  B6 해석 범위(외부용)",
    { x:M, y:4.30, w:CW-0.6, h:0.4, fontFace:F, fontSize:12, color:INKSUB, isTextBox:true, margin:0 });
  foot(s, true);
}

/* ───────────────────────── B1 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B1", "반변 15 에서 순위는 회전에 무관하다");
  table(s, [
    hrow(["반변","회전 대조군","시드 대조군","격차","판정"]),
    drow(["15.000","0.905","0.927","−0.022","방향 무관"],
      [{bold:true},null,null,{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["18.000","0.897","0.979","−0.082","방향 의존"], [null,null,null,{color:ACCENT},{color:ACCENT}])
  ], { x:M, y:1.95, w:7.3, colW:[1.3,1.6,1.6,1.3,1.5], rowH:0.42, size:11 });
  plain(s, "med, 위치내 · 수용체를 그리드 중심 기준으로 회전(zyx 37/23/51°)시켜 재도킹",
    { x:M, y:3.42, w:7.3, size:10.5, italic:true, color:GRAY });
  bullets(s, [
    [{text:"시드 대조군이 반드시 필요하다. ", bold:true}, {text:"회전하면 탐색 경로 자체가 달라지므로, 동등한 조건에서 다른 탐색을 했을 때의 차이를 재는 자가 없으면 해석이 불가능하다."}],
    [{text:"최대차 0.520 은 무시할 수 없다. ", bold:true, color:ACCENT}, {text:"288개 표준편차 0.298 대비 1.7σ. 평균차가 0.005 라 편향은 없고 소수가 크게 튄다 — 순위는 견디고 절대값은 못 견딘다.", color:ACCENT}]
  ], { x:M, y:3.90, w:7.3, h:2.2, size:12, gap:10 });
  const bx=8.25, bw=4.4;
  card(s, bx, 1.95, bw, 4.4, TINT2, LINE);
  plain(s, "역설", { x:bx+0.3, y:2.15, w:bw-0.6, size:12, bold:true, color:ACCENT });
  plain(s, "재현성이 높다고 옳은 것이 아니다", { x:bx+0.3, y:2.48, w:bw-0.6, size:14, bold:true, color:TEXT, ls:1.15 });
  plain(s, "시드 대조군에서 반변 18(0.979)이 반변 15(0.927)보다 재현성이 높다. r18 은 포즈의 79% 가 한 군집에 박혀 있어 시드를 바꿔도 중앙값이 안 움직이기 때문이다. r15 는 여러 군집이 거의 같은 에너지로 경쟁하니 조금씩 흔들린다.",
    { x:bx+0.3, y:3.10, w:bw-0.6, size:11, color:TEXT, ls:1.2 });
  card(s, bx+0.3, 4.90, bw-0.6, 1.25, "FBEEE2", "E8C9A8");
  plain(s, "안정적으로 엉뚱한 것을 재는 것보다 조금 흔들리며 맞는 것을 재는 편이 낫다.",
    { x:bx+0.52, y:5.12, w:bw-1.04, size:11.5, bold:true, color:ACCENT, ls:1.2 });
  foot(s);
}

/* ───────────────────────── B2 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B2", "주 지표를 med 로 정한 근거는 박스 섭동 안정성이다");
  table(s, [
    hrow(["지표","ρ 위치내 (반변 15 → 18, CE2 중심)"]),
    drow(["med  (50런 중앙값)","0.918  /  0.943"], [{bold:true,color:DEEP},{bold:true,color:DEEP}]),
    drow(["best  (극단값 하나)","0.674  /  0.691"], [null,{color:ACCENT}])
  ], { x:M, y:2.05, w:7.0, colW:[2.8,4.2], rowH:0.45, size:12 });
  bullets(s, [
    [{text:"best 는 새 최소점이 생기면 즉시 갈아탄다. ", bold:true}, {text:"med 는 50런 앙상블의 강건 통계라 지형이 조금 바뀌어도 순위를 유지한다."}],
    [{text:"n = 27 로 지표를 고르면 과적합이다. ", bold:true}, {text:"실측 Km 상관으로 지표를 선택하는 대신 박스 섭동 안정성이라는 독립 기준으로 답했다."}],
    [{text:"med = −median(50런), IQR 병기. ", bold:true, color:DEEP}, {text:"모든 보고에서 산포를 함께 적는다.", color:DEEP}]
  ], { x:M, y:3.85, w:7.0, h:2.4, size:12.5, gap:12 });
  const bx=8.05, bw=4.6;
  card(s, bx, 2.05, bw, 2.35, TINT2, LINE);
  plain(s, "50런 에너지 범위", { x:bx+0.3, y:2.28, w:bw-0.6, size:11.5, bold:true, color:GRAY, align:"center" });
  plain(s, "0.22", { x:bx+0.3, y:2.62, w:bw-0.6, size:40, bold:true, color:DEEP, align:"center" });
  plain(s, "kcal/mol  (반변 15)", { x:bx+0.3, y:3.32, w:bw-0.6, size:11, color:GRAY, align:"center" });
  plain(s, "지형이 평평하다 — 극단값 하나를 믿을 근거가 없다.",
    { x:bx+0.3, y:3.72, w:bw-0.6, size:10.5, italic:true, color:GRAY, align:"center", ls:1.15 });
  card(s, bx, 4.60, bw, 1.65, TINT, LINE);
  plain(s, "실행 조건 (재현성 필수)", { x:bx+0.3, y:4.80, w:bw-0.6, size:11.5, bold:true, color:DEEP });
  plain(s, "--nrun 50 --seed 1\n--heuristics 0 --autostop 0", { x:bx+0.3, y:5.15, w:bw-0.6, size:12, color:TEXT, ls:1.25 });
  plain(s, "뒤 둘이 빠지면 런 수가 가변이 되어 재현성이 붕괴한다.", { x:bx+0.3, y:5.85, w:bw-0.6, size:10, color:GRAY });
  foot(s);
}

/* ───────────────────────── B3 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B3", "도킹이 단순 기술자를 앞서지만 n = 21 에서 차이는 확정적이지 않다");
  table(s, [
    hrow(["기술자","ρ (실측 Km)","상동체 적용"]),
    drow(["ad_med_lig  (도킹)","0.903","가능"], [{bold:true,color:DEEP},{bold:true,color:DEEP},null]),
    drow(["ad_best_ds  (도킹)","0.876","가능"]),
    drow(["ad_med_ds  (도킹)","0.856","가능"]),
    drow(["dchg  (변이 전하 변화)","0.815","불가"], [{bold:true,color:ACCENT},{bold:true,color:ACCENT},{color:ACCENT}]),
    drow(["cdock  (C-Docker)","0.719","가능"]),
    drow(["chg","0.594","불가"]),
    drow(["ev_dist","0.465","가능"]),
    drow(["dvol / vol","0.458 / 0.439","불가"])
  ], { x:M, y:1.95, w:6.8, colW:[3.0,2.0,1.8], rowH:0.365, size:10.5 });
  const bx=7.75, bw=4.9;
  bullets(s, [
    [{text:"dchg 는 상동체에 적용할 수 없다. ", bold:true}, {text:"변이체에만 정의되는 기술자이므로 족 전체로 확장하려면 도킹이 필요하다."}],
    [{text:"도킹이 전하 이상의 정보를 담는지는 미확인. ", bold:true, color:ACCENT}, {text:"dchg 를 통제한 뒤 도킹의 증분 설명력을 288 전수로 확인해야 한다 — 계산만 하면 되며 아직 수행하지 않았다.", color:ACCENT}]
  ], { x:bx, y:1.95, w:bw, h:2.6, size:12, gap:12 });
  card(s, bx, 4.55, bw, 1.85, TINT2, LINE);
  plain(s, "8X9F 기하 (참고값)", { x:bx+0.3, y:4.75, w:bw-0.6, size:11.5, bold:true, color:DEEP });
  plain(s, "체인 A 의 EV 무게중심이 그리드 중심에서 12.65 Å, 체인 A-2 는 16.89 Å. 두 EV 사이 28.03 Å 로 그리드 중심이 두 자리 사이에 놓인다. 1SU6 protomer 1 에 중첩 시 CA RMSD 0.28 Å.",
    { x:bx+0.3, y:5.10, w:bw-0.6, size:10.5, color:TEXT, ls:1.2 });
  foot(s);
}

/* ───────────────────────── B4 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B4", "고정 상수 — 세션 간 인터페이스");
  table(s, [
    hrow(["항목","값","근거"]),
    drow(["클러스터 임계","0.70","스윕 0.30–0.90 에서 선택. 초기 분석 219개 재현"]),
    drow(["박스 중심","각 구조의 41번 CA","CE2 는 상동체 11.7%(Ile 등)에 없음"], [null,{bold:true,color:DEEP},null]),
    drow(["박스 크기","npts 80 × 0.375 = 반변 15.000 Å","반변 18 이상에서 위치내 상관 0.93 → 0.78 붕괴"], [null,{bold:true,color:DEEP},null]),
    drow(["리간드","ev_ds.pdbqt","원자 16 / 전하합 1.998 / TORSDOF 3"]),
    drow(["실행","--nrun 50 --seed 1 --heuristics 0 --autostop 0","뒤 둘이 빠지면 런 수 가변 → 재현성 붕괴"]),
    drow(["주 지표","med (50런 중앙값)","박스 섭동에 best 보다 강건"]),
    drow(["기준 점수","Wildtype ad_med_ds = −4.38, IQR 0.09","engines_288_CA_r15.csv"]),
    drow(["시드","20260916","선정 · 분할 · 실험셋 재현"])
  ], { x:M, y:1.95, w:CW, colW:[1.9,4.3,5.8], rowH:0.42, size:10.5 });
  card(s, M, 5.95, CW, 1.05, "FBEEE2", "E8C9A8");
  plain(s, "리간드를 obabel / prepare_ligand4 로 재생성하면 전하합이 +2 → 0.000 이 되어 신호가 통째로 사라진다.",
    { x:M+0.32, y:6.12, w:CW-0.64, size:12, bold:true, color:ACCENT });
  plain(s, "출력 숫자는 정상 범위라 알아채기 어렵다. GPF 생성기 사용도 금지 — parameter_file None 줄이 autogrid4 SIGSEGV 를 일으킨다.",
    { x:M+0.32, y:6.48, w:CW-0.64, size:10.5, color:TEXT });
  foot(s);
}

/* ───────────────────────── B5 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B5", "명시할 한계 네 가지");
  const L = [
    ["Fe-S 클러스터 부재","FES(2Fe-2S)가 41번 CA 로부터 5.13 Å, 즉 결합 부위 바로 옆에 있으나 receptor_types 에 금속이 없어 두 세트 모두에서 빠져 있다. 비교는 내부적으로 일관되지만 정전기 환경은 불완전하다. AlphaFold 가 보조인자를 예측할 수 없어 대안이 없으며, 이 상태로도 실측 Km 상관 0.88 이 나온다는 사실은 순위 예측에 이 클러스터가 필수가 아님을 시사한다."],
    ["고유사도 상동체의 부재","90% 이상이 3~4개, 62.5 – 72.5% 는 공백이다. 모델은 처음부터 30 – 50% 영역에서 일해야 한다."],
    ["저유사도 구간의 구조 신뢰도","< 30 구간의 조립 RMSD 가 가장 부정확하고 클러스터도 10개뿐이라, 학습 곡선의 마지막 증분이 7개에 그친다."],
    ["구조 예측 방법의 혼재","321개 중 23개(7%)가 ColabFold 다. struct_source 로 기록했으며 도킹 후 mode · 임계충돌 · pLDDT 와 함께 회귀로 확인한다."]
  ];
  let y = 1.85;
  L.forEach((l,i) => {
    const h = i===0 ? 1.40 : 0.95;
    card(s, M, y, CW, h, i===0?"FBEEE2":TINT2, i===0?"E8C9A8":LINE);
    badge(s, M+0.28, y+0.22, 0.40, String(i+1), i===0?ACCENT:MID);
    plain(s, l[0], { x:M+0.82, y:y+0.20, w:CW-1.2, size:13, bold:true, color:i===0?ACCENT:TEXT });
    plain(s, l[1], { x:M+0.82, y:y+0.55, w:CW-1.2, size:11, color:TEXT, ls:1.2 });
    y += h + 0.15;
  });
  foot(s);
}

/* ───────────────────────── B6 ───────────────────────── */
{
  const s = newSlide();
  head(s, "BACKUP  B6  ·  외부 발표 전환용", "외부로 옮길 때는 도킹 점수의 해석 범위를 명시해야 한다");
  card(s, M, 1.95, 6.1, 0.95, TINT, LINE);
  plain(s, "내부 발표 — 직접 표현 사용", { x:M+0.28, y:2.10, w:5.54, size:11, bold:true, color:GRAY });
  plain(s, "“EV 결합부위에 도킹하여 결합친화도를 예측한다”",
    { x:M+0.28, y:2.42, w:5.54, size:12, bold:true, color:DEEP });
  card(s, 7.05, 1.95, 5.6, 0.95, "FBEEE2", "E8C9A8");
  plain(s, "외부 원고 · 발표 — 조정 필요", { x:7.33, y:2.10, w:5.04, size:11, bold:true, color:GRAY });
  plain(s, "“41번 주변 변이 영역을 훑는 상자에서 얻은 점수가 실측 Km 과 상관한다”",
    { x:7.33, y:2.42, w:5.04, size:11, bold:true, color:ACCENT, ls:1.15 });
  plain(s, "근거 — 이것은 자리 특이적 결합친화도의 직접 측정이 아니다", { x:M, y:3.20, w:CW, size:13, bold:true, color:TEXT });
  bullets(s, [
    [{text:"50런 에너지 범위가 반변 15 에서 0.22 kcal/mol 이다. ", bold:true}, {text:"에너지 지형이 평평하다."}],
    [{text:"박스를 키우면 항상 더 먼 곳에 더 좋은 자리가 나타난다. ", bold:true}, {text:"17.23 → 21.42 → 23.62 Å — 수렴하는 결합 자리가 없다."}],
    [{text:"8X9F 결정 EV 무게중심이 그리드 중심에서 12.65 Å 떨어져 있다.", bold:true}],
    [{text:"MDW 와 GLW 변이체 결정구조에서 EV 위치가 서로 크게 다르다. ", bold:true}, {text:"결정학적 EV 자리는 보존되지 않으므로 박스 설계의 기준이 될 수 없다."}]
  ], { x:M, y:3.60, w:CW, h:2.2, size:12, gap:11 });
  card(s, M, 5.95, CW, 0.95, TINT2, LINE);
  plain(s, "논문 제목의 enhanced electron mediator affinity 는 목표 서술이며, Results 에서는 순위 예측(ranking)으로 일관되게 기술한다.",
    { x:M+0.32, y:6.20, w:CW-0.64, size:11.5, bold:true, color:DEEP });
  foot(s);
}

const out = __dirname + "/20260916_progress_ko.pptx";
pres.writeFile({ fileName: out }).then(() => console.log("WROTE " + out + "  slides=" + pageNo));
