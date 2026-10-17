# Ni · 이량체 ablation — 점수가 올라가는가

**묻는 것 둘**

```
Δ(Ni)     = holo − apo            Ni 을 넣으면 신뢰도가 오르는가
Δ(dimer)  = 이량체 − 단량체        이량체가 필요한가
```

**절대 점수는 쓰지 않는다.** Ni 을 넣으면 토큰 수가 달라져 `confidence_score`,
`complex_plddt`, 전체 `ptm` 은 apo·holo 비교가 불가능하다. 양쪽에 똑같이 존재하는
**단백질 잔기** 값만 쓴다.

| 쓰는 값 | 왜 |
|---|---|
| 배위 잔기 pLDDT (CIF B-factor) | apo·holo 양쪽에 같은 잔기로 존재 |
| 배위 잔기끼리의 PAE 블록 | 금속 자리 = 네 잔기의 상대 기하. 가장 직접적 |
| `chains_ptm` 중 해당 사슬 | 사슬 단위라 비교 가능 |
| `pair_chains_iptm` A–B | **이량체 조건끼리만** |

**ChCODH2 는 넣지 않는다** (`WITH_CODH = False`). 후보의 Cys4 자리는 네 잔기 전부
후보 자신의 사슬에 있고, 이종 간 쌍이라 paired MSA 가 성립하지 않아 계면 점수가
교란된다. 원 도킹과 맞추고 싶으면 플래그만 켜면 된다.

---

## CELL A1 — 입력 생성

```python
# =============================================================================
# CELL A1 | ablation 입력 생성
#
#   조건 2 × 2 (단량체/이량체 × apo/holo) + 대조군
#   대조군이 핵심이다. CooC1 은 결정구조(3kji)에서 금속이 이량체 계면에 걸쳐
#   있는 것을 안다. 거기서 Δ(dimer) 가 안 나오면 이 측정 자체를 못 쓴다.
# =============================================================================
from pathlib import Path
import json, re, shutil, datetime

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
FAA   = {"BL21":   TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"bl21_db_match_qjz.faa",
         "MG1655": TOOLS/"database"/"protein_list"/"mg1655_protein.faa",
         "Y19":    TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"y19_db_match.faa"}

ABL   = BASE/"result"/"ablation"; ABL.mkdir(parents=True, exist_ok=True)
INP   = ABL/"inputs"
NREP  = 5                      # 조건당 반복. 1 이면 Δ 를 읽을 수 없다
WITH_CODH = False              # True 면 ChCODH2 를 사슬로 추가

# ── 대상 정의 ────────────────────────────────────────────────────────────
#   site = 배위 잔기 번호 (1-based). CIF 의 label_seq_id 와 맞아야 한다.
TARGETS = {
    "QJZ12568.1": {"faa": "BL21", "site": [99, 102, 112, 115], "role": "candidate"},
    # 양성대조군 — 서열은 아래 SEQ_OVERRIDE 로 직접 넣어도 된다
    "CooC1":      {"faa": None,   "site": [112, 114],           "role": "positive"},
    # 음성대조군 — Cys 없는 비슷한 길이. 아래에서 자동으로 고른다
    "NEGCTRL":    {"faa": "BL21", "site": [],                   "role": "negative"},
}

SEQ_OVERRIDE = {
    # "CooC1": "MSK...",     ← ChCooC1 서열을 여기에 붙여넣는다
    # "ChCODH2": "M...",
}

# ── FASTA 에서 서열 읽기 ─────────────────────────────────────────────────
def read_faa(path):
    seqs, k, buf = {}, None, []
    for ln in open(path):
        if ln.startswith(">"):
            if k: seqs[k] = "".join(buf)
            k, buf = ln[1:].split()[0], []
        else:
            buf.append(ln.strip())
    if k: seqs[k] = "".join(buf)
    return seqs

FAA_CACHE = {}
def get_seq(name, which):
    if name in SEQ_OVERRIDE: return SEQ_OVERRIDE[name]
    if which is None:
        raise SystemExit(f"[중단] {name} 서열이 없다. SEQ_OVERRIDE 에 넣을 것")
    if which not in FAA_CACHE: FAA_CACHE[which] = read_faa(FAA[which])
    d = FAA_CACHE[which]
    for k, v in d.items():
        if k.split("|")[0] == name or name in k: return v
    raise SystemExit(f"[중단] {FAA[which].name} 에서 {name} 을 못 찾았다")

# 음성대조군 자동 선정 — Cys 0개, 길이 ±20%
def pick_negctrl(ref_len):
    if "BL21" not in FAA_CACHE: FAA_CACHE["BL21"] = read_faa(FAA["BL21"])
    lo, hi = int(ref_len * 0.8), int(ref_len * 1.2)
    for k, v in FAA_CACHE["BL21"].items():
        if lo <= len(v) <= hi and v.count("C") == 0:
            return k.split("|")[0], v
    raise SystemExit("[중단] Cys 없는 동일 길이 단백질을 못 찾았다")

# ── 서열 확정 ────────────────────────────────────────────────────────────
SEQ = {}
for name, cfg in TARGETS.items():
    if name == "NEGCTRL": continue
    SEQ[name] = get_seq(name, cfg["faa"])
    print(f"  {name:14s} {len(SEQ[name]):4d} aa  Cys {SEQ[name].count('C')}")
    for r in cfg["site"]:
        if r > len(SEQ[name]):
            raise SystemExit(f"[중단] {name} 길이 {len(SEQ[name])} < 배위 잔기 {r}")
        if SEQ[name][r-1] not in "CH":
            print(f"    ⚠ {name} {r} 번이 {SEQ[name][r-1]} 다 (Cys/His 가 아님) — 번호 확인")

_nm, _sq = pick_negctrl(len(SEQ["QJZ12568.1"]))
SEQ["NEGCTRL"] = _sq
print(f"  NEGCTRL      = {_nm}  {len(_sq)} aa  Cys 0")
TARGETS["NEGCTRL"]["real_id"] = _nm

CODH = SEQ_OVERRIDE.get("ChCODH2") if WITH_CODH else None
if WITH_CODH and not CODH:
    raise SystemExit("[중단] WITH_CODH=True 인데 ChCODH2 서열이 없다")

# ── YAML 쓰기 ────────────────────────────────────────────────────────────
def write_yaml(path, chains, with_ni):
    """chains = [(id, seq), ...]. 리간드는 항상 맨 뒤 — PAE 토큰 순서가 예측 가능해진다."""
    L = ["version: 1", "sequences:"]
    for cid, s in chains:
        L += [f"  - protein:", f"      id: {cid}", f"      sequence: {s}"]
    if with_ni:
        L += ["  - ligand:", "      id: L", "      ccd: NI"]
    path.write_text("\n".join(L) + "\n")

# 입력 폴더는 통째로 새로 만든다. 남은 파일이 섞이면 결과가 오염된다.
if INP.exists():
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.move(str(INP), str(ABL/f"inputs_old_{ts}"))
    print(f"  기존 inputs → inputs_old_{ts}")
INP.mkdir(parents=True)

PLAN = {}
for name, cfg in TARGETS.items():
    s = SEQ[name]
    for asm in ("mono", "dimer"):
        for lig in ("apo", "holo"):
            chains = [("A", s)] + ([("B", s)] if asm == "dimer" else [])
            if WITH_CODH: chains.append(("Z", CODH))
            tag = f"{name}__{asm}__{lig}"
            write_yaml(INP/f"{tag}.yaml", chains, lig == "holo")
            PLAN[tag] = {"protein": name, "assembly": asm, "ligand": lig,
                         "site": cfg["site"], "role": cfg["role"],
                         "n_chain": len(chains), "len": len(s),
                         "chain_ids": [c for c, _ in chains]}

# manifest 는 입력 폴더 **밖**에 쓴다. 안에 두면 boltz 가 파싱하려다 죽는다.
(ABL/"ablation_plan.json").write_text(json.dumps(PLAN, indent=2, ensure_ascii=False))

print(f"\n  입력 {len(PLAN)} 개 → {INP}")
print(f"  계획  → {ABL/'ablation_plan.json'}")
print(f"  반복 {NREP} 회/조건  (diffusion_samples)")
print("\n  조건:")
for t in sorted(PLAN): print("   ", t)
```

---

## CELL A2 — 실행

```python
# =============================================================================
# CELL A2 | 백그라운드 실행
#   GPU 가 비었는지 확인하고, 다른 boltz 가 돌고 있으면 물러난다.
# =============================================================================
import subprocess, os, time, datetime
from pathlib import Path

OUT   = ABL/"out"; OUT.mkdir(exist_ok=True)
LOG   = ABL/"ablation.log"
GPU   = 0
MIN_FREE = 12000          # MiB

def gpu_free(i):
    q = subprocess.run(["nvidia-smi", f"--id={i}", "--query-gpu=memory.free",
                        "--format=csv,noheader,nounits"], capture_output=True, text=True)
    return int(q.stdout.strip().split("\n")[0])

def other_boltz():
    q = subprocess.run(["pgrep", "-af", "boltz predict"], capture_output=True, text=True)
    return [l for l in q.stdout.strip().split("\n") if l]

busy = other_boltz()
if busy:
    print("  다른 boltz 가 돌고 있다. 끝나고 다시 실행할 것:")
    for b in busy: print("   ", b[:110])
elif gpu_free(GPU) < MIN_FREE:
    print(f"  GPU {GPU} 여유 {gpu_free(GPU)} MiB < {MIN_FREE} — 기다릴 것")
else:
    # 로그는 매번 새로 연다. append 로 열면 지난 실행의 오류를 현재 것으로 읽게 된다.
    if LOG.exists(): LOG.rename(LOG.with_suffix(
        f".{datetime.datetime.now():%Y%m%d_%H%M%S}.log"))
    cmd = (f'boltz predict "{INP}" --out_dir "{OUT}" --use_msa_server '
           f'--diffusion_samples {NREP} --recycling_steps 3 --override')
    env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(GPU),
               PYTORCH_CUDA_ALLOC_CONF="max_split_size_mb:512")
    p = subprocess.Popen(["setsid", "bash", "-lc", f'{cmd} > "{LOG}" 2>&1'],
                         env=env, start_new_session=True)
    print(f"  시작 PID {p.pid}  GPU {GPU}  샘플 {NREP}")
    print(f"  로그 {LOG}")

# ── 진행 확인 ────────────────────────────────────────────────────────────
def running():
    """setsid 쉘이 아니라 실제 boltz 자식이 살아 있는지 본다."""
    return bool(other_boltz())

def progress():
    pred = next(iter(OUT.glob("boltz_results_*/predictions")), None)
    done = len(list(pred.iterdir())) if pred else 0
    print(f"  {'실행 중' if running() else '종료'}  완료 {done}/{len(PLAN)}")
    if LOG.exists():
        for ln in LOG.read_text(errors="ignore").strip().split("\n")[-4:]:
            print("   ", ln[:120])

progress()
```

---

## CELL A3 — 집계

```python
# =============================================================================
# CELL A3 | Δ 계산
#
#   pLDDT 는 npz 가 아니라 **CIF 의 B-factor** 에서 읽는다. 리간드가 있으면
#   토큰 수가 달라져 npz 인덱스가 밀리는데, CIF 는 사슬·잔기번호로 직접 고를
#   수 있어 그 문제가 없다.
# =============================================================================
import json, math, numpy as np, pandas as pd
from pathlib import Path
from collections import defaultdict

PLAN = json.loads((ABL/"ablation_plan.json").read_text())
PRED = next(iter(ABL.glob("out/boltz_results_*/predictions")), None)
if PRED is None: raise SystemExit("[중단] 예측 폴더가 없다. A2 가 끝났는지 확인")

def cif_rows(p):
    """(chain, resnum, atom, element, x, y, z, bfactor) — label_asym_id 기준."""
    rows, cols, inloop = [], {}, False
    for ln in open(p):
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols[s.split(".")[1]] = len(cols); inloop = True; continue
        if inloop and (s.startswith("ATOM") or s.startswith("HETATM")):
            f = s.split()
            g = lambda k: f[cols[k]] if k in cols else "."
            try: rn = int(g("label_seq_id"))
            except ValueError: rn = -1
            rows.append((g("label_asym_id"), rn, g("label_atom_id"),
                         g("type_symbol").upper(),
                         float(g("Cartn_x")), float(g("Cartn_y")), float(g("Cartn_z")),
                         float(g("B_iso_or_equiv"))))
        elif inloop and s.startswith("#"):
            inloop = False
    return rows

def residue_plddt(rows, chain):
    """잔기당 pLDDT — CA 의 B-factor."""
    return {rn: b for ch, rn, at, el, x, y, z, b in rows
            if ch == chain and at == "CA"}

def pae_block(npz_path, idx, ntok_expect):
    if not npz_path.exists(): return None
    m = np.load(npz_path)
    k = "pae" if "pae" in m else m.files[0]
    P = m[k]
    if P.ndim == 3: P = P[0]
    if P.shape[0] < ntok_expect: return None          # 토큰 수가 안 맞으면 포기
    sub = P[np.ix_(idx, idx)]
    off = ~np.eye(len(idx), dtype=bool)
    return float(sub[off].mean()) if len(idx) > 1 else None

rows_out = []
for tag, cfg in PLAN.items():
    d = PRED/tag
    if not d.is_dir():
        print(f"  ⚠ {tag} 결과 없음"); continue
    site = cfg["site"]
    for cif in sorted(d.glob(f"{tag}_model_*.cif")):
        rep = int(cif.stem.rsplit("_", 1)[1])
        R   = cif_rows(cif)
        pl  = residue_plddt(R, "A")
        if not pl:
            print(f"  ⚠ {cif.name} 사슬 A 가 비었다"); continue

        rec = {"tag": tag, "rep": rep, "protein": cfg["protein"],
               "assembly": cfg["assembly"], "ligand": cfg["ligand"],
               "role": cfg["role"], "len": cfg["len"]}

        if site:
            miss = [r for r in site if r not in pl]
            if miss:
                print(f"  ⚠ {cif.name} 잔기 {miss} 없음 — 번호 확인"); continue
            rest = [r for r in pl if r not in site]
            rec["plddt_site"] = float(np.mean([pl[r] for r in site]))
            rec["plddt_rest"] = float(np.mean([pl[r] for r in rest]))
            # PAE: 단백질 사슬을 YAML 순서대로 썼으므로 토큰 = 잔기 순서
            idx = [r - 1 for r in site]
            rec["pae_site"] = pae_block(d/f"pae_{tag}_model_{rep}.npz", idx, max(idx)+1)
        else:
            rec["plddt_site"] = rec["plddt_rest"] = rec["pae_site"] = np.nan
            rec["plddt_rest"] = float(np.mean(list(pl.values())))

        # 동형이량체 대칭 — 같은 서열 두 사슬이 같은 값을 받아야 한다
        if cfg["assembly"] == "dimer" and site:
            plB = residue_plddt(R, "B")
            if all(r in plB for r in site):
                b = float(np.mean([plB[r] for r in site]))
                rec["plddt_site_B"] = b
                rec["chain_asym"]   = abs(rec["plddt_site"] - b)

        # Ni 배위 — holo 에서만
        if cfg["ligand"] == "holo":
            ni = [r for r in R if r[3] == "NI"]
            if ni:
                _, _, _, _, nx, ny, nz, _ = ni[0]
                ds = [(math.dist((nx, ny, nz), (x, y, z)), ch)
                      for ch, rn, at, el, x, y, z, bb in R
                      if el == "S" and math.dist((nx, ny, nz), (x, y, z)) <= 3.0]
                rec["ni_min_dist"]   = round(min(d_ for d_, _ in ds), 2) if ds else np.nan
                rec["ni_n_S"]        = len(ds)
                rec["ni_both_chain"] = len({c for _, c in ds}) > 1

        # confidence json
        cj = d/f"confidence_{tag}_model_{rep}.json"
        if cj.exists():
            c = json.loads(cj.read_text())
            cp = c.get("chains_ptm") or {}
            rec["chain_ptm"] = cp.get("A", cp.get("0"))
            pc = c.get("pair_chains_iptm") or {}
            rec["pair_iptm"] = (pc.get("A") or pc.get("0") or {}).get("B") \
                            or (pc.get("A") or pc.get("0") or {}).get("1")
        rows_out.append(rec)

D = pd.DataFrame(rows_out)
if D.empty: raise SystemExit("[중단] 읽은 결과가 없다")
D.to_csv(ABL/"ablation_raw.csv", index=False)
print(f"  원자료 {len(D)} 행 → {ABL/'ablation_raw.csv'}\n")

# ── 조건별 평균 ± 표준편차 ───────────────────────────────────────────────
MET = [c for c in ["plddt_site", "plddt_rest", "pae_site", "chain_ptm",
                   "pair_iptm", "chain_asym", "ni_min_dist"] if c in D]
G = D.groupby(["protein", "assembly", "ligand"])[MET].agg(["mean", "std", "count"])
print("=" * 100); print("조건별"); print("=" * 100)
print(G.round(3).to_string())

# ── Δ 계산 ───────────────────────────────────────────────────────────────
def cell(p, a, l, m):
    try:  return float(G.loc[(p, a, l), (m, "mean")])
    except KeyError: return np.nan
def sd(p, a, l, m):
    try:  return float(G.loc[(p, a, l), (m, "std")])
    except KeyError: return np.nan

out = []
for p in D.protein.unique():
    for a in ("mono", "dimer"):
        dS = cell(p, a, "holo", "plddt_site") - cell(p, a, "apo", "plddt_site")
        dR = cell(p, a, "holo", "plddt_rest") - cell(p, a, "apo", "plddt_rest")
        pooled = math.hypot(sd(p, a, "holo", "plddt_site") or 0,
                            sd(p, a, "apo",  "plddt_site") or 0)
        out.append({"protein": p, "비교": f"Δ(Ni) · {a}",
                    "Δ자리": round(dS, 2), "Δ나머지": round(dR, 2),
                    "DiD": round(dS - dR, 2),
                    "반복분산": round(pooled, 2),
                    "판정": "신호" if (dS - dR) > pooled else "잡음 수준"})
    for l in ("holo", "apo"):
        dS = cell(p, "dimer", l, "plddt_site") - cell(p, "mono", l, "plddt_site")
        dR = cell(p, "dimer", l, "plddt_rest") - cell(p, "mono", l, "plddt_rest")
        pooled = math.hypot(sd(p, "dimer", l, "plddt_site") or 0,
                            sd(p, "mono",  l, "plddt_site") or 0)
        out.append({"protein": p, "비교": f"Δ(dimer) · {l}",
                    "Δ자리": round(dS, 2), "Δ나머지": round(dR, 2),
                    "DiD": round(dS - dR, 2),
                    "반복분산": round(pooled, 2),
                    "판정": "신호" if (dS - dR) > pooled else "잡음 수준"})

S = pd.DataFrame(out)
S.to_csv(ABL/"ablation_delta.csv", index=False)
print("\n" + "=" * 100); print("Δ  (DiD = Δ자리 − Δ나머지)"); print("=" * 100)
print(S.to_string(index=False))

# ── 대조군 판정 ──────────────────────────────────────────────────────────
print("\n" + "=" * 100); print("대조군"); print("=" * 100)
ctl = S[(S.protein == "CooC1") & (S["비교"] == "Δ(dimer) · holo")]
if len(ctl):
    r = ctl.iloc[0]
    if r["DiD"] > r["반복분산"]:
        print(f"  ✅ CooC1 Δ(dimer) DiD = {r['DiD']} > 분산 {r['반복분산']}")
        print("     지표가 작동한다. 후보 결과를 해석해도 된다.")
    else:
        print(f"  ❌ CooC1 Δ(dimer) DiD = {r['DiD']} ≤ 분산 {r['반복분산']}")
        print("     정답을 아는 케이스에서 신호가 안 난다.")
        print("     → 이 측정으로는 이량체 필요 여부를 판정할 수 없다. 후보 결과를 쓰지 말 것.")
else:
    print("  CooC1 결과가 없다 — 서열을 SEQ_OVERRIDE 에 넣고 A1 부터 다시")

neg = S[(S.protein == "NEGCTRL")]
if len(neg) and (neg.DiD.abs() > neg["반복분산"]).any():
    print("  ⚠ 음성대조군(Cys 없음)에서도 신호가 난다 — 지표가 Ni 과 무관한 것을 재고 있다")

# 동형이량체 대칭
if "chain_asym" in D:
    A = D[D.assembly == "dimer"].groupby("protein").chain_asym.mean().dropna()
    if len(A):
        print("\n  동형이량체 사슬 간 pLDDT 차 (작을수록 좋다):")
        for p, v in A.items():
            print(f"    {p:14s} {v:6.2f}" + ("   ⚠ 두 사슬이 비대칭" if v > 5 else ""))
```

---

## 판정

| 결과 | 뜻 |
|---|---|
| CooC1 `Δ(dimer)` DiD > 반복분산 | **지표 작동.** 후보 해석 가능 |
| CooC1 에서 신호 없음 | **여기서 멈춘다.** 후보를 돌려도 해석할 수 없다 |
| 후보 `Δ(Ni)` DiD > 0 | 모델이 Ni 을 원한다 — 자리가 실재할 근거 |
| 후보 `Δ(dimer)` DiD > 0 | **이량체가 필요하다** — 재예측을 진행할 근거 |
| 후보 `Δ(dimer)` ≈ 0 | 이량체 재예측을 **안 해도 된다** |
| `ni_both_chain = True` | Ni 이 두 사슬에 걸쳐 배위 — CooC1 형 자리 |
| 음성대조군에서도 신호 | 지표가 Ni 과 무관한 것을 재고 있다 |

**반복분산보다 작은 Δ 는 읽지 않는다.** 같은 입력에 ipTM 이 0.102~0.309 벌어진
계다. `NREP = 1` 로 돌리면 이 표의 모든 판정이 무의미하다.

## 실행 전 채울 것

1. `SEQ_OVERRIDE["CooC1"]` — ChCooC1 서열. **이게 없으면 대조군이 없고, 대조군이
   없으면 이 실험 전체가 무의미하다**
2. `TARGETS["CooC1"]["site"]` — 기본 `[112, 114]`. 실제 서열 번호와 맞는지 A1 의
   경고 줄을 볼 것
3. `WITH_CODH` — 원 도킹 구성과 맞추려면 True, 그러면 `SEQ_OVERRIDE["ChCODH2"]` 도 필요
