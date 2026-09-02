"""
CODH/EV docking verification pipeline  (288 variants x 2 engines x 2 ligands)

Reproduces the run executed on 2026-09-02 on the ChCODH2 A559W variant set.
Each CELL block below was run as one Jupyter cell, in order.

Fixed conditions (recovered from the produced artifacts, not from prose):
    center      (3.266, -25.350, -3.469)   = F41 CE2 of the first protomer
    AutoDock    npts 80^3 x spacing 0.375  -> 30 A cube, identical for all 288
                autodock_gpu_128wi (AutoDock-GPU, AD4 scoring function)
                --nrun 50 --seed 1 --heuristics 0 --autostop 0
                ** --heuristics 0 --autostop 0 must not be dropped: without them
                   the run count is chosen adaptively and seed reproducibility
                   and the fixed 50-run sample are both lost. **
    RxDock      SphereSiteMapper radius 15 / small 1.5 / large 6.0, grid-step 0.5,
                min-cav-vol 50, max-cav 1;  rxcmd dock -p dock.json -n 50 -s 1
                The 15 A here is a cavity-detection sphere, NOT the search box:
                the search space is the detected cavity, which varies per
                structure (371-1723 A^3 across the 288).  Always record cav_vol.
    receptor    ATOM only (all HETATM stripped: SF4/FES/NFS/HOH)
                -> pdb2pqr --ff=AMBER --with-ph=8.0 --titration-state-method=propka
                -> RxDock:   obabel {pqr} -O {mol2} --partialcharge gasteiger
                -> AutoDock: prepare_receptor -U nphs -A None  (flags immaterial
                             here; all candidate flag sets gave identical output)
    ligand      ** Meeko only.  obabel and prepare_ligand4 both yield a PDBQT with
                   charge sum 0.000 instead of +2, silently destroying the signal
                   that explains most of the score.  Never regenerate the ligand
                   with them. **
                ev_ds  = oxEV_prepared.sd  (Discovery Studio, the C-Docker ligand)
                ev_lig = crystal reference conformer (ev_ref.sd, planar, no H)

Every stage was validated by reproducing an existing artifact byte-for-byte or
score-for-score before the batch was launched (CELL 5, CELL 7).
"""

# ============================================================ CELL 1 : settings
import os, re, csv, json, glob, shutil, hashlib, subprocess, tempfile, time
import numpy as np
from pathlib import Path

BASE  = Path("/mnt/af2results/mingyu/workspace/docking_simulation")
IN    = BASE / "input"
WORK  = IN / "veri_engines"          # RxDock workspace (flat; rxcmd resolves from cwd)
PREP  = WORK / "prep"                # receptor prep + legacy AutoDock output
VAR   = IN / "1su6_varient"          # 288 structures
EVSRC = IN / "ev_source"
RES   = BASE / "result" / "veri_engines"
KMX   = IN / "variant_km_real.xlsx"  # measured Km / kcat, 21 variants
RES.mkdir(parents=True, exist_ok=True)

SOFT   = Path("/mnt/af2results/mingyu/soft/docking")
ADFR   = SOFT / "ADFRsuite/bin"
PREPR  = str(ADFR / "prepare_receptor")
AUTOGR = str(ADFR / "autogrid4")
ADGPU  = str(SOFT / "AutoDock-GPU/bin/autodock_gpu_128wi")
MK     = shutil.which("mk_prepare_ligand.py")
TMP    = Path(tempfile.mkdtemp(prefix="dock_"))

CENTER     = (3.266, -25.350, -3.469)
NPTS, SPC  = 80, 0.375
NRUN, SEED = 50, 1
REC_TYPES  = "A C HD N NA OA SA"
LIG_TYPES  = "A C N"
DEVFLAG    = ["--devnum"]                          # 2x RTX 3090
LIGS       = {"ds": WORK / "ev_ds.pdbqt", "lig": WORK / "ev_lig.pdbqt"}
LIGS_SD    = {"ds": "ev_ds.sd",           "lig": "ev_lig.sd"}
RXOUT      = {"ds": "RXDS",               "lig": "RXLG"}

GPF = """npts {n} {n} {n}
gridfld {tag}.maps.fld
spacing {s}
receptor_types {rt}
ligand_types {lt}
receptor {tag}.pdbqt
gridcenter {cx:.3f} {cy:.3f} {cz:.3f}
smooth 0.5
map {tag}.A.map
map {tag}.C.map
map {tag}.N.map
elecmap {tag}.e.map
dsolvmap {tag}.d.map
dielectric -0.1465
"""

_ERE  = re.compile(r"Estimated Free Energy of Binding\s*=\s*([-\d.]+)")
_SPAT = re.compile(r"^>\s*<rxdock\.score>")      # NOT "SCORE"; the field is rxdock.score


def dlg_energies(p):
    return [float(x) for x in _ERE.findall(open(p, errors="ignore").read())]


def sd_scores(p):
    t, out = open(p, errors="ignore").read().splitlines(), []
    for i, l in enumerate(t):
        if _SPAT.match(l.strip()):
            try:
                out.append(float(t[i + 1].strip()))
            except Exception:
                pass
    return out


def cav_vol(tag):
    f = WORK / ("RX_%s-docking-site.json" % tag)
    if not f.exists():
        return None
    try:
        c = json.load(open(f))["docking-site"]["cavities"][0]["coordinates"]
        return len(c) * 0.5 ** 3
    except Exception:
        return None


# ============================================================ CELL 2 : targets
MUT = re.compile(r"^([A-Z])(\d+)([A-Z])_A559W$")
rows = []
for p in sorted(VAR.rglob("*.pdb")):
    m = MUT.match(p.stem)
    rows.append(dict(tag=p.stem, folder=p.parent.name, path=str(p),
                     pos=(m.group(1) + m.group(2)) if m else "BG",
                     wt=m.group(1) if m else "", mut=m.group(3) if m else ""))
TAGS = [r["tag"] for r in rows]
with open(RES / "manifest.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
# 288 = 15 positions x 19 (+ N59N) + A559W background + Wildtype


# ============================================================ CELL 3 : ligands
def build_ligands():
    """Meeko is mandatory. .sd must be handed over as .sdf (extension check)."""
    DS = EVSRC / "oxEV_prepared.sd"
    shutil.copy(str(DS), str(WORK / "ev_ds.sd"))          # legacy charge column already present
    sdf = TMP / "ev_ds.sdf"; shutil.copy(str(DS), str(sdf))
    out = TMP / "ev_ds.pdbqt"
    subprocess.run([MK, "-i", str(sdf), "-o", str(out)], capture_output=True, text=True)
    shutil.copy(str(out), str(WORK / "ev_ds.pdbqt"))


def pq_sig(p):
    at = [l for l in open(p) if l.startswith(("ATOM", "HETATM"))]
    td = [l.split()[1] for l in open(p) if l.startswith("TORSDOF")]
    return (len(at), round(sum(float(l[70:76]) for l in at), 3),
            sorted({l[77:79].strip() for l in at}), td[0] if td else "-")
# both ligands: (16, 1.998, ['A','C','N'], '3')  -- identical except the conformer


# ============================================================ CELL 4 : receptors
def atom_md5(path):
    h = hashlib.md5()
    for l in open(path, "rb"):
        if l[:4] == b"ATOM":
            h.update(l[12:54])
    return h.hexdigest()


def prep_one(rec, force=False):
    """strip -> pdb2pqr -> {mol2 for RxDock, pdbqt for AutoDock}.
    Staleness is decided by comparing the strip file against the current source,
    so a receptor prepared from an older copy of the structure is caught."""
    tag, src = rec["tag"], rec["path"]
    strip, pqr = PREP / (tag + "_strip.pdb"), PREP / (tag + ".pqr")
    ph8, pdbqt = PREP / (tag + "_ph8.pdb"), PREP / (tag + ".pdbqt")
    mol2 = WORK / (tag + ".mol2")
    if not force and all(p.exists() for p in (strip, pqr, ph8, pdbqt, mol2)) \
       and atom_md5(strip) == atom_md5(src):
        return tag, "ok(reuse)", ""
    with open(src) as fi, open(strip, "w") as fo:
        for l in fi:
            if l.startswith("ATOM"):
                fo.write(l)
        fo.write("END\n")
    r = subprocess.run(["pdb2pqr", "--ff=AMBER", "--with-ph=8.0",
                        "--titration-state-method=propka",
                        "--pdb-output", str(ph8), str(strip), str(pqr)],
                       capture_output=True, text=True)
    if not pqr.exists():
        return tag, "fail_pdb2pqr", (r.stderr or r.stdout)[-200:]
    subprocess.run(["obabel", str(pqr), "-O", str(mol2),
                    "--partialcharge", "gasteiger"], capture_output=True)
    if not mol2.exists():
        return tag, "fail_obabel", ""
    r = subprocess.run([PREPR, "-r", ph8.name, "-o", pdbqt.name, "-U", "nphs", "-A", "None"],
                       cwd=str(PREP), capture_output=True, text=True)
    if not pdbqt.exists():
        return tag, "fail_prepare_receptor", (r.stderr or r.stdout)[-200:]
    return tag, "regenerated", ""


# ==================================================== CELL 6 : AutoDock batch
def dock_one(tag, devnum):
    """One grid per structure, shared by both ligands (ligand_types are identical).
    Maps are deleted immediately so disk stays flat."""
    outs = {k: PREP / k / (tag + ".dlg") for k in LIGS}
    if all(p.exists() and len(dlg_energies(p)) == NRUN for p in outs.values()):
        return tag, "skip", ""
    sc = PREP / "_grid" / tag
    shutil.rmtree(sc, ignore_errors=True); sc.mkdir(parents=True)
    try:
        shutil.copy(str(PREP / (tag + ".pdbqt")), str(sc / (tag + ".pdbqt")))
        (sc / (tag + ".gpf")).write_text(GPF.format(
            n=NPTS, s=SPC, rt=REC_TYPES, lt=LIG_TYPES, tag=tag,
            cx=CENTER[0], cy=CENTER[1], cz=CENTER[2]))
        r = subprocess.run([AUTOGR, "-p", tag + ".gpf", "-l", tag + ".glg"],
                           cwd=str(sc), capture_output=True, text=True)
        if not (sc / (tag + ".maps.fld")).exists():
            return tag, "fail_autogrid", (r.stderr or r.stdout)[-200:]
        dev = ["--devnum", str(devnum)] if DEVFLAG else []
        for k, lp in LIGS.items():
            r = subprocess.run([ADGPU, "--ffile", tag + ".maps.fld", "--lfile", str(lp),
                                "--resnam", k, "--nrun", str(NRUN), "--seed", str(SEED),
                                "--heuristics", "0", "--autostop", "0"] + dev,
                               cwd=str(sc), capture_output=True, text=True)
            if not (sc / (k + ".dlg")).exists():
                return tag, "fail_dock_" + k, (r.stderr or r.stdout)[-200:]
            shutil.move(str(sc / (k + ".dlg")), str(outs[k]))
        return tag, "ok", ""
    finally:
        shutil.rmtree(sc, ignore_errors=True)


# ====================================================== CELL 8 : RxDock batch
RXJSON = """{
  "media-type": "application/vnd.rxdock.parameters",
  "title": "RX_%(tag)s",
  "version": "0.1.0",
  "sections": ["receptor", "mapper", "cavity"],
  "receptor": {"file": "%(tag)s.mol2"},
  "mapper": {
    "site-mapper": "SphereSiteMapper",
    "center": "(%(cx).3f,%(cy).3f,%(cz).3f)",
    "radius": 15.0,
    "small-sphere-radius": 1.5,
    "large-sphere-radius": 6.0,
    "excluded-volume-radius-increment": 0,
    "grid-step": 0.5,
    "minimum-cavity-volume": 50,
    "maximum-cavities": 1
  },
  "cavity": {"scoring-function": "CavityGridSF", "weight": 1}
}
"""


def rx_one(tag):
    """The cavity does not depend on the ligand, so it is detected once and
    reused by both docking runs.  No -H: RxDock reads polar hydrogens only by
    default, which puts the 18-H DS ligand and the 0-H reference ligand on the
    same heavy-atom footing, leaving the conformer as the only difference."""
    prm  = WORK / ("RX_%s.json" % tag)
    site = WORK / ("RX_%s-docking-site.json" % tag)
    outs = {k: WORK / ("%s_%s_out.sd" % (RXOUT[k], tag)) for k in LIGS_SD}
    if all(p.exists() and len(sd_scores(p)) == NRUN for p in outs.values()):
        return tag, "skip", ""
    if not prm.exists():
        prm.write_text(RXJSON % dict(tag=tag, cx=CENTER[0], cy=CENTER[1], cz=CENTER[2]))
    if not site.exists():
        r = subprocess.run(["rxcmd", "cavity-search", "-r", prm.name, "-W"],
                           cwd=str(WORK), capture_output=True, text=True)
        if not site.exists():
            return tag, "fail_cavity", (r.stderr or r.stdout)[-200:]
    for k, sdname in LIGS_SD.items():
        o = outs[k]
        if o.exists() and len(sd_scores(o)) == NRUN:
            continue
        r = subprocess.run(["rxcmd", "dock", "-i", sdname, "-o", o.name,
                            "-r", prm.name, "-p", "dock.json",
                            "-n", str(NRUN), "-s", str(SEED)],
                           cwd=str(WORK), capture_output=True, text=True)
        if not o.exists() or len(sd_scores(o)) != NRUN:
            return tag, "fail_dock_" + k, (r.stderr or r.stdout)[-200:]
    return tag, "ok", ""


# ==================================================== CELL 9-10 : aggregation
CHG = {"D": -1, "E": -1, "K": 1, "R": 1}      # pH 8; His treated as neutral
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
       "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
       "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
       "Y": 193.6, "V": 140.0}


def stat(v):
    if not v:
        return dict(n=0, best=np.nan, med=np.nan, iqr=np.nan)
    a = np.array(v)
    return dict(n=len(a), best=float(a.min()), med=float(np.median(a)),
                iqr=float(np.percentile(a, 75) - np.percentile(a, 25)))


def best_pose(dlg, nat=16):
    txt = open(dlg, errors="ignore").read()
    E = [float(x) for x in _ERE.findall(txt)]
    at = [l[8:] for l in txt.splitlines() if l.startswith("DOCKED: ATOM")]
    if not E or len(at) < nat:
        return None
    k = int(np.argmin(E)); blk = at[k * nat:(k + 1) * nat]
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in blk])


# C-Docker scores: discovery_mutant_training.fasta carries them in the headers
# (">G40Y_A559W 23.4125"), 287 entries, verified identical to mechanism.csv on
# the 172 overlapping variants.  Only N59N_A559W is absent.
