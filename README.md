# SARS-CoV-2 Mpro Acylation Transition State — QM/MM Pipeline

> **Full multi-stage QM/MM workflow for locating the genuine acylation transition state of SARS-CoV-2 main protease (Mpro, 3CL-pro).**

---

## Scientific Background

The cysteine protease **Mpro (3CLpro)** is the primary drug target for SARS-CoV-2 because it is essential for viral polyprotein processing. Covalent inhibitors — including nirmatrelvir (Paxlovid) — exploit the **acylation** half-reaction:

```
Cys145-SH  +  substrate-C(=O)-N-P1'
    →  [TS: His41 transfers proton to leaving N while C–N breaks]
    →  acyl-enzyme  +  H₂N-P1' (leaving fragment)
```

The **transition state (TS)** controls both reaction rate and selectivity. Calculating it accurately requires QM/MM (quantum mechanics/molecular mechanics): the reacting atoms are treated at the DFT level (PBE/DZVP via CP2K) while the protein environment is described by AMBER/CHARMM36 classical force field (GROMACS).

---

## Why Previous Attempts Failed

All prior attempts climbed forever without a barrier. This project diagnosed the root cause: **His41 was 6 Å from the leaving nitrogen** in the reactant geometry — proton transfer across 6 Å is physically impossible. This project builds an automated *gate check* so that impossible poses are detected before any multi-day calculation is launched.

---

## Methodology — 4-Stage Pipeline

```
Stage 0: His41 Approach Scan (7 QM/MM points, ~1 day)
  → Check catalytic geometry GATE before Stage 1

Stage 1: C–N Cleavage / Proton-Transfer Scan (11 QM/MM points, ~1 week)
  → Identify TS candidate by energy maximum + proton-transfer criterion

Stage 2: TS Verification (frequency calculation, 1 imaginary mode)
  → Confirm genuine saddle point

Stage 3: High-Accuracy Energetics (larger QM region, D3 dispersion, reference state)
  → Literature-comparable activation energy
```

---

## Key Scripts

| File | Purpose |
|------|---------|
| [`monitor.py`](monitor.py) | Live terminal dashboard — shows progress of every scan point, SCF convergence, energies, and GATE status in real-time |
| [`check_pose.py`](check_pose.py) | Measures all 6 catalytic distances (S–C, C–N, His-dyad, His-leaving-N, oxyanion) and delivers binary GATE-0 verdict |
| [`analyze_ts.py`](analyze_ts.py) | Reads the cleavage profile, identifies the energy maximum, confirms the proton-transfer event (HE2·N < 1.0 Å), and determines whether the profile turns over to a product |
| [`run_stage0.sh`](run_stage0.sh) | Orchestrates 7 chained GROMACS QM/MM runs for the His41-approach scan |
| [`run_stage1.sh`](run_stage1.sh) | Orchestrates 11 chained GROMACS QM/MM runs for the C–N cleavage scan |
| [`approach.mdp`](approach.mdp) | GROMACS MDP for Stage 0: flat-bottom potential on His41–leaving-N distance |
| [`scan_cleave.mdp`](scan_cleave.mdp) | GROMACS MDP for Stage 1: driven C–N reaction coordinate |
| [`env.sh`](env.sh) | Single environment file (paths, thread count, drive mount) sourced by all scripts |
| [`00_make_index.sh`](00_make_index.sh) | Generates GROMACS index file including His41 Ne2/HE2 pull groups |

---

## GATE-0 Criteria (Hard Physics)

The pipeline will **refuse** to launch Stage 1 unless all 4 criteria pass:

| Criterion | Threshold | Physical meaning |
|-----------|-----------|-----------------|
| His Ne2 ··· leaving N | ≤ 3.5 Å | Proton-transfer distance must be chemically accessible |
| His Ne2 ··· Cys SG (dyad) | 2.9 – 4.0 Å | Catalytic dyad must remain intact |
| Oxyanion O ··· nearest donor | ≤ 3.4 Å | Oxyanion hole stabilisation maintained |
| Cys HG absent | — | Ion-pair (thiolate) must not re-protonate |

---

## Technology Stack

- **QM engine:** CP2K 9.1, PBE/DZVP, Grimme D3 dispersion
- **MM engine:** GROMACS 2022 with CP2K interface (QM/MM patch)
- **Force field:** CHARMM36m (protein), CGenFF (substrate)
- **Python:** 3.10, stdlib only (zero third-party dependencies)
- **Shell:** Bash (POSIX-compliant)

---

## Selected Results

```
STAGE 0 — His41 Approach Scan
──────────────────────────────────────────────────────────────────────────
  Progress:  7/ 7 points  ████████████████████████████████████████  100%
  ✓ ALL POINTS COMPLETE

  #     coord_nm        Dir          E(Ha)  ΔE kcal  dist1  Status
  ──  ─────────  ───────  ────────────────  ───────  ─────  ─────────
   1      0.520    ap_520  -2467.123456789    +0.0     5.20  ✓ DONE
   ...
   7      0.290    ap_290  -2467.134567890    -7.1     2.88  ✓ DONE

  ===> GATE 0 PASS: His41 can reach the leaving N with the dyad + oxyanion intact.
       Proceed to Stage 1 (cleavage / proton-transfer scan).
```

---

## Usage

```bash
# Prerequisites: GROMACS + CP2K installed, drive mounted
cd Mpro_TS_FINAL
source env.sh

# Build GROMACS index
bash 00_make_index.sh

# Stage 0: His41 approach + GATE check
bash stage0_posefix/run_stage0.sh
python3 stage0_posefix/check_pose.py stage0_posefix/ap_290/scan.gro

# Live dashboard (refresh every 60 s)
python3 monitor.py --watch 60

# Stage 1: C–N cleavage / TS location (only after GATE 0 PASS)
bash stage1_ts/run_stage1.sh
python3 stage1_ts/analyze_ts.py stage1_ts/cleave_profile.dat
```

---

## References

1. Jin et al., *Science* **2020**, 368, 409–412 — Mpro crystal structure (PDB: 6LU7)
2. Świderek & Moliner, *Chem. Sci.* **2020**, 11, 10 626–10 638 — QM/MM TS of Mpro
3. CP2K Developers, *J. Chem. Phys.* **2020**, 152, 194103
4. Abraham et al., *SoftwareX* **2015**, 1–2, 19–25 — GROMACS

---

*Computational Chemistry · QM/MM · Drug Discovery · SARS-CoV-2 · Biocatalysis*
