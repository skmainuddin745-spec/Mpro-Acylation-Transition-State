#!/bin/bash
# =============================================================================
# env.sh  —  central environment + inputs for the Mpro acylation TS search.
# EVERY script sources this. Edit ONLY this file if paths/mount change.
# =============================================================================
# --- 1. The external drive mount (change here if it remounts as rgrc1, rgrc2 ...) ---
DRIVE="/media/tasnim/rgrc11"
[ -d "$DRIVE/cp2k+gromacs" ] || DRIVE="$(dirname "$(find /media/tasnim -maxdepth 2 -name 'cp2k+gromacs' -type d 2>/dev/null | head -1)")"

# --- 2. CP2K-linked GROMACS (double precision, QM/MM) ---
export PATH="$DRIVE/cp2k+gromacs/gromacs-cp2k/bin:$DRIVE/cp2k+gromacs/cp2k-2025.1/exe/local:$PATH"
export LD_LIBRARY_PATH="$DRIVE/cp2k+gromacs/gromacs-cp2k/lib:${LD_LIBRARY_PATH:-}"
export CP2K_DATA_DIR="$DRIVE/cp2k+gromacs/cp2k-2025.1/data"
XSMM=$(find "$DRIVE/cp2k+gromacs" "$DRIVE/cp2k_2025.1" -maxdepth 6 -name 'libxsmm*' -type d 2>/dev/null | head -1)
[ -n "$XSMM" ] && export LD_LIBRARY_PATH="$XSMM/lib:$LD_LIBRARY_PATH"
GMX=gmx_mpi_d
command -v "$GMX" >/dev/null 2>&1 || { echo "FATAL: gmx_mpi_d (CP2K build) not on PATH. Fix DRIVE in env.sh."; return 1 2>/dev/null; exit 1; }

# --- 3. Threads (set to your physical core count) ---
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12}
export OMP_PROC_BIND=close OMP_PLACES=cores

# --- 4. Inputs (absolute; topology includes resolve from its own directory) ---
MPRO="$DRIVE/Mpro/Scenario_B_IonPair"
export TOP="$MPRO/topol_B.top"                       # master topology (do not move; has includes)
export CP2K_TMPL="$MPRO/tslocate_B/02_cp2k_scan_B.inp"
export ORIG_NDX="$MPRO/tslocate_B/index_tsB.ndx"     # original index (we augment it -> index_TS.ndx)
export TI_GRO="$MPRO/tslocate_B/CORRECTED/a_180/scan.gro"   # tetrahedral intermediate (Phase-1 endpoint)

# --- 5. This folder + the augmented index it builds ---
export TSROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export NDX="$TSROOT/index_TS.ndx"

# --- 6. Fixed atom indices (VERIFIED in em_B.gro) ---
export SG=2243        # Cys145 SG (thiolate / thioether)
export Cc=9454        # P1-Gln carbonyl C (scissile)
export Oo=9455        # P1-Gln O (oxyanion)
export NL=9456        # P1'-Ser N (leaving nitrogen / proton acceptor)
export HNE2=617       # His41 Ne2 (proton donor heavy atom)
export HHE2=618       # His41 HE2 (the transferring proton)

echo "[env] GMX=$(command -v gmx_mpi_d)  DRIVE=$DRIVE  OMP=$OMP_NUM_THREADS"
echo "[env] TOP=$TOP"
echo "[env] TI_GRO=$TI_GRO"
