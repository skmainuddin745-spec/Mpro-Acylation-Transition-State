#!/bin/bash
# =============================================================================
# run_stage1.sh — STAGE 1: cleavage / proton-transfer scan -> the acylation TS.
# START ONLY after Stage 0 GATE 0 passed. Seed = the competent TI from Stage 0.
# Drives C-N 0.15 -> 0.28 nm (chained). Records FORCE_EVAL energy + geometry.
# Output: cleave_profile.dat ; run analyze_ts.py at the end.
#
# CRASH RECOVERY: Safe to restart after loadshedding. Each point uses done.tag
# checkpointing. If power cut mid-point, trjconv extracts the last saved frame
# (nstxout=50) and restarts with reduced nsteps. WFN backup is also restored.
#
# RESUME AFTER LOADSHEDDING: just re-run this exact script.
#   bash stage1_ts/run_stage1.sh
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/../env.sh"
cd "$HERE"
MDP="$HERE/scan_cleave.mdp"
OUT="$HERE/cleave_profile.dat"
NSTEPS_FULL=800

# Seed = CLOSEST-approach structure from Stage 0 (smallest His-N target = competent TI).
S0="$HERE/../stage0_posefix"
CLOSEST=$(ls "$S0" 2>/dev/null | grep -E '^ap_[0-9]+$' | sort -t_ -k2 -n | head -1)
SEED="$S0/$CLOSEST/scan.gro"
if [ -z "$CLOSEST" ] || [ ! -f "$SEED" ]; then
  echo "FATAL: Stage-0 seed not found ($SEED). Run Stage 0 first and pass GATE 0."
  exit 1
fi
echo "seed (competent TI) = $SEED"

CP2KH="$HERE/cp2k_hard.inp"
sed -e 's/MAX_SCF 30/MAX_SCF 50/g' -e 's/MAX_SCF 20/MAX_SCF 40/g' "$CP2K_TMPL" > "$CP2KH"

# C-N breaking-bond schedule (nm): TI 0.15 -> broken 0.28, fine through the expected TS (~0.19-0.23)
GRID="0.150 0.160 0.170 0.180 0.190 0.200 0.210 0.220 0.240 0.260 0.280"

[ ! -f "$OUT" ] && echo "# CN_nm  E_Ha(FORCE_EVAL,min-last-100)  CN_A  HE2-N_A  NE2-N_A" > "$OUT"
prev="$SEED"; prev_wfn=""

for cn in $GRID; do
  tag="cn_$(python3 -c "print(f'{int(round($cn*1000)):03d}')")"
  mkdir -p "$HERE/$tag"; cd "$HERE/$tag"
  echo "[$(date +%T)] cleavage C-N=$cn nm  ($tag)"

  # ── already done: skip ──────────────────────────────────────────────────────
  if [ -f done.tag ]; then
    echo "   [skip] done.tag present"
    prev="$HERE/$tag/scan.gro"
    prev_wfn="$HERE/$tag/GROMACS-RESTART.wfn"
    cd "$HERE"; continue
  fi

  # ── CRASH RECOVERY: if a previous trr exists, extract last geometry ─────────
  NSTEPS=$NSTEPS_FULL
  START_GRO="$(readlink -f "$prev")"

  if [ -f scan.trr ] && [ -f scan.tpr ]; then
    echo "   [!] Found interrupted scan.trr — attempting crash recovery..."
    echo "0" | $GMX trjconv -f scan.trr -s scan.tpr -o crash_recover.gro \
      -dump 9999999 -pbc mol >> run.log 2>&1 || true

    if [ -f crash_recover.gro ] && [ -s crash_recover.gro ]; then
      LAST_STEP=$(grep -E "^[ ]+[0-9]+ +[0-9]+\." scan.log 2>/dev/null | \
                  awk '{print $1}' | tail -1 || echo "0")
      LAST_STEP=${LAST_STEP:-0}
      if [ "$LAST_STEP" -gt 0 ] 2>/dev/null; then
        NSTEPS=$(( NSTEPS_FULL - LAST_STEP ))
        if [ "$NSTEPS" -lt 10 ]; then NSTEPS=10; fi
        START_GRO="$HERE/$tag/crash_recover.gro"
        echo "   [!] Last completed step: $LAST_STEP. Resuming with $NSTEPS remaining steps."
      else
        START_GRO="$HERE/$tag/crash_recover.gro"
        echo "   [!] Using crash_recover.gro as start (step count unreadable)."
      fi
    else
      echo "   [!] trjconv failed. Restarting point from beginning."
    fi
  fi

  # ── Write MDP with correct CN target and nsteps ─────────────────────────────
  sed -e "s/CN_DIST/$cn/" -e "s/^nsteps.*=.*/nsteps = $NSTEPS/" "$MDP" > scan.mdp

  # ── WFN: prefer local backup, then neighbour ────────────────────────────────
  if [ -f GROMACS-RESTART.wfn.bak-1 ] && [ -s GROMACS-RESTART.wfn.bak-1 ]; then
    echo "   [!] Restoring WFN from local backup"
    cp GROMACS-RESTART.wfn.bak-1 GROMACS-RESTART.wfn
  elif [ ! -f GROMACS-RESTART.wfn ] || [ ! -s GROMACS-RESTART.wfn ]; then
    if [ -n "$prev_wfn" ] && [ -f "$prev_wfn" ]; then
      echo "   [!] Using WFN from previous point: $prev_wfn"
      cp "$prev_wfn" GROMACS-RESTART.wfn
    fi
  fi

  # ── grompp ──────────────────────────────────────────────────────────────────
  $GMX grompp -f scan.mdp -c "$START_GRO" -p "$TOP" -n "$NDX" -o scan.tpr \
    -maxwarn 2 > grompp.log 2>&1
  if [ ! -f scan.tpr ]; then
    echo "  !! grompp FAILED:"; tail -15 grompp.log | sed 's/^/     /'
    exit 1
  fi

  # ── mdrun ───────────────────────────────────────────────────────────────────
  $GMX mdrun -s scan.tpr -deffnm scan -ntmpi 1 -ntomp "$OMP_NUM_THREADS" \
    -pin on >> run.log 2>&1 || \
  $GMX mdrun -s scan.tpr -deffnm scan -ntomp "$OMP_NUM_THREADS" \
    -pin on >> run.log 2>&1 || true

  if ! ls scan_cp2k*.out > /dev/null 2>&1; then
    echo "  !! mdrun no CP2K output:"; tail -15 run.log | sed 's/^/     /'
    exit 1
  fi

  # ── Extract energy ───────────────────────────────────────────────────────────
  E=$(grep "ENERGY| Total FORCE_EVAL" scan_cp2k*.out | awk '{print $NF}' | \
      tail -100 | sort -g | head -1)

  # ── Measure distances ────────────────────────────────────────────────────────
  read CNA HHN NEN < <(python3 - scan.gro "$Cc" "$NL" "$HHE2" "$HNE2" <<'PY'
import sys, math
at=open(sys.argv[1]).read().splitlines()[2:-1]; c,n,he2,ne2=[int(x) for x in sys.argv[2:]]
p=lambda i:(float(at[i-1][20:28])*10,float(at[i-1][28:36])*10,float(at[i-1][36:44])*10)
d=lambda a,b:math.dist(p(a),p(b)); print(f"{d(c,n):.2f} {d(he2,n):.2f} {d(ne2,n):.2f}")
PY
)

  # ── Record (no duplicates) and mark done ────────────────────────────────────
  grep -q "^$cn " "$OUT" 2>/dev/null || echo "$cn $E $CNA $HHN $NEN" >> "$OUT"
  touch done.tag
  echo "     -> E=$E Ha   C-N=$CNA A   HE2..N=$HHN A (proton transfer if ->1.0)   NE2..N=$NEN A"

  prev="$HERE/$tag/scan.gro"
  prev_wfn="$HERE/$tag/GROMACS-RESTART.wfn"
  cd "$HERE"
done

echo ""; python3 "$HERE/analyze_ts.py" "$OUT"
