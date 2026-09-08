#!/bin/bash
# =============================================================================
# run_stage0.sh — STAGE 0: chained "approach" scan. Gently draws His41 Ne2 to
# the leaving N (0.55 -> 0.29 nm) from the tetrahedral intermediate, recording
# the ENERGY COST at each step. Ends with the GATE-0 verdict.
#
# CRASH RECOVERY: Safe to restart after loadshedding. Each point uses done.tag
# checkpointing. If power cut mid-point, trjconv extracts the last saved frame
# (nstxout=50) and restarts from there with reduced nsteps. WFN is also saved.
#
# Result files: approach_profile.dat  and  ap_*/scan.gro (relaxed structures)
# The final ap_029/scan.gro is the competent TI seed for Stage 1 (if GATE 0 passes).
#
# RESUME AFTER LOADSHEDDING: just re-run this exact script.
#   bash stage0_posefix/run_stage0.sh
# It will skip completed points (done.tag present) and resume the interrupted one.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
source "$HERE/../env.sh"
cd "$HERE"
MDP="$HERE/approach.mdp"
OUT="$HERE/approach_profile.dat"
NSTEPS_FULL=800

# harden CP2K SCF a touch (more outer/inner cycles) for the reorienting steps
CP2KH="$HERE/cp2k_hard.inp"
sed -e 's/MAX_SCF 30/MAX_SCF 50/g' -e 's/MAX_SCF 20/MAX_SCF 40/g' "$CP2K_TMPL" > "$CP2KH"

# His-Ne2 .. leaving-N target schedule (nm): from ~0.58 (current) down to 0.29 (H-bond poised)
GRID="0.520 0.470 0.420 0.380 0.340 0.310 0.290"

[ ! -f "$OUT" ] && echo "# HisN_target_nm  E_Ha(FORCE_EVAL)  HisNe2-N_actual_A  HisHE2-N_A" > "$OUT"
prev="$TI_GRO"; prev_wfn=""

for tgt in $GRID; do
  tag="ap_$(python3 -c "print(f'{int(round($tgt*1000)):03d}')")"
  mkdir -p "$HERE/$tag"; cd "$HERE/$tag"
  echo "[$(date +%T)] approach His-N target=$tgt nm  ($tag)"

  # ── already done: skip ──────────────────────────────────────────────────────
  if [ -f done.tag ]; then
    echo "   [skip] done.tag present"
    prev="$HERE/$tag/scan.gro"
    prev_wfn="$HERE/$tag/GROMACS-RESTART.wfn"
    cd "$HERE"; continue
  fi

  # ── MDP: substitute target distance ────────────────────────────────────────
  sed -e "s/APPROACH_DIST/$tgt/" "$MDP" > scan_full.mdp

  # ── CRASH RECOVERY: if a previous trr exists, extract last geometry ─────────
  NSTEPS=$NSTEPS_FULL
  START_GRO="$(readlink -f "$prev")"

  if [ -f scan.trr ] && [ -f scan.tpr ]; then
    echo "   [!] Found interrupted scan.trr — attempting crash recovery..."
    # Extract last frame from trajectory
    echo "0" | $GMX trjconv -f scan.trr -s scan.tpr -o crash_recover.gro \
      -dump 9999999 -pbc mol >> run.log 2>&1 || true

    if [ -f crash_recover.gro ] && [ -s crash_recover.gro ]; then
      # Count frames to estimate last step (nstxout=50)
      NFRAMES=$(grep -c "^[ ]*[0-9]" scan.trr 2>/dev/null || echo "0")
      # More reliable: check scan.log for last step line
      LAST_STEP=$(grep -E "^[ ]+[0-9]+ +[0-9]+\." scan.log 2>/dev/null | \
                  awk '{print $1}' | tail -1 || echo "0")
      LAST_STEP=${LAST_STEP:-0}
      if [ "$LAST_STEP" -gt 0 ] 2>/dev/null; then
        NSTEPS=$(( NSTEPS_FULL - LAST_STEP ))
        if [ "$NSTEPS" -lt 10 ]; then NSTEPS=10; fi
        START_GRO="$HERE/$tag/crash_recover.gro"
        echo "   [!] Last completed step: $LAST_STEP. Resuming with $NSTEPS remaining steps."
      else
        echo "   [!] Could not read last step from log. Using crash_recover.gro as start."
        START_GRO="$HERE/$tag/crash_recover.gro"
      fi
    else
      echo "   [!] trjconv failed. Restarting from beginning of this point."
    fi
  fi

  # ── Write the final mdp with correct nsteps ─────────────────────────────────
  sed -e "s/APPROACH_DIST/$tgt/" -e "s/^nsteps.*=.*/nsteps = $NSTEPS/" "$MDP" > scan.mdp

  # ── WFN: prefer local backup, then neighbour, then none ─────────────────────
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

  # ── mdrun (try ntmpi 1 first, fallback to no ntmpi flag) ────────────────────
  $GMX mdrun -s scan.tpr -deffnm scan -ntmpi 1 -ntomp "$OMP_NUM_THREADS" \
    -pin on >> run.log 2>&1 || \
  $GMX mdrun -s scan.tpr -deffnm scan -ntomp "$OMP_NUM_THREADS" \
    -pin on >> run.log 2>&1 || true

  if ! ls scan_cp2k*.out > /dev/null 2>&1; then
    echo "  !! mdrun produced no CP2K output:"; tail -15 run.log | sed 's/^/     /'
    exit 1
  fi

  # ── Extract energy (min of last 100 FORCE_EVAL readings) ────────────────────
  E=$(grep "ENERGY| Total FORCE_EVAL" scan_cp2k*.out | awk '{print $NF}' | \
      tail -100 | sort -g | head -1)

  # ── Measure distances from the relaxed structure ─────────────────────────────
  read HN HH < <(python3 - scan.gro "$HNE2" "$HHE2" "$NL" <<'PY'
import sys, math
at=open(sys.argv[1]).read().splitlines()[2:-1]; ne2,he2,n=[int(x) for x in sys.argv[2:]]
p=lambda i:(float(at[i-1][20:28])*10,float(at[i-1][28:36])*10,float(at[i-1][36:44])*10)
d=lambda a,b:math.dist(p(a),p(b)); print(f"{d(ne2,n):.2f} {d(he2,n):.2f}")
PY
)

  # ── Record and mark done ────────────────────────────────────────────────────
  # Only append if this point is not already in the profile
  grep -q "^$tgt " "$OUT" 2>/dev/null || echo "$tgt $E $HN $HH" >> "$OUT"
  touch done.tag
  echo "     -> E=$E Ha   His-Ne2..N=$HN A   His-HE2..N=$HH A"

  prev="$HERE/$tag/scan.gro"
  prev_wfn="$HERE/$tag/GROMACS-RESTART.wfn"
  cd "$HERE"
done

echo ""; echo "=== STAGE 0 ENERGY COST of the approach (relative to first point) ==="
python3 - "$OUT" <<'PY'
import sys; H=627.5094740631
r=[l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
if not r: print("no data"); raise SystemExit
e0=float(r[0][1])
print(f"  {'HisN_tgt':>9}{'dE_kcal':>9}{'HisNe2-N_A':>12}{'HisHE2-N_A':>12}")
for a in r:
    print(f"  {float(a[0]):>9.3f}{(float(a[1])-e0)*H:>+9.1f}{float(a[2]):>12.2f}{float(a[3]):>12.2f}")
print("\n  If the last dE is small (<~10 kcal/mol) and His-Ne2..N reached <=3.5 A -> pose is ACCESSIBLE.")
print("  If dE is large / His-Ne2..N stayed far -> pose inaccessible; use a crystal-structure complex.")
PY

echo ""; echo "=== GATE 0 check on the CLOSEST-approach structure ==="
CLOSEST_TAG=$(ls "$HERE" 2>/dev/null | grep -E '^ap_[0-9]+$' | sort -t_ -k2 -n | head -1)
CLOSEST_GRO="$HERE/$CLOSEST_TAG/scan.gro"
if [ -f "$CLOSEST_GRO" ]; then
  python3 "$HERE/check_pose.py" "$CLOSEST_GRO"
else
  echo "  (Stage 0 not yet complete — no closest-approach structure yet)"
fi
