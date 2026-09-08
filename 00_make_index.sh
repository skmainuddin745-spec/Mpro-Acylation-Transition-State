#!/bin/bash
# =============================================================================
# 00_make_index.sh — build index_TS.ndx = original index + His41 Ne2/HE2 groups.
# Run ONCE, first. Safe to re-run.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/env.sh"

cp "$ORIG_NDX" "$NDX"
# Append the two single-atom pull groups we need (His41 Ne2 = 617, HE2 = 618).
# (Cys145_SG_pull, P1Gln_C, P1prime_N already exist in the original index.)
grep -q '\[ His41_NE2 \]' "$NDX" || printf '\n[ His41_NE2 ]\n%d\n' "$HNE2" >> "$NDX"
grep -q '\[ His41_HE2 \]' "$NDX" || printf '\n[ His41_HE2 ]\n%d\n' "$HHE2" >> "$NDX"

echo "Built $NDX"
echo "Groups now present (the ones we use):"
grep -E '\[ (Cys145_SG_pull|P1Gln_C|P1prime_N|His41_NE2|His41_HE2) \]' "$NDX"
echo ""
echo "Sanity — atom identities at our indices (from the TI structure):"
python3 - "$TI_GRO" "$SG" "$Cc" "$Oo" "$NL" "$HNE2" "$HHE2" <<'PY'
import sys
gro=sys.argv[1]; idx=[int(x) for x in sys.argv[2:]]
at=open(gro).read().splitlines()[2:-1]
lab={2243:'Cys145 SG',9454:'Gln C',9455:'Gln O',9456:'Ser N(leaving)',617:'His41 Ne2',618:'His41 HE2'}
for i in idx:
    l=at[i-1]; print(f"  {i}: {l[5:10].strip():6}{l[10:15].strip():5}  <- {lab.get(i,'')}")
PY
