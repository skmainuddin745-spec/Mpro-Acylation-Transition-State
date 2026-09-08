#!/usr/bin/env python3
# analyze_ts.py <cleave_profile.dat> — find the acylation TS on the cleavage path.
import sys
H = 627.5094740631
rows = []
for ln in open(sys.argv[1]):
    if ln.startswith('#') or not ln.strip(): continue
    a = ln.split()
    try: rows.append((float(a[0]), float(a[1]), float(a[2]), float(a[3]), float(a[4])))
    except: pass
if not rows:
    print("no data yet"); sys.exit()
e0 = rows[0][1]                       # tetrahedral intermediate = reference
print(f"  {'C-N_nm':>6} {'dE_kcal':>8} {'C-N_A':>6} {'HE2-N_A':>8} {'NE2-N_A':>8}  note")
print("  " + "-"*58)
peak = (None, -1e9)
for cn, e, cna, hhn, nen in rows:
    dE = (e - e0) * H
    if dE > peak[1]: peak = (cn, dE, cna, hhn, nen)
    note = ""
    if hhn < 1.3: note = "proton TRANSFERRED"
    elif hhn < 1.8: note = "proton transferring"
    if cna > 2.2: note += " | C-N broken"
    print(f"  {cn:>6.3f} {dE:>+8.1f} {cna:>6.2f} {hhn:>8.2f} {nen:>8.2f}  {note}")
print()
if peak[0] is not None:
    cn, dE, cna, hhn, nen = peak
    print(f"  >>> Cleavage/PT BARRIER (PBE, potential energy) = {dE:+.1f} kcal/mol")
    print(f"      TS candidate at C-N = {cn:.3f} nm  (C-N {cna} A, HE2..N {hhn} A)")
    # is it a real TS (a maximum with product below) or still rising?
    last = (rows[-1][1]-e0)*H
    if peak[0] != rows[-1][0] and last < dE - 1.0:
        print(f"      -> profile TURNS OVER (product/acyl-enzyme at {last:+.1f}) = genuine barrier. GOOD.")
        print(f"      Next: verify this point with a frequency calc (stage2_verify).")
    else:
        print(f"      -> still rising at the end (product not reached). Extend the C-N grid to 0.30-0.34,")
        print(f"         and confirm the proton actually transferred (HE2..N should reach ~1.0 A).")
        print(f"         If HE2..N never drops, His is still not close enough -> revisit Stage 0 / pose.")
