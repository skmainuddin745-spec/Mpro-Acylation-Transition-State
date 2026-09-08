#!/usr/bin/env python3
# check_pose.py <structure.gro> — measure all catalytic distances + GATE 0 verdict.
import sys, math
IDX = dict(SG=2243, C=9454, O=9455, N=9456, HNE2=617, HE2=618,
           G143N=2218, S144N=2225, C145N=2236)
g = sys.argv[1]
at = open(g).read().splitlines()[2:-1]
def p(i):
    l = at[i-1]; return (float(l[20:28])*10, float(l[28:36])*10, float(l[36:44])*10)
def d(a, b): return round(math.dist(p(a), p(b)), 2)
def nm(i):
    l = at[i-1]; return l[10:15].strip()

sc   = d(IDX['SG'], IDX['C'])
cn   = d(IDX['C'],  IDX['N'])
dyad = d(IDX['HNE2'], IDX['SG'])
hisN = d(IDX['HNE2'], IDX['N'])
hN   = d(IDX['HE2'],  IDX['N'])
oxy  = min(d(IDX['O'], IDX['G143N']), d(IDX['O'], IDX['C145N']))
hg   = 'HG' in [nm(i) for i in range(2236, 2246)]

print(f"structure: {g}")
print(f"  S-C (thioether/attack)      = {sc} A")
print(f"  C-N (scissile)              = {cn} A")
print(f"  His Ne2 - Cys SG (dyad)     = {dyad} A")
print(f"  His Ne2 - leaving N          = {hisN} A   <-- proton-donor / acceptor gap")
print(f"  His HE2 - leaving N          = {hN} A   <-- proton to acceptor")
print(f"  oxyanion O - nearest donor  = {oxy} A")
print(f"  ion pair intact (no Cys HG) = {not hg}")
print()
ok = (hisN <= 3.5) and (2.9 <= dyad <= 4.0) and (oxy <= 3.4) and (not hg)
if ok:
    print("  ===> GATE 0 PASS: His41 can reach the leaving N with the dyad + oxyanion intact.")
    print("       Proceed to Stage 1 (cleavage / proton-transfer scan).")
else:
    print("  ===> GATE 0 not yet met:")
    if hisN > 3.5:  print(f"       - His Ne2..leaving N = {hisN} A (need <=3.5). Run more approach steps, or the pose is inaccessible.")
    if not (2.9 <= dyad <= 4.0): print(f"       - dyad His Ne2..SG = {dyad} A (His pulled out of the dyad -> lower CV3 k or use a crystal structure).")
    if oxy > 3.4:   print(f"       - oxyanion lost ({oxy} A).")
    if hg:          print(f"       - Cys re-protonated (ion pair broken).")
