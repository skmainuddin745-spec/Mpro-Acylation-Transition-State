#!/usr/bin/env python3
# =============================================================================
# monitor.py — Live dashboard for Mpro_TS_FINAL Stage 0 + Stage 1
# Usage:  python3 monitor.py            (single shot)
#         python3 monitor.py --watch 60 (refresh every 60s)
# =============================================================================
import os, sys, glob, time, math
from datetime import datetime, timedelta

ROOT  = os.path.dirname(os.path.abspath(__file__))
S0    = os.path.join(ROOT, "stage0_posefix")
S1    = os.path.join(ROOT, "stage1_ts")
HARTREE = 627.5094740631

RST="\033[0m"; BOLD="\033[1m"; DIM="\033[2m"
G="\033[32m"; Y="\033[33m"; C="\033[36m"; R="\033[31m"; M="\033[35m"

S0_GRID = [0.520, 0.470, 0.420, 0.380, 0.340, 0.310, 0.290]
S1_GRID = [0.150, 0.160, 0.170, 0.180, 0.190, 0.200, 0.210, 0.220, 0.240, 0.260, 0.280]

def ptdir(base, prefix, val):
    return os.path.join(base, f"{prefix}{int(round(val*1000)):03d}")

def tail(path, nbytes=80000):
    try:
        with open(path,'rb') as f:
            f.seek(0,2); sz=f.tell(); f.seek(max(0,sz-nbytes))
            return f.read().decode('utf-8','ignore')
    except: return ""

def last_fe(d):
    for cp in glob.glob(os.path.join(d,"scan_cp2k*.out")):
        vals=[float(l.split()[-1]) for l in tail(cp).splitlines()
              if "ENERGY| Total FORCE_EVAL" in l]
        if vals: return vals[-1]
    return None

def last_scf(d):
    for cp in glob.glob(os.path.join(d,"scan_cp2k*.out")):
        for ln in reversed(tail(cp).splitlines()):
            if "OT DIIS" in ln:
                try: return int(ln.split()[0])
                except: pass
    return 0

def gmx_step(d):
    log=os.path.join(d,"scan.log")
    for ln in reversed(tail(log).splitlines()):
        parts=ln.strip().split()
        if len(parts)==2:
            try:
                s,t=int(parts[0]),float(parts[1])
                if abs(s-t)<1e-6*max(1,abs(t)): return s
            except: pass
    return 0

def read_profile(path):
    pts={}
    if not os.path.isfile(path): return pts
    for ln in open(path):
        if ln.startswith('#') or not ln.strip(): continue
        a=ln.split()
        try: pts[float(a[0])]=tuple(float(x) for x in a[1:])
        except: pass
    return pts

def bar(frac,n=40):
    k=int(frac*n)
    return G+"█"*k+DIM+"░"*(n-k)+RST

def section(title, grid, base, prefix, profile_path, nsteps):
    pts = read_profile(profile_path)
    done=[]; active=None
    for v in grid:
        d=ptdir(base, prefix, v)
        if os.path.isfile(os.path.join(d,"done.tag")):
            done.append(v)
        elif os.path.isdir(d) and os.path.isfile(os.path.join(d,"scan.log")) and active is None:
            active=v

    ref=pts.get(grid[0],None)
    ref_e=ref[0] if ref else None

    print(f"\n{'═'*74}")
    print(f"  {BOLD}{C}{title}{RST}")
    print(f"{'═'*74}")
    frac=len(done)/len(grid)
    print(f"  Progress: {len(done):>2}/{len(grid)} points  {bar(frac)}  {100*frac:.0f}%")

    if active is not None:
        d=ptdir(base,prefix,active)
        step=gmx_step(d); scf=last_scf(d); lE=last_fe(d)
        pct=f"{100*step/nsteps:.0f}%" if step else "0%"
        print(f"  {Y}Active:{RST} {prefix}{int(round(active*1000)):03d}  step {step}/{nsteps} ({pct})  SCF {scf}")
        if lE: print(f"  {M}Live E:{RST} {lE:.9f} Ha")
    elif len(done)==len(grid):
        print(f"  {G}✔ ALL POINTS COMPLETE{RST}")
    else:
        print(f"  {DIM}  Waiting...{RST}")

    # Table
    hdr="coord_nm" if "ap_" in prefix else "C-N_nm"
    print(f"\n  {BOLD}{'#':>3}  {hdr:>9}  {'Dir':>7}  {'E(Ha)':>16}  {'ΔE kcal':>9}  {'dist1':>7}  Status{RST}")
    print(f"  {DIM}{'─'*3}  {'─'*9}  {'─'*7}  {'─'*16}  {'─'*9}  {'─'*7}  {'─'*9}{RST}")
    for i,v in enumerate(grid,1):
        tag=f"{prefix}{int(round(v*1000)):03d}"
        d=ptdir(base,prefix,v)
        done_f=os.path.isfile(os.path.join(d,"done.tag"))
        running=os.path.isdir(d) and os.path.isfile(os.path.join(d,"scan.log")) and not done_f

        row=pts.get(v)
        e_str=f"{row[0]:.9f}" if row else "─"*14
        dE_str=f"{(row[0]-ref_e)*HARTREE:+8.1f}" if (row and ref_e) else "    N/A"
        d1_str=f"{row[1]:.2f}" if (row and len(row)>1) else "  ─"

        status=f"{G}✔ DONE{RST}" if done_f else (f"{Y}⟳ RUNNING{RST}" if running else f"{DIM}waiting{RST}")
        print(f"  {i:>3}  {v:>9.3f}  {tag:>7}  {e_str:>16}  {dE_str:>9}  {d1_str:>7}  {status}")

def run():
    os.system('clear')
    now=datetime.now()
    print(f"{'═'*74}")
    print(f"  {BOLD}Mpro Acylation TS Monitor — {now.strftime('%Y-%m-%d %H:%M:%S')}{RST}")
    print(f"  Root: {ROOT}")
    print(f"{'═'*74}")

    section("STAGE 0 — His41 Approach Scan (ap_NNN, 7 points)",
            S0_GRID, S0, "ap_", os.path.join(S0,"approach_profile.dat"), 800)

    # Check GATE 0
    closest_tag = sorted([x for x in os.listdir(S0) if x.startswith("ap_") and
                          os.path.isfile(os.path.join(S0,x,"done.tag"))],
                         key=lambda x: int(x.split("_")[1]))
    if closest_tag:
        closest_gro = os.path.join(S0, closest_tag[0], "scan.gro")
        if os.path.isfile(closest_gro):
            import subprocess
            r=subprocess.run(["python3", os.path.join(S0,"check_pose.py"), closest_gro],
                             capture_output=True, text=True)
            gate_line=[l for l in r.stdout.splitlines() if "GATE" in l]
            if gate_line:
                col=G if "PASS" in gate_line[0] else R
                print(f"\n  {BOLD}GATE 0:{RST} {col}{gate_line[0].strip()}{RST}")

    s1_profile=os.path.join(S1,"cleave_profile.dat")
    if os.path.isfile(s1_profile) or any(os.path.isdir(os.path.join(S1,f"cn_{int(round(v*1000)):03d}")) for v in S1_GRID):
        section("STAGE 1 — C-N Cleavage / Proton Transfer Scan (cn_NNN, 11 points)",
                S1_GRID, S1, "cn_", s1_profile, 800)

    print(f"\n{'═'*74}")
    print(f"  {DIM}RESUME COMMAND (after loadshedding):")
    print(f"  Stage 0:  cd {ROOT} && bash stage0_posefix/run_stage0.sh")
    print(f"  Stage 1:  cd {ROOT} && bash stage1_ts/run_stage1.sh{RST}")
    print(f"{'═'*74}\n")

if __name__=="__main__":
    delay=60
    if "--watch" in sys.argv:
        try: delay=int(sys.argv[sys.argv.index("--watch")+1])
        except: pass
        try:
            while True:
                run()
                print(f"  {DIM}Refreshing every {delay}s — Ctrl-C to stop{RST}\n")
                time.sleep(delay)
        except KeyboardInterrupt:
            print("\nMonitor stopped.")
    else:
        run()
