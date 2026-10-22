#!/usr/bin/env python3
"""Compare the drop velocity in postProcessing/drop.dat with the
Young-Goldstein-Block terminal velocity
    V = 2/((2 + k')(2 + 3 mu')) * (-sigmaTemp) |grad T| R / mu
(mu' = k' = 1, sigmaTemp = -1, grad T = 1, R = 0.5, mu = 1).
The drop velocity is averaged over the last 20 % of the run.
Usage: evaluate.py [caseDir]"""
import os
import sys

case = sys.argv[1] if len(sys.argv) > 1 else "."
V_ygb = 2.0/((2.0 + 1.0)*(2.0 + 3.0))*1.0*1.0*0.5/1.0

rows = []
with open(os.path.join(case, "postProcessing", "drop.dat")) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])

t_end = rows[-1][0]
tail = [r for r in rows if r[0] >= 0.8*t_end]
w = sum(r[3] for r in tail)/len(tail)
print("t_end = %.3f, volume change = %.2e"
      % (t_end, (rows[-1][1] - rows[0][1])/rows[0][1]))
for frac in (0.2, 0.4, 0.6, 0.8, 1.0):
    r = min(rows, key=lambda r: abs(r[0] - frac*t_end))
    print("  t = %6.3f  w_c/V_YGB = %.4f" % (r[0], r[3]/V_ygb))
print("drop velocity %.5e, V_YGB %.5e, ratio %.4f (error %+.2f %%)"
      % (w, V_ygb, w/V_ygb, 100.0*(w/V_ygb - 1.0)))
