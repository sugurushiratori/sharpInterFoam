#!/usr/bin/env python3
"""Compare postProcessing/bubble.dat with the reference of Hysing et al.
(2009), test case 1, group TP2D at h = 1/320:
    minimum circularity 0.9013 at t = 1.9000
    maximum rise velocity 0.2417 at t = 0.9213
    centroid height at t = 3: 1.0813
Usage: evaluate.py [caseDir]"""
import os
import sys

case = sys.argv[1] if len(sys.argv) > 1 else "."
rows = []
with open(os.path.join(case, "postProcessing", "bubble.dat")) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        rows.append([float(x) for x in line.split()])

t = [r[0] for r in rows]
yc = [r[3] for r in rows]
vc = [r[5] for r in rows]
circ = [r[6] for r in rows]

i_c = min(range(len(circ)), key=lambda i: circ[i])
i_v = max(range(len(vc)), key=lambda i: vc[i])

ref = {
    "circularity_min": (circ[i_c], 0.9013),
    "t(circularity_min)": (t[i_c], 1.9000),
    "rise_velocity_max": (vc[i_v], 0.2417),
    "t(rise_velocity_max)": (t[i_v], 0.9213),
    "y_c(t=3)": (yc[-1], 1.0813),
}
print("t_end = %.4f" % t[-1])
for k, (x, r) in ref.items():
    print("%-22s %8.4f   reference %8.4f   rel. error %+6.2f %%"
          % (k, x, r, 100.0*(x - r)/r))
