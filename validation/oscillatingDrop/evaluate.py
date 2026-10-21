#!/usr/bin/env python3
"""Oscillation frequency and damping of the l = 2 mode from
postProcessing/shape.dat (D = a_z - a_r).

Inviscid (Lamb): omega0^2 = l(l-1)(l+2) sigma/(rho R^3) = 8 (sigma = rho = R = 1,
the gas density 1e-3 is neglected). The viscosity is recovered from the
measured damping rate Gamma and frequency omega with the finite-viscosity
correction (zeta = Gamma/omega):
    nu = R^2 Gamma/5 [1 + 3 sqrt(10)/25 zeta^(1/2) + 27/125 zeta
                      + 2817 sqrt(10)/12500 zeta^(3/2) + 2049/3125 zeta^2]
and compared with the input nu (Lamb's relation nu = R^2 Gamma/5 is shown
for comparison).
Usage: evaluate.py [caseDir] [nu]"""
import math
import os
import sys

case = sys.argv[1] if len(sys.argv) > 1 else "."
nu_in = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
omega0 = math.sqrt(8.0)

t, d = [], []
with open(os.path.join(case, "postProcessing", "shape.dat")) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        v = [float(x) for x in line.split()]
        t.append(v[0])
        d.append(v[1] - v[2])

# extrema of |D| by parabolic interpolation
ext = []
for i in range(1, len(d) - 1):
    if (d[i] - d[i - 1])*(d[i + 1] - d[i]) < 0:
        t0, t1, t2 = t[i - 1], t[i], t[i + 1]
        y0, y1, y2 = d[i - 1], d[i], d[i + 1]
        den = (t0 - t1)*(t0 - t2)*(t1 - t2)
        A = (t2*(y1 - y0) + t1*(y0 - y2) + t0*(y2 - y1))/den
        B = (t2*t2*(y0 - y1) + t1*t1*(y2 - y0) + t0*t0*(y1 - y2))/den
        C = (t1*t2*(t1 - t2)*y0 + t2*t0*(t2 - t0)*y1 + t0*t1*(t0 - t1)*y2)/den
        tm = -B/(2*A)
        ext.append((tm, A*tm*tm + B*tm + C))
ext = [e for e in ext if e[0] > 0.2]          # skip the start-up
if len(ext) < 3:
    sys.exit("not enough extrema")

# omega from the spacing of the extrema (half periods), Gamma from a
# least-squares fit of ln|D| at the extrema
n = len(ext)
half = (ext[-1][0] - ext[0][0])/(n - 1)
omega_d = math.pi/half
ts = [e[0] for e in ext]
ls = [math.log(abs(e[1])) for e in ext]
tm, lm = sum(ts)/n, sum(ls)/n
Gamma = -sum((a - tm)*(b - lm) for a, b in zip(ts, ls))/sum((a - tm)**2 for a in ts)
# undamped frequency of the fitted damped oscillation
omega = math.sqrt(omega_d**2 + Gamma**2)
zeta = Gamma/omega
s = math.sqrt(zeta)
nu_lamb = Gamma/5.0
nu_corr = nu_lamb*(1 + 3*math.sqrt(10)/25*s + 27/125*zeta
                   + 2817*math.sqrt(10)/12500*zeta*s + 2049/3125*zeta**2)
print("extrema: %d, t = %.3f ... %.3f" % (n, ext[0][0], ext[-1][0]))
print("omega_d = %.4f, Gamma = %.5f, zeta = %.4f" % (omega_d, Gamma, zeta))
print("omega_d/omega0(Lamb) = %.4f" % (omega_d/omega0))
print("nu (Lamb, R^2 Gamma/5)        = %.5f  (%+.2f %%)" % (nu_lamb, 100*(nu_lamb/nu_in - 1)))
print("nu (finite-viscosity corr.)   = %.5f  (%+.2f %%)" % (nu_corr, 100*(nu_corr/nu_in - 1)))
