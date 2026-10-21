#!/usr/bin/env python3
"""Compare the amplitude in postProcessing/amplitude.dat with the exact
initial-value solution of Prosperetti (1981) for equal densities and
viscosities (as used by Popinet 2009). Prints the RMS error normalised by
the initial amplitude for 0 <= omega0 t <= 25 and writes amplitude.csv.
Needs numpy and scipy.  Usage: evaluate.py [caseDir]"""
import os
import sys
import numpy as np
from scipy.special import erfc

case = sys.argv[1] if len(sys.argv) > 1 else "."
sigma, rho, lam, La = 1.0, 1.0, 1.0, 3000.0
nu = np.sqrt(sigma*rho*lam/La)/rho
k = 2*np.pi/lam
beta = 0.25                                   # rho1 rho2/(rho1 + rho2)^2
w0 = np.sqrt(sigma*k**3/(2*rho))
nk = nu*k*k
z = np.roots([1, -4*beta*np.sqrt(nk), 2*(1 - 6*beta)*nk,
              4*(1 - 3*beta)*nk**1.5, (1 - 4*beta)*nk**2 + w0**2])


def exact(t):
    t = np.asarray(t, dtype=complex)
    r = 4*(1 - 4*beta)*nk**2/(8*(1 - 4*beta)*nk**2 + w0**2)*erfc(np.sqrt(nk*t))
    for i in range(4):
        Z = np.prod([z[j] - z[i] for j in range(4) if j != i])
        r = r + z[i]/Z*w0**2/(z[i]**2 - nk)*np.exp((z[i]**2 - nk)*t)*erfc(z[i]*np.sqrt(t))
    return r.real


d = np.loadtxt(os.path.join(case, "postProcessing", "amplitude.dat"))
t, a = d[:, 0], d[:, 1]
a0 = a[0]
ae = exact(t)
mask = w0*t <= 25.0
err = np.sqrt(np.mean((a[mask]/a0 - ae[mask])**2))
np.savetxt(os.path.join(case, "amplitude.csv"),
           np.column_stack([w0*t, a/a0, ae]), delimiter=",",
           header="omega0_t,a_over_a0,exact", comments="")
print("a0 = %.6e (set 0.01), omega0 t_end = %.2f" % (a0, w0*t[-1]))
print("RMS error of a/a0 for omega0 t <= 25: %.4e" % err)
