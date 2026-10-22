#!/usr/bin/env python3
"""Compute the regression metrics of one finished case.

    metrics.py <caseName> <caseDir>

Prints a JSON dictionary {metric: value}.
"""

import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from foamfield import read_internal, latest, mag  # noqa: E402


def centres(case):
    return read_internal(os.path.join(case, "0", "C"))


def zalesak(case):
    """Interface advection by rigid rotation (one revolution).
    E_alpha = sum|alpha - alpha0| / sum alpha0 (uniform mesh)."""
    a0 = read_internal(os.path.join(case, "0", "alpha.water"))
    a = read_internal(os.path.join(case, latest(case), "alpha.water"))
    s0 = sum(a0)
    return {
        "E_alpha": sum(abs(x - y) for x, y in zip(a, a0)) / s0,
        "volume_change": abs(sum(a) - s0) / s0,
        "alpha_min": min(a),
        "alpha_max": max(a),
    }


def static_droplet(case):
    """Static droplet (radius 1, centre at the origin, x-z plane, sigma = 1):
    spurious velocity and pressure error (Laplace pressure 1)."""
    n = len(centres(case))
    C = centres(case)
    t = latest(case)
    U = read_internal(os.path.join(case, t, "U"), n)
    p = read_internal(os.path.join(case, t, "p"), n)
    pe = [1.0 if math.hypot(c[0], c[2]) < 1.0 else 0.0 for c in C]
    um = [mag(u) for u in U]
    return {
        "U_max": max(um),
        "E_p": sum(abs(x - y) for x, y in zip(p, pe)) / n,
    }


def translating_drop(case):
    """Drop translating with the uniform velocity (1 0 0), density ratio
    1000, no surface tension: the velocity must stay uniform."""
    n = len(centres(case))
    U = read_internal(os.path.join(case, latest(case), "U"), n)
    return {"U_error_max": max(mag((u[0] - 1.0, u[1], u[2])) for u in U)}


def two_layer_conduction(case):
    """Steady conduction through two layers (k = 1 below z = 0.5,
    k = 0.1 above), T = 1 at z = 0 and T = 0 at z = 1.
    Exact heat flux q = 1/(0.5/1 + 0.5/0.1)."""
    C = centres(case)
    n = len(C)
    T = read_internal(os.path.join(case, latest(case), "T"), n)
    order = sorted(range(n), key=lambda i: C[i][2])
    dz = 1.0 / n
    q_exact = 1.0 / (0.5 / 1.0 + 0.5 / 0.1)
    q_bottom = 1.0 * (1.0 - T[order[0]]) / (0.5 * dz)
    q_top = 0.1 * (T[order[-1]] - 0.0) / (0.5 * dz)
    return {
        "q_bottom_error": abs(q_bottom - q_exact) / q_exact,
        "q_top_error": abs(q_top - q_exact) / q_exact,
    }


def thermocapillary_cavity(case):
    """Thermocapillary cavity (x-y plane), short run: maximum velocity and
    the difference between the serial and the 2-process runs."""
    n = len(centres(case))
    t = latest(case)
    U = read_internal(os.path.join(case, t, "U"), n)
    a = read_internal(os.path.join(case, t, "alpha.water"), n)
    par = os.path.join(case, "parallel")
    Up = read_internal(os.path.join(par, latest(par), "U"), n)
    ap = read_internal(os.path.join(par, latest(par), "alpha.water"), n)
    umax = max(mag(u) for u in U)
    return {
        "U_max": umax,
        "parallel_U_diff": max(
            mag((u[0] - v[0], u[1] - v[1], u[2] - v[2])) for u, v in zip(U, Up)
        ) / umax,
        "parallel_alpha_diff": max(abs(x - y) for x, y in zip(a, ap)),
    }


def moving_mesh(case):
    """Moving wall (mesh motion): liquid volume sum(alpha V) is conserved."""
    t = latest(case)
    V0 = read_internal(os.path.join(case, "0", "V"))
    n = len(V0)
    a0 = read_internal(os.path.join(case, "0", "alpha.water"), n)
    V = read_internal(os.path.join(case, t, "V"))
    a = read_internal(os.path.join(case, t, "alpha.water"), n)
    vol0 = sum(x * y for x, y in zip(a0, V0))
    vol = sum(x * y for x, y in zip(a, V))
    return {
        "volume_change": abs(vol - vol0) / vol0,
        "alpha_max": max(a),
        "alpha_min": min(a),
    }


CASES = {
    "zalesak": zalesak,
    "staticDroplet": static_droplet,
    "translatingDrop": translating_drop,
    "twoLayerConduction": two_layer_conduction,
    "thermocapillaryCavity": thermocapillary_cavity,
    "movingMesh": moving_mesh,
}


if __name__ == "__main__":
    name, case = sys.argv[1], sys.argv[2]
    print(json.dumps(CASES[name](case), indent=2))
