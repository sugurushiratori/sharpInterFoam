#!/usr/bin/env python3
"""Compare computed metrics with tests/regression/reference.json.

    check.py <caseName> <metrics.json>

Reference entries per metric:
    {"value": v, "rtol": r}   |x - v| <= r*|v|
    {"value": v, "atol": a}   |x - v| <= a
    {"max": m}                x <= m
    {"min": m}                x >= m
Exit status 1 if any check fails.
"""

import json
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
name, mfile = sys.argv[1], sys.argv[2]
with open(os.path.join(here, "reference.json")) as f:
    ref = json.load(f)[name]
with open(mfile) as f:
    got = json.load(f)

ok = True
for key, r in ref.items():
    if key.startswith("_"):
        continue
    x = got.get(key)
    if x is None:
        print("  %-22s MISSING" % key)
        ok = False
        continue
    if "value" in r:
        v = r["value"]
        tol = r["rtol"] * abs(v) if "rtol" in r else r["atol"]
        passed = abs(x - v) <= tol
        detail = "ref %.6g, tol %.2g" % (v, tol)
    elif "max" in r:
        passed = x <= r["max"]
        detail = "<= %.3g" % r["max"]
    else:
        passed = x >= r["min"]
        detail = ">= %.3g" % r["min"]
    print("  %-22s %.6g  (%s)  %s" % (key, x, detail, "ok" if passed else "FAILED"))
    ok = ok and passed

sys.exit(0 if ok else 1)
