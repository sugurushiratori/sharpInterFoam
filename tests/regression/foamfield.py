"""Minimal reader for ASCII OpenFOAM volume fields (internalField only).

Pure Python (no numpy) so that it runs in plain OpenFOAM containers.
"""

import os
import re


def _number_list(text):
    return [float(x) for x in text.split()]


def read_internal(path, ncells=None):
    """Return the internalField as a list of floats (scalar) or of
    tuples (vector, tensor). A uniform field needs ncells."""
    with open(path) as f:
        s = f.read()
    i = s.index("internalField")
    s = s[i + len("internalField"):]
    s = s.lstrip()
    if s.startswith("uniform"):
        s = s[len("uniform"):].lstrip()
        if s.startswith("("):
            j = s.index(")")
            val = tuple(_number_list(s[1:j]))
        else:
            val = float(re.match(r"[^;\s]+", s).group(0))
        if ncells is None:
            raise ValueError("uniform field %s needs ncells" % path)
        return [val] * ncells
    # nonuniform List<...> N ( ... )
    m = re.match(r"nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    if not m:
        raise ValueError("cannot parse internalField of " + path)
    kind, n = m.group(1), int(m.group(2))
    body = s[m.end():]
    if kind == "scalar":
        vals = body.split(None, n)[:n]
        return [float(v) for v in vals]
    out = []
    pos = 0
    for _ in range(n):
        a = body.index("(", pos)
        b = body.index(")", a)
        out.append(tuple(_number_list(body[a + 1:b])))
        pos = b + 1
    return out


def time_dirs(case):
    """Sorted numerical time directories (excluding 0)."""
    ts = []
    for d in os.listdir(case):
        try:
            t = float(d)
        except ValueError:
            continue
        if t > 0 and os.path.isdir(os.path.join(case, d)):
            ts.append((t, d))
    return [d for _, d in sorted(ts)]


def latest(case):
    return time_dirs(case)[-1]


def mag(v):
    return sum(c * c for c in v) ** 0.5
