#!/bin/sh
#------------------------------------------------------------------------------
# Regression tests of sharpInterFoam / sharpInterThermoFoam
#
#   tests/regression/run.sh [case ...]
#
# Requires the OpenFOAM environment and the built solvers (./Allwmake).
# Each case is copied to $REGRESSION_WORK (default: ./regression-work),
# run, evaluated (metrics.py) and compared with reference.json (check.py).
# Exit status: number of failed cases.
#------------------------------------------------------------------------------
here=$(cd "${0%/*}" && pwd)
work=${REGRESSION_WORK:-$PWD/regression-work}
python=${PYTHON:-python3}

cases="$*"
[ -n "$cases" ] || cases="zalesak staticDroplet translatingDrop twoLayerConduction thermocapillaryCavity movingMesh"

mkdir -p "$work"
nFailed=0
for c in $cases
do
    echo "=== $c"
    rm -rf "$work/$c"
    cp -r "$here/cases/$c" "$work/$c"
    start=$(date +%s)
    if ! "$work/$c/Allrun" > "$work/$c/log.Allrun" 2>&1
    then
        echo "  Allrun FAILED (see $work/$c)"
        nFailed=$((nFailed + 1))
        continue
    fi
    end=$(date +%s)
    if ! $python "$here/metrics.py" "$c" "$work/$c" > "$work/$c/metrics.json"
    then
        echo "  metrics FAILED (see $work/$c)"
        nFailed=$((nFailed + 1))
        continue
    fi
    echo "  run time $((end - start)) s"
    $python "$here/check.py" "$c" "$work/$c/metrics.json" || nFailed=$((nFailed + 1))
done

echo "=== $nFailed failed case(s)"
exit $nFailed
