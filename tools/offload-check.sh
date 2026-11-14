#!/usr/bin/env bash
# Check a GPU (offload) build of GIMIC: run the benzene 3D current-density
# test case on the CPU path and on the GPU, compare jvec.vti and report the
# wall times. Run from the top of the source tree after building with
# ./setup --offload; pass the build directory (default: build).
set -euo pipefail
build=${1:-build}
here=$(cd "$(dirname "$0")/.." && pwd)
gimic="$here/$build/bin/gimic"
[ -x "$gimic" ] || { echo "no gimic in $build/bin; build with ./setup --offload first" >&2; exit 1; }
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cp -rL "$here/test/benzene/3d" "$work/cpu"
cp -rL "$here/test/benzene/3d" "$work/gpu"

run() { # dir mode
    ( cd "$1" && s=$(date +%s.%N) && GIMIC_OFFLOAD=$2 "$gimic" > gimic.out 2>&1 && e=$(date +%s.%N) && python3 -c "print(f'{$e-$s:.3f}')" )
}
echo "threads: OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
t_cpu=$(run "$work/cpu" off)
echo "CPU path:  ${t_cpu}s"
t_gpu=$(run "$work/gpu" auto)
echo "GPU path:  ${t_gpu}s   ($(grep -m1 'Offload:' "$work/gpu/gimic.out" || echo 'no offload message: no device found, CPU path was used'))"

python3 - "$work/cpu/jvec.vti" "$work/gpu/jvec.vti" <<'PY'
import sys
def load(f):
    body=open(f).read().split('Format="ascii">')[1].split('</DataArray>')[0]
    return [float(x) for x in body.split()]
a=load(sys.argv[1]); b=load(sys.argv[2])
assert len(a)==len(b), "different sizes"
dmax=max(abs(x-y) for x,y in zip(a,b))
amax=max(abs(x) for x in a)
ndiff=sum(1 for x,y in zip(a,b) if x!=y)
print(f"jvec.vti: {len(a)} values, max |diff| = {dmax:.3e} (max |value| = {amax:.3e}), {ndiff} values differ")
print("OK" if dmax <= 1e-9*max(amax,1e-30) else "MISMATCH")
PY
