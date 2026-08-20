#!/bin/bash
# Run BOSS vs. OpenMM vs. GROMACS single-point energy comparison for one
# molecule, given a real Zmatrix. See ../README.md for the full
# methodology, setup, and known gotchas -- read that first.
#
# Usage:
#   ./compare.sh <resname> <path-to-zmatrix> [charge]
#
# <resname> MUST be exactly 3 characters -- LigParGen's PDB writer uses a
# fixed 3-column residue-name field; anything else silently overflows into
# the coordinate columns and produces garbage (or a hard parse error from
# a strict reader like OpenMM's).
#
# Example:
#   ./compare.sh PHN ~/Codes/WLJ/boss/molecules/small/phenol.z 0
#
# To validate a molecule from the live production Space instead of a BOSS
# reference file: submit it via the Space's API, download the output zip,
# extract the <resid>.z file from it, and pass that path here instead --
# see ../README.md's "Testing the current pipeline, not just legacy
# reference files" section for the exact curl commands.
set -euo pipefail
cd "$(dirname "$0")"

resname="$1"
zmat="$2"
charge="${3:-0}"

if [ ${#resname} -ne 3 ]; then
    echo "ERROR: resname must be exactly 3 characters (got '$resname', ${#resname} chars)." >&2
    echo "LigParGen's PDB writer uses a fixed 3-column residue-name field." >&2
    exit 1
fi

for img in ligpargen-openmm:dev ligpargen-gmx:dev; do
    if ! docker image inspect "$img" >/dev/null 2>&1; then
        echo "ERROR: $img not found. Run ./build.sh first." >&2
        exit 1
    fi
done

WORKDIR="$(mktemp -d)"
trap 'rm -rf "$WORKDIR"' EXIT
cp "$zmat" "$WORKDIR/"
zmat_basename="$(basename "$zmat")"
LP="$(cd .. && cd .. && pwd)/LigParGen"

echo "=== $resname ($zmat_basename, charge=$charge) ==="

gen_log="$WORKDIR/gen.log"
docker run --rm -v "$LP":/app/LigParGen -v "$(pwd)":/tools -v "$WORKDIR":/tmp \
    --entrypoint python3 ligpargen-openmm:dev /tools/gen_and_boss_energy.py "$resname" "$zmat_basename" "$charge" \
    > "$gen_log" 2>&1
if ! grep -q WROTE_FILES "$gen_log"; then
    echo "Generation/BOSS energy step failed:"
    cat "$gen_log"
    exit 1
fi
grep -E 'BOSS_ENERGY|BOSS_TERMS' "$gen_log"

omm_log="$WORKDIR/omm.log"
docker run --rm -v "$(pwd)":/tools -v "$WORKDIR":/tmp \
    --entrypoint python3 ligpargen-openmm:dev /tools/eval_openmm_energy.py "$resname" \
    > "$omm_log" 2>&1 || true
if grep -q OPENMM_ENERGY "$omm_log"; then
    grep -E 'OPENMM_ENERGY|OPENMM_TERMS' "$omm_log"
else
    echo "OpenMM evaluation failed:"
    cat "$omm_log"
fi

gmx_log="$WORKDIR/gmx.log"
docker run --rm -v "$(pwd)":/tools -v "$WORKDIR":/tmp \
    --entrypoint bash ligpargen-gmx:dev -c "cd /tmp && /tools/eval_gromacs_energy.sh $resname" \
    > "$gmx_log" 2>&1 || true
if grep -q GROMACS_ENERGY "$gmx_log"; then
    grep -E 'GROMACS_ENERGY|GROMACS_TERMS' "$gmx_log"
else
    echo "GROMACS evaluation failed:"
    cat "$gmx_log"
fi
