#!/usr/bin/env bash
# Builds the three derived validation images on top of the main repo's own
# ligpargen:dev (build that one first with the main repo's ./build.sh --
# see ../README.md).
set -euo pipefail
cd "$(dirname "$0")"

if ! docker image inspect ligpargen:dev >/dev/null 2>&1; then
    echo "ERROR: ligpargen:dev not found. Build it first from the repo root:" >&2
    echo "  ./build.sh /path/to/your/boss/install" >&2
    exit 1
fi

# ligpargen:dev is amd64-only (BOSS is a 32-bit x86 binary -- see the main
# Dockerfile and docs/adr/0001), and on a non-amd64 host (e.g. Apple
# Silicon) buildkit won't resolve a same-name local image's platform for a
# derived build without being told explicitly, or it'll error with
# "no match for platform in manifest" -- request it explicitly here so
# this script works the same way on any host.
echo "Building ligpargen-openmm:dev..."
docker build -f Dockerfile.openmm --platform linux/amd64 -t ligpargen-openmm:dev .

echo "Building ligpargen-gmx:dev..."
docker build -f Dockerfile.gmx --platform linux/amd64 -t ligpargen-gmx:dev .

echo "Building ligpargen-lammps:dev..."
docker build -f Dockerfile.lammps --platform linux/amd64 -t ligpargen-lammps:dev .

echo "Building ligpargen-tinker:dev (compiles TINKER from source -- slower, a few minutes)..."
docker build -f Dockerfile.tinker --platform linux/amd64 -t ligpargen-tinker:dev .

echo "Done. All four images are ready -- see ../README.md for how to run a comparison."
