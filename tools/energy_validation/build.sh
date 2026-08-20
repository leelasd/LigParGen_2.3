#!/usr/bin/env bash
# Builds the two derived validation images on top of the main repo's own
# ligpargen:dev (build that one first with the main repo's ./build.sh --
# see ../README.md).
set -euo pipefail
cd "$(dirname "$0")"

if ! docker image inspect ligpargen:dev >/dev/null 2>&1; then
    echo "ERROR: ligpargen:dev not found. Build it first from the repo root:" >&2
    echo "  ./build.sh /path/to/your/boss/install" >&2
    exit 1
fi

echo "Building ligpargen-openmm:dev..."
docker build -f Dockerfile.openmm -t ligpargen-openmm:dev .

echo "Building ligpargen-gmx:dev..."
docker build -f Dockerfile.gmx -t ligpargen-gmx:dev .

echo "Done. Both images are ready -- see ../README.md for how to run a comparison."
