#!/usr/bin/env bash
# Builds the LigParGen Docker image against your own licensed BOSS install.
#
# BOSS is proprietary, "all rights reserved" academic software (see
# docs/adr/0001-boss-binary-supplied-locally-never-published.md) -- it can't
# live in this repo or be baked into any published image. This script stages
# a trimmed copy of your BOSS install into docker/boss/ (gitignored) so the
# Dockerfile can COPY it into the build context, then builds the image.
#
# Usage: ./build.sh /path/to/your/boss/install [image-tag]

set -euo pipefail

BOSS_HOST_DIR="${1:-}"
IMAGE_TAG="${2:-ligpargen:dev}"
STAGE_DIR="docker/boss"

if [ -z "$BOSS_HOST_DIR" ]; then
    echo "Usage: ./build.sh /path/to/your/boss/install [image-tag]" >&2
    echo "  e.g. ./build.sh ~/Codes/WLJ/boss" >&2
    exit 1
fi

if [ ! -x "$BOSS_HOST_DIR/BOSS" ]; then
    echo "ERROR: $BOSS_HOST_DIR/BOSS not found or not executable." >&2
    echo "This image requires your own licensed BOSS install -- see" >&2
    echo "docs/adr/0001-boss-binary-supplied-locally-never-published.md" >&2
    exit 1
fi

echo "Staging BOSS from $BOSS_HOST_DIR into $STAGE_DIR (excluding molecules/ and testjobs/, ~381MB of BOSS's own example/test data LigParGen never touches)..."
mkdir -p "$STAGE_DIR"
rsync -a --delete --exclude 'molecules/' --exclude 'testjobs/' "$BOSS_HOST_DIR"/ "$STAGE_DIR"/

echo "Building $IMAGE_TAG..."
docker build --platform linux/amd64 -t "$IMAGE_TAG" .

echo "Done. Example run:"
echo "  docker run --rm -v \$(pwd):/work -w /work $IMAGE_TAG -s 'c1ccc(cc1)O' -r PHN -c 0 -o 0 -l"
