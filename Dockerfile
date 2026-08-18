# Runs LigParGen against a real BOSS install. BOSS is proprietary, "all rights
# reserved" academic software -- it is never committed to this repo and this
# image is never published anywhere (see docs/adr/0001). Build with ./build.sh
# pointing at your own licensed BOSS install; see that script or the README
# for how BOSS gets staged into the build context.
#
# BOSS is a 32-bit x86 binary built for a "GNU/Linux 2.6.18"-era glibc (see
# docs/research/docker-boss-base-image.md, issue #5) -- this image is built
# for linux/amd64 with i386 multiarch libraries, and needs to run under
# emulation on non-amd64 hosts (e.g. `docker build --platform linux/amd64`).
FROM --platform=linux/amd64 python:3.12-slim-bookworm

# --- BOSS runtime dependencies: csh/tcsh for BOSS's own driver scripts, and
# i386 multiarch libs for the 32-bit BOSS binary. libgl1 covers rdkit/openbabel's
# optional GL-dependent bits defensively (see #5's findings).
RUN dpkg --add-architecture i386 \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        csh \
        tcsh \
        libc6:i386 \
        libstdc++6:i386 \
        libgcc-s1:i386 \
        libgl1 \
    && rm -rf /var/lib/apt/lists/*

# --- BOSS, supplied locally at build time only (see docs/adr/0001). Never
# baked into any image that leaves this machine. Defaults to docker/boss/,
# staged there by ./build.sh; override with --build-arg BOSS_SRC_DIR=<path>
# if you've already staged a trimmed BOSS copy somewhere else in the build
# context yourself.
ARG BOSS_SRC_DIR=docker/boss
COPY ${BOSS_SRC_DIR} /opt/boss
RUN test -x /opt/boss/BOSS || (echo "ERROR: /opt/boss/BOSS is missing from the build context." \
        "This image requires your own licensed BOSS install -- run ./build.sh /path/to/your/boss" \
        "instead of building this Dockerfile directly. See docs/adr/0001-boss-binary-supplied-locally-never-published.md." \
        >&2 && exit 1)

ENV BOSSdir=/opt/boss
# MCPRO is dropped entirely (see docs/adr/0002) -- MCPROdir is deliberately
# never set, so LigParGen's existing MCPRO-absent fallback path is what runs.

# --- LigParGen itself. Pinned dependencies (numpy/pandas/rdkit/networkx/
# openbabel) come from setup.py's install_requires (see #6's research) and
# install from prebuilt wheels -- no system Open Babel package is needed
# once the babel-CLI call sites are on the openbabel Python API (see #8).
WORKDIR /app
COPY setup.py .
COPY LigParGen/ LigParGen/
RUN pip install --no-cache-dir -e .

ENTRYPOINT ["LigParGen"]
