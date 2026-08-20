# Adds Q (the Aqvist-lab MD/EVB/FEP package), built from source (freely
# available, no license needed -- see qusers/Q6 on GitHub) on top of the
# OpenMM image, for single-point energy evaluation of generated .lib/
# .Q.prm files. See ../README.md.
#
# Build (after building ligpargen-openmm:dev via Dockerfile.openmm):
#   docker build --platform linux/amd64 -f Dockerfile.q -t ligpargen-q:dev .
FROM ligpargen-openmm:dev
USER root
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
        gfortran make git ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 https://github.com/qusers/Q6.git /opt/q6-src
WORKDIR /opt/q6-src/src
RUN make all COMP=gcc
RUN find /opt/q6-src -maxdepth 3 -iname "q*6" -type f
ENV PATH="/opt/q6-src/bin:${PATH}"
ENTRYPOINT ["/bin/bash"]
