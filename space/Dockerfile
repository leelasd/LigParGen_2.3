# LigParGen Hugging Face Space. BOSS is fetched at container startup from a
# private HF Dataset repo (see docs/adr/0003 in the main repo) -- never
# committed here, never baked into this image.
FROM --platform=linux/amd64 python:3.12-slim-bookworm

RUN dpkg --add-architecture i386 \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        git \
        csh \
        tcsh \
        libc6:i386 \
        libstdc++6:i386 \
        libgcc-s1:i386 \
        libgl1 \
    && rm -rf /var/lib/apt/lists/*

# Hugging Face Docker Spaces run as UID 1000 -- create that user up front.
RUN useradd -m -u 1000 user
USER user
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH
WORKDIR $HOME/app

RUN pip install --no-cache-dir --upgrade pip

# LigParGen itself, from the public GitHub repo -- the Space repo only holds
# this space/ subdirectory's contents, not the whole LigParGen_2.3 checkout.
RUN pip install --no-cache-dir "git+https://github.com/leelasd/LigParGen_2.3.git@dev"

# Space-specific extras on top of LigParGen's own pinned deps.
COPY --chown=user requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=user app.py .

ENV GRADIO_SERVER_NAME=0.0.0.0 \
    PYTHONUNBUFFERED=1
EXPOSE 7860

CMD ["python", "app.py"]
