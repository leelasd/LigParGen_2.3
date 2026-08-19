---
title: LigParGen
emoji: 🧪
colorFrom: blue
colorTo: green
sdk: docker
app_port: 7860
pinned: false
tags:
  - mcp-server
---

# LigParGen

A Hugging Face Space mirroring the core of [zarbi.chem.yale.edu/ligpargen](https://zarbi.chem.yale.edu/ligpargen/) — an OPLS-AA/CM1A force-field parameter generator for organic ligands, built on [LigParGen](https://github.com/leelasd/LigParGen_2.3).

Submit a molecule by SMILES, by drawing a structure, or by uploading a PDB/MOL file, and get back parameter/topology files for OpenMM, GROMACS, CHARMM, LAMMPS, TINKER, CNS/X-PLOR, Q, DESMOND, and BOSS/MCPRO — plus a 3D preview of the optimized geometry.

BOSS is proprietary software supplied to this Space at runtime from a private, access-controlled store — see `docs/adr/0001` and `docs/adr/0003` in the [main repo](https://github.com/leelasd/LigParGen_2.3) for how and why.

## Using this Space as an MCP tool

Besides the web form, this Space runs as an [MCP](https://modelcontextprotocol.io) server (Gradio's built-in `mcp_server=True`), so an MCP-aware agent can call `run_ligpargen` directly instead of a human filling out the form. Same 200-atom limit, job timeout, and rate limit as the UI apply.

SSE endpoint:

```
https://lsdodda-ligpargen.hf.space/gradio_api/mcp/sse
```

For an SSE-capable client (e.g. Claude Desktop's `mcpServers` config):

```json
{
  "mcpServers": {
    "ligpargen": {
      "url": "https://lsdodda-ligpargen.hf.space/gradio_api/mcp/sse"
    }
  }
}
```

For a stdio-only client, bridge via [`mcp-remote`](https://www.npmjs.com/package/mcp-remote):

```json
{
  "mcpServers": {
    "ligpargen": {
      "command": "npx",
      "args": ["mcp-remote", "https://lsdodda-ligpargen.hf.space/gradio_api/mcp/sse"]
    }
  }
}
```

The running Space also has its own "MCP" tab (bottom of the page) with the same endpoint and a live tool schema.
