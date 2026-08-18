"""
LigParGen Hugging Face Space.

Mirrors the core of https://zarbi.chem.yale.edu/ligpargen/: SMILES input
(typed or drawn via Ketcher), PDB/MOL file upload, optimization level,
charge model, and a results zip with every output format the CLI already
produces -- plus a 3D preview of the optimized geometry the original site
doesn't have.

BOSS is fetched into this container at startup from a private HF Dataset
repo (see docs/adr/0003 in the main repo) -- never committed here, never
baked into the image.
"""
import glob
import os
import shutil
import tempfile
import uuid

import gradio as gr
from gradio_molecule3d import Molecule3D
from huggingface_hub import snapshot_download
from rdkit import Chem
from rdkit.Chem import AllChem

from LigParGen.Converter import convert

MAX_ATOMS = 200
BOSS_DIR = "/home/user/boss"
BOSS_REPO = os.environ.get("BOSS_ASSET_REPO", "lsdodda/ligpargen-boss-assets")


def fetch_boss():
    """Pull BOSS from the private asset store into BOSS_DIR and fix permissions.

    HF's storage does not preserve the executable bit (confirmed empirically
    while provisioning the asset store), so every file needs chmod +x after
    the fetch -- a plain read-only volume mount can't do this, which is why
    this is a snapshot_download to a writable directory instead.
    """
    if os.path.isfile(os.path.join(BOSS_DIR, "BOSS")) and os.access(
        os.path.join(BOSS_DIR, "BOSS"), os.X_OK
    ):
        return  # already fetched (e.g. hot-reload during dev)

    token = os.environ.get("HF_TOKEN")
    if not token:
        raise RuntimeError(
            "HF_TOKEN secret is not set on this Space -- required to fetch "
            "BOSS from the private asset store. Set it in Space Settings."
        )

    os.environ.setdefault("HF_HUB_DOWNLOAD_TIMEOUT", "30")
    try:
        snapshot_download(
            repo_id=BOSS_REPO,
            repo_type="dataset",
            token=token,
            local_dir=BOSS_DIR,
        )
    except Exception as exc:  # noqa: BLE001 -- surface any failure loudly at startup
        raise RuntimeError(
            f"Failed to fetch BOSS from {BOSS_REPO}: {exc}"
        ) from exc

    for root, _dirs, files in os.walk(BOSS_DIR):
        for name in files:
            path = os.path.join(root, name)
            os.chmod(path, os.stat(path).st_mode | 0o111)

    boss_bin = os.path.join(BOSS_DIR, "BOSS")
    if not os.access(boss_bin, os.X_OK):
        raise RuntimeError(f"{boss_bin} still isn't executable after fetch+chmod")


def count_heavy_and_h_atoms(smiles=None, file_path=None):
    """Best-effort atom count (heavy + explicit/implicit H) for the 200-atom guard."""
    mol = None
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
    elif file_path:
        if file_path.lower().endswith(".pdb"):
            mol = Chem.MolFromPDBFile(file_path, sanitize=False)
        else:
            mol = Chem.MolFromMolFile(file_path, sanitize=False)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    return mol.GetNumAtoms()


def run_ligpargen(smiles_text, upload_file, opt_iters, charge_model, charge, progress=gr.Progress()):
    smiles_text = (smiles_text or "").strip()
    upload_path = upload_file if upload_file is not None else None

    if not smiles_text and not upload_path:
        raise gr.Error("Provide a SMILES string (typed or drawn) or upload a PDB/MOL file.")
    if smiles_text and upload_path:
        raise gr.Error("Submit either SMILES or a PDB/MOL file, not both.")

    n_atoms = count_heavy_and_h_atoms(smiles=smiles_text or None, file_path=upload_path)
    if n_atoms is not None and n_atoms > MAX_ATOMS:
        raise gr.Error(f"Molecule has {n_atoms} atoms; maximum allowed is {MAX_ATOMS}.")

    lbcc = charge_model == "1.14*CM1A-LBCC (neutral molecules)"
    resolved_charge = 0 if lbcc else int(charge)

    resname = "LPG" + uuid.uuid4().hex[:5].upper()
    job_dir = tempfile.mkdtemp(prefix="ligpargen_")

    kwargs = dict(opt=int(opt_iters), charge=resolved_charge, lbcc=lbcc, resname=resname)
    if upload_path:
        staged = os.path.join(job_dir, os.path.basename(upload_path))
        shutil.copyfile(upload_path, staged)
        if upload_path.lower().endswith(".pdb"):
            kwargs["pdb"] = os.path.basename(staged)
        else:
            kwargs["mol"] = os.path.basename(staged)
    else:
        kwargs["smiles"] = smiles_text

    progress(0.1, desc="Running BOSS + LigParGen...")
    starting_dir = os.getcwd()
    try:
        os.chdir(job_dir)
        convert(**kwargs)
    finally:
        os.chdir(starting_dir)

    zip_path = os.path.join(job_dir, f"{resname}.zip")
    if not os.path.isfile(zip_path):
        raise gr.Error("LigParGen did not produce output -- check the molecule is valid.")

    preview_pdb_candidates = glob.glob(f"/tmp/{resname}.pdb")
    preview_pdb = preview_pdb_candidates[0] if preview_pdb_candidates else None

    progress(1.0, desc="Done")
    return zip_path, preview_pdb, f"Done -- {resname}"


KETCHER_HTML = """
<div id="loading" style="display:flex;justify-content:center;align-items:center">
<p style="padding:0.2rem 1rem 0 0;color:#888; font-size:1rem">loading structure editor</p>
</div>
<div id="root" style="height:420px"></div>
"""

# Ketcher 2.7.2, pinned per docs/research/ketcher-and-molecule3d-integration.md --
# newer versions have a known getSmiles() regression on some structures.
KETCHER_LOAD_JS = """
async () => {
  let url = "https://huggingface.co/datasets/simonduerr/ketcher-2.7.2/raw/main/static/css/main.6a646761.css";
  fetch(url).then(r => r.text()).then(text => {
    const style = document.createElement('style');
    style.textContent = text;
    document.head.appendChild(style);
  });
  url = "https://huggingface.co/datasets/simonduerr/ketcher-2.7.2/resolve/main/static/js/main.5445f351.js";
  fetch(url).then(r => r.text()).then(text => {
    const script = document.createElement('script');
    script.src = URL.createObjectURL(new Blob([text], { type: 'application/javascript' }));
    document.head.appendChild(script);
    document.getElementById('loading').style.display = 'none';
  });
}
"""

KETCHER_GET_SMILES_JS = """
async () => {
  return ketcher.getSmiles().then(function(smi){ return smi; });
}
"""


def build_ui():
    with gr.Blocks(title="LigParGen") as demo:
        gr.Markdown("# LigParGen\nOPLS-AA/CM1A force-field parameter generator for organic ligands.")

        with gr.Row():
            with gr.Column():
                gr.Markdown("### Step 1: Input structure")
                smiles_box = gr.Textbox(label="SMILES", placeholder="Enter SMILES, e.g. c1ccccc1")
                gr.Button("Sample: Benzene").click(
                    lambda: "c1ccccc1", inputs=None, outputs=smiles_box
                )

                gr.Markdown("**Or draw a structure:**")
                ketcher_html = gr.HTML(KETCHER_HTML)
                ketcher_hidden = gr.Textbox(visible=False)
                # Two separate event bindings, not a .click().then() chain -- matches
                # the proven-working simonduerr/gradio-2dmoleculeeditor pattern. A
                # chained .then() after a JS-only (fn=None) step does not reliably
                # forward the JS-computed value to the next step (reproduced live:
                # the hidden textbox never left its default " " placeholder value).
                gr.Button("Use drawn structure").click(
                    fn=None, inputs=[], outputs=[ketcher_hidden], js=KETCHER_GET_SMILES_JS
                )
                ketcher_hidden.change(fn=lambda s: s, inputs=[ketcher_hidden], outputs=[smiles_box])

                gr.Markdown("**Or upload a MOL/PDB file** (must include all hydrogens):")
                upload = gr.File(label="MOL or PDB file", file_types=[".mol", ".pdb"])

                gr.Markdown("### Step 2: Options")
                opt_iters = gr.Dropdown(["0", "1", "2", "3"], value="0", label="Molecule optimization iterations")
                charge_model = gr.Radio(
                    ["1.14*CM1A-LBCC (neutral molecules)", "1.14*CM1A (neutral or charged)"],
                    value="1.14*CM1A-LBCC (neutral molecules)",
                    label="Charge model",
                )
                charge = gr.Dropdown(["0", "-1", "-2", "1", "2"], value="0", label="Molecule charge")

                submit = gr.Button("Submit Molecule", variant="primary")

            with gr.Column():
                gr.Markdown("### Results")
                status = gr.Textbox(label="Status", interactive=False)
                output_zip = gr.File(label="Download all output formats (.zip)")
                preview = Molecule3D(label="3D preview (optimized geometry)", reps=[{"style": "stick"}])

        submit.click(
            run_ligpargen,
            inputs=[smiles_box, upload, opt_iters, charge_model, charge],
            outputs=[output_zip, preview, status],
        )

        demo.load(fn=None, inputs=None, outputs=None, js=KETCHER_LOAD_JS)

    demo.queue(default_concurrency_limit=1)
    return demo


if __name__ == "__main__":
    fetch_boss()
    os.environ["BOSSdir"] = BOSS_DIR
    build_ui().launch(server_name="0.0.0.0", server_port=7860)
