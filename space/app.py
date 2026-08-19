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
import random
import shutil
import string
import tempfile

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
    if smiles_text and upload_path and not upload_path.lower().endswith(".pdb"):
        raise gr.Error(
            "SMILES + upload is only supported for PDB uploads (used to fix bond "
            "orders/missing Hs) -- submit either SMILES or a MOL file, not both."
        )

    # Prefer the SMILES for the atom count when both are given: PDB uploads are
    # often missing hydrogens, which would undercount against MAX_ATOMS.
    n_atoms = count_heavy_and_h_atoms(smiles=smiles_text or None, file_path=None if smiles_text else upload_path)
    if n_atoms is not None and n_atoms > MAX_ATOMS:
        raise gr.Error(f"Molecule has {n_atoms} atoms; maximum allowed is {MAX_ATOMS}.")

    lbcc = charge_model == "1.14*CM1A-LBCC (neutral molecules)"
    resolved_charge = 0 if lbcc else int(charge)

    # PDB's resName field is a strict 3-character fixed-width column (spec cols
    # 18-20) -- LigParGen's own PDB/GRO/etc. writers use fixed-width formatting
    # too, so anything longer silently overflows into the chainID/resSeq/coordinate
    # columns downstream. Confirmed live: an 8-char resname ("LPG"+5 hex chars)
    # shifted every coordinate column and produced a non-numeric resSeq, which
    # RDKit's own PDB parser (MolFromPDBFile) outright rejected and PyMOL misread.
    # 3 random letters -- matches the CLI's own convention ("-r, should be a 3
    # LETTER WORD") and gives 26**3 = 17,576 combinations, plenty for this app's
    # single-worker, serialized-queue usage.
    resname = "".join(random.choices(string.ascii_uppercase, k=3))
    job_dir = tempfile.mkdtemp(prefix="ligpargen_")

    kwargs = dict(opt=int(opt_iters), charge=resolved_charge, lbcc=lbcc, resname=resname)
    if upload_path:
        staged = os.path.join(job_dir, os.path.basename(upload_path))
        shutil.copyfile(upload_path, staged)
        if upload_path.lower().endswith(".pdb"):
            kwargs["pdb"] = os.path.basename(staged)
            # A SMILES supplied alongside the PDB is used as a trusted template
            # to fix connectivity/bond orders/missing Hs -- see
            # LigParGen.mol_boss.convert_pdb2mol_with_smiles.
            if smiles_text:
                kwargs["smiles"] = smiles_text
        else:
            kwargs["mol"] = os.path.basename(staged)
    else:
        kwargs["smiles"] = smiles_text

    progress(0.1, desc="Running BOSS + LigParGen...")
    starting_dir = os.getcwd()
    try:
        os.chdir(job_dir)
        try:
            convert(**kwargs)
        except ValueError as exc:
            raise gr.Error(str(exc))
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
<button id="ketcher-use-btn" type="button"
        style="width:100%;margin-top:8px;padding:8px;border-radius:6px;border:none;
               background:#4b5563;color:white;cursor:pointer;font-size:1rem">
  Use drawn structure
</button>
"""

# Ketcher 2.7.2, pinned per docs/research/ketcher-and-molecule3d-integration.md --
# newer versions have a known getSmiles() regression on some structures.
#
# The "Use drawn structure" button lives inside this HTML block as a plain
# <button>, wired here via direct DOM manipulation rather than Gradio's
# js=/fn=None event-return mechanism -- that mechanism did not reliably fire
# on real user clicks in this Gradio version (reproduced live: getSmiles()
# worked when called from the console, but neither a real click nor a
# programmatic .click() on a Gradio-wired button ever invoked it, per
# instrumentation that counted calls to a wrapped getSmiles()). Writing
# straight to the target textbox's <textarea> and dispatching a native
# "input" event is the standard, robust way to feed a Gradio-bound
# component from custom JS.
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

  document.getElementById('ketcher-use-btn').addEventListener('click', async () => {
    const smi = await ketcher.getSmiles();
    const textarea = document.querySelector('#smiles_box textarea');
    const setter = Object.getOwnPropertyDescriptor(window.HTMLTextAreaElement.prototype, 'value').set;
    setter.call(textarea, smi);
    textarea.dispatchEvent(new Event('input', { bubbles: true }));
  });
}
"""

def build_ui():
    with gr.Blocks(title="LigParGen") as demo:
        gr.Markdown("# LigParGen\nOPLS-AA/CM1A force-field parameter generator for organic ligands.")

        with gr.Row():
            with gr.Column():
                gr.Markdown("### Step 1: Input structure")
                smiles_box = gr.Textbox(
                    label="SMILES", placeholder="Enter SMILES, e.g. c1ccccc1", elem_id="smiles_box"
                )
                gr.Button("Sample: Benzene").click(
                    lambda: "c1ccccc1", inputs=None, outputs=smiles_box
                )

                gr.Markdown("**Or draw a structure:**")
                ketcher_html = gr.HTML(KETCHER_HTML)

                gr.Markdown(
                    "**Or upload a MOL/PDB file:** MOL files must include all "
                    "hydrogens. For a PDB upload, also supply the SMILES above "
                    "-- it's used to fix bond orders and add any missing "
                    "hydrogens (PDB files don't encode bond order and are "
                    "often missing Hs)."
                )
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
