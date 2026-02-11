# Re-run the provided parsing code carefully on the uploaded file,
# extract energies directly from REMARK lines, and regenerate outputs.

from pathlib import Path
import os
import re
import csv
from variables import address_dict, subfolders


def parse_swissdock_ligand_file(pdb_path, out_dir):
    """
    Parse swissdock output file "result_dock4" containing all the ligand coordinates, to get the top pose for the top 5 clusters
    """

    # regex
    re_int = {
        "CLUSTER_NUM": re.compile(r"\bCLUSTER_NUM\b\s*[:=]?\s*(-?\d+)\b"),
        "CLUSTER_MEMBER": re.compile(r"\bCLUSTER_MEMBER\b\s*[:=]?\s*(-?\d+)\b"),
    }
    re_float = {
        "MEMBER_ENERGY": re.compile(r"\bMEMBER_ENERGY\b\s*[:=]?\s*(-?\d+(?:\.\d+)?)"),
        "SP-dG": re.compile(r"\bSP-dG\b\s*[:=]?\s*(-?\d+(?:\.\d+)?)"),
    }

    # split on TER
    lines = pdb_path.read_text().splitlines(True)
    blocks, cur = [], []
    for l in lines:
        cur.append(l)
        if l.startswith("TER"):
            blocks.append(cur)
            cur = []
    if cur:
        blocks.append(cur)

    records = []
    selected_blocks = {}

    for i, block in enumerate(blocks):
        cn = cm = None
        me = sdg = None
        for l in block:
            if l.startswith("REMARK"):
                if cn is None:
                    m = re_int["CLUSTER_NUM"].search(l)
                    if m: cn = int(m.group(1))
                if cm is None:
                    m = re_int["CLUSTER_MEMBER"].search(l)
                    if m: cm = int(m.group(1))
                if me is None:
                    m = re_float["MEMBER_ENERGY"].search(l)
                    if m: me = float(m.group(1))
                if sdg is None:
                    m = re_float["SP-dG"].search(l)
                    if m: sdg = float(m.group(1))

        records.append((i, cn, cm, me, sdg))

        if cn in {0,1,2,3,4} and cm == 1 and me is not None and sdg is not None:
            if me < 0 and sdg < 0:
                selected_blocks[cn] = (block, me, sdg)

    # write CSV
    csv_path = out_dir / "dock_parsed.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["pose_index","CLUSTER_NUM","CLUSTER_MEMBER","MEMBER_ENERGY","SP-dG"])
        for r in records:
            w.writerow(r)

    # write pdbs
    written = {}
    for cn, (block, me, sdg) in selected_blocks.items():
        p = out_dir / f"lig_{cn}-1.pdb"
        with open(p, "w") as f:
            f.writelines(block)
        written[cn] = (p, me, sdg)

if __name__=='__main__':
    # --- Configuration ---
    os.chdir('../')
    data_folder = address_dict['PIPS2'] # address_dict['PON-Sol2']
    struct_name = 'CviUPO'
    data_subfolder = f'UPOs_peroxygenation_analysis/docked/swissdock/{struct_name}_S82_swissdock'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    ligand_result_fname = 'result_dock4.pdb'
    pdb_path = Path(f'{pdb_dir}{ligand_result_fname}')
    out_dir = Path(pdb_dir)
    out_dir.mkdir(exist_ok=True)

    # parse docked poses
    parse_swissdock_ligand_file(pdb_path, out_dir)