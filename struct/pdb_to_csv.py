import os
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List
from variables import address_dict, subfolders, aaList, mapping_rev, \
    kyte_doolittle_hydrophobicity_index, hopp_woods_polarity_index, \
    aa_sidechain_volume, aa_polarity_mapping


def parse_pdb_atom_line(line: str) -> Dict[str, Optional[object]]:
    """
    Parse a single ATOM/HETATM line from a PDB file using fixed-width columns.
    Returns a dict of fields. Missing numeric fields become None.
    """
    line = line.rstrip("\n")
    padded = line + " " * (80 - len(line)) if len(line) < 80 else line

    def to_int(s: str) -> Optional[int]:
        s = s.strip()
        return int(s) if s else None

    def to_float(s: str) -> Optional[float]:
        s = s.strip()
        return float(s) if s else None

    # get residue
    res_name = padded[17:20].strip()
    res_symbol = mapping_rev[res_name] if res_name in mapping_rev else None

    # compose row_dict
    row_dict = {
        "serial": to_int(padded[6:11]),           # 7–11
        "res_num": to_int(padded[22:26]),         # 23–26
        "res_name": res_name,  # 18–20
        "res": res_symbol,
        "atom_name": padded[12:16].strip(),       # 13–16
        "chain_id": padded[21:22].strip() or None,# 22
        "x": to_float(padded[30:38]),             # 31–38
        "y": to_float(padded[38:46]),             # 39–46
        "z": to_float(padded[46:54]),             # 47–54
        "occupancy": to_float(padded[54:60]),     # 55–60
        "temp_factor": to_float(padded[60:66]),   # 61–66
        "element": padded[76:78].strip() or None, # 77–78

    }
    # update aa properties
    if res_symbol in aaList:
        row_dict.update({
            "aa_polarity": aa_polarity_mapping[res_symbol],
            "kd_hydro": kyte_doolittle_hydrophobicity_index[res_symbol],
            "hw_polarity": hopp_woods_polarity_index[res_symbol],
            "aa_vol": aa_sidechain_volume[res_symbol],
        })

    return row_dict


def pdb_to_dataframe(pdb_path: Path) -> pd.DataFrame:
    """
    Parse ATOM and HETATM records from a PDB file into a DataFrame.
    """
    rows: List[Dict[str, Optional[object]]] = []

    with pdb_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            rec = line[0:6].strip()
            if rec in {"ATOM", "HETATM"}:
                rows.append(parse_pdb_atom_line(line))

    return pd.DataFrame(rows)


if __name__ == "__main__":
    os.chdir('../')

    # ---- user input ----
    data_folder = address_dict['PIPS2']
    data_subfolder = 'UPOs_peroxygenation_analysis/structure_csv/' # 'CARs' # 'sidestream_cocktail' #
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder
    pdb_fname_list = [f for f in os.listdir(pdb_dir) if f.find('.pdb')>-1]

    # process files
    for pdb_fname in pdb_fname_list:
        print(pdb_fname)

        # get filepaths
        pdb_fpath = Path(pdb_dir + pdb_fname)
        out_csv = str(pdb_fpath).replace('.pdb', '.csv')

        # process pdb
        df = pdb_to_dataframe(pdb_fpath)
        df.to_csv(out_csv, index=False)

        # get backbone of protein only
        df_backbone = df[df['atom_name']=='CA']
        df_backbone.to_csv(out_csv.replace('.csv','_backbone.csv'), index=False)

        print(f"Parsed {len(df)} atoms")
        print(f"Saved CSV to: {out_csv}")