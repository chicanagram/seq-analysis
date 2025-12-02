import subprocess
import os
from variables import address_dict, subfolders
chimerax_dir = "/Applications/ChimeraX-1.10.1.app/Contents/MacOS/"
def main():
    # --- Configuration ---
    os.chdir('../')
    data_folder = address_dict['PON-Sol2']
    data_subfolder = ''
    cif_dir = data_folder + 'cif/' + data_subfolder + '/'
    pdb_dir = cif_dir.replace('cif/', 'pdb/')
    cif_fname_list = [f for f in os.listdir(cif_dir) if f.find('.cif')>-1]
    print(cif_fname_list)

    # --- Run ChimeraX ---
    print("Converting CIF to PDB...")

    for cif_fname in cif_fname_list:
        cif_fpath = cif_dir + cif_fname
        pdb_fpath = pdb_dir + cif_fname.replace('.cif', '.pdb').replace('_model_0','')
        cmd = [
            os.path.join(chimerax_dir, "ChimeraX"),
            "--nogui",
            "--cmd", "pwd",
            "--cmd", f"open {cif_fpath}",
            "--cmd", f"save {pdb_fpath}",
            "--cmd", "exit"
        ]
        try:
            subprocess.run(cmd, check=True)
            print(f"CIF > PDB complete: {pdb_fpath.replace('.pdb', '')}.")
        except subprocess.CalledProcessError as e:
            print(f"ChimeraX conversion failed: {e}")
        except FileNotFoundError:
            print("Could not find ChimeraX executable. Check your `chimerax_dir` path.")

if __name__ == "__main__":
    main()