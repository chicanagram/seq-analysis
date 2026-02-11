import os
import yasara
from variables import address_dict, subfolders
chimerax_dir = "/Applications/ChimeraX-1.10.1.app/Contents/MacOS/"

def initialize_yasara():
    yasara.info.mode = 'txt'
    yasara.Console('Off')

def load_struct(struct_fpath, delete_water=True):
    yasara.Clear()
    if struct_fpath.find('.sce')>-1:
        yasara.LoadSce(struct_fpath)
    elif struct_fpath.find('.pdb')>-1:
        yasara.LoadPDB(struct_fpath)
    if delete_water:
        yasara.DelWater()

def main():
    # --- Configuration ---
    os.chdir('../')
    data_folder = address_dict['PIPS2'] # address_dict['PON-Sol2']
    data_subfolder = f'UPOs_peroxygenation_analysis/docked/ALL_aligned/'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    fname_list = [f for f in os.listdir(pdb_dir) if (f.find('.pdb')>-1 & f.find('swissdock')>-1)]

    # --- Run YASARA ---
    initialize_yasara()

    print("Combining struct files...")
    for fname in fname_list:
        struct_name = fname.replace('.pdb','')
        print(struct_name)
        yasara.Clear()
        yasara.LoadPDB(pdb_dir + fname, transfer=True)
        yasara.NameMol('protein', 'A')
        yasara.NameMol('not protein and Res HEM', 'B')
        yasara.NameMol('not protein and Res Mg', 'C')
        yasara.NameMol('not A and not B and not C', 'D')
        yasara.SavePDB(1, pdb_dir + fname, transform=False)
        print(f'Saved: {pdb_dir + fname}')
if __name__ == "__main__":
    main()

    # SavePDB All, TEST.pdb,Format=PDB,Transform=Yes,UseCIF=No - OK