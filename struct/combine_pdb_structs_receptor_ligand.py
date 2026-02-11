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
    struct_name = 'TE314'
    data_subfolder = f'UPOs_peroxygenation_analysis/docked/swissdock/{struct_name}_S82_swissdock'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    struct_1_list = [f for f in os.listdir(pdb_dir) if (f.find('.pdb')>-1 & f.find('lig')==-1 & f.find('dock4')==-1 & f.find('TEST')==-1)]
    struct_1_list.sort()
    struct_2_list = [f for f in os.listdir(pdb_dir) if (f.find('.pdb') > -1 & f.find('lig') > -1)]
    struct_2_list.sort()

    # get struct pairs
    struct_pairs = [(struct_1, struct_2) for struct_1 in struct_1_list for struct_2 in struct_2_list]

    # --- Run YASARA ---
    initialize_yasara()

    print("Combining struct files...")
    for (struct_1_fname, struct_2_fname) in struct_pairs:
        struct_1 = struct_1_fname.replace('.pdb','')
        struct_2 = struct_2_fname.replace('.pdb', '')
        print(struct_1, struct_2)
        yasara.Clear()
        yasara.LoadPDB(pdb_dir + struct_1_fname, transfer=True)
        yasara.LoadPDB(pdb_dir + struct_2_fname, transfer=True)
        yasara.JoinObj(2,1, center=False)
        save_fname = f'{struct_1}_S82_swissdock_{struct_2.split("_")[1].split("-")[0]}.pdb'
        yasara.SavePDB(1, pdb_dir + save_fname, transform=False)
        print(f'Saved: {pdb_dir + save_fname}')
if __name__ == "__main__":
    main()

    # SavePDB All, TEST.pdb,Format=PDB,Transform=Yes,UseCIF=No - OK