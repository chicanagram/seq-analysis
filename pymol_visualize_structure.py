import os
import pymol
from pymol import cmd
from variables import address_dict, subfolders

def launch_pymol():
    pymol.finish_launching()
def load_pdb_pymol(pdb_abs_fpath):
    cmd.load(pdb_abs_fpath)
    cmd.zoom()
    cmd.png(pdb_abs_fpath.replace('.cif', '.png'), 300, 200)

def rename_res(sel_str_old, sel_str_new, sel_name):
    cmd.select(sel_name, f'resn {sel_str_old}')
    cmd.alter('lig_to_rename', f"resn = '{sel_str_new}'")

def save_pdb_pymol(pdb_fpath):
    cmd.save(pdb_fpath, selection='(all)')

if __name__=='__main__':

    data_folder = address_dict['influenza-resistance']
    data_subfolder = 'PA-A37T-I38T/PA-H3N2-Chai1'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    pdb_fname = 'pred.rank_0.cif'
    pdb_fpath = pdb_dir + pdb_fname
    pdb_abs_fpath = os.path.abspath(pdb_fpath)
    print(pdb_abs_fpath)

    # launch_pymol()
    load_pdb_pymol(pdb_fpath)
    rename_res("LIG2", 'E4Z', 'lig_to_rename')
    rename_res("LIG3", ' MN', 'lig_to_rename')
    rename_res("LIG4", ' MN', 'lig_to_rename')
    save_pdb_pymol(pdb_abs_fpath.replace('.cif', '.pdb'))
