import os
import shutil
import yasara
from yasara import yasaradir
from variables import address_dict, subfolders


if __name__=='__main__':
    # Path to the FASTA file
    data_folder = address_dict['influenza-resistance']
    data_subfolder = 'PA_benchmark_mod'
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'

    folder_list = [
        # 'PA-H1N1-E23K_BXR_8t81',
        # 'PA-H1N1-E23K_cpd23_8t81',
        # 'PA-H1N1-I38T_BXR_8t5w',
        # 'PA-H1N1-WT_BXR_6fs6-A',
        # 'PA-H1N1-WT_BXR_6fs6-B',
        # 'PA-H1N1-WT_BXR_6fs6-C',
        # 'PA-H1N1-WT_BXR_6fs6-D',
        # 'PA-H1N1-WT_BXR_6fs6-E',
        # 'PA-H1N1-WT_BXR_6fs6-F',
        # 'PA-H1N1-I38T_cpd23_8t6z',
        # 'PA-H1N1-WT_cpd23_8t94'
    ]

    yasara.info.mode = 'txt'
    yasara.Console('Off')

    for folder in folder_list:
        target = [f for f in os.listdir(pdb_dir+folder) if f.find('.fasta')>-1][0]
        print(target)
        # set target
        yasara.Clear()
        yasara.MacroTarget(filename=f'{pdb_dir}{folder}/{target}')
        yasara.PlayMacro(f'{yasaradir}mcr/hm_build_1template.mcr', label='')


    # clean up results folder
    sce_flist = [f for f in os.listdir(pdb_dir) if f.find('.sce')>-1]
    # save all SCE files as PDB in current folder
    for sce_fname in sce_flist:
        print(sce_fname)
        sce_fpath = pdb_dir + sce_fname
        yasara.Clear()
        yasara.LoadSce(sce_fpath)
        yasara.SavePDB(1, sce_fpath.replace('sce','pdb'))

        # convert all sce files to right format
        yasara.SplitObj(1)
        yasara.JoinObj(2,1)
        yasara.SwapObj(3,2)
        yasara.SaveSce(sce_fpath.replace('pdb','sce'))

