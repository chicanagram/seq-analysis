import subprocess
import time
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta
import yasara

if __name__ == '__main__':
    # Define input files and output file
    data_folder = address_dict['ECOHARVEST'] # address_dict['PIPS2'] #
    data_subfolder = 'lipases' # 'ET096' #
    pdb_fname = 'DockedRMLOleicAcidpreOpt.pdb' # 'DockedET096S82.pdb' #
    pdb_fpath = data_folder + subfolders['pdb'] + data_subfolder + '/' + pdb_fname

    ligand_res = '0A'.split(' ') # '0A 244A'.split(' ') #
    res_to_mutate = '82A 83A 88A 91A 92A 94A 108A 204A 205A 208A 254A 258A 265A 266A 267A'.split(' ') # '60A 64A 70A 73A 74A 75A 76A 77A 78A 80A 81A 178A 218A 220A 223A'.split(' ') #
    res_to_fix = '28A 81A 111A 143A 144A 145A 177A 178A 209A 257A'.split(' ') # '36A'.split(' ') #
    res_all = res_to_fix + res_to_mutate + ligand_res

    # initialize yasara
    yasara.info.mode = 'txt'
    yasara.Console('Off')
    yasara.Clear()

    # load PDB
    yasara.LoadPDB(pdb_fpath)

    for res_chain in res_all:
        res = int(res_chain[:-1])
        chain = res_chain[-1]
        aa = yasara.SequenceRes(f'Mol {chain} & Res {res}')
        # aa = yasara.ListRes(f'Mol {chain} & Res {res}')
        print(res_chain, aa)


