import os
import pandas as pd
import numpy as np
pd.set_option('display.max_columns',None)
import yasara
from variables import address_dict, subfolders
from utils import opsys, run_msa, fetch_sequences_from_fasta


if __name__ == '__main__':
    data_folder = address_dict['PIPS2']
    pdb_dir = data_folder + subfolders['pdb']
    data_fbase = 'UPOpanel_selected'
    struct_fnames = [
        'ET096.pdb',
        'CviUPO.pdb',
    ]
    pos_to_hybridize = [38, 39, 60, 64, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 103, 171, 174, 175, 176, 177, 178, 179, 182, 218, 219, 220, 223]
    struct_ali_csv = f'{pdb_dir}{data_fbase}/{struct_fnames[0].split('.')[0]}_{struct_fnames[1].split('.')[0]}_yasaraStructAli.csv'
    pdb_hybridized_fpath = struct_ali_csv.replace('_yasaraStructAli.csv', '_hybrid.pdb')

    # read CSV of aligned structures
    df = pd.read_csv(struct_ali_csv, index_col=0)

    # hybridize the two structures at the selected positions
    df_hybridized = df.copy()
    # update sequence AA
    df_hybridized.loc[df_hybridized['pos_0'].isin(pos_to_hybridize), 'seq_0'] = df_hybridized.loc[df_hybridized['pos_0'].isin(pos_to_hybridize), 'seq_1']
    # update PDB segments
    for pos_0 in pos_to_hybridize:
        mask = df_hybridized['pos_0']==pos_0
        pos_0_str = ' ' + str(pos_0) + ' '
        pos_1_str = ' ' + str(df_hybridized.loc[mask, 'pos_1'].iloc[0]) + ' '
        print(pos_0_str, pos_1_str)
        pdb_str_residue = df_hybridized.loc[mask, 'pdb_by_residue_1'].iloc[0]
        pdb_residue_str_mod = pdb_str_residue.replace(pos_1_str, pos_0_str)
        df_hybridized.loc[mask, 'pdb_by_residue_0'] = pdb_residue_str_mod
    df_hybridized_list = df_hybridized.loc[df_hybridized['pos_0'] != -1, ['pdb_by_residue_0']]['pdb_by_residue_0'].tolist()

    # update atom numbering
    atom_count = 0
    atom_str_list = []
    atom_attribute_array = []
    spaces = [(6, 'L', 0), (7, 'R', 2), (4, 'L', 0), (4, 'L', 0), (2, 'L', 0), (7, 'R', 4), (8, 'R', 0), (8, 'R', 0),
              (8, 'R', 0), (6, 'R', 0), (17, 'R', 11), (1, 'L', 0)] # (length, alignment, offset)
    # 'ATOM      1  N   MET A   1      73.504 -80.983   2.168  1.00114.13           N'
    # iterate through atoms in protein
    for i, residue_str in enumerate(df_hybridized_list):
        atom_str_list_byresidue = residue_str.split('\n')
        residue_str_updated = residue_str
        for atom_str in atom_str_list_byresidue:
            atom_count += 1
            # detect atom number
            atom_str_parsed = atom_str.split()
            atom_num_original = atom_str_parsed[1]
            # update atom number in parsed list and string
            atom_str_parsed_updated = atom_str_parsed.copy()
            atom_str_parsed_updated[1] = str(atom_count)

            # fix spaces in string
            atom_str_updated = ''
            for j, atom_property in enumerate(atom_str_parsed_updated):
                (num_spaces, alignment, offset) = spaces[j]
                atom_str_property = [' ']*num_spaces
                entry_aslist = list(atom_property)
                entry_len = len(entry_aslist)
                if alignment=='L':
                    start_idx = offset
                elif alignment=='R':
                    start_idx = num_spaces-offset-entry_len
                atom_str_property[start_idx:start_idx+entry_len] = entry_aslist
                atom_str_updated += ''.join(atom_str_property)
            atom_str_updated = atom_str_updated.replace('           ', '     ')
            # print(atom_str_updated)

            # update atom_str_list
            atom_str_list.append(atom_str_updated)

            # append to atom_attribute_array
            if len(atom_str_parsed) >= 9:
                atom_str_toappend = atom_str_parsed[:9]
            else:
                atom_str_toappend = atom_str_parsed + ['']*(9-len(atom_str_parsed))
            atom_attribute_array.append(atom_str_toappend)

    atom_attribute_array = np.array(atom_attribute_array)
    atoms_df = pd.DataFrame(atom_attribute_array, columns=['ATOM', 'atom_num', 'element_w_idx', 'amino_acid', 'chain', 'residue_num', 'x', 'y', 'z'])

    # save the hybridized PDB
    pdb_hybridized = '\n'.join(atom_str_list)
    with open(pdb_hybridized_fpath, 'w') as f:
        f.write(pdb_hybridized)
    print('Written hybrid PDB file.')

    # run YASARA energy minimization
    if opsys == 'Windows':
        process_name = 'YASARA.exe'
    else:
        process_name = 'yasara'

    # set parameters
    move = 'all' # '!backbone
    ff = 'YASARA2' # 'AMBER15FB'
    # initialize YASARA
    yasara.info.mode = 'txt'
    yasara.info.licenseshown = 0
    time = yasara.SystemTime()
    yasara.info.pid = int(time)+5
    yasara.Console('Off')
    yasara.Processors(1)
    yasara.EnergyUnit('kcal/mol')
    yasara.ForceField(ff, setpar='yes')
    yasara.Clear()

    # load PDB
    yasara.LoadPDB(pdb_hybridized_fpath)
    print('Loaded PDB in YASARA')
    yasara.CleanAll()
    yasara.CellAuto(extension='10')
    yasara.FixAll()

    # perform minimization
    print('Performing energy minimization')
    yasara.ExperimentMinimization(convergence=0.01)
    yasara.Experiment('On')
    yasara.Wait('ExpEnd')

    # save minimized PDB
    yasara.SavePDB(1, pdb_hybridized_fpath.replace('.pdb', '_minimized.pdb'))
    print('Completed minimization and saved updated PDB.')
