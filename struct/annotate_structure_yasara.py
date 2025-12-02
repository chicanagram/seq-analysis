import os

import pandas as pd

import yasara
from yasara import yasaradir
from variables import address_dict, subfolders


if __name__=='__main__':
    os.chdir('../')

    # Path to the FASTA file
    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'lipases' # 'CARs' #
    struct_fname = 'Docked_RMLmut_SucroseOleate_Chai.pdb'  # 'Boltz2_MpCAR-A_Oleoyl-AMP.pdb' #
    csv_fname = 'RML-propeptide-mature_ss_nondistal_singlesite_SELECTED.csv' # 'MpCAR-A_nondistal_singlesite_SELECTED.csv' #
    struct_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    csv_dir = data_folder + subfolders['mutagenesis_proposal'] + data_subfolder + '/'
    res_col = 'mutations' # Position
    residues_to_annotate = None # [75, 214, 155, 509, 261, 130, 389, 394, 63, 106, 152, 133, 307, 565, 465, 280, 154, 285, 172, 100, 388, 213, 574, 201, 303, 23, 384, 561, 207, 224, 18, 628]
    show_only_annnotated_residues = True
    pos_offset_from_pdb = 71 # 0 #

    # get residues to annotate from CSV
    if csv_fname is not None and residues_to_annotate is None:
        csv_fpath = csv_dir + csv_fname
        residues_to_annotate = pd.read_csv(csv_fpath)[res_col].tolist()
    print('residues_to_annotate:', residues_to_annotate)

    # initialize yasara
    yasara.info.mode = 'txt'
    yasara.Console('Off')
    yasara.Clear()

    # load structure
    if struct_fname.find('.pdb') > -1:
        yasara.LoadPDB(filename=f'{struct_dir}{struct_fname}')
    elif struct_fname.find('.sce') > -1:
        yasara.LoadPDB(filename=f'{struct_dir}{struct_fname}')

    # hide all residues
    if show_only_annnotated_residues:
        yasara.HideRes('protein')

    # iterate through residues
    for res in residues_to_annotate:
        # check if residue or mutation
        if isinstance(res, int):
            pos = res
        elif isinstance(res, str):
            pos = int(res[1:-1])
        # adjust for position offset
        pos_in_pdb = pos - pos_offset_from_pdb
        yasara.ShowRes(f'protein and {pos_in_pdb}')
        if isinstance(res, int):
            yasara.LabelRes(f'protein and {pos_in_pdb}', pos)
        elif isinstance(res, str):
            yasara.LabelRes(f'protein and {pos_in_pdb}', res)
        print(f'Annotated {res} (Res {pos_in_pdb} in PDB).')

    # save structure as scene
    yasara.SaveSce(f'{struct_dir}{struct_fname}'.replace('pdb','sce').replace('.sce', '_annotated.sce'))



