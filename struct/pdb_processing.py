import os
from variables import address_dict, subfolders, hetatm_non_metal_ion, hetatm_metal_ion, hetatm_anion
from utils import fetch_sequences_from_fasta
import yasara


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

def get_chains():
    chain_list = list(set(yasara.NameMol('all')))
    chain_list.sort()
    return chain_list

def parse_objs_format_scene(res_to_keep=hetatm_non_metal_ion+hetatm_metal_ion, res_to_delete=hetatm_anion, sce_fpath=None):
    # split all molecules into separate objects
    yasara.SplitObj(1)
    # get a dict of object number to res list
    obj_list = yasara.ListObj('all')
    obj_dict = {}
    for obj_num in obj_list:
        non_protein_res = yasara.NameRes(f'Obj {obj_num} & not Protein')
        if len(non_protein_res) > 0:
            obj_dict.update({obj_num:non_protein_res})
    print(obj_dict)
    obj_to_delete = []
    obj_with_ligand = []
    for obj_num, res_list in obj_dict.items():
        for res in res_list:
            if res not in res_to_keep:
                if res in res_to_delete:
                    obj_to_delete.append(obj_num)
                else:
                    obj_with_ligand.append(obj_num)
    obj_to_delete = list(set(obj_to_delete))
    obj_with_ligand = list(set(obj_with_ligand))
    obj_to_join = [obj_num for obj_num in list(obj_dict.keys()) if obj_num not in obj_to_delete+obj_with_ligand]
    print('Object # to delete:', *obj_to_delete)
    print('Object # with ligand:', *obj_with_ligand)
    print('Object # to join with protein:', *obj_to_join)

    # delete objects to delete
    for obj_num in obj_to_delete:
        yasara.DelObj(obj_num)
    # join all except ligand object to protein object number (1)
    for obj_num in obj_to_join:
        yasara.JoinObj(obj_num,1)
    # swapres for ligand object #2
    yasara.SwapObj(obj_with_ligand[0], 2)

    # print result
    print('After processing objects:')
    for obj_num in yasara.ListObj('all'):
        print(f"Obj # {obj_num}: {yasara.NameRes(f'Obj {obj_num}')}")

    # save final Scene
    if sce_fpath is not None:
        yasara.SaveSce(sce_fpath)


def split_into_separate_chains(struct_fpath, chain_list=None, pdb_dir=None, sce_dir=None, struct_suffix=''):
    struct_name = os.path.basename(struct_fpath).split('.')[0]
    if chain_list is None:
        chain_list = get_chains()
    print('Chain list:', *chain_list)
    # iterate through chains and save individually
    for chain in chain_list:
        print(chain)
        struct_suffix_wchain = struct_suffix
        if len(chain_list)>1:
            struct_suffix_wchain = f'-{chain}{struct_suffix}'
        yasara.Clear()
        load_struct(struct_fpath, delete_water=True)
        # delete all other chains but chain of interest
        yasara.DelMol(f'not {chain}')
        # save as PDB
        if pdb_dir is not None:
            yasara.SavePDB(1, f'{pdb_dir}{struct_name}{struct_suffix_wchain}.pdb')
        # save as SCE
        if sce_dir is not None:
            parse_objs_format_scene(
                res_to_keep=hetatm_non_metal_ion+hetatm_metal_ion,
                res_to_delete=hetatm_anion,
                sce_fpath=f'{sce_dir}{struct_name}{struct_suffix_wchain}.sce'
            )


if __name__ == '__main__':
    # Define input files and output file
    data_folder = address_dict['influenza-resistance'] # address_dict['PIPS2'] #
    data_subfolder = 'PA_benchmark' # 'ET096' #
    pdb_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    sce_dir = data_folder + subfolders['sce'] + data_subfolder + '/'
    pdb_fname_list = [
        'PA-H1N1-I38T_BXR_6fs7.pdb',
        'PA-B-WT_BXR_6fs8.pdb',
        'PA-B-I38T_BXR_6fs9.pdb'
        ]
    pdb_fpath_list = [f'{pdb_dir}/{pdb_fname}' for pdb_fname in pdb_fname_list]

    # initialize yasara
    initialize_yasara()

    # load PDB
    for pdb_fpath in pdb_fpath_list:
        print(pdb_fpath)

        # load structure
        load_struct(pdb_fpath, delete_water=True)

        # split into separate chains
        chain_list = get_chains()
        split_into_separate_chains(pdb_fpath, chain_list, pdb_dir, sce_dir, struct_suffix='_loopdel')






