import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from variables import address_dict, subfolders, mapping_rev
from get_surface_residues import yasara_get_surface_residues

class MolSelect(Select):
    def __init__(self, mol_type):
        self.mol_type = mol_type
    def accept_residue(self, residue):
        if self.mol_type=='Protein':
            return is_aa(residue, standard=True)
        elif self.mol_type=='Heme':
            # return (not is_aa(residue, standard=True)) and residue.id[0]=='H_HEM'
            return (not is_aa(residue, standard=True)) and residue.get_resname() == 'HEM'
        elif self.mol_type=='MetalIon':
            return (not is_aa(residue, standard=True)) and residue.get_resname() in ['MG', 'ZN', 'CA', 'FE', 'MN', 'CO', 'CU']
        elif self.mol_type=='Ligand':
            return (not is_aa(residue, standard=True)) and residue.id[0]=='UNK'

def parse_protein_hetatom(pdb_fpath, data_fbase=''):
    from Bio import PDB
    repo = PDB.PDBList()
    repo.retrieve_pdb_file(data_fbase, file_format="pdb")
    pdb = PDBParser().get_structure(data_fbase, pdb_fpath)
    io = PDBIO()
    io.set_structure(pdb)
    for mol_type in ['Protein', 'Ligand']: # ['Protein', 'Heme', 'MetalIon', 'Ligand']:
        io.save(pdb_fpath.replace('.pdb', f'_{mol_type}.pdb'), MolSelect(mol_type))
        print(f'Saved {mol_type} pdb.')

def get_dssp_properties(pdb_fpath, data_fbase):
    from Bio.PDB.DSSP import DSSP
    dssp_property_names = [
        'RealPos',
        'aa',
        'secondary_structure', # {H,B,E,G,I,P,T,S} H = α-helix; B = residue in isolated β-bridge; E = extended strand, participates in β ladder; G = 310-helix; I = π-helix; P = κ-helix (poly-proline II helix); T = hydrogen-bonded turn; S = bend
        'relative ASA',  # relative accessible solvent area
        'phi',  # peptide backbone torsion angle phi
        'psi',  # peptide backbone torsion angle psi
        'NH_O_1_relidx',  # relative index of H-bond 1 (between N-H group of this residue with O of another residue)
        'NH_O_1_energy',  # energy of H-bond 1 (between N-H group of this residue with O of another residue)
        'O_NH_1_relidx',  # relative index of H-bond 1 (between O group of this residue with N-H of another residue)
        'O_NH_1_energy',  # energy of H-bond 1 (between O group of this residue with N-H of another residue)
        'NH_O_2_relidx',  # relative index of H-bond 2 (between N-H group of this residue with O of another residue)
        'NH_O_2_energy',  # energy of H-bond 2 (between N-H group of this residue with O of another residue)
        'O_NH_2_relidx',  # relative index of H-bond 2 (between O group of this residue with N-H of another residue)
        'O_NH_2_energy'  # energy of H-bond 2 (between O group of this residue with N-H of another residue)
    ]
    # get model with protein chain
    p = PDBParser()
    structure = p.get_structure(data_fbase, pdb_fpath)
    model = structure[0]
    chain_ids = [chain.get_id() for chain in model.get_chains()]
    print(chain_ids)
    protein_chain_id = chain_ids[0]
    chain_ids_to_remove = [id for id in chain_ids if id!=protein_chain_id]
    for id in chain_ids_to_remove:
        model.detach_child(id)
        print(f'Removed chain id {id}.')
    num_res = len(list(model[protein_chain_id].get_residues()))

    # get DSSP properties for each residue
    dssp = DSSP(model, pdb_fpath, dssp='mkdssp')
    dssp_res = []
    for res_num in range(1,num_res+1):
        res_key = (protein_chain_id, (' ', res_num, ' '))
        res_vals = dssp[res_key]
        res_dict = {property_name: val for property_name, val in zip(dssp_property_names, res_vals)}
        dssp_res.append(res_dict)
    dssp_res = pd.DataFrame(dssp_res)

    return dssp_res
def calc_res_to_lig_distances(pdb_fpath, chain_id='A', ligand_resname='UNK'):

    # Parse PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_fpath)
    model = structure[0]

    # Select ligand atoms
    ligand_atoms = []
    model = structure[0]
    for chain in model:
        for residue in chain:
            # select target (ligand, heme, metal ion) atoms only, excluding H
            if residue.get_resname().strip() == ligand_resname:
                ligand_atoms.extend([atom for atom in residue if atom.element != 'H'])

    # Compute distances between ligand atoms and heavy atoms of each residue
    distance_results = []
    for chain in model:
        # match the
        if chain.id == chain_id:
            for residue in chain:
                # select only the protein residues
                if is_aa(residue, standard=True):
                    res_id = (residue.get_resname(), residue.get_id()[1])
                    aa = mapping_rev[res_id[0]]
                    min_distances = []
                    avg_distances = []
                    # iterate through atoms in the protein residue, excluding H atoms
                    for atom in residue:
                        dist_res_atom_to_all_target_atoms = []
                        if atom.element != 'H':  # exclude hydrogen atoms
                            # calculate distances between all residue atoms and target atoms
                            for lig_atom in ligand_atoms:
                                dist_atom2atom = np.linalg.norm(atom.coord - lig_atom.coord)
                                dist_res_atom_to_all_target_atoms.append(dist_atom2atom)
                            avg_dist_res_atom_to_target = np.mean(np.array(dist_res_atom_to_all_target_atoms))
                            avg_distances.append(avg_dist_res_atom_to_target)
                            min_distances.append(min(dist_res_atom_to_all_target_atoms))
                    distance_results.append({'RealPos':res_id[1], 'aa': aa, f'min_distance_to_{ligand_resname}':min(min_distances), f'avg_distance_{ligand_resname}':min(avg_distances)})

    distance_results = pd.DataFrame(distance_results)
    return distance_results

def get_remaining_positions_residues(struct_df):
    filt_pos = struct_df["RealPos"].tolist()
    filt_aa = struct_df["aa"].tolist()
    filt_res = [str(pos)+aa for pos,aa in zip(filt_pos,filt_aa)]
    print(*filt_res)
    return filt_pos, filt_res

if __name__ == '__main__':
    # input settings
    data_folder = address_dict['PIPS2'] # address_dict['ECOHARVEST'] #
    data_fbase = 'ET096' # 'lipases' #
    pdb_fname =  'Docked_ET096_S82.pdb'  # 'Docked_RML_OleicAcid_preOpt.pdb' # # 'ETS83096.pdb' #
    ligand_resname_list =  ['UNK'] # ['HEM', 'UNK'] #
    conservation_analysis_fname = 'ET096_UPO_aligned_clustalo_sift_selected.csv' # 'RML_blastp_nr_E1e-05_mafft_sift_selected.csv'
    pdb_fpath = f"{data_folder}{subfolders['pdb']}{data_fbase}/{pdb_fname}"
    conservation_analysis_fpath = f"{data_folder}{subfolders['conservation_analysis']}{data_fbase}/{conservation_analysis_fname}"

    # get residue to ligand distances
    for ligand_resname in ligand_resname_list:
        struct_df = calc_res_to_lig_distances(pdb_fpath, chain_id='A', ligand_resname=ligand_resname)

    # # get residue to surface distances, and surface patches
    # surfdist = yasara_get_surface_residues(pdb_fpath, surft=2.55, distt=12, minepisize=13)
    # struct_df = struct_df.merge(surfdist, on='RealPos', how='left')

    # # parse structure into protein vs. non protein
    # parse_protein_hetatom(pdb_fpath, data_fbase)
    # # get DSSP properties
    # dssp_df = get_dssp_properties(pdb_fpath.replace('.pdb', '_Protein.pdb'), data_fbase)
    # struct_df = struct_df.merge(dssp_df, on='RealPos', how='left')

    # combine with conservation analysis scores and shortlist residues
    if conservation_analysis_fpath is not None:
        conservation_df = pd.read_csv(conservation_analysis_fpath)[['RealPos', 'entropy', 'sift_avg']]
        struct_df = struct_df.merge(conservation_df, on='RealPos', how='left')

    struct_df = struct_df.round(3)
    print(struct_df)
    struct_df.to_csv(f"{data_folder}{subfolders['pdb']}{data_fbase}/{data_fbase}_StructProperties.csv")

    # filter
    conservation_score_label = 'sift_avg'
    frac_positions_to_exclude = 0.15
    dist_to_ligand_thres = 6.65
    pos_suffix = 'A'
    conservation_scores_sorted = np.sort(struct_df[conservation_score_label].to_numpy())
    num_positions_to_exclude = int(frac_positions_to_exclude*len(conservation_scores_sorted))
    conservation_score_thres = conservation_scores_sorted[num_positions_to_exclude]
    print('\n# of positions to exclude by conservation:', f'{num_positions_to_exclude}/{len(conservation_scores_sorted)}', '; conservation score threshold:', conservation_score_thres)
    struct_df_filt = struct_df.copy()
    struct_df_filt = struct_df_filt.loc[(struct_df_filt['min_distance_to_UNK']<dist_to_ligand_thres)]
    filt_pos_1, filt_res_1 = get_remaining_positions_residues(struct_df_filt)
    print('\n', f'{len(filt_pos_1)} positions after filtering for ligand distance.', *filt_pos_1)
    struct_df_filt = struct_df_filt.loc[(struct_df_filt[conservation_score_label] > conservation_score_thres)]
    filt_pos_2, filt_res_2 = get_remaining_positions_residues(struct_df_filt)
    print('\n', f'{len(filt_pos_2)} positions after filtering for ligand distance and conservation.', *[str(pos)+pos_suffix for pos in filt_pos_2])
    print('\n', f'{len(filt_pos_1)-len(filt_pos_2)} conserved positions removed:', *[pos for pos in filt_pos_1 if pos not in filt_pos_2])
    print(*[res for res in filt_res_1 if res not in filt_res_2])



