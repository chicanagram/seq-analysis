import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from variables import address_dict, subfolders
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
    for mol_type in ['Protein', 'Heme', 'MetalIon', 'Ligand']:
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
                    distance_results.append({'RealPos':res_id[1], f'min_distance_to_{ligand_resname}':min(min_distances), f'avg_distance_{ligand_resname}':min(avg_distances)})

    distance_results = pd.DataFrame(distance_results)
    return distance_results

if __name__ == '__main__':
    # input settings
    data_folder = address_dict['PIPS2']
    data_fbase = 'ET096'
    pdb_fname = 'Docked_ETS83096_B99990314_S82.pdb' # 'ETS83096.pdb' #
    pdb_fpath = f"{data_folder}{subfolders['pdb']}{data_fbase}/{pdb_fname}"

    # parse structure into protein vs. non protein
    parse_protein_hetatom(pdb_fpath, data_fbase)

    # get DSSP properties
    struct_df = get_dssp_properties(pdb_fpath.replace('.pdb', '_Protein.pdb'), data_fbase)

    # get residue to ligand distances
    for ligand_resname in ['HEM', 'UNK']:
        dist_df = calc_res_to_lig_distances(pdb_fpath, chain_id='A', ligand_resname=ligand_resname)
        struct_df = struct_df.merge(dist_df, on='RealPos', how='left')

    # get residue to surface distances, and surface patches
    surfdist = yasara_get_surface_residues(pdb_fpath, surft=2.55, distt=12, minepisize=13)
    struct_df = struct_df.merge(surfdist, on='RealPos', how='left')

    struct_df = struct_df.round(3)
    print(struct_df)
    struct_df.to_csv(f"{data_folder}{subfolders['pdb']}{data_fbase}/{data_fbase}_StructProperties.csv")

