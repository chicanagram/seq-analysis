import os
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select, is_aa
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta, run_msa
from variables import address_dict, subfolders

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
def calc_res_to_lig_distances(pdb_fpath, chain_id='A', ligand_resname='H_UNK'):

    # Parse PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_fpath)

    # Select ligand atoms
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_resname:
                    ligand_atoms.extend([atom for atom in residue if atom.element != 'H'])
    print(ligand_atoms)

    # Compute distances between ligand atoms and heavy atoms of each residue
    distance_results = {}
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.get_resname().strip() != ligand_resname:
                        res_id = (residue.get_resname(), residue.get_id()[1])
                        distances = []
                        for atom in residue:
                            if atom.element != 'H':  # exclude hydrogen atoms
                                for lig_atom in ligand_atoms:
                                    distance = np.linalg.norm(atom.coord - lig_atom.coord)
                                    distances.append(distance)
                        distance_results[res_id] = distances

    # Example output
    for res_id, distances in distance_results.items():
        print(f"Residue {res_id} - Min distance: {min(distances):.2f} Å, Max distance: {max(distances):.2f} Å")



if __name__ == '__main__':
    # input settings
    data_folder = address_dict['PIPS2']
    data_fbase = 'ET096'
    pdb_fname = 'Docked_ETS83096_B99990314_S82.pdb' # 'ETS83096.pdb' #
    pdb_fpath = f"{data_folder}{subfolders['pdb']}{data_fbase}/{pdb_fname}"

    # parse structure into protein vs. non protein
    parse_protein_hetatom(pdb_fpath, data_fbase)

    # get DSSP properties
    dssp_res = get_dssp_properties(pdb_fpath.replace('.pdb', '_Protein.pdb'), data_fbase)
    # print(dssp_res)

    # get residue to ligand distances
    calc_res_to_lig_distances(pdb_fpath, chain_id='A', ligand_resname='HEM')