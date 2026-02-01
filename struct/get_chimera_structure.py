import os
import pandas as pd
import numpy as np
pd.set_option('display.max_columns',None)
from variables import address_dict, subfolders


class ChimerizePDB:
    def __init__(
            self,
            data_folder,
            data_fbase='',
            pdb_subfolder='pdb/',
            msa_subfolder='msa/',
            sequences_subfolder='sequences/',
    ):
        self.data_folder = data_folder
        self.pdb_dir = data_folder + pdb_subfolder
        self.seq_dir = data_folder + sequences_subfolder
        self.msa_dir = data_folder + msa_subfolder
        self.data_fbase = data_fbase

    def get_aligned_structural_info(self, struct_ali_csv, struct_names=None):
        # check if alignment already exists
        if struct_ali_csv is not None and os.path.exists(struct_ali_csv):
            df = pd.read_csv(struct_ali_csv, index_col=0)
        # run structural alignment and extract aligned residue positions as CSV
        else:
            print('Performing structural alignment of', *struct_names, '...')
            from align.run_structure_alignment import AlignStruct
            struct_fnames = [f+'.pdb' for f in struct_names]
            seq_align_fpath = f'{self.msa_dir}{self.data_fbase}/{struct_names[0]}_{struct_names[1]}_yasaraStructAli.fasta'
            # perform alignment and parsing of aligned PDB residue info
            align_struct = AlignStruct(self.data_folder, subfolders['pdb'], self.data_fbase)
            df = align_struct.run_pipeline(struct_fnames, seq_align_fpath)
        return df

    def hybridize_pdbs(self, df_hybridized, pos_to_hybridize, pdb_hybridized_fpath):
        """hybridize the two structures at the selected positions"""
        # update sequence AA
        df_hybridized.loc[df_hybridized['pos_0'].isin(pos_to_hybridize), 'seq_0'] = df_hybridized.loc[
            df_hybridized['pos_0'].isin(pos_to_hybridize), 'seq_1']

        # get positions to hybridize wrt 2nd sequence
        pos_to_hybridize_wrt_seq2 = df_hybridized.loc[df_hybridized['pos_0'].isin(pos_to_hybridize), 'pos_1'].tolist()

        # update PDB segments residue by residue
        for pos_0 in pos_to_hybridize:
            mask = df_hybridized['pos_0'] == pos_0
            pos_0_str = ' ' + str(pos_0) + ' '
            pos_1_str = ' ' + str(df_hybridized.loc[mask, 'pos_1'].iloc[0]) + ' '
            print(pos_0_str, pos_1_str)
            pdb_str_residue = df_hybridized.loc[mask, 'pdb_by_residue_1'].iloc[0]
            pdb_residue_str_mod = pdb_str_residue.replace(pos_1_str, pos_0_str)
            df_hybridized.loc[mask, 'pdb_by_residue_0'] = pdb_residue_str_mod
        df_hybridized_list = df_hybridized.loc[df_hybridized['pos_0'] != -1, ['pdb_by_residue_0']][
            'pdb_by_residue_0'].tolist()

        # update atom numbering
        atom_count = 0
        atom_str_list = []
        atom_attribute_array = []
        spaces = [(6, 'L', 0), (7, 'R', 2), (4, 'L', 0), (4, 'L', 0), (2, 'L', 0), (7, 'R', 4), (8, 'R', 0), (8, 'R', 0), (8, 'R', 0), (6, 'R', 0), (17, 'R', 11), (1, 'L', 0)]  # (length, alignment, offset)
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
                    atom_str_property = [' '] * num_spaces
                    entry_aslist = list(atom_property)
                    entry_len = len(entry_aslist)
                    if alignment == 'L':
                        start_idx = offset
                    elif alignment == 'R':
                        start_idx = num_spaces - offset - entry_len
                    atom_str_property[start_idx:start_idx + entry_len] = entry_aslist
                    atom_str_updated += ''.join(atom_str_property)
                atom_str_updated = atom_str_updated.replace('           ', '     ')

                # update atom_str_list
                atom_str_list.append(atom_str_updated)

                # append to atom_attribute_array
                if len(atom_str_parsed) >= 9:
                    atom_str_toappend = atom_str_parsed[:9]
                else:
                    atom_str_toappend = atom_str_parsed + [''] * (9 - len(atom_str_parsed))
                atom_attribute_array.append(atom_str_toappend)

        atom_attribute_array = np.array(atom_attribute_array)
        atoms_df = pd.DataFrame(atom_attribute_array,
                                columns=['ATOM', 'atom_num', 'element_w_idx', 'amino_acid', 'chain', 'residue_num', 'x',
                                         'y', 'z'])

        # save the hybridized PDB
        pdb_hybridized = '\n'.join(atom_str_list)
        with open(pdb_hybridized_fpath, 'w') as f:
            f.write(pdb_hybridized)
        print('Written hybrid PDB file to:', pdb_hybridized_fpath)
        return pos_to_hybridize_wrt_seq2, df_hybridized, atoms_df

    def run_structural_energy_minimization(self, pdb_hybridized_fpath):
        import yasara
        # run YASARA energy minimization
        print('Running energy minimization of hybridized structure...')
        yasara.Console('Off')
        # Load PDB file into YASARA
        yasara.Clear()
        yasara.LoadPDB(pdb_hybridized_fpath)  # or yasara.Load('your_structure.pdb')
        # Apply the macro
        # yasara.ApplyMacro('/Applications/YASARA.app/Contents/yasara/mcr/em_run.mcr', targets='')
        yasara.PlayMacro('/Applications/YASARA.app/Contents/yasara/mcr/em_run.mcr', label='')
        # Save the minimized structure
        pdb_minimized_fpath = pdb_hybridized_fpath.replace('.pdb', '_minimized.pdb')
        yasara.SavePDB(1, pdb_minimized_fpath)
        print('Saved energy minimized PDB file to:', pdb_minimized_fpath)
        return pdb_minimized_fpath

    def parse_atom_coords_from_structure(self, pdb_fpath, selected_residues=None):
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_fpath)
        if selected_residues is None:
            selected_residues = [residue.get_id()[1] for residue in structure[0][0]]
        # get selected residue coordinates
        coords_byresidue = []
        all_atom_coords = []
        # parse structure
        for model in structure:
            for chain in model:
                for residue in chain:
                    res_id = residue.get_id()
                    # res_id is a tuple like (' ', 150, ' ')
                    resnum = res_id[1]
                    if resnum in selected_residues:
                        coords_atoms_in_residue = []
                        for atom in residue:
                            if atom.element != 'H':  # Skip hydrogen atoms
                                coords_atoms_in_residue.append(atom.coord)
                                all_atom_coords.append(atom.coord)
                        coords_byresidue.append(np.array(coords_atoms_in_residue))

        alpha_carbon_coords = np.array([list(coords[1]) for coords in coords_byresidue])
        residue_centroid_coords = np.array([list(np.mean(coords, axis=0)) for coords in coords_byresidue])
        all_atom_coords = np.array(all_atom_coords)
        return alpha_carbon_coords, residue_centroid_coords, all_atom_coords, coords_byresidue

    def run_pipeline(self,
                     struct_names,
                     pos_to_hybridize,
                     align_structures=False,
                     hybridize_structures=True,
                     ):
        if hybridize_structures:
            # get aligned structural info
            struct_ali_csv = f'{self.pdb_dir}{self.data_fbase}/{struct_names[0]}_{struct_names[1]}_yasaraStructAli.csv'
            pdb_hybridized_fpath = struct_ali_csv.replace('_yasaraStructAli.csv', '_hybrid.pdb')
            if align_structures or not os.path.exists(struct_ali_csv):
                struct_ali_csv = None
            df_hybridized = self.get_aligned_structural_info(struct_ali_csv, struct_names)

            # hybridize PDBs
            pos_to_hybridize_wrt_seq2, df_hybridized, atoms_df = self.hybridize_pdbs(df_hybridized, pos_to_hybridize, pdb_hybridized_fpath)

            # minimize energy of hybridized structure
            pdb_minimized_fpath = self.run_structural_energy_minimization(pdb_hybridized_fpath)

        # calculate average distance of selected residues to centroid for each structure:
        pdb_fnames = [
            f'{struct_names[0]}_aligned_{struct_names[1]}.pdb', # backbone structure
            f'{struct_names[0]}_{struct_names[1]}_hybrid.pdb', # raw hybridized structure
            f'{struct_names[0]}_{struct_names[1]}_hybrid_minimized.pdb', # minimized hybridized structure
        ]
        f_list = []
        pdb_fpaths = [f'{self.pdb_dir}{self.data_fbase}/' + pdb_fname for pdb_fname in pdb_fnames]
        coords_dict = {}
        for i, pdb_fpath in enumerate(pdb_fpaths):
            f = os.path.basename(pdb_fpath).split('.')[0]
            f_list.append(f)
            print(f)
            alpha_carbon_coords, residue_centroid_coords, all_atom_coords, coords_byresidue = self.parse_atom_coords_from_structure(pdb_fpath, pos_to_hybridize)
            avg_distance_from_centroid = round(np.linalg.norm(all_atom_coords - all_atom_coords.mean(axis=0), axis=1).mean(),4)
            print('avg_distance_from_centroid', avg_distance_from_centroid)
            # update coords_dict
            coords_dict[i] = {
                'alpha_carbon_coords':alpha_carbon_coords,
                'residue_centroid_coords':residue_centroid_coords,
                'all_atom_coords':all_atom_coords,
                'coords_byresidue':coords_byresidue
            }
            if i==2:
                # backbone discrepancy
                backbone_diff = {}
                residue_centroid_diff = {}
                for j in [0,1]:
                    backbone_diff[j] = round(np.linalg.norm(coords_dict[2]['alpha_carbon_coords'] - coords_dict[j]['alpha_carbon_coords']).mean(),4)
                    print(f'backbone diff vs [{j}] {f_list[j]}: {backbone_diff[j]}')
                    residue_centroid_diff[j] = round(np.linalg.norm(coords_dict[2]['residue_centroid_coords'] - coords_dict[j]['residue_centroid_coords']).mean(), 4)
                    print(f'residue centroid diff vs [{j}] {f_list[j]}: {residue_centroid_diff[j]}')

if __name__ == '__main__':
    data_folder = address_dict['PIPS2']
    pdb_subfolder = subfolders['pdb']
    data_fbase = 'UPOpanel_selected'
    struct_names = ['ET096', 'CviUPO']
    align_structures = True
    hybridize_structures = False
    pos_to_hybridize = [38, 39, 60, 64, 70, 71, 73, 74, 75, 76, 77, 78, 79, 80, 81, 103, 171, 174, 175, 176, 177, 178, 179, 182, 218, 219, 220, 223]

    chimerize_pdbs = ChimerizePDB(data_folder, data_fbase)
    chimerize_pdbs.run_pipeline(struct_names, pos_to_hybridize, align_structures, hybridize_structures)


