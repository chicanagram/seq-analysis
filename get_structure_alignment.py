import os
import pandas as pd
import yasara
from variables import address_dict, subfolders
from utils import run_msa, fetch_sequences_from_fasta, get_mutations_on_sk_wrt_s0

class AlignStruct:
    def __init__(self, data_folder, pdb_subfolder, data_fbase):
        self.data_folder = data_folder
        self.pdb_subfolder = pdb_subfolder
        self.data_fbase = data_fbase

    def save_aligned_structures(self, struct_fpaths, output_sce_fpath=None, join_obj_in_pdb=False, save_indiv_aligned_structs=False):
        print('struct_fpaths:', len(struct_fpaths), struct_fpaths)
        # save SCE
        if output_sce_fpath is not None:
            yasara.SaveSce(output_sce_fpath)
            print('Saved SCE output:', output_sce_fpath)
        # save PDB of joined objects
        num_obj = yasara.CountObj('All')
        if join_obj_in_pdb:
            for obj_num in range(2, num_obj + 1):
                yasara.JoinObj(obj_num, 1, center='No')
        # save individual aligned PDBs
        if save_indiv_aligned_structs:
            for i in range(num_obj):
                obj_num = i + 1
                struct_fpath_unaligned = struct_fpaths[i]
                other_struct_fpaths = [struct_fpaths[k] for k in range(len(struct_fpaths)) if i != k]
                other_structs_str = '-'.join([os.path.basename(f).replace('.pdb', '') for f in other_struct_fpaths])
                output_pdb_fpath = struct_fpath_unaligned.replace('.pdb', f'_aligned_{other_structs_str}.pdb')
                yasara.SavePDB(obj_num, output_pdb_fpath)
                print('Saved PDB output:', output_pdb_fpath)

    def check_and_color(self, residues, calist):
        groups = {
            'hypho': (['ALA', 'ILE', 'LEU', 'MET', 'VAL', 'PHE', 'TRP', 'TYR'], 'yellow'),
            'polar': (['ASN', 'GLN', 'SER', 'THR', 'GLY', 'PRO'], 'green'),
            'posit': (['LYS', 'ARG'], 'blue'),
            'negat': (['GLU', 'ASP'], 'red'),
            'speci': (['CYS', 'HIS'], 'cyan'),
        }
        num_iden = 0
        num_sim = 0
        num_iden_hypho = 0
        num_sim_hypho = 0

        # Process aligned residues
        for j in range(residues):
            atom1 = calist[j * 2]
            atom2 = calist[j * 2 + 1]
            k = yasara.NameRes(f'Atom {atom1}')[0]
            l = yasara.NameRes(f'Atom {atom2}')[0]
            # print('\n', j, k, l)
            # iterate through aa groups
            for aagroup_type, (aalist_group, aagroup_color) in groups.items():
                if k in aalist_group and l in aalist_group:
                    # print(aagroup_type, 'similar', end=' ')
                    num_sim += 1
                    if aagroup_type == 'hypho':
                        num_sim_hypho += 1
                    yasara.ColorRes(f'Atom {atom1} or Atom {atom2}', aagroup_color)
                    if k == l:
                        # print('+ identical', end=' ')
                        num_iden += 1
                        if aagroup_type == 'hypho':
                            num_iden_hypho += 1

        return num_iden, num_iden_hypho, num_sim, num_sim_hypho

    def yasara_align_structures(self, struct_fpaths, seq_align_fpath, output_sce_fpath, join_obj_in_pdb=False, save_indiv_aligned_structs=False):

        # Initialize YASARA
        yasara.info.mode = 'txt'
        yasara.Console('off')
        yasara.FormatRes('RESName')

        # Load PDB files for all molecules
        for i, struct_fpath in enumerate(struct_fpaths):
            yasara.LoadPDB(struct_fpath)
            yasara.ColorObj(i + 1, 'magenta')
        # remove unnecessary components
        yasara.DelWater()

        # Align all other objects to object 1 using MUSTANGPP
        obj_ref = yasara.NameObj(1)
        for obj_num in range(2, len(struct_fpaths) + 1):
            res = yasara.AlignObj(f'{obj_num} and Protein', '1 and Protein', method='MUSTANGPP', results=4)
            rmsd, percent_identity, residues = res[0], res[1], res[2]
            calist = res[3:]
            print('RMSD:', rmsd)
            print('% identity:', percent_identity)
            print('residues:', residues)
            print('calist:', calist)

            # Process aligned residues
            num_iden, num_iden_hypho, num_sim, num_sim_hypho = self.check_and_color(residues, calist)

            # print stats
            obj = yasara.NameObj(obj_num)
            print(obj_ref, obj, '# of aligned residues:', residues, 'RMSD:', f'{rmsd:.2f}', '% identity:', round(percent_identity, 3))
            print(f'[IDENTICAL] all: {num_iden} ({100 * num_iden / residues:.2f}%); hydrophobic: {num_iden_hypho} ({100 * num_iden_hypho / residues:.2f}%); non-hydrophobic: {num_iden - num_iden_hypho} ({100 * (num_iden - num_iden_hypho) / residues:.2f}%)')
            print(f'[SIMILAR] all: {num_sim} ({100 * num_sim / residues:.2f}%); hydrophobic: {num_sim_hypho} ({100 * num_sim_hypho / residues:.2f}%); non-hydrophobic: {num_sim - num_sim_hypho} ({100 * (num_sim - num_sim_hypho) / residues:.2f}%)')

        # Save the alignment as FASTA
        yasara.SaveAli('!1', '1', filename=seq_align_fpath, format='FASTA')

        # save aligned structures
        self.save_aligned_structures(struct_fpaths, output_sce_fpath, join_obj_in_pdb, save_indiv_aligned_structs)

    def parse_aligned_struct_info(self, seq_align_fpath, struct_fpaths, csv_fpath):
        # get aligned sequences
        seqs_aligned, seq_names, _ = fetch_sequences_from_fasta(seq_align_fpath)
        # get aligned PDB structures
        pdbs_as_string = []
        pdbs_as_string_byresidue = []
        for i, struct_fpath_unaligned in enumerate(struct_fpaths):
            other_struct_fpaths = [struct_fpaths[k] for k in range(len(struct_fpaths)) if i != k]
            other_structs_str = '-'.join([os.path.basename(f).replace('.pdb', '') for f in other_struct_fpaths])
            aligned_pdb_fpath = struct_fpath_unaligned.replace('.pdb', f'_aligned_{other_structs_str}.pdb')
            # load PDB data
            with open(aligned_pdb_fpath, 'r') as f:
                pdb_string = f.read()
                pdbs_as_string.append(pdb_string)
                pdb_string_list = [l for l in pdb_string.split('\n') if
                                   (l[:4] == 'ATOM' or (l[:6] == 'HETATM' and l.find('HIP') > -1) or l[:3] == 'TER')]
                pdb_as_string_byresidue = {}
                # keep only lines with residue information
                for l in pdb_string_list:
                    l_list = ' '.join(l.split()).split(' ')
                    if l_list[0] == 'ATOM' or (l_list[0] == 'HETATM' and l_list[3] in ['HIP']):
                        res_num = int(l_list[5])
                    elif l_list[0] == 'TER':
                        res_num = int(l_list[4])
                    # get residue number
                    if res_num not in pdb_as_string_byresidue:
                        pdb_as_string_byresidue[res_num] = []
                    # combine all lines belonging to the same residue
                    pdb_as_string_byresidue[res_num].append(l)

            # join lines as a string for each residue
            for res_num, pdb_string_list_res in pdb_as_string_byresidue.items():
                pdb_as_string_byresidue[res_num] = '\n'.join(pdb_string_list_res)
                print(res_num, pdb_as_string_byresidue[res_num], '\n')

            # update overall list
            pdbs_as_string_byresidue.append(pdb_as_string_byresidue)

        df = pd.DataFrame()
        for i, (seq_ali, seq_name, pdb_as_string_byresidue) in enumerate(
                zip(seqs_aligned, seq_names, pdbs_as_string_byresidue)):
            seq_nogaps = seq_ali.replace('-', '')
            pos_list = list(pdb_as_string_byresidue.keys())
            pos_list_ali = []
            pos_list.sort()
            print(len(seq_nogaps), len(pos_list), [pos for pos in range(1, pos_list[-1] + 1) if pos not in pos_list])
            pdb_str_byresidue = []
            pos_idx = 0
            for aa in seq_ali:
                if aa != '-':
                    print(pos_idx, end=' ')
                    pos = pos_list[pos_idx]
                    pos_list_ali.append(int(pos))
                    pdb_str_byresidue.append(pdb_as_string_byresidue[pos])
                    pos_idx += 1
                else:
                    pos_list_ali.append(-1)
                    pdb_str_byresidue.append('')
            df[f'seq_{i}'] = list(seq_ali)
            df[f'pos_{i}'] = pos_list_ali
            df[f'pdb_by_residue_{i}'] = pdb_str_byresidue
            print()
        print(df)
        df.to_csv(csv_fpath)
        return df

    def run_pipeline(
            self,
            struct_names,
            seq_align_fpath,
            run_structure_alignment=True,
            parse_seq_struct_alignments=True,
            save_sce=False
    ):
        struct_fpaths = [self.data_folder + self.pdb_subfolder + self.data_fbase + '/' + f for f in struct_names + '.pdb']
        csv_fpath = seq_align_fpath.replace('msa/', 'pdb/').replace('.fasta', '.csv')
        output_sce_fpath = None
        if save_sce:
            output_sce_fpath = seq_align_fpath.replace('msa/', 'sce/').replace('.fasta', '.sce')

        # align structures
        if run_structure_alignment:
            self.yasara_align_structures(struct_fpaths, seq_align_fpath, output_sce_fpath)

        # parse aligned struct info
        if parse_seq_struct_alignments:
            df = self.parse_aligned_struct_info(seq_align_fpath, struct_fpaths, csv_fpath)
            return df

if __name__ == '__main__':
    print(os.getcwd())
    data_folder = address_dict['PIPS2']
    data_fbase = ''
    seq_dir = data_folder + subfolders['sequences'] + data_fbase + '/'
    msa_dir = data_folder + subfolders['msa'] + data_fbase + '/'
    run_structure_alignment = False
    join_obj_in_pdb = False
    save_indiv_aligned_structs = False
    parse_seq_struct_alignments = True
    struct_names = [
        'MpCAR-A',
        'NiCAR-A',
        'MmCAR-A',
        # 'CviUPO',
        # 'ET096',
        # 'CmaUPO',
    ]
    seq_align_fpath = f'{msa_dir}{"_".join(struct_names)}_StructAli.fasta'
    struct_fpaths = [f'{data_folder}{subfolders["pdb"]}{data_fbase}/{f}.pdb' for f in struct_names]
    print(seq_align_fpath)
    print(struct_fpaths)

    # perform alignment and parsing of aligned PDB residue info
    align_struct = AlignStruct(data_folder, subfolders['pdb'], data_fbase)
    align_struct.yasara_align_structures(
        struct_fpaths,
        seq_align_fpath,
        seq_align_fpath.replace('/msa/','/sce/').replace('.fasta','.sce'),
        join_obj_in_pdb,
        save_indiv_aligned_structs,
    )
    # df = align_struct.run_pipeline(struct_fnames, seq_align_fpath)

    # convert mutations from one basis (indexed wrt a 1st sequence) to other bases
    mutations_ref_s0 = None # ['F88A', 'T157A'] # ['T125A', 'A129G', 'A247H', 'T244A', 'F243G']
    reorder_seqs = None # [2,0,1] #
    if mutations_ref_s0 is not None:
        mutations_conversion = get_mutations_on_sk_wrt_s0(seq_align_fpath, mutations_ref_s0, reorder_seqs)

    ##############################################
    # use alignment as seed align more sequences #
    ##############################################
    seq_fname = None
    output_msa_fpath = None
    if seq_fname is not None and output_msa_fpath is not None:
        seq_fname = 'MpCAR-A.fasta'  # None
        run_msa(seq_fname, output_msa_fpath, 'mafft', seq_dir, msa_dir, fmt='fasta', seed_ali=seq_align_fpath)