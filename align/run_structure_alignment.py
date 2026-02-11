import os
import pandas as pd
import numpy as np
import yasara
from variables import address_dict, subfolders
from utils import run_msa, fetch_sequences_from_fasta, get_mutations_on_sk_wrt_s0, run_tmalign

class AlignStruct:
    def __init__(self, pdb_dir, sce_dir, msa_dir, delete_residue_str=None, superpose_method='resnum'):
        self.pdb_dir = pdb_dir
        self.sce_dir = sce_dir
        self.msa_dir = msa_dir
        self.delete_residue_str = delete_residue_str
        self.superpose_method = superpose_method

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
                # other_struct_fpaths = [struct_fpaths[k] for k in range(len(struct_fpaths)) if i != k]
                # other_structs_str = '-'.join([os.path.basename(f).replace('.pdb', '') for f in other_struct_fpaths])
                # output_pdb_fpath = struct_fpath_unaligned.replace('.pdb', f'_aligned_{other_structs_str}.pdb')
                output_pdb_fpath = struct_fpath_unaligned
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

    def yasara_align_structures(
            self,
            struct_fpaths,
            seq_align_fpath,
            output_sce_fpath,
            join_obj_in_pdb=False,
            save_indiv_aligned_structs=False,
            delete_not_protein=False
    ):

        # Initialize YASARA
        yasara.info.mode = 'txt'
        yasara.Console('off')
        yasara.FormatRes('RESName')

        # Load PDB files for all molecules
        for i, struct_fpath in enumerate(struct_fpaths):
            yasara.LoadPDB(struct_fpath)
            yasara.ColorObj(i + 1, 'magenta')

            # remove unnecessary components
            if delete_not_protein:
                yasara.DelWater()
                yasara.DelRes('not Protein')

            # delete selected residues and save temporary structures
            if self.delete_residue_str is not None:
                if isinstance(self.delete_residue_str,str):
                    del_target = self.delete_residue_str
                elif isinstance(self.delete_residue_str, list):
                    del_target =self.delete_residue_str[i]
                yasara.DelRes(f'Obj {i+1} and {del_target}')
                print(f'Deleted selected residues for Obj {i+1}:', del_target)

            # save PDB
            yasara.SavePDB(f'Obj {i+1}', struct_fpath.replace('.pdb', '_TEMP.pdb'))

        # Align all other objects to object 1 using MUSTANGPP
        obj_ref = os.path.basename(struct_fpaths[0])
        num_res_ref = yasara.CountRes(f'Obj 1 and Protein')
        resnum_list_ref = yasara.ListRes('Obj 1 and Protein', format='RESNUM')
        coords_ref = np.reshape(np.array(yasara.PosRes('Obj 1 and Protein')), shape=(-1, 3))
        seq_ref = yasara.SequenceRes('Obj 1 and Protein')[0]
        print(f'Ref Object: {obj_ref} ({num_res_ref} residues)', '\n')

        # iterate through objects 2 and up
        for obj_num in range(2, len(struct_fpaths) + 1):
            # get struct_fpath
            struct_fpath = struct_fpaths[obj_num-1]

            # # run TM-align (Zhang Yang group 2022 version)
            if delete_not_protein:
                tm_res = round(run_tmalign(struct_fpaths[0].replace('.pdb','_TEMP.pdb'), struct_fpaths[obj_num-1].replace('.pdb','_TEMP.pdb')),3)
            else:
                tm_res = None

            # align structure to reference
            res = yasara.AlignObj(f'{obj_num} and Protein', '1 and Protein', method='MUSTANGPP', results=4)
            rmsd_ali_residues, percent_identity_ali_residues, num_ali_residues = res[0], res[1], res[2]
            pairwise_ali_fpath = struct_fpaths[obj_num-1].replace('.pdb', '_ALI.fasta')
            yasara.SaveAli(obj_num, 1, method='Structure', filename=pairwise_ali_fpath, format='FASTA')
            ali_seqs, _, _ = fetch_sequences_from_fasta(pairwise_ali_fpath)
            os.remove(pairwise_ali_fpath)

            if obj_num==3:
                ali_seqs_TEMPLATE = ali_seqs

            # get superposed RMSD
            if obj_num<3:
                rmsd_sup = 0
            else:
                resnum_list = yasara.ListRes(f'Obj {obj_num} and Protein', format='RESNUM')
                residx_ref = 0
                residx = 0
                pos_ref_ali = []
                pos_target_ali = []
                if self.superpose_method=='resnum':
                    for resnum_ref, resnum in zip(resnum_list_ref, resnum_list):
                        pos_ref = yasara.PosRes(f'Obj 1 and Protein and Res {resnum_ref}')
                        pos_target = yasara.PosRes(f'Obj {obj_num} and Protein and Res {resnum}')
                        pos_ref_ali.append(pos_ref)
                        pos_target_ali.append(pos_target)
                elif self.superpose_method=='struct':
                    for res_ref, res_target in zip(ali_seqs_TEMPLATE[0], ali_seqs_TEMPLATE[1]):
                        pos_ref, pos_target = [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]
                        # advance residue counters
                        if res_ref !='-':
                            residx_ref += 1
                            if residx_ref - 1 < len(resnum_list_ref):
                                resnum_ref = resnum_list_ref[residx_ref - 1]
                                pos_ref = yasara.PosRes(f'Obj 1 and Protein and Res {resnum_ref}')
                        if res_target !='-':
                            residx += 1
                            if residx-1<len(resnum_list):
                                resnum = resnum_list[residx - 1]
                                pos_target = yasara.PosRes(f'Obj {obj_num} and Protein and Res {resnum}')
                                if len(pos_target) != 3: print(obj_num, resnum, pos_target, residx, resnum_list[residx - 1], resnum_list)
                        pos_ref_ali.append(pos_ref)
                        pos_target_ali.append(pos_target)
                # get rmsd
                pos_ref_ali = np.array(pos_ref_ali)
                pos_target_ali = np.array(pos_target_ali)
                try:
                    rmsd_sup_byres = np.sqrt(np.sum((pos_ref_ali - pos_target_ali) ** 2, axis=1))
                    rmsd_sup = round(np.nanmean(rmsd_sup_byres), 2)
                except:
                    rmsd_sup = None

            # parse alignment results
            calist = res[3:]
            obj = os.path.basename(struct_fpaths[obj_num-1])
            num_res = yasara.CountRes(f'Obj {obj_num} and Protein')
            percent_ali_res = round(num_ali_residues/num_res*100,1)
            print(f'{obj} ({num_res} residues)')
            print('TM-align score:', tm_res)
            print('RMSD (superposed residues):', rmsd_sup)
            print('RMSD (aligned residues):', round(rmsd_ali_residues,2))
            print(f'# of aligned residues: {num_ali_residues}/{num_res} ({percent_ali_res}%)')

            # Process aligned residues
            num_iden, num_iden_hypho, num_sim, num_sim_hypho = self.check_and_color(num_ali_residues, calist)

            # print stats
            print(f'[IDENTICAL] all: {num_iden} ({100 * num_iden / num_ali_residues:.1f}%); hydrophobic: {num_iden_hypho} ({100 * num_iden_hypho / num_ali_residues:.1f}%); non-hydrophobic: {num_iden - num_iden_hypho} ({100 * (num_iden - num_iden_hypho) / num_ali_residues:.1f}%)')
            print(f'[SIMILAR] all: {num_sim} ({100 * num_sim / num_ali_residues:.1f}%); hydrophobic: {num_sim_hypho} ({100 * num_sim_hypho / num_ali_residues:.1f}%); non-hydrophobic: {num_sim - num_sim_hypho} ({100 * (num_sim - num_sim_hypho) / num_ali_residues:.1f}%)')
            print()

            # delete TEMP file
            os.remove(struct_fpath.replace('.pdb', '_TEMP.pdb'))

        # Save the alignment as FASTA
        if seq_align_fpath is not None:
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
            seq_align_fname,
            run_structure_alignment=True,
            parse_seq_struct_alignments=False,
            save_sce=False,
            join_obj_in_pdb=False,
            save_indiv_aligned_structs=False,
            delete_not_protein=False
    ):
        # get inputs
        struct_fpaths = [f'{self.pdb_dir}{f}.pdb' for f in struct_names]
        seq_align_fpath = self.msa_dir + seq_align_fname if seq_align_fname is not None else None
        output_sce_fpath = None
        if save_sce:
            output_sce_fpath = self.sce_dir + seq_align_fname.replace('.fasta', '.sce')

        # align structures
        if run_structure_alignment:
            self.yasara_align_structures(struct_fpaths, seq_align_fpath, output_sce_fpath, join_obj_in_pdb, save_indiv_aligned_structs, delete_not_protein)

        # parse aligned struct info
        if parse_seq_struct_alignments:
            csv_fpath = seq_align_fpath.replace('msa/', 'pdb/').replace('.fasta', '.csv')
            df = self.parse_aligned_struct_info(seq_align_fpath, struct_fpaths, csv_fpath)
            return df

    def use_seed_alignment_to_get_msa(self, seq_fname, output_msa_fpath, seq_align_fpath, seq_dir, msa_dir):
        run_msa(seq_fname, output_msa_fpath, 'mafft', seq_dir, msa_dir, fmt='fasta', seed_ali=seq_align_fpath)


if __name__ == '__main__':
    os.chdir('../')
    print('CWD:', os.getcwd())
    # directories
    data_folder = address_dict['PIPS2'] # address_dict['ECOHARVEST']
    data_fbase = 'UPOs_peroxygenation_analysis/docked/ALL_aligned' # 'CARs/PoET2_NiCAR-MpCAR_Hybrids'
    pdb_dir = data_folder + subfolders['pdb'] + data_fbase + '/'
    sce_dir = pdb_dir.replace('pdb','sce')
    seq_dir = data_folder + subfolders['sequences'] + data_fbase + '/'
    msa_dir = None # data_folder + subfolders['msa'] + data_fbase + '/'
    # seq_align_fname = None
    # seq_align_fname = f'{"_".join(struct_names)}_StructAli.fasta'
    # seq_align_fname = f'PoET2_NiCARseq_MpCARstruct_Res383hybrid_06-10_StructAli.fasta'
    seq_align_fname = 'UPOs_peroxygenation.fasta'

    # parameters
    run_structure_alignment = True
    save_sce = False
    join_obj_in_pdb = False
    save_indiv_aligned_structs = True
    parse_seq_struct_alignments = False
    delete_not_protein = False
    struct_base = 'OA167'
    struct_name_ref = f'{struct_base}_S82_swissdock_0'
    struct_names = [struct_name_ref] + \
                   [f.replace('.pdb','') for f in os.listdir(pdb_dir) if (f.find('.pdb')>-1 and f.replace('.pdb','')!=struct_name_ref and f.find(struct_base)>-1)]
    # struct_names = [
        # # 'Boltz2_NiCAR-A_WT_Oleoyl-AMP',
        # 'Boltz2_MpCAR-A_WT_Oleoyl-AMP',
        # 'Boltz2_MpCAR-A_WT_Oleoyl-AMP',
        # 'Boltz2_NiCAR-A_WT_Oleoyl-AMP',
        #
        # 'PoET2_NiCAR-MpCAR_Hybrid01',
        # 'PoET2_NiCAR-MpCAR_Hybrid02',
        # 'PoET2_NiCAR-MpCAR_Hybrid03',
        # 'PoET2_NiCAR-MpCAR_Hybrid04',
        # 'PoET2_NiCAR-MpCAR_Hybrid05',

        # 'PoET2_NiCAR-MpCAR_Hybrid06',
        # 'PoET2_NiCAR-MpCAR_Hybrid07',
        # 'PoET2_NiCAR-MpCAR_Hybrid08',
        # 'PoET2_NiCAR-MpCAR_Hybrid09',
        # 'PoET2_NiCAR-MpCAR_Hybrid10',
    # ]
    print(struct_names)
    delete_residue_str = None
    # delete_residue_str = ['Protein and 1-369']*2 + ['Protein and 1-383']*6 # vs MpCAR ref, hybrids #1-5
    # delete_residue_str = ['Protein and 1-383']*2 + ['Protein and 1-369']*1 + ['Protein and 1-383']*5  # vs NiCAR ref, hybrids #1-5
    # delete_residue_str = ['Protein and not 1-369']*2 + ['Protein and not 1-383']*6 # vs MpCAR ref, hybrids #6-10
    # delete_residue_str = ['Protein and not 1-383']*2 + ['Protein and not 1-369']*1 + ['Protein and not 1-383']*5  # vs NiCAR ref, hybrids #6-10
    # delete_residue_str = ['Protein and 1-383']*2 + ['Protein and 1-369']*1 + ['Protein and 1-383']*5
    superpose_method = 'struct' # 'resnum' #

    # perform alignment and parsing of aligned PDB residue info
    align_struct = AlignStruct(
        pdb_dir,
        sce_dir,
        msa_dir,
        delete_residue_str,
        superpose_method
    )

    # run alignment pipeline
    df = align_struct.run_pipeline(
        struct_names,
        seq_align_fname,
        run_structure_alignment,
        parse_seq_struct_alignments,
        save_sce,
        join_obj_in_pdb,
        save_indiv_aligned_structs,
        delete_not_protein
    )

    # convert mutations from one basis (indexed wrt a 1st sequence) to other bases
    mutations_ref_s0 = None # ['F88A', 'T157A'] # ['T125A', 'A129G', 'A247H', 'T244A', 'F243G']
    reorder_seqs = None # [2,0,1] #
    if mutations_ref_s0 is not None:
        seq_align_fpath = msa_dir + seq_align_fname
        mutations_conversion = get_mutations_on_sk_wrt_s0(seq_align_fpath, mutations_ref_s0, reorder_seqs)
