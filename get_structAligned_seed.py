import yasara
from variables import address_dict, subfolders
from utils import run_msa

def save_aligned_structures(output_sce_fpath, save_pdb=False):
    # save SCE
    yasara.SaveSce(output_sce_fpath)
    print('Saved SCE output:', output_sce_fpath)
    num_obj = yasara.CountObj('All')
    # save PDB
    if save_pdb:
        for obj_num in range(2,num_obj+1):
            yasara.JoinObj(obj_num, 1, center='No')
        output_pdb_fpath = output_sce_fpath.replace('.sce', '.pdb').replace('sce/', 'pdb/')
        yasara.SavePDB(1, output_pdb_fpath)
        print('Saved PDB output:', output_pdb_fpath)

def check_and_color(residues, calist):
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
        print('\n', j, k, l)
        # iterate through aa groups
        for aagroup_type, (aalist_group, aagroup_color) in groups.items():
            if k in aalist_group and l in aalist_group:
                print(aagroup_type, 'similar', end=' ')
                num_sim += 1
                if aagroup_type == 'hypho':
                    num_sim_hypho += 1
                yasara.ColorRes(f'Atom {atom1} or Atom {atom2}', aagroup_color)
                if k==l:
                    print('+ identical', end=' ')
                    num_iden += 1
                    if aagroup_type == 'hypho':
                        num_iden_hypho += 1

    return num_iden, num_iden_hypho, num_sim, num_sim_hypho

def yasara_align_structures(struct_fpaths, alignment_fpath, output_sce_fpath):

    # Initialize YASARA
    yasara.Console('off')
    yasara.FormatRes('RESName')

    # Load PDB files for all molecules
    for i, struct_fpath in enumerate(struct_fpaths):
        yasara.LoadPDB(struct_fpath)
        yasara.ColorObj(i+1, 'magenta')
    # remove unnecessary components
    yasara.DelWater()
    yasara.DelMol('not Protein')
    for obj_num in range(1, len(struct_fpaths)+1):
        num_mol = yasara.CountMol('Obj '+str(obj_num))
        if num_mol > 1:
            yasara.DelMol('Obj '+str(obj_num)+' and not Mol A')

    # Align all other objects to object 1 using MUSTANGPP
    obj_ref = yasara.NameObj(1)
    for obj_num in range(2,len(struct_fpaths)+1):
        num_iden = num_iden_hypho = num_sim = num_sim_hypho = 0
        res = yasara.AlignObj(obj_num, 1, method='MUSTANGPP', results=4)
        rmsd, percent_identity, residues = res[0], res[1], res[2]
        calist = res[3:]

        # Process aligned residues
        num_iden, num_iden_hypho, num_sim, num_sim_hypho = check_and_color(residues, calist)

        # print stats
        obj = yasara.NameObj(obj_num)
        print(obj_ref, obj, '# of aligned residues:', residues, 'RMSD:', f'{rmsd:.2f}', '% identity:', round(percent_identity,3))
        print(f'[IDENTICAL] all: {num_iden} ({100*num_iden/residues:.2f}%); hydrophobic: {num_iden_hypho} ({100*num_iden_hypho/residues:.2f}%); non-hydrophobic: {num_iden-num_iden_hypho} ({100*(num_iden-num_iden_hypho)/residues:.2f}%)')
        print(f'[SIMILAR] all: {num_sim} ({100*num_sim/residues:.2f}%); hydrophobic: {num_sim_hypho} ({100*num_sim_hypho/residues:.2f}%); non-hydrophobic: {num_sim-num_sim_hypho} ({100*(num_sim-num_sim_hypho)/residues:.2f}%)')

    # Save the alignment as FASTA
    yasara.SaveAli('!1', '1', filename=alignment_fpath, format='FASTA')

    # save aligned structures
    save_aligned_structures(output_sce_fpath)


if __name__ == '__main__':
    data_folder = address_dict['ECOHARVEST']
    seq_dir = data_folder + subfolders['sequences']
    msa_dir = data_folder + subfolders['msa']
    struct_fnames = [
        'CALB_1tca.pdb',
        'CALA_2veo.pdb'
        # 'CALBonly.pdb',
        # 'CALAonly.pdb',
    ]
    struct_ali_fpath = msa_dir + 'CALB-CALA_yasaraStructAli.fasta'
    seq_fname = None
    output_msa_fpath = None

    # perform structural alignment
    output_sce_fpath = struct_ali_fpath.replace('msa/','sce/').replace('.fasta', '.sce')
    struct_fpaths = [data_folder+subfolders['pdb']+f for f in struct_fnames]
    yasara_align_structures(struct_fpaths, struct_ali_fpath, output_sce_fpath)

    # use alignment to get MSA alignment
    if seq_fname is not None and output_msa_fpath is not None:
        run_msa(seq_fname, output_msa_fpath, 'mafft', seq_dir, msa_dir, fmt='fasta', seed_ali=struct_ali_fpath)