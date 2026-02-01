import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns',None)
import matplotlib.pyplot as plt
import yasara
from yasara import yasaradir
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta


if __name__=='__main__':
    os.chdir('../')

    # Path to the FASTA file
    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'CARs'
    struct_fname_dict = {
        'NiCAR-A': [
            'NiCAR-A/Boltz2_NiCAR-A_WT_ATP/NiCAR-A_WT_ATP_model_0.pdb',
            'NiCAR-A/Boltz2_NiCAR-A_WT_ATP/NiCAR-A_WT_ATP_model_0_out/NiCAR-A_WT_ATP_model_0_out.pdb',
            'NiCAR-A/Boltz2_NiCAR-A_WT_Oleoyl-AMP_Sucrose/NiCAR-A_WT_Oleoyl-AMP_Sucrose_model_0.pdb',
            'NiCAR-A/Boltz2_NiCAR-A_WT_Oleoyl-AMP_Sucrose/NiCAR-A_WT_Oleoyl-AMP_Sucrose_model_0_out/NiCAR-A_WT_Oleoyl-AMP_Sucrose_model_0_out.pdb',
        ],
        'MpCAR-A': [
            'MpCAR-A/Boltz2_MpCAR-A_WT_ATP/MpCAR-A_WT_ATP_model_0.pdb',
            'MpCAR-A/Boltz2_MpCAR-A_WT_ATP/MpCAR-A_WT_ATP_model_0_out/MpCAR-A_WT_ATP_model_0_out.pdb',
            'MpCAR-A/Boltz2_MpCAR-A_WT_Oleoyl-AMP_Sucrose/MpCAR-A_WT_Oleoyl-AMP_Sucrose_model_0.pdb',
            'MpCAR-A/Boltz2_MpCAR-A_WT_Oleoyl-AMP_Sucrose/MpCAR-A_WT_Oleoyl-AMP_Sucrose_model_0_out/MpCAR-A_WT_Oleoyl-AMP_Sucrose_model_0_out.pdb',
        ]
    }
    struct_dir = data_folder + subfolders['pdb'] + data_subfolder + '/'
    msa_dir = data_folder + subfolders['msa'] + data_subfolder + '/'
    pocket_ligand_dist_min_thres = 5
    pocket_residue_dist_min_thres = 6.5
    ligname_list = ['ATP', 'Oleoyl-AMP']
    struct_grp_list = list(struct_fname_dict.keys())
    dist_structAli_fpath = f"{struct_dir}{'_'.join(struct_grp_list)}_ligands_pockets_DIST.csv"
    align_calculate_distances = False
    plot_results = True

    # align structures and calculate distances from residues to ligand
    if align_calculate_distances:

        # initialize yasara
        yasara.info.mode = 'txt'
        yasara.Console('Off')
        yasara.Clear()

        # --- Align structures and get MSA ---
        for struct_num, (struct_grp, struct_fnames) in enumerate(struct_fname_dict.items()):
            # load pdb
            yasara.LoadPDB(filename=f'{struct_dir}{struct_fnames[0]}')
            # align structure
            if struct_num > 0:
                yasara.AlignMol(f'Obj {struct_num + 1} and protein', 'Obj 1 and protein', method='MUSTANG')
                print(f'Aligned Obj {struct_num + 1} to Obj 1.')
        # save alignment
        seq_align_fpath = f"{msa_dir}{'_'.join(list(struct_fname_dict.keys()))}_StructAli.fasta"
        yasara.SaveAli('!1', '1', filename=seq_align_fpath, format='FASTA')
        seqs_ali, seq_names, _ = fetch_sequences_from_fasta(seq_align_fpath)
        df_structAli = {}
        for seq_ali, seq_name in zip(seqs_ali, list(struct_fname_dict.keys())):
            pos_ali = []
            pos = 1
            for aa in seq_ali:
                if aa=='-':
                    pos_ali.append(-1)
                else:
                    pos_ali.append(pos)
                    pos += 1
            df_structAli[f'position_{seq_name}'] = pos_ali
            df_structAli[f'aa_{seq_name}'] = list(seq_ali)
        df_structAli = pd.DataFrame(df_structAli)
        n = len(df_structAli)
        df_structAli.insert(0, 'ali_index', list(range(n)))
        print(df_structAli)

        # --- Get distances between residues and ligands / pockets ---
        df_structgrp_combined_list = []
        # load structure
        for struct_grp, struct_fnames in struct_fname_dict.items():
            # clear scene
            yasara.Clear()
            df_list_bystructgrp = []
            for i, struct_fname in enumerate(struct_fnames):
                is_pocket_struct = False
                yasara.LoadPDB(filename=f'{struct_dir}{struct_fname}')
                if i>0:
                    yasara.AlignMol(f'Obj {i+1} and protein', 'Obj 1 and protein', method='MUSTANG')
                    print(f'Aligned Obj {i+1} to Obj 1.')

                # list residue names
                res_names = yasara.NameRes(f'Obj {i+1} and not protein')

                # if structure with ligands but without pockets
                if res_names[0]!='STP':
                    # delete 2nd ligand (sucrose)
                    yasara.DelMol(f'Obj {i+1} and C')

                # if structure with pockets
                else:
                    # update is_pocket_struct bool
                    is_pocket_struct = True
                    # delete protein
                    yasara.DelMol(f'Obj {i+1} and protein')

                    # get res_ids
                    res_list = yasara.ListRes(f'Obj {i + 1} and not protein')

                    # rename pocket residues
                    res_names_updated = [f'P{k}' if k>9 else f'P0{k}' for k in range(1,len(res_names)+1)]
                    for res_id, res_name_updated in zip(res_list, res_names_updated):
                        yasara.NameRes(res_id, res_name_updated)
                    res_names = yasara.NameRes(f'Obj {i + 1} and not protein')
                    print(i + 1, len(res_names), res_names)

                    # get minimum distances between pockets and ligands
                    pocket_names_to_keep = []
                    pocket_names_to_remove = []
                    for res_id, res_name in zip(res_list, res_names):
                        dist_min = round(yasara.GroupDistance(f'Obj {i+1} and {res_id}', f'Obj {i} and Mol B',center='Closest'),3)
                        if dist_min < pocket_ligand_dist_min_thres:
                            pocket_names_to_keep.append(res_name)
                            print(res_name, dist_min)
                        else:
                            pocket_names_to_remove.append(res_name)

                    # delete all pockets not within 6.5 angstrom of any atom of ligand
                    for pocket_name in pocket_names_to_remove:
                            yasara.DelRes(f'Obj {i+1} and {pocket_name}')

                # get final list of residue names
                res_names = yasara.NameRes(f'Obj {i + 1} and not protein')
                res_ids = yasara.ListRes(f'Obj {i + 1} and not protein')
                print(i + 1, len(res_names), res_names)

                # get distances between protein residues and ligands & binding pocket
                if i%2==1:

                    # get ligname
                    ligname = ligname_list[int(np.floor(i / 2))]

                    # get protein residues
                    protein_res_id_list = yasara.ListRes(f'Obj {i} and protein', format='RESNUM')
                    protein_res_aa_list = list(yasara.SequenceRes(f'Obj {i} and protein')[0])
                    residue_dist_df = []
                    pocket_residue_list = []
                    for protein_res_id, protein_res_aa in zip(protein_res_id_list, protein_res_aa_list):

                        # initialize res_dist_dict
                        res_dist_dict = {'position': protein_res_id, 'aa': protein_res_aa}

                        # distance to pocket
                        for pocket_res_name in res_names:
                            dist_min = round(yasara.GroupDistance(f'Obj {i} and Mol A and Res {protein_res_id}', f'Obj {i+1} and Res {pocket_res_name}', center='Closest'), 3)
                            dist_avg = round(yasara.GroupDistance(f'Obj {i} and Mol A and Res {protein_res_id}', f'Obj {i + 1} and Res {pocket_res_name}', center='Geometric'), 3)
                            res_dist_dict.update({f'dist_min_{pocket_res_name}_{ligname}': dist_min, f'dist_avg_{pocket_res_name}_{ligname}': dist_avg})

                        # distance to ligand
                        dist_min = round(yasara.GroupDistance(f'Obj {i} and Mol A and Res {protein_res_id}', f'Obj {i} and Mol B', center='Closest'), 3)
                        dist_avg = round(yasara.GroupDistance(f'Obj {i} and Mol A and Res {protein_res_id}', f'Obj {i} and Mol B {pocket_res_name}', center='Geometric'), 3)
                        res_dist_dict.update({f'dist_min_{ligname}': dist_min, f'dist_avg_{ligname}': dist_avg})
                        residue_dist_df.append(res_dist_dict)

                    # get distance dataframe
                    residue_dist_df = pd.DataFrame(residue_dist_df)
                    df_list_bystructgrp.append(residue_dist_df)

            # combine dataframes in df_list_bystructgrp
            df_structgrp_combined = df_list_bystructgrp[0]
            for df in df_list_bystructgrp[1:]:
                df_structgrp_combined = df_structgrp_combined.merge(df, on=['position', 'aa'])
            df_structgrp_combined = df_structgrp_combined[['position','aa'] + [c for c in df_structgrp_combined.columns if c.find('dist_min')>-1]]
            df_structgrp_combined_list.append(df_structgrp_combined)
            # save dataframe as CSV
            df_structgrp_combined.to_csv(f"{struct_dir}{struct_grp}_ligands_pockets_DIST.csv")

            # save aligned structures for structgroup as scene
            yasara.SaveSce(f"{struct_dir.replace('pdb/','sce/')}{struct_grp}_ligands_pockets.sce")

        # --- merge into overall struct aligned dataframe ---
        for struct_num, struct_grp in enumerate(struct_fname_dict.keys()):
            df_structgrp_combined = df_structgrp_combined_list[struct_num]
            df_structgrp_combined = df_structgrp_combined.rename(columns={c:f'{c}_{struct_grp}' for c in df_structgrp_combined.columns.tolist()})
            df_structAli = df_structAli.merge(df_structgrp_combined, on=[f'position_{struct_grp}', f'aa_{struct_grp}'], how='outer')
        df_structAli = df_structAli.sort_values(by='ali_index').reset_index(drop=True)
        df_structAli.to_csv(dist_structAli_fpath)


    if plot_results:
        # --- filter for pockets close to ligands ---
        df_structAli = pd.read_csv(dist_structAli_fpath, index_col=0)
        df_structAli_ligandpocket = df_structAli[
            (df_structAli['dist_min_ATP_NiCAR-A']<pocket_residue_dist_min_thres) | \
            (df_structAli['dist_min_Oleoyl-AMP_NiCAR-A'] < pocket_residue_dist_min_thres) | \
            (df_structAli['dist_min_ATP_MpCAR-A'] < pocket_residue_dist_min_thres) | \
            (df_structAli['dist_min_Oleoyl-AMP_MpCAR-A'] < pocket_residue_dist_min_thres)
            ]

        print(df_structAli_ligandpocket)

        # --- plot distances of pocket residues to ligands ---
        x = np.arange(len(df_structAli_ligandpocket))# df_structAli_ligandpocket['ali_index'].to_numpy()
        color_list = ['purple', 'orange', 'b', 'r', ]
        marker_dict = {'ATP': '^', 'Oleoyl-AMP':'o'}
        legend = []
        plot_num = 0
        fig, ax = plt.subplots(1,1, figsize=(20,8))
        for ligand in ['ATP', 'Oleoyl-AMP', ]: # ['Oleoyl-AMP']: #
            for protein in ['NiCAR-A', 'MpCAR-A']:
                ax.plot(x, df_structAli_ligandpocket[f'dist_min_{ligand}_{protein}'].to_numpy(), c=color_list[plot_num], marker=marker_dict[ligand], markersize=4, alpha=0.7, linewidth=0.7)
                legend.append(f'{protein} <> {ligand}')
                plot_num += 1
        ax.legend(legend)
        xmin, xmax = ax.get_xlim()
        dists_NiCAR_OAMP = df_structAli_ligandpocket[f'dist_min_Oleoyl-AMP_NiCAR-A'].to_numpy()
        dists_MpCAR_OAMP = df_structAli_ligandpocket[f'dist_min_Oleoyl-AMP_MpCAR-A'].to_numpy()
        for xpoint in x:
            res_NiCAR = str(df_structAli_ligandpocket.iloc[xpoint]['position_NiCAR-A'])
            aa_NiCAR = df_structAli_ligandpocket.iloc[xpoint]['aa_NiCAR-A']
            res_MpCAR = str(df_structAli_ligandpocket.iloc[xpoint]['position_MpCAR-A'])
            aa_MpCAR = df_structAli_ligandpocket.iloc[xpoint]['aa_MpCAR-A']
            res_aa_NiCAR = res_NiCAR + aa_NiCAR
            res_aa_MpCAR = res_MpCAR + aa_MpCAR
            if aa_NiCAR=='-' or aa_MpCAR=='-' or aa_NiCAR==aa_MpCAR:
                fontweight='light'
            else:
                fontweight = 'bold'
            dist_NiCAR = dists_NiCAR_OAMP[xpoint]
            dist_MpCAR = dists_MpCAR_OAMP[xpoint]
            ni_gt_mp = (dist_NiCAR>dist_MpCAR)*1
            marg = 0.25
            ax.text(xpoint-0.6, dist_NiCAR-((-1)**ni_gt_mp)*marg, res_aa_NiCAR, fontsize=5, c='b', fontweight=fontweight, verticalalignment='center')
            ax.text(xpoint-0.6, dist_MpCAR+((-1)**ni_gt_mp)*marg, res_aa_MpCAR, fontsize=5, c='r', fontweight=fontweight, verticalalignment='center')

        ax.hlines(pocket_residue_dist_min_thres, xmin, xmax, linestyle='--', color='k', linewidth=0.8)
        ax.set_ylabel('Distance (angstrom)')
        ax.set_xlabel('Binding pocket residue index')
        plt.title('Distance of binding pocket residues to docked ligands (ATP, Oleoyl-AMP)')
        plt.show()




