import Bio.PDB
import os
import re
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from variables import address_dict, subfolders
from utils import sort_list, fetch_sequences_from_fasta

def dist(p1, p2):
    dx = p1[0] - p2[0]
    dy = p1[1] - p2[1]
    dz = p1[2] - p2[2]
    return math.sqrt(dx**2 + dy**2 + dz**2)
def read_atoms(file, chain=".", model=1):
    pattern = re.compile(chain)
    current_model = model
    atoms = []
    for line in file:
        line = line.strip()
        if line.startswith("ATOM") and current_model == model:
            type = line[12:16].strip()
            chain = line[21:22]
            if type == "CA" and re.match(pattern, chain):
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append((x, y, z))
        elif line.startswith("MODEL"):
            current_model = int(line[10:14].strip())
    return atoms
def compute_distances(atoms, threshold):
    distances = np.zeros((len(atoms), len(atoms)))
    contacts = []
    for i in range(len(atoms)-2):
        for j in range(i+2, len(atoms)):
            pairdist = dist(atoms[i], atoms[j])
            distances[i,j] = pairdist
            distances[j, i] = pairdist
            if pairdist < threshold:
                contacts.append([i+1, j+1])
    return distances, contacts
def write_output(contacts, file):
    for c in contacts:
        file.write("\t".join(map(str, c))+"\n")

def contacts_to_map(contacts):
    all_res = []
    for pairs in contacts:
        all_res += pairs
    all_res = list(set(all_res))
    all_res.sort()
    final_res = all_res[-1]
    res_list = list(range(1,final_res+1))
    contact_map = pd.DataFrame(np.zeros((final_res,final_res)), index=res_list, columns=res_list)
    for (res1, res2) in contacts:
        contact_map.loc[res1, res2] = 1
        contact_map.loc[res2, res1] = 1
    return contact_map

def pdb_to_cm(file, threshold, chain=".", model=1):
    atoms = read_atoms(open(file, 'r'), chain, model)
    return compute_distances(atoms, threshold)

def map_to_contacts(fpath, sheet_name, upper_filt_thres, filter_out_positions=(50, 200)):

    # get map
    scores = pd.read_excel(pd.ExcelFile(fpath), sheet_name=sheet_name, index_col=0).to_numpy()
    # zero-out scores on upper triangle
    scores = np.tril((scores>upper_filt_thres)*1)

    # get pairs where value is above a threshold
    idxs = np.argwhere(scores==1) + 1
    idxs_filt = idxs.copy()

    # filter out certain positions
    if filter_out_positions is not None:
        idxs_filt = idxs_filt[(~(idxs_filt[:,0]<=filter_out_positions[0]) & ~(idxs_filt[:,1]<=filter_out_positions[0]))]
        idxs_filt = idxs_filt[(~(idxs_filt[:, 0] >= filter_out_positions[1]) & ~(idxs_filt[:, 1] >= filter_out_positions[1]))]
        print(f'# of residue pairs >> BEFORE / AFTER filtering out positions: {len(idxs)} / {len(idxs_filt)}')

    # reorder pairs
    idxs_filt_sorted = []
    for idxpair in idxs_filt:
        idxs_filt_sorted.append(sort_list(list(idxpair)))
    idxs_filt_sorted = sort_list(idxs_filt_sorted)
    return idxs_filt_sorted

def print_residue_pairs_in_contact_map(scores):

    # zero-out scores on upper triangle
    scores = np.tril(scores)

    # get pairs where value is above a threshold
    idxs = np.argwhere(scores == 1) + 1
    idxs_filt = idxs.copy()
    #
    # # filter out certain positions
    # idxs_filt = idxs_filt[(~(idxs_filt[:, 0] <= 50) & ~(idxs_filt[:, 1] <= 50))]
    # idxs_filt = idxs_filt[(~(idxs_filt[:, 0] >= 200) & ~(idxs_filt[:, 1] >= 200))]
    # print(f'# of residue pairs >> BEFORE / AFTER filtering out positions: {len(idxs)} / {len(idxs_filt)}')

    # reorder pairs
    idxs_filt_sorted = []
    for idxpair in idxs_filt:
        idxs_filt_sorted.append(sort_list(list(idxpair)))
    idxs_filt_sorted = sort_list(idxs_filt_sorted)
    print(idxs_filt_sorted)
    return idxs_filt_sorted


def main():
    data_folder = address_dict['PIPS2']
    data_subfolder = 'UPO_batch1/'
    msa_fname = f'UPO_aligned_clustalo.fasta'
    target_seqs_fpath = f'{data_folder}{subfolders["sequences"]}UPO_batch1_shortlisted2.fasta'
    res_fpath = f'{data_folder}{subfolders["conservation_analysis"]}{data_subfolder}'
    EC_analysis_fname_suffix = '_prefiltered'
    EC_distance_thres = 5 # 12
    EC_upper_filt_thres = 0.2
    filter_out_positions = None  # (50, 200)
    print_residue_pairs_in_contact = True
    plot_raw_contact_map_only = True

    target_seqs, target_seqs_name_list, _ = fetch_sequences_from_fasta(target_seqs_fpath)
    for seq_name, seq in zip(target_seqs_name_list, target_seqs):
        # get contact map
        pdb_filename = seq_name + '.pdb'
        pdb_fpath = data_folder + subfolders['pdb'] + data_subfolder[:-1] + '_enzymeOnly/' + pdb_filename
        print(seq_name, pdb_fpath)
        distances, contacts = pdb_to_cm(pdb_fpath, EC_distance_thres)
        print('Seq len:', len(seq))
        print('Distance map size:', distances.shape)
        contact_map_thres = contacts_to_map(contacts)
        # write_output(contacts, open(data_folder + subfolders['pdb'] + data_subfolder + seq_name + f'_contacts_thres={EC_distance_thres}.txt', "w"))

        if print_residue_pairs_in_contact:
            contact_map_idx_pairs = print_residue_pairs_in_contact_map(contact_map_thres)

        # plot contact maps with EC coupled residues overlaid
        fig, ax = plt.subplots(2,1, figsize=(6,10))
        # plot distance map
        ax[0].imshow(distances, cmap='viridis_r')
        # plot contact map
        ax[1].imshow(contact_map_thres)
        # turn off border
        for ax_idx in [0,1]:
            for edge in ['bottom', 'top', 'right', 'left']:
                ax[ax_idx].spines[edge].set_color('None')
        ax[0].set_title('Distance map')
        ax[1].set_title('Contact map (thresholded distance map)')
        plt.suptitle(seq_name, y=0.93)

        if not plot_raw_contact_map_only:
            # draw lines for conservation analysis positions
            conservation_analysis_fpath = f'{res_fpath}{msa_fname[:-6]}_sift-ShanEntropy_selected_RESULTS.csv'
            conservation_analysis_res = pd.read_csv(conservation_analysis_fpath, index_col=0)
            seq_conservation_labels = conservation_analysis_res.loc[seq_name].to_numpy()
            idxs_unconserved = np.argwhere(seq_conservation_labels==1).reshape(-1,)
            idxs_conserved = np.argwhere(seq_conservation_labels==4).reshape(-1,)

            # overlay with EC coupled residues from Potts model
            EC_analysis_fpath = f'{data_folder}{subfolders["conservation_analysis"]}{data_subfolder}{msa_fname[:-6]}{EC_analysis_fname_suffix}_PottsEC-APC_selected.xlsx'
            ec_pairs = np.array(map_to_contacts(EC_analysis_fpath, seq_name, upper_filt_thres=EC_upper_filt_thres, filter_out_positions=filter_out_positions))

            ## annotate plot distance map ##
            # draw lines for conserved / unconserved positions
            ax[0].scatter(idxs_unconserved, np.ones(len(idxs_unconserved))*-len(distances)/100, color='orange', s=0.1)
            ax[0].scatter(idxs_conserved, np.ones(len(idxs_conserved))*-len(distances)/100,color='purple', s=0.1)
            ax[0].scatter(np.ones(len(idxs_unconserved))*len(distances)*1.01, idxs_unconserved, color='orange', s=0.1)
            ax[0].scatter(np.ones(len(idxs_conserved))*len(distances)*1.01, idxs_conserved, color='purple', s=0.1)
            # scatter EC positions
            ax[0].scatter(ec_pairs[:,0]-1, ec_pairs[:,1]-1, color='r', s=0.4, alpha=0.7)
            ax[0].scatter(ec_pairs[:,1]-1, ec_pairs[:,0]-1, color='r', s=0.4, alpha=0.7)

            # draw lines for conserved / unconserved positions
            ax[1].scatter(idxs_unconserved, np.ones(len(idxs_unconserved))*-len(distances)/100,  color='orange', s=0.1)
            ax[1].scatter(idxs_conserved, np.ones(len(idxs_conserved))*-len(distances)/100, color='purple', s=0.1)
            ax[1].scatter(np.ones(len(idxs_unconserved))*len(distances)*1.01, idxs_unconserved, color='orange', s=0.1)
            ax[1].scatter(np.ones(len(idxs_conserved))*len(distances)*1.01, idxs_conserved, color='purple', s=0.1)
            # scatter EC positions
            ax[1].scatter(ec_pairs[:,0]-1, ec_pairs[:,1]-1, color='r', s=0.4, alpha=0.7)
            ax[1].scatter(ec_pairs[:,1]-1, ec_pairs[:,0]-1, color='r', s=0.4, alpha=0.7)

        plt.savefig(f'{res_fpath}{data_subfolder[:-1]}_proteinMap_{seq_name}.png', bbox_inches='tight')
        plt.show()



if __name__ == "__main__":
    main()