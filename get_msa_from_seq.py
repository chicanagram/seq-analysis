#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
from variables import address_dict, subfolders
from utils import run_msa
from plot_utils import visualize_msa

def main():
    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'lipases' # 'sidestream_cocktail' #
    seq_dir = data_folder + subfolders['sequences'] + data_subfolder + '/'
    msa_dir = data_folder + subfolders['msa'] + data_subfolder + '/'
    seq_fname = 'RML_ali.fasta' # 'RML-propeptide-mature_blastp_nr_E1e-05.fasta' # 'exoglucanase_TrichodermaHarzianum.fasta' # 'CARs_litsearch.fasta' #
    msa_method = 'mafft' # 'clustalo' #
    fname_suffix = '' # '_filt_trimmed' #
    msa_fname = seq_fname[:seq_fname.find('.fa')] + f'_{msa_method}' + fname_suffix + '.fasta'
    get_msa = True # False
    plot_msa = 'seaborn' # 'pymsaviz' # None #
    plot_msa_pos_range = None # [200,236] #
    wrap_length = 100 # 400 # 600
    filter_by_refseq = None # 0 #
    label_residues = 'ref' # None # 'consensus' #
    ytick_interval = 1 # 100 #
    show_all_sequences = True # False #
    if ytick_interval==1:
        show_seq_names = True
        pos_int_to_label = 5
    else:
        show_seq_names = False
    savefig = None
    if plot_msa_pos_range is None:
        fontsize = 8
        xtick_interval = 5 # 20
    else:
        fontsize = 20
        xtick_interval = 5

    # get MSA
    if get_msa:
        run_msa(seq_fname, msa_fname, method=msa_method, seq_dir=seq_dir, msa_dir=msa_dir)

    # visualize MSA
    if plot_msa is not None:
        msa_fpath = msa_dir + msa_fname
        visualize_msa(msa_fpath, how=plot_msa, color_scheme='Taylor', plot_msa_pos_range=plot_msa_pos_range,
                      wrap_length=wrap_length, xtick_interval=xtick_interval, ytick_interval=ytick_interval, pos_int_to_label=pos_int_to_label,
                      show_seq_names=show_seq_names, label_residues=label_residues, show_all_sequences=show_all_sequences, fontsize=fontsize, filter_by_refseq=filter_by_refseq, savefig=savefig)

if __name__ == "__main__":
    main()