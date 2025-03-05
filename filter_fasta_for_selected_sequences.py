#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import os
import numpy as np
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta
from variables import address_dict, subfolders
data_folder = address_dict['PIPS2']
sequence_subfolder = subfolders['sequences']
msa_subfolder = subfolders['msa']

fasta_fname = 'UPO.fasta'
fasta_fpath = f'{data_folder}{sequence_subfolder}{fasta_fname}'
fasta_fname_selected = 'UPO_filtered'

## GET SEQUENCES TO KEEP ##
fasta_fpath_w_seq_to_keep = f'{data_folder}{msa_subfolder}UPO_cdhit-cluster1_aligned.fasta'
sequence_list, sequence_names, sequence_descriptions = fetch_sequences_from_fasta(fasta_fpath_w_seq_to_keep)
sequence_names_to_keep = sequence_names + ['PBK73312', 'KAE8306975']

## GET SEQUENCES FROM FILE TO FILTER FROM ##
sequence_list, sequence_names, sequence_descriptions = fetch_sequences_from_fasta(fasta_fpath)
sequence_list_selected, sequence_names_selected, sequence_descriptions_selected = [], [], []
sequence_list_unselected, sequence_names_unselected, sequence_descriptions_unselected = [], [], []
print(f'Original # of sequences {len(sequence_names)}')

# get sequences to keep first
for seq, name, desc in zip(sequence_list, sequence_names, sequence_descriptions):
    if name in sequence_names_to_keep:
        sequence_list_selected.append(seq)
        sequence_names_selected.append(name)
        sequence_descriptions_selected.append(desc)
    else:
        sequence_list_unselected.append(seq)
        sequence_names_unselected.append(name)
        sequence_descriptions_unselected.append(desc)
print(f'Selected {len(sequence_names_selected)} sequences')

# further randomly select sequences from those that remain at a given rate
proportion_of_unselected_seq_to_include = 1 # 0.65
num_of_unselected_seq_to_include = int(proportion_of_unselected_seq_to_include*len(sequence_names_unselected))
idxs_to_add = np.random.choice(np.arange(len(sequence_names_unselected)), size=num_of_unselected_seq_to_include, replace=False)
seq_to_add = [sequence_list_unselected[i] for i in idxs_to_add]
seqname_to_add = [sequence_names_unselected[i] for i in idxs_to_add]
seqdesc_to_add = [sequence_descriptions_unselected[i] for i in idxs_to_add]
print(f'Adding {len(seqname_to_add)} sequences')

# append to selected lists
sequence_list_selected += seq_to_add
sequence_names_selected += seqname_to_add
sequence_descriptions_selected += seqdesc_to_add
print(f'Final # of sequences:  {len(sequence_names_selected)}')

# create sequence file
write_sequence_to_fasta(sequence_list_selected, sequence_names_selected, fasta_fname_selected, f'{data_folder}{sequence_subfolder}')