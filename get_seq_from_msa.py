#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta
from variables import address_dict, subfolders
data_folder = address_dict['PIPS2']
sequence_subfolder = subfolders['sequences']
msa_subfolder = subfolders['msa']

msa_fname = 'UPO_cdhit-cluster1_aligned.fasta' # 'UPO_aligned.fasta'
msa_fpath = f'{data_folder}{msa_subfolder}{msa_fname}'
seq_fname = msa_fname[:-6]
strip_dashes = True

# get sequences and names
sequence_list, sequence_names, sequence_descriptions = fetch_sequences_from_fasta(msa_fpath)
# strip dashes from sequence
if strip_dashes:
    sequence_list_stripped = [s.replace('-','') for s in sequence_list]
else:
    sequence_list_stripped = sequence_list

# create sequence file -- strip dashes from sequence
write_sequence_to_fasta(sequence_list_stripped, sequence_names, seq_fname, f'{data_folder}{sequence_subfolder}')