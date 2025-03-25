#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import os
import numpy as np
import pandas as pd
pd.set_option('display.max_columns', None)
from Bio import SeqIO
from variables import address_dict, subfolders, aaList
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta, split_fasta, deduplicate_fasta_sequences, get_sequences_in_folder, run_pairwise_alignment, get_mutagenesis_sequences, parse_genbank_to_df

data_folder = address_dict['PIPS2'] # address_dict['ECOHARVEST']
sequence_subfolder = subfolders['sequences']
sequence_fname = 'ET096.fasta'
fbase = sequence_fname.split('.')[0] # 'GOh1052mut'
get_all_ss_mutants = True
mutagenesis_dataset_fpath = data_folder + subfolders["ml_prediction"] + 'Input/' + fbase + '/' + fbase + '.csv'

# get all sequences in directory
if sequence_fname is None:
    seq_fnames = [f for f in os.listdir(f'{data_folder}{sequence_subfolder}') if f.endswith('.fasta')]
    sequence_list = get_sequences_in_folder(seq_fnames, seq_dir=f'{data_folder}{sequence_subfolder}')
# get sequence backbone from specific fasta file
else:
    sequences, seq_names, _ = fetch_sequences_from_fasta(data_folder + sequence_subfolder + sequence_fname)
    seq_backbone = sequences[0]
    seq_name = seq_names[0]

# get mutation list by mutating every single position to every possible amino acid
if get_all_ss_mutants:
    mutations_list = []
    for i,wt_aa in enumerate(seq_backbone):
        for mt_aa in aaList:
            if mt_aa != wt_aa:
                mutations_list.append(wt_aa + str(i+1) + mt_aa)
# get mutation list from CSV file
else:
    df = pd.read_csv(mutagenesis_dataset_fpath, index_col=0)
    mutations_list = df['Mutation'].tolist()

# get mutants sequences
seq_mutants_list = get_mutagenesis_sequences(mutations_list, seq_backbone, split_on='_')
write_sequence_to_fasta(seq_mutants_list, [f'{seq_name}_{mut}' for mut in mutations_list], f'{fbase}_mutagenesis', f'{data_folder}{subfolders["sequences"]}')

