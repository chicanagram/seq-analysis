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
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta, split_fasta, deduplicate_fasta_sequences, get_sequences_in_folder, run_pairwise_alignment, get_mutagenesis_sequences, parse_genbank_to_df
data_folder = address_dict['ECOHARVEST']
sequence_subfolder = subfolders['sequences']

# print sequences in dir
seq_fnames = [f for f in os.listdir(f'{data_folder}{sequence_subfolder}') if f.endswith('.fasta')]
sequence_list = get_sequences_in_folder(seq_fnames, seq_dir=f'{data_folder}{sequence_subfolder}')

# get mutation list
fbase = 'GOh1052mut'
mutagenesis_dataset_fpath = f'{data_folder}{subfolders["ml_prediction"]}Input/{fbase}/{fbase}.csv'
df = pd.read_csv(mutagenesis_dataset_fpath, index_col=0)
mutations_list = df['Mutation'].tolist()

# get backbone
seq_fpath = f'{data_folder}{subfolders["sequences"]}GOh1052.fasta'
seqs, _, _ = fetch_sequences_from_fasta(seq_fpath)
seq_backbone = seqs[0]

# get mutants
seq_mutants_list = get_mutagenesis_sequences(mutations_list, seq_backbone, split_on='_')
write_sequence_to_fasta(seq_mutants_list, [f'GOh1052_{mut}' for mut in mutations_list], 'GOh1052_mutagenesis', f'{data_folder}{subfolders["sequences"]}')

