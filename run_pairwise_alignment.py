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
data_folder = address_dict['pips-insilico'] # address_dict['ECOHARVEST']
sequence_subfolder = subfolders['sequences']

# compare 2 sequences
seq_fname_1 = 'GOh1001b.fasta' # 'CALB.fasta'
seq_fname_2 = 'FGGALOX.fasta'

target_seq = [record for record in SeqIO.parse(f'{data_folder}{sequence_subfolder}{seq_fname_1}', "fasta")][0].seq
query_seq = [record for record in SeqIO.parse(f'{data_folder}{sequence_subfolder}{seq_fname_2}', "fasta")][0].seq

# run pairwise alignment
alignments = run_pairwise_alignment(target_seq, query_seq, mode='global', match_score=2, mismatch_score=-1,
                                    open_gap_score=-0.5, extend_gap_score=-0.1, target_end_gap_score=0.0,
                                    query_end_gap_score=0.0, print_alignments=False)

# get first alignment
al = alignments[0]
target_al = al[0]
query_al = al[1]
# get mutations
res_idx = 0
mut_list = []
actual_wt_aa_list = []
for target_aa, query_aa in zip(target_al, query_al):
    if query_aa != '-':
        res_idx += 1
        if target_aa != query_aa and target_aa != '-':
            mut = query_aa + str(res_idx) + target_aa
            mut_list.append(mut)
            actual_wt_aa_list.append(query_seq[res_idx - 1] + str(res_idx))
print(target_al)
print(query_al)
print(mut_list)