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
from variables import address_dict, subfolders
from utils import split_fasta, deduplicate_fasta_sequences

######################################
## SPLIT FASTA INTO SUB-FASTA FILES ##
######################################
data_folder = address_dict['PIPS2']
seq_fasta = 'UPO_batch2.fasta'
fasta_dir = data_folder + subfolders['sequences']
max_res_per_fasta = 8000 # 200000
sort_by_seq_name = True

# deduplicate sequences
seq_fasta_deduped, seqs_deduped, seq_names_deduped = deduplicate_fasta_sequences(seq_fasta, sort_by_seq_name, fasta_fname_out=seq_fasta.replace('.fasta','_deduped'), fasta_dir=fasta_dir)
print(seq_fasta_deduped)

# split fasta into sub-fasta files
split_fasta(max_res_per_fasta, os.path.basename(seq_fasta_deduped), fasta_dir)
