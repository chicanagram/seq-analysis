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
data_folder = address_dict['SoluProtMut']
seq_fasta = 'soluprotmutdb_2level_LGK.fasta'
fasta_dir = data_folder + subfolders['sequences']
max_res_per_fasta = 200000

# deduplicate sequences
seq_fasta_deduped, seqs_deduped, seq_names_deduped = deduplicate_fasta_sequences(seq_fasta, seq_fasta.replace('.fasta','_deduped'), fasta_dir)
# split fasta into sub-fasta files
split_fasta(max_res_per_fasta, seq_fasta_deduped, fasta_dir)
