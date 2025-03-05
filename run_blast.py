#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from variables import address_dict, subfolders
sequence_subfolder = subfolders['sequences']
blast_subfolder = subfolders['sequences']

def run_blastp(protein_sequence, blast_fname, data_folder, blast_subfolder=blast_subfolder, expect=10.0, hitlist_size=50):
    """
    Run a blast search given an input sequence
    """
    from Bio.Blast import NCBIWWW
    print('Starting blast search...')
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence, expect=expect, hitlist_size=hitlist_size)
    print(f'Completed blast search.')
    # save result
    if blast_fname is not None:
        with open(f"{data_folder}{blast_subfolder}{blast_fname}.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        print()
    return result_handle

def run_psiblast():
    from Bio.Blast.Applications import NcbipsiblastCommandline
    cline = NcbipsiblastCommandline(help=True)
    NcbipsiblastCommandline(cmd='psiblast', help=True)

def parse_blast(blast_fname, data_folder, blast_subfolder=blast_subfolder, sequence_subfolder=sequence_subfolder):
    """
    Parse XML BLAST output with BioPython
    """
    from Bio.Blast import NCBIXML
    results = []
    fasta = ''
    blast_records = NCBIXML.parse(open(f"{data_folder}{blast_subfolder}{blast_fname}.xml"))
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print('score:', hsp.score)
                print('identities:', hsp.identities)
                print('positives:', hsp.positives)
                print('query:', hsp.query)
                print('match:', hsp.match)
                print('subject:', hsp.sbjct)

                # append to dataframe output
                results.append({
                    "sequence_title": alignment.title,
                    "length": alignment.length,
                    "e_value": hsp.expect,
                    "score": hsp.score,
                    "identity": hsp.identities,
                    "query": hsp.query,
                    "match": hsp.match,
                    "subject": hsp.sbjct
                })

                # append to fasta output
                fasta += '>' + alignment.title + '\n'
                fasta += hsp.sbjct.replace('-','') + '\n'

    # Convert list to DataFrame
    df = pd.DataFrame(results)
    df.to_csv(f"{data_folder}{blast_subfolder}{blast_fname}.csv")

    # save fasta file
    with open(f"{data_folder}{sequence_subfolder}{blast_fname}.fasta", 'w') as f:
        f.write(fasta)

    return df, fasta


def main():

    data_folder = address_dict['examples']
    # seq_fnames = os.listdir(f'{data_folder}{sequence_subfolder}')
    seq_fname = 'GOh1001b'
    record_to_blast = SeqIO.read(f'{data_folder}{sequence_subfolder}{seq_fname}.fasta', format="fasta")
    blast_fname = f'{seq_fname}_blastp'

    # get protein sequence
    protein_sequence = record_to_blast.format("fasta")

    # run blast search
    expect = 10.0
    hitlist_size = 50
    result_handle = run_blastp(protein_sequence, blast_fname=blast_fname, data_folder=data_folder, blast_subfolder=blast_subfolder, expect=expect, hitlist_size=hitlist_size)

    # read blast output
    df, fasta = parse_blast(blast_fname, data_folder, blast_subfolder, sequence_subfolder)

if __name__ == "__main__":
    main()