###
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import pyhmmer
from pyhmmer.hmmer import hmmscan, hmmsearch
from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile
import pandas as pd
import pickle
import os
from variables import address_dict, subfolders
database_folder = address_dict['databases']
sequence_subfolder = subfolders['sequences']
hmm_subfolder = subfolders['hmm']

def invert_result_dict(res_byquery):
    """
    Convert res_byquery to res_bytarget
    """
    res_bytarget = {}
    for i, (query, match_res_list) in enumerate(res_byquery.items()):
        for match in match_res_list:
            target_id = (match['target_accession'], match['target_name'])
            if target_id not in res_bytarget:
                res_bytarget.update({target_id: []})
            res_bytarget[target_id].append({
                'target_accession': match['query_accession'],
                'target_name': match['query_name'],
                'target_description': match['query_description'],
                'target_sequence': match['query_sequence'],
                'evalue': match['evalue'],
                'query_accession': match['target_accession'],
                'query_name': match['target_name'],
                'query_description': match['target_description'],
                'query_sequence': match['target_sequence'],
            })
    return res_bytarget
def domtbl2csv(domtbl_fpath, delete_domtbl=True):
    """
    loading function adapted from SearchIO
    https://github.com/biopython/biopython/blob/master/Bio/SearchIO/HmmerIO/hmmer3_domtab.py
    """
    all_lines = []
    names = ['target_name', 'target_accession', 'target_len', 'query_name', 'query_accession',
             'query_len', 'full_e-val', 'full_score', 'full_bias', 'dom_number',
             'dom_of', 'dom_c-Evalue', 'dom_i-Evalue', 'dom_score', 'dom_bias',
             'hmmcoord_from', 'hmmcoord_to', 'alncoord_from', 'alncoord_to', 'envcoord_from',
             'envcoord_to', 'acc', 'target_description']

    with open(domtbl_fpath, 'r') as f:
        lines = (line for idx, line in enumerate(f) if idx > 2)
        for line in lines:
            cols = [x for x in line.strip().split(' ') if x]
            # if len(cols) > 23, we have extra description columns
            # combine them all into one string in the 19th column
            if len(cols) > 23:
                cols[22] = ' '.join(cols[22:])
                del (cols[23:])
                # assert len(cols) == 23
            elif len(cols) < 23:
                cols.append('')
                assert len(cols) == 23
            all_lines.append(cols)
        df = pd.DataFrame(all_lines, columns=names)
        df.to_csv(domtbl_fpath[:-7]+'.csv')
    # remove domtbl after saving CSV
    if delete_domtbl:
        os.remove(domtbl_fpath)
    return df
def combine_domtbl_csv(res_byquery, domtbl_fileprefix, domtbl_dir, delete_csv=True):
    # combine individual domtbl results
    domtbl_flist = [f for f in list(os.listdir(domtbl_dir)) if f.find(domtbl_fileprefix) > -1]
    df_all = []
    for f in domtbl_flist:
        d = pd.read_csv(domtbl_dir+f, index_col=0).to_dict(orient='records')
        df_all += d
    df_all = pd.DataFrame.from_dict(df_all)

    # add columns from res_byquery: query_description, query_sequence, target_sequence
    query_description_list = []
    query_sequence_list = []
    target_sequence_list = []
    for query_accession, query_name, target_accession, target_name in zip(df_all['query_accession'].tolist(), df_all['query_name'].tolist(), df_all['target_accession'].tolist(), df_all['target_name'].tolist()):
        res_byquery_df = pd.DataFrame.from_dict(res_byquery[(query_accession, query_name)])
        match = res_byquery_df[(res_byquery_df.target_accession==target_accession) & (res_byquery_df.target_name==target_name)].iloc[0]
        query_description_list.append(match['query_description'])
        query_sequence_list.append(match['query_sequence'])
        target_sequence_list.append(match['target_sequence'])
    df_all['query_description'] = query_description_list
    df_all['query_sequence'] = query_sequence_list
    df_all['target_sequence'] = target_sequence_list

    # save as CSV
    df_all.to_csv(f'{domtbl_dir}{domtbl_fileprefix}_hits.csv')
    print(f'Saved results to {domtbl_dir}{domtbl_fileprefix}_hits.csv')

    # remove individual CSVs after combining them
    if delete_csv:
        for f in domtbl_flist:
            os.remove(domtbl_dir+f)
    return df_all
def get_hmm_seq_attributes(idx, hmm_or_seq_list, hmm_or_seq):
    if hmm_or_seq == 'hmm':
        sequence_idx = hmm_or_seq_list[idx].consensus
        description_idx = hmm_or_seq_list[idx].description.decode('utf-8')
    elif hmm_or_seq == 'seq':
        sequence_idx = hmm_or_seq_list[idx].sequence
        description_idx = hmm_or_seq_list[idx].description
    return sequence_idx, description_idx
def perform_hmmsearch_or_hmmscan(seq_fname, hmm_fname, hmm_dir, res_fname, data_folder, get_matching_hmm_or_seq, hmmsearch_or_hmmscan='hmmsearch', print_res=False):
    query_target_labels_dict = {
        'hmmsearch': {'query': 'hmm', 'target': 'seq'},
        'hmmscan': {'query': 'seq', 'target': 'hmm'}
    }
    target_type = query_target_labels_dict[hmmsearch_or_hmmscan]['target']
    query_type = query_target_labels_dict[hmmsearch_or_hmmscan]['query']
    alphabet = pyhmmer.easel.Alphabet.amino()
    domtbl_dir = f'{data_folder}{hmm_subfolder}'

    # iterate over profiles
    with HMMFile(f"{hmm_dir}{hmm_fname}.hmm") as hmms:
        # iterate over sequences
        with SequenceFile(f'{data_folder}{sequence_subfolder}{seq_fname}.fasta', digital=True) as seqs:
            # search for matching sequences given a HMM profile
            if hmmsearch_or_hmmscan=='hmmsearch':
                all_hits = list(hmmsearch(hmms, seqs))
                total = sum(len(hits) for hits in all_hits)
            # scan a given sequence against a database of HMM profiles
            elif hmmsearch_or_hmmscan=='hmmscan':
                all_hits = list(hmmscan(seqs, hmms))
                total = sum(len(hits) for hits in all_hits)
            print(f'{total} hits found.')

    # get hmms
    with HMMFile(f"{hmm_dir}{hmm_fname}.hmm") as hmms:
        hmm_list = list(hmms)
        hmm_dict = {(hmm.accession.decode('utf-8'), hmm.name.decode('utf-8')): i for i, hmm in enumerate(hmm_list)}
        print(len(hmm_list))
    # get sequences
    with SequenceFile(f'{data_folder}{sequence_subfolder}{seq_fname}.fasta') as seqs:
        seq_list = seqs.read_block()
        seq_dict = {(seq.accession.decode('utf-8'), seq.name.decode('utf-8')): i for i, seq in enumerate(seq_list)}
        print(len(seq_list))

    # aggregate and print results
    res_byquery = {}
    for i, hits_byquery in enumerate(all_hits):
        match_res_list = []
        if len(hits_byquery)>0:
            # get query id & description
            query_accession = hits_byquery.query_accession.decode('utf-8') if hits_byquery.query_accession is not None else ''
            query_name = hits_byquery.query_name.decode('utf-8') if hits_byquery.query_name is not None else ''
            # get query sequence and description
            if query_type=='hmm':
                query_sequence, query_description = get_hmm_seq_attributes(i, hmm_list, query_type)
            elif query_type=='seq':
                query_sequence, query_description = get_hmm_seq_attributes(i, seq_list, query_type)

            # save hit results for each query
            domtbl_fname = f'{res_fname}_{query_type}query{i}_hits.domtbl'
            with open(f"{domtbl_dir}{domtbl_fname}", "wb") as f:
                hits_byquery.write(f, format='domains')
            df = domtbl2csv(f'{data_folder}{hmm_subfolder}{domtbl_fname}', delete_domtbl=True)

            # iterate over matches
            for j, match in enumerate(hits_byquery):
                target_accession = match.accession.decode('utf-8') if match.accession is not None else ''
                target_name = match.name.decode('utf-8') if match.name is not None else ''
                query_accession = query_accession if query_accession!='' else '-'
                query_name = query_name if query_name != '' else '-'
                # get target sequence
                if target_type=='hmm':
                    target_sequence = hmm_list[hmm_dict[(target_accession, target_name)]].consensus
                elif target_type=='seq':
                    target_sequence = seq_list[seq_dict[(target_accession, target_name)]].sequence

                match_res = {
                    'target_accession': target_accession if target_accession!='' else '-',
                    'target_name': target_name if target_name!='' else '-',
                    'target_description': match.description.decode('utf-8') if match.description is not None else None,
                    'target_sequence': target_sequence,
                    'target_length': match.length,
                    'evalue': match.evalue,
                    'query_accession': query_accession if query_accession!='' else '-',
                    'query_name': query_name if query_name!='' else '-',
                    'query_description': query_description,
                    'query_sequence': query_sequence,
                    'domains': match.domains,
                }
                match_res_list.append(match_res)
                if print_res:
                    print(f'query {i} ({query_accession}, {query_name}) <> hit {j}: {match_res}')
                    # print matching domain alignments
                    for k, domain in enumerate(match_res['domains']):
                        ali = domain.alignment
                        print(f'Alignment of hit {j} domain {k} with query hmm profile')
                        print(ali)
            res_byquery.update({(query_accession, query_name): match_res_list})

    # combine all relevant results into one file, update with sequence information extracted from res_byquery
    df_all = combine_domtbl_csv(res_byquery, f'{res_fname}_{query_type}query', domtbl_dir, delete_csv=True)

    # get results in terms of desired output type (matching HMMs or sequences)
    if get_matching_hmm_or_seq != target_type:
        res = invert_result_dict(res_byquery)
        output_col_prefix = 'query'
    else:
        res = res_byquery.copy()
        output_col_prefix = 'target'

    if get_matching_hmm_or_seq=='hmm':
        # get HMM accessions & seed sequences
        output_accession_all = df_all[f'{output_col_prefix}_accession'].tolist()
        output_name_all = df_all[f'{output_col_prefix}_name'].tolist()
        output_seq_all = df_all[f'{output_col_prefix}_sequence'].tolist()
        # # get matching HMMs given reference sequences
        # if get_matching_hmm_or_seq == 'hmm':
        #     ## extract relevant HMM profile and save separately
        # # get matching sequences given reference HMM profile
        # if get_matching_hmm_or_seq=='seq':
        #     ## get MSA for sequences
        #
        #     ## get HMM from MSA

    print('\nRESULTS')
    for query, hits in res.items():
        print(f'{query} hits:')
        for hit in hits:
            print(hit)

    return res, df_all

def main():

    # get sequences to query from sequence file
    data_folder = address_dict['examples']
    seq_fname = 'GOh1001b'
    hmm_fname = 'Pfam-A'
    if hmm_fname in ['Pfam-A']:
        hmm_dir = database_folder
    else:
        hmm_dir = data_folder + hmm_subfolder
    res_fname = f'{seq_fname}_vs_{hmm_fname}'
    get_matching_hmm_or_seq = 'hmm'
    hmmsearch_or_hmmscan = 'hmmsearch'  # 'hmmscan'#
    print_res = False

    res, df_all = perform_hmmsearch_or_hmmscan(seq_fname, hmm_fname, hmm_dir, res_fname, data_folder, get_matching_hmm_or_seq,
                                               hmmsearch_or_hmmscan, print_res)
    return res, df_all

if __name__ == "__main__":
    res, df_all = main()
