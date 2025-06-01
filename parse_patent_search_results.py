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
from utils import sort_list, split_mutation, fetch_sequences_from_fasta, write_sequence_to_fasta, split_fasta, deduplicate_fasta_sequences, get_sequences_in_folder, run_pairwise_alignment, get_mutagenesis_sequences, parse_genbank_to_df
from scrape_sequence_patents import extract_patent_data
data_folder = address_dict['ECOHARVEST']
sequence_subfolder = subfolders['sequences']
patents_subfolder = subfolders['patents']
seq_name_base_list = ['CALB', 'CALA', 'RML', 'TLL', 'BSL', 'UML', 'TCL'] # ['MmCAR'] #
overall_results_fname_prefix = 'Lipase'
parse_patent_search_results = True
scrape_patent_details = True
analyse_patent_search_results = True
perform_pairwise_alignment = False

if parse_patent_search_results:
    query_seq_count = 0
    for k, seq_name_base in enumerate(seq_name_base_list):
        print(seq_name_base)
        query_seq_fname = seq_name_base + '.fasta'
        target_seq_fname = seq_name_base+'_patent_hits.fasta'
        target_seq_aligned_fname = seq_name_base+ '_patent_hits_aligned.fasta'
        blastp_description = seq_name_base+'_patent_hits_Descriptions.csv'
        db_fname = seq_name_base+'_patent_hits_GenBank.txt'

        ###########################
        # get pairwise alignments #
        ###########################
        # get sequences
        seqs, seq_names, seq_descriptions = fetch_sequences_from_fasta(data_folder+patents_subfolder+target_seq_fname)
        seqs_aligned, seq_names_aligned, _ = fetch_sequences_from_fasta(data_folder + patents_subfolder + target_seq_aligned_fname)
        # get alignment start and end positions
        ali_starts = []
        ali_ends = []
        ali_strs = [s[s.find(':')+1:] for s in seq_names_aligned]
        for ali_str in ali_strs:
            ali_startend = [int(v) for v in ali_str.split('-')]
            ali_starts.append(ali_startend[0])
            ali_ends.append(ali_startend[1])
        print('# of target sequences:', len(seqs))

        # get blast results description
        df = pd.read_csv(data_folder+patents_subfolder+blastp_description)
        # get db info
        with open(data_folder+patents_subfolder+db_fname, 'r') as f:
            db_txt = f.read().split('//')

        db = []
        for entry in db_txt:
            if entry.find('LOCUS')>-1:
                entry_dict = parse_genbank_to_df(entry).iloc[0].to_dict()
                db.append(entry_dict)
        db = pd.DataFrame(db)
        # db.to_csv(data_folder+patents_subfolder+db_fname.replace('.txt','.csv'))
        print(db)

        print('# of seqs:', len(seqs))
        query_seq = [record for record in SeqIO.parse(f'{data_folder}{sequence_subfolder}{query_seq_fname}', "fasta")][0].seq
        num_mut_thres = 25 #50 #
        query_cover_thres = 30 # 20
        res_parsed = []
        for i, (target_seq, target_seq_aligned, seq_name, seq_desc, ali_start, ali_end) in enumerate(zip(seqs, seqs_aligned, seq_names, seq_descriptions, ali_starts, ali_ends)):
            seq_desc = seq_desc.replace(seq_name+' ', '')
            hit = df.loc[df['Description']==seq_desc].iloc[0].to_dict()
            seq_len = len(target_seq)
            max_score = int(hit['Max Score'])
            total_score = int(hit['Total Score'])
            query_cover = float(hit['Query Cover'][:-1])
            percent_identity = float(hit['Per. ident'])
            db_entry = db[db['Accession']==seq_name].iloc[0]
            title = db_entry['Title']
            journal = db_entry['Journal']
            authors = db_entry['Authors']
            date = db_entry['Date'][-4:]
            num_mut = int(np.ceil(len(query_seq)*query_cover/100*(1-percent_identity/100)))
            print(i, seq_name, seq_desc, '; seq len:', seq_len, '; query cover:', query_cover, '; percent identity:', percent_identity, '; # of mutations:', num_mut)

            ### CONDITIONS FOR ADDING ENTRY ###
            if num_mut > num_mut_thres or query_cover < query_cover_thres:
                print(f'Alignment not performed: # mut = {num_mut}; query cover: {query_cover}')
                continue

            # perform pairwise alignment
            if perform_pairwise_alignment:
                target_seq_to_align = target_seq
            else:
                target_seq_to_align = target_seq_aligned
            alignments = run_pairwise_alignment(target_seq, query_seq,
                                   mode='global', match_score=2, mismatch_score=-1,
                                   open_gap_score=-0.5, extend_gap_score=-0.1,
                                   target_end_gap_score=0.0, query_end_gap_score=0.0, print_alignments=False)
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
                    if target_aa != query_aa and target_aa!='-':
                        mut = query_aa+str(res_idx)+target_aa
                        mut_list.append(mut)
                        actual_wt_aa_list.append(query_seq[res_idx-1]+str(res_idx))
            # print(target_al)
            # print(query_al)

            seq_desc = seq_desc.replace('Sequence','Seq')
            patent_id = seq_desc[seq_desc.find('patent ')+7:].replace(' ','')
            print(mut_list)
            res_parsed.append({
                'query_seq': seq_name_base,
                'accession': seq_name,
                'description': seq_desc,
                'patent_id': patent_id,
                'title': title,
                'journal': journal,
                'date': date,
                'seq_len': seq_len,
                'max_score': max_score,
                'total_score': total_score,
                'query_cover': query_cover,
                'percent_identity': percent_identity,
                'mutations': ', '.join(mut_list),
                'ali_start': ali_start,
                'ali_end': ali_end,
                'sequence': target_seq,
                'sequence_aligned': target_seq_aligned,
                'authors':authors,
                'abstract': '',
                'claims': ''
                })

        if len(res_parsed)>0:
            res = pd.DataFrame(res_parsed)
            res = res.sort_values(by=['patent_id', 'percent_identity'], ascending=[True, False], ignore_index=True)

            # update with scraped details of patent
            if scrape_patent_details:
                patent_id_list = list(set(res['patent_id'].tolist()))
                print('patent_id_list:', patent_id_list)
                for patent_id in patent_id_list:
                    data = extract_patent_data(patent_id)
                    if data:
                        # get index of first occurrence of patent_id
                        idx = int(res.loc[res['patent_id']==patent_id].index[0])
                        res.loc[idx, 'abstract'] = data['abstract']
                        res.loc[idx, 'claims'] = data['claims']
            # save results for given query sequence
            res.to_csv(data_folder+patents_subfolder+target_seq_fname.replace('.fasta','.csv'))
            print(res)

            # aggregate with overall results
            if query_seq_count==0: res_all = res.copy()
            else: res_all = pd.concat([res_all, res], ignore_index=True)
            query_seq_count += 1

    print(res_all)
    res_all.to_csv(data_folder+patents_subfolder+overall_results_fname_prefix+'_patent_hits_all.csv')

###################
# ANALYSE RESULTS #
###################
if analyse_patent_search_results:
    patent_analysis = []
    df = pd.read_csv(data_folder+patents_subfolder+overall_results_fname_prefix+'_patent_hits_all.csv')
    for k, seq_name_base in enumerate(seq_name_base_list):
        print(seq_name_base)
        df_seqbase = df[df['query_seq']==seq_name_base]
        num_hits = len(df_seqbase)
        num_patents = len(list(set(df_seqbase['title'].tolist())))
        num_hits_nomut = len(df_seqbase[df_seqbase['mutations'].isnull()])
        unique_mutants = list(set([mut for mut in df_seqbase['mutations'] if isinstance(mut,str)]))
        all_mutations = []
        all_residues = []
        mut_dict = {}
        mut_count_dict = {}
        # get unique mutations, mutated positions
        for mutstr in unique_mutants:
            muts = mutstr.split(', ')
            for mut in muts:
                if mut not in all_mutations:
                    all_mutations.append(mut)
                    WT_aa, res, MT_aa = split_mutation(mut, aa_letter_representation=True)
                    WT_res = WT_aa + str(res)
                    if res not in all_residues:
                        all_residues.append(res)
                        mut_dict[res] = (WT_res, [])
                    mut_dict[res][1].append(MT_aa)
        unique_mutants.sort()
        all_mutations.sort()
        all_residues.sort()
        mut_dict_str = []
        for res in all_residues:
            res_row = f'{mut_dict[res][0]}:{','.join(sort_list(mut_dict[res][1]))}'
            mut_dict_str.append(res_row)
        print(f'Total # of entries: {num_hits}')
        print(f'Total # unique patents: {num_patents}')
        print(f'# of entries with no mutations: {num_hits_nomut}')
        print(f'# of unique mutants: {len(unique_mutants)}')
        print(unique_mutants)
        print(f'# of unique mutations: {len(all_mutations)}')
        print(all_mutations)
        print(f'# of positions mutated: {len(all_residues)}')
        print(all_residues)
        print(mut_dict)
        print()

        patent_analysis_seqbase = {
            'query_seq': seq_name_base,
            '# patents': num_patents,
            '# seq': num_hits,
            '# seq (no mut)': num_hits_nomut,
            '# mutants': len(unique_mutants),
            '# mutations': len(all_mutations),
            '# residues': len(all_residues),
            'mutations': '; '.join(mut_dict_str)
        }
        patent_analysis.append(patent_analysis_seqbase)
    patent_analysis = pd.DataFrame(patent_analysis)
    patent_analysis.to_csv(data_folder+patents_subfolder+overall_results_fname_prefix+'_patent_analysis_all.csv')
    print(patent_analysis)



