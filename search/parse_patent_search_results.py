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
from utils import sort_list, split_mutation, fetch_sequences_from_fasta, run_msa, write_sequence_to_fasta, split_fasta, deduplicate_fasta_sequences, get_sequences_in_folder, run_pairwise_alignment, get_mutagenesis_sequences, parse_genbank_to_df
from scrape_sequence_patents import extract_patent_data
data_folder = address_dict['ECOHARVEST']
sequence_subfolder = subfolders['sequences']
msa_subfolder = subfolders['msa']
patents_subfolder = subfolders['patents']
seq_name_base_list = ['CALB', 'CALA', 'RML', 'TLL', 'BSL', 'UML', 'TCL'] # ['MmCAR'] #
overall_results_fname_prefix = 'Lipase' # 'CAR' #
parse_patent_search_results = False
scrape_patent_details = False
analyse_patent_search_results = True
filter_hits = {'percent_identity':50, 'query_cover':30, 'seq_len':50}


if parse_patent_search_results:
    query_seq_count = 0
    for k, seq_name_base in enumerate(seq_name_base_list):
        print(seq_name_base)
        query_seq_fname = seq_name_base + '.fasta'
        target_seq_fname = seq_name_base+'_patent_hits.fasta'
        blastp_description = seq_name_base+'_patent_hits_Descriptions.csv'
        db_fname = seq_name_base+'_patent_hits_GenBank.txt'

        ###########################
        # get pairwise alignments #
        ###########################
        # get query seq
        query_seq = str([record for record in SeqIO.parse(f'{data_folder}{sequence_subfolder}{query_seq_fname}', "fasta")][0].seq)

        # get sequences
        seqs, seq_names, seq_descriptions = fetch_sequences_from_fasta(data_folder+patents_subfolder+target_seq_fname)

        # write temporary fasta with query seq at the top
        temp_seq_fpath = write_sequence_to_fasta([query_seq]+seqs, [seq_name_base]+seq_names, seq_name_base+'_patent_hits_wquery', data_folder+sequence_subfolder)

        # perform alignment for all sequences
        fname = os.path.basename(temp_seq_fpath)
        run_msa(fname, fname, method='mafft', seq_dir=data_folder+sequence_subfolder, msa_dir=data_folder+msa_subfolder)
        os.remove(temp_seq_fpath)

        # get alignment for each sequence
        seqs_ali_wquery, seq_names_wquery, _ = fetch_sequences_from_fasta(data_folder + msa_subfolder + fname)

        # get residue numbering for each sequence (excluding dashes)
        seqs_pos_wquery = []
        for seq in seqs_ali_wquery:
            seq_pos_list = []
            count = 0
            for aa in seq:
                if aa != '-':
                    count+=1
                    seq_pos_list.append(count)
                else:
                    seq_pos_list.append('-')
            seqs_pos_wquery.append(seq_pos_list)

        # parse aligned sequences and positions, split into query & target sequences
        query_al = seqs_ali_wquery[0]
        query_al_pos = seqs_pos_wquery[0]
        target_seqs_aligned = seqs_ali_wquery[1:]
        target_seqs_aligned_pos = seqs_pos_wquery[1:]

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
        # print(db)
        print('# of seqs:', len(seqs))

        res_parsed = []
        for i, (target_seq, target_al, target_al_pos, seq_name, seq_desc) in enumerate(zip(seqs, target_seqs_aligned, target_seqs_aligned_pos, seq_names, seq_descriptions)):
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

            query_al_trimmed = ''
            target_al_trimmed = ''
            insertion_list = []
            deletion_list = []
            mutation_list = []
            latest_nongap_query_pos = 0
            latest_nongap_target_pos = 0
            for query_pos, query_aa, target_pos, target_aa in zip(query_al_pos, list(query_al), target_al_pos, list(target_al)):
                # update latest_nongap pos
                if query_pos != '-':
                    latest_nongap_query_pos = query_pos
                if target_pos != '-':
                    latest_nongap_target_pos = target_pos
                # handle gaps '-'
                if query_aa != '-' and target_aa != '-':
                    query_al_trimmed += query_aa
                    target_al_trimmed += target_aa
                # get insertions
                if query_aa == '-' and target_aa != '-':
                    insertion = f'{latest_nongap_query_pos}>{target_aa}'
                    insertion_list.append(insertion)
                # get deletions
                elif query_aa != '-' and target_aa == '-':
                    deletion = f'{query_aa}{query_pos}'
                    deletion_list.append(deletion)
                # get mutations
                elif query_aa != '-' and target_aa != '-' and query_aa!=target_aa:
                    mutation = f'{query_aa}{query_pos}{target_aa}'
                    mutation_list.append(mutation)
            print('Q:', query_al_trimmed)
            print('T:', target_al_trimmed)

            seq_desc = seq_desc.replace('Sequence','Seq')
            patent_id = seq_desc[seq_desc.find('patent ')+7:].replace(' ','')
            res_parsed.append({
                'query_seq': seq_name_base,
                'accession': seq_name,
                'description': seq_desc,
                'patent_id': patent_id,
                'title': title,
                'abstract': '',
                'claims': '',
                'journal': journal,
                'authors': authors,
                'date': date,
                'seq_len': seq_len,
                'max_score': max_score,
                'total_score': total_score,
                'query_cover': query_cover,
                'percent_identity': percent_identity,
                'num_mutations': len(mutation_list),
                'num_insertions': len(insertion_list),
                'num_deletions': len(deletion_list),
                'mutations': ', '.join(mutation_list),
                'insertions': ', '.join(insertion_list),
                'deletions': ', '.join(deletion_list),
                'sequence': target_seq,
                'sequence_aligned': target_al_trimmed,
                'query_sequence_aligned': query_al_trimmed,
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
    # load results
    df = pd.read_csv(data_folder+patents_subfolder+overall_results_fname_prefix+'_patent_hits_all.csv', index_col=0)
    print(df)
    # filter results on sequence identity, length, query cover conditions
    df_filt = df.copy()
    for col_to_filter, thres_to_filter in filter_hits.items():
        df_filt = df_filt[df_filt[col_to_filter]>thres_to_filter]
    # filter by num_mutations
    df_filt = df_filt[df_filt['num_mutations']<50]
    # remove duplicates sequences
    df_filt = df_filt.drop_duplicates(subset=['sequence'])

    # iterate through individual seq_name_base in
    seq_name_base_list_filt = [seq_name_base for seq_name_base in seq_name_base_list if seq_name_base in list(set(df_filt['query_seq'].tolist()))]
    for k, seq_name_base in enumerate(seq_name_base_list_filt):
        df_seqbase_orig = df[df['query_seq']==seq_name_base]
        df_seqbase = df_filt[df_filt['query_seq']==seq_name_base]
        print(f'[{seq_name_base}] Pre-filter: n={len(df_seqbase_orig)}; Post-filter: n={len(df_seqbase)}')

        # update scraped patent info, if filtered out
        patent_id_list = set(list(df_seqbase['patent_id'].tolist()))
        for patent_id in patent_id_list:
            # get index of first row
            df_patentid_match = df_seqbase_orig.loc[df_seqbase_orig['patent_id']==patent_id].iloc[0]
            abstract = df_patentid_match['abstract']
            claims = df_patentid_match['claims']
            df_seqbase.loc[df_seqbase['patent_id']==patent_id, 'abstract'] = abstract
            df_filt.loc[df_filt['patent_id']==patent_id, 'abstract'] = abstract
            df_seqbase.loc[df_seqbase['patent_id']==patent_id, 'claims'] = claims
            df_filt.loc[df_filt['patent_id']==patent_id, 'claims'] = claims
        # save filtered hits for seqbase
        df_seqbase.sort_values(by='percent_identity', ascending=False).reset_index(drop=True).to_csv(data_folder+patents_subfolder+seq_name_base+'_patent_hits_filtered.csv')

        # perform analysis
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
            res_row = f'{mut_dict[res][0]}:{".".join(sort_list(mut_dict[res][1]))}'
            mut_dict_str.append(res_row)
        print(f'Total # of entries: {num_hits}')
        print(f'Total # unique patents: {num_patents}')
        print(f'# of entries with no mutations: {num_hits_nomut}')
        print(f'# of positions mutated: {len(all_residues)}')
        print(f'# of unique mutants: {len(unique_mutants)}')
        # print(unique_mutants)
        print(f'# of unique mutations: {len(all_mutations)}')
        # print(all_mutations)
        print('all_residues:', len(all_residues), all_residues)
        print('mut_dict:', len(mut_dict), mut_dict)
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

    # save overall filtered results
    df_filt = df_filt.sort_values(by=['query_seq', 'percent_identity'], ascending=[True,False])
    df_filt.reset_index(drop=True).to_csv(data_folder + patents_subfolder + overall_results_fname_prefix+'_patent_hits_all_filtered.csv')

    # save patent analysis csv
    patent_analysis = pd.DataFrame(patent_analysis)
    print(patent_analysis)
    patent_analysis = patent_analysis.sort_values(by='query_seq')
    patent_analysis.reset_index(drop=True).to_csv(data_folder+patents_subfolder+overall_results_fname_prefix+'_patent_analysis_all_filtered.csv')
