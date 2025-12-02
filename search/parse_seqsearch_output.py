import os
import pandas as pd
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta
from Bio import SeqIO

def get_seq_organism(description, search_type='qblast'):
    # get organism and strain
    organism_strain = description[description.find('[')+1:description.find(']')]
    # get organism only
    organism = ' '.join(organism_strain.split(' ')[:2])
    strain = organism_strain.replace(organism, '')
    strain = strain.replace('(','').replace(')','').strip()
    return organism, strain

def parse_hmmer_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir, max_target_seqs=None, parse_description=False):

    cols_init = [
        'E-value (Full Seq)', 'Score (Full Seq)', 'Bias (Full Seq)', 'E-value (Best 1 Domain)', 'Score (Best 1 Domain)', 'Bias (Best 1 Domain)','#dom exp','N', 'Seq Name', 'Description',
        'HMM From', 'HMM To', 'Ali From', 'Ali To', 'Env From', 'Env To', 'Acc', 'Target Seq Ali', 'Query Seq Ali', 'Target Seq'
        ]
    columns = ['Seq Name', 'Description', 'E-value (Full Seq)', 'Score (Full Seq)', 'Bias (Full Seq)',
            'E-value (Best 1 Domain)', 'Score (Best 1 Domain)', 'Bias (Best 1 Domain)',
            'HMM From', 'HMM To', 'Ali From', 'Ali To', 'Target Seq Ali', 'Query Seq Ali']

    ### UPDATE CODE TO HANDLE JACKHMMER OUTPUTS WITH MULTIPLE BLOCKS OF SEARCHES ###

    with open(seqsearch_dir + input_fname, 'r') as f:
        data = f.read()
        table = data[data.find('Scores for complete sequences'):data.find('>>')]
        entries = data[data.find('>>'):].split('>>')
        entries = entries[1:]
        print('# of entries:', len(entries))

    # parse table portion
    rows = table.split('\n')
    print('# of rows:', len(rows))
    for i, row in enumerate(rows):
        print(i, row)
    print()
    data_rowstart = index = next(i for i, row in enumerate(rows) if 'E-value' in row) + 2
    data_rows = rows[data_rowstart:-3]
    df_rows = []
    for i, row in enumerate(data_rows):
        print(i, row)
        if row.find('------')==-1 and len(row)>0:
            if row[0]=='+':
                row = row[1:]
            rowsplit = row.strip().split()
            rowsplit = rowsplit[:9] + [''.join(rowsplit[9:])]
            rowsplit = [float(v) if i in range(8) else v for i,v in enumerate(rowsplit)]
            df_rows.append(rowsplit[:-1])
            print(rowsplit)

    for i, entry in enumerate(entries):
        print(i, df_rows[i])
        lines = entry.split('\n')
        # get description
        headerline = lines[0].strip()
        seq_name = headerline[:headerline.find(' ')]
        print(seq_name)
        seq_desc = headerline[headerline.find(' ')+1:]
        print(seq_desc)
        df_rows[i].append(seq_desc)
        # get stats
        statsrow = lines[3]
        statsrow = statsrow.strip().split()
        statsrow = [float(v) for i,v in enumerate(statsrow) if i in [6,7,9,10,12,13,15]]
        statsrow[:-1] = [int(v) for v in statsrow[:-1]]
        df_rows[i] += statsrow
        # get aligned sequences
        seq_lines = [l.strip() for l in lines[7:]]
        target_seq_row_idxs = [l_idx for l_idx, l in enumerate(seq_lines) if l.find(seq_name)>-1]
        query_seq_row_idxs = [l_idx-2 for l_idx in target_seq_row_idxs]
        target_seq_ali = ''.join([seq_lines[l_idx].strip().split()[2] for l_idx in target_seq_row_idxs]).upper()
        query_seq_ali = ''.join([seq_lines[l_idx].strip().split()[2] for l_idx in query_seq_row_idxs]).upper()
        target_seq = target_seq_ali.replace('-','')
        df_rows[i] += [target_seq_ali, query_seq_ali, target_seq]
    df = pd.DataFrame(df_rows, columns=cols_init)
    df = df[columns]

    # filter to max target # of sequences
    if max_target_seqs is not None and max_target_seqs<len(df):
        df = df.iloc[:max_target_seqs,:]
        print(f'Filtered top {max_target_seqs} hits.')
    df.to_csv(seqsearch_dir + output_fname)

    # get fasta of sequences
    seq_names = df['Seq Name'].tolist()
    sequences = df['Target Seq'].tolist()
    sequences = [s.replace('U','X') if isinstance(s, str) else None for s in sequences]

    # append query sequence
    if isinstance(query_seq_input, str):
        query_seq, query_seq_name, _ = fetch_sequences_from_fasta(sequences_dir + query_seq_input)
    elif isinstance(query_seq_input, tuple):
        (query_seq, query_seq_name) = query_seq_input
    _ = write_sequence_to_fasta(query_seq+sequences, query_seq_name+seq_names, output_fasta_fname, sequences_dir)
    return df

def parse_blastp_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir, max_target_seqs=None, parse_description=False):

    # get query sequence
    if isinstance(query_seq_input, str):
        query_seq, query_seq_name, _ = fetch_sequences_from_fasta(sequences_dir + query_seq_input)
    elif isinstance(query_seq_input, tuple):
        (query_seq, query_seq_name) = query_seq_input
        query_seq, query_seq_name = [query_seq], [query_seq_name]

    # parse out file
    with open(seqsearch_dir + input_fname, 'r') as f:
        print(seqsearch_dir + input_fname)
        data = f.read()
        data = data[data.find('Sequences producing significant alignments'):]
        data = data.replace('->', 'to').replace('[[', '[')
        entries = data[data.find('>'):].split('>')
        entries = entries[1:]
        print('# of entries:', len(entries))

    df_rows = []
    for i, entry in enumerate(entries):
        # get metadata
        metadata_txt = entry[:entry.find('Query')]
        seqname_description = metadata_txt[:metadata_txt.find('Length')].replace('\n', '')
        seqname = seqname_description[:seqname_description.find(' ')]
        description = seqname_description[seqname_description.find(' ')+1:]

        # Extract metadata
        row = {}
        if parse_description:
            organism, strain = get_seq_organism(description, search_type='qblast')
            row['Target Organism'] = organism
            row['Target Strain'] = strain
        row['Query Seq Name'] = query_seq_name[0]
        row['Seq Name'] = seqname
        row['Description'] = description


        # Update Score and Expect
        meta_lines = metadata_txt.split("\n")
        score_line = next(line for line in meta_lines if "Score =" in line)
        score_parts = score_line.split("Score =")[1].split(",")
        row['Score'] = float(score_parts[0].split()[0])
        row['Expect'] = float(score_parts[1].split("=")[1].strip())

        # Extract Identities, Positives, and Gaps
        alignment_line = next(line for line in meta_lines if "Identities =" in line)
        alignment_parts = alignment_line.split(",")
        row['Identities'] = float(alignment_parts[0].split("(")[1].split("%")[0])
        row['Positives'] = float(alignment_parts[1].split("(")[1].split("%")[0])
        row['Gaps'] = float(alignment_parts[2].split("(")[1].split("%")[0])
        # print(i, row)

        # get aligned sequences
        entry = entry[entry.find('Query'):]
        seq_lines = [l.strip() for l in entry.split('\n')]
        target_seq_row_idxs = [l_idx for l_idx, l in enumerate(seq_lines) if l.find('Sbjct')>-1]
        target_seq_lines = [seq_lines[l_idx].strip().split()[2] for l_idx in target_seq_row_idxs]
        target_seq_ali = ''.join(target_seq_lines).upper()
        query_seq_row_idxs = [l_idx-2 for l_idx in target_seq_row_idxs]
        query_seq_lines = [seq_lines[l_idx].strip().split()[2] for l_idx in query_seq_row_idxs]
        query_seq_ali = ''.join(query_seq_lines).upper()
        row['Target Seq Ali'] = target_seq_ali
        row['Query Seq Ali'] = query_seq_ali
        row['Target Seq'] = target_seq_ali.replace('-','')
        df_rows.append(row)
    df = pd.DataFrame(df_rows)
    df = df.rename(columns={'Expect':'E-value'})

    # filter to max target # of sequences
    if max_target_seqs is not None and max_target_seqs<len(df):
        df = df.iloc[:max_target_seqs,:]
        print(f'Filtered top {max_target_seqs} hits.')
    df.to_csv(seqsearch_dir + output_fname)

    # get fasta of sequences
    seq_names = df['Seq Name'].tolist()
    sequences = df['Target Seq'].tolist()
    sequences = [s.replace('U','X') if isinstance(s, str) else None for s in sequences]

    # write to fasta
    _ = write_sequence_to_fasta(query_seq+sequences, query_seq_name+seq_names, output_fasta_fname, sequences_dir)

    return df

def match_query_seq(input_fname, query_seqname_list, query_seq_list):
    match_bool = [query_seqname in input_fname for query_seqname in query_seqname_list]
    match_index = match_bool.index(True)
    query_seqname_match = query_seqname_list[match_index]
    query_seq_match = query_seq_list[match_index]
    return (query_seq_match, query_seqname_match)

if __name__=='__main__':
    # Input and output file paths
    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'sidestream_cocktail' # 'lipases' #
    seqsearch_dir = data_folder + subfolders['seqsearch'] + data_subfolder + '/'
    sequences_dir = data_folder + subfolders['sequences'] + data_subfolder + '/'
    query_seq_fname = 'queryseqs_enzcocktail.fasta' #
    input_fname_list = [f for f in os.listdir(seqsearch_dir) if f.find('_nr.out')>-1]
    parse_description = True

    if isinstance(input_fname_list, str):
        input_fname_list = [input_fname_list]
    # get query seqs
    query_seqs, query_seqnames, _ = fetch_sequences_from_fasta(sequences_dir + query_seq_fname)
    print(len(query_seqnames), query_seqnames)

    # iterate through inputs
    for i, input_fname in enumerate(input_fname_list):
        output_fname = input_fname.replace('.out', '.csv')
        output_fasta_fname = input_fname.replace('.out', '')
        query_seq_input = match_query_seq(input_fname, query_seqnames, query_seqs)
        print('Processing', i, input_fname, query_seq_input[1])

        # parse results
        if input_fname.find('hmmer')>-1:
            df = parse_hmmer_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir, parse_description=parse_description)
        elif input_fname.find('blastp')>-1:
            df = parse_blastp_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir, parse_description=parse_description)

        # append to full results
        if len(input_fname_list)>1:
            if i==0:
                df_all = df.copy()
            else:
                df_all = pd.concat([df_all, df], axis=0, ignore_index=True)
            df_all.to_csv(seqsearch_dir + data_subfolder + '.csv')


