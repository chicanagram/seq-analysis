import pandas as pd
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta

def parse_hmmer_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir):

    cols_init = [
        'E-value (Full Seq)', 'Score (Full Seq)', 'Bias (Full Seq)', 'E-value (Best 1 Domain)', 'Score (Best 1 Domain)', 'Bias (Best 1 Domain)','#dom exp','N', 'Seq Name', 'Description',
        'HMM From', 'HMM To', 'Ali From', 'Ali To', 'Env From', 'Env To', 'Acc', 'Target Seq Ali', 'Query Seq Ali'
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
        df_rows[i] += [target_seq_ali, query_seq_ali]

    df = pd.DataFrame(df_rows, columns=cols_init)
    df = df[columns]
    print(df)
    df.to_csv(seqsearch_dir + output_fname)

    # get fasta of sequences
    seq_names = df['Seq Name'].tolist()
    sequences = df['Target Seq Ali'].tolist()
    sequences = [s.replace('U','X') if isinstance(s, str) else None for s in sequences]

    # append query sequence
    if isinstance(query_seq_input, str):
        query_seq, query_seq_name, _ = fetch_sequences_from_fasta(sequences_dir + query_seq_input)
    elif isinstance(query_seq_input, tuple):
        (query_seq, query_seq_name) = query_seq_input
    _ = write_sequence_to_fasta(query_seq+sequences, query_seq_name+seq_names, output_fasta_fname, sequences_dir)

def parse_blastp_output(input_fname, output_fname, output_fasta_fname, query_seq_input, sequences_dir, seqsearch_dir):
    with open(seqsearch_dir + input_fname, 'r') as f:
        data = f.read()
        entries = data[data.find('>'):].split('>')
        entries = entries[1:]
        print('# of entries:', len(entries))

    df_rows = []
    for i, entry in enumerate(entries):
        # get metadata
        metadata_txt = entry[:entry.find('Query')]
        # Extract metadata
        row = {}
        # Extract sequence name and description
        seqname_description = metadata_txt[:metadata_txt.find('Length')].replace('\n', '')
        row['Seq Name'] = seqname_description[:seqname_description.find(' ')]
        row['Description'] = seqname_description[seqname_description.find(' ')+1:]
        # Extract Score and Expect
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
        lines = entry.split('\n')
        seq_lines = [l.strip() for l in lines[7:]]
        target_seq_row_idxs = [l_idx for l_idx, l in enumerate(seq_lines) if l.find('Sbjct')>-1]
        query_seq_row_idxs = [l_idx-2 for l_idx in target_seq_row_idxs]
        target_seq_ali = ''.join([seq_lines[l_idx].strip().split()[2] for l_idx in target_seq_row_idxs]).upper()
        query_seq_ali = ''.join([seq_lines[l_idx].strip().split()[2] for l_idx in query_seq_row_idxs]).upper()
        row['Target Seq Ali'] = target_seq_ali
        row['Query Seq Ali'] = query_seq_ali
        df_rows.append(row)
    df = pd.DataFrame(df_rows)
    df = df.rename(columns={'Expect':'E-value'})
    df.to_csv(seqsearch_dir + output_fname)

    # get fasta of sequences
    seq_names = df['Seq Name'].tolist()
    sequences = df['Target Seq Ali'].tolist()
    sequences = [s.replace('U','X') if isinstance(s, str) else None for s in sequences]

    # append query sequence
    if isinstance(query_seq_input, str):
        query_seq, query_seq_name, _ = fetch_sequences_from_fasta(sequences_dir + query_seq_input)
    elif isinstance(query_seq_input, tuple):
        (query_seq, query_seq_name) = query_seq_input
        query_seq, query_seq_name = [query_seq], [query_seq_name]
    _ = write_sequence_to_fasta(query_seq+sequences, query_seq_name+seq_names, output_fasta_fname, sequences_dir)


if __name__=='__main__':
    # Input and output file paths
    data_folder = address_dict['ECOHARVEST']
    seqsearch_dir = data_folder + subfolders['seqsearch']
    sequences_dir = data_folder + subfolders['sequences']
    query_seq_fname = 'CALA.fasta'
    input_fname = 'CALA_phmmer_uniprot_trembl_incE1e-03_E1e-03'  # 'CALB_jackhmmer_Pfam-A_incE1e-05_E1e-03.out' #
    output_fname = input_fname.replace('.out', '.csv')
    output_fasta_fname = input_fname.replace('.out', '')
    # parse results
    if input_fname.find('hmmer')>-1:
        parse_hmmer_output(input_fname, output_fname, output_fasta_fname, query_seq_fname, sequences_dir, seqsearch_dir)
    elif input_fname.find('blastp')>-1:
        parse_blastp_output(input_fname, output_fname, output_fasta_fname, query_seq_fname, sequences_dir, seqsearch_dir)



