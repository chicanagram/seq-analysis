import re
import csv
import pandas as pd
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta

# Input and output file paths
data_folder = address_dict['ECOHARVEST']
input_fname = 'CALB_phmmer_uniprot_trembl_incE=1e-03_E=1e-03.out' # 'CALB_jackhmmer_Pfam-A_incE=1e-05_E=1e-03.out' #
output_fname = input_fname.replace('.out', '.csv')
output_fasta_fname = input_fname.replace('.out', '')
input_fpath = data_folder + subfolders['seqsearch'] + input_fname
output_fpath = data_folder + subfolders['seqsearch'] + output_fname
query_seq_fpath = data_folder + subfolders['sequences'] + 'CALB.fasta'

cols_init = [
    'E-value (Full Seq)', 'Score (Full Seq)', 'Bias (Full Seq)', 'E-value (Best 1 Domain)', 'Score (Best 1 Domain)', 'Bias (Best 1 Domain)','#dom exp','N', 'Seq Name', 'Description',
    'HMM From', 'HMM To', 'Ali From', 'Ali To', 'Env From', 'Env To', 'Acc', 'Target Seq Ali', 'Query Seq Ali'
    ]
columns = ['Seq Name', 'Description', 'E-value (Full Seq)', 'Score (Full Seq)', 'Bias (Full Seq)',
        'E-value (Best 1 Domain)', 'Score (Best 1 Domain)', 'Bias (Best 1 Domain)',
        'HMM From', 'HMM To', 'Ali From', 'Ali To', 'Target Seq Ali', 'Query Seq Ali']

### UPDATE CODE TO HANDLE JACKHMMER OUTPUTS WITH MULTIPLE BLOCKS OF SEARCHES ###

with open(input_fpath, 'r') as f:
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
df.to_csv(output_fpath)

# get fasta of sequences
seq_names = df['Seq Name'].tolist()
sequences = df['Target Seq Ali'].tolist()
# append query sequence
query_seq, query_seq_name, _ = fetch_sequences_from_fasta(query_seq_fpath)
_ = write_sequence_to_fasta(query_seq+sequences, query_seq_name+seq_names, output_fasta_fname, data_folder + subfolders['sequences'])
