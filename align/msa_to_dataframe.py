import os
import pandas as pd
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta

def convert_msa_to_dataframe(seqs, seq_names, pos_offset={}):
    """
    Build alignment dataframe
    """
    aln_len = len(seqs[0])
    if any(len(seq) != aln_len for seq in seqs):
        raise ValueError("Alignment length mismatch among sequences")

    res_counters = {seq_name: 0 if seq_name not in pos_offset else pos_offset[seq_name] for seq_name in seq_names}  # start from 0; increment before recording
    rows = []
    for i in range(aln_len):
        row = {"index": i+1}
        for seq, seq_name in zip(seqs, seq_names):
            aa = seq[i]
            if aa != "-":
                res_counters[seq_name] += 1
                row[f"{seq_name}_res_num"] = res_counters[seq_name]
            else:
                row[f"{seq_name}_res_num"] = ""
            row[f"{seq_name}_res_aa"] = aa
        rows.append(row)

    df_cols = ['index']
    for seq_name in seq_names:
        df_cols += [f'{seq_name}_res_num', f'{seq_name}_res_aa']
    df = pd.DataFrame(rows, columns=df_cols)
    return df

def convert_dataframe_to_msa(df):
    cols = df.columns.tolist()
    seq_names = [c.replace('_res_num', '') for c in cols if c.find('_res_num')>-1]
    print(seq_names)
    seqs = []
    for seq_name in seq_names:
        seq = ''.join(df[f'{seq_name}_res_aa'].tolist())
        seqs.append(seq)
    return seqs, seq_names


if __name__ == '__main__':
    os.chdir('../')
    print('CWD:', os.getcwd())
    data_folder = address_dict['PIPS2']
    data_fbase = ''
    seq_dir = data_folder + subfolders['sequences'] + data_fbase + '/'
    msa_dir = data_folder + subfolders['msa']
    msa_fname = 'UPOs_peroxygenation.fasta'
    msa_fpath = f'{msa_dir}{msa_fname}'
    seqs, seq_names, _ = fetch_sequences_from_fasta(msa_fpath)

    df = convert_msa_to_dataframe(seqs, seq_names)
    df.to_csv(msa_fpath.replace('.fasta', '.csv'))