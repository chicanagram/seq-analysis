import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
pd.set_option('display.max_columns', None)
from variables import address_dict, subfolders, aaList, aa2idx
from utils import fetch_sequences_from_fasta

def one_hot_encode_sequence(sequence, max_length=None):
    """
    Converts a single protein sequence into a one-hot encoded matrix.
    Args:
        sequence (str): Protein sequence.
        max_length (int, optional): Fixed length for padding/truncation.
    Returns: np.ndarray: One-hot encoded array of shape (max_length or len(sequence), 20).
    """
    sequence = sequence.upper()  # Ensure uppercase amino acids
    seq_len = len(sequence)
    if max_length is None:
        max_length = seq_len  # Use original sequence length if not specified
    # Initialize a zero matrix (max_length x 20)
    one_hot = np.zeros((max_length, len(aaList)), dtype=np.float32)
    # Fill the matrix with one-hot encoding
    for i, aa in enumerate(sequence[:max_length]):  # Truncate if needed
        if aa in aa2idx:  # Ignore unknown characters
            one_hot[i, aa2idx[aa]] = 1
    return one_hot

def one_hot_encode_sequences(sequences, max_length=None):
    """
    Converts a list of protein sequences into one-hot encoded matrices.
    Args:
        sequences (list of str): List of protein sequences.
        max_length (int, optional): Fixed length for all sequences.
    Returns: np.ndarray: Array of shape (num_sequences, max_length, 20).
    """
    if max_length is None:
        max_length = max(len(seq) for seq in sequences)  # Use max sequence length
    encoded_seqs = [one_hot_encode_sequence(seq, max_length) for seq in sequences]
    return np.array(encoded_seqs, dtype=np.float32)

def main():
    data_folder = address_dict['PIPS2'] # address_dict['PIPS'] # address_dict['PON-Sol2']
    data_subfolder = 'ET096_mutagenesis_Lib3' # 'ET096_mutagenesis_Purified_2025-12-05' # 'GOh1052_mutagenesis' # '''ET096_mutagenesis_Round3-batch1' # 'ET096_mutagenesis_purified_activity' #
    input_fname = 'ET096_mutagenesis_Lib3.csv' # 'ET096_mutagenesis_Purified_2025-12-05.csv' # 'GOh1052_mutagenesis.csv' # 'ET096_mutagenesis_Round3-batch1.csv' # 'ET096_mutagenesis_purified_activity.csv'
    output_fname = 'ohe.npz' # 'ohe_rev.npz'
    ohe_dir = data_folder + subfolders['ohe'] + data_subfolder + '/'
    df = pd.read_csv(data_folder + subfolders['expdata'] + data_subfolder + '/' + input_fname)
    sequence_base_list = df['sequence_base'].tolist()
    sequence_base = sequence_base_list[0]
    sequence_list = df['sequence'].tolist()
    seq_base_ohe = one_hot_encode_sequences(sequence_base_list)
    seq_ohe = one_hot_encode_sequences(sequence_list)
    ohe_3D = np.concatenate((seq_base_ohe,seq_ohe), axis=1)
    ohe = ohe_3D.reshape(ohe_3D.shape[0], ohe_3D.shape[1]*ohe_3D.shape[2])
    W = seq_ohe.shape[1]
    cols = [f'{i+1}_{aa}_{WT_or_MT}' for WT_or_MT in ['WT','MT'] for i,wt_aa in enumerate(sequence_base) for aa in aaList]
    colsMT = [f'{wt_aa}{i+1}{aa}' for i,wt_aa in enumerate(sequence_base) for aa in aaList]
    ohe_df = pd.DataFrame(ohe, columns=cols, index=df['name'].tolist())
    print(colsMT)

    # Convert to sparse format
    # WT and MT
    ohe_sparse = sp.csr_matrix(ohe)
    sp.save_npz(ohe_dir + output_fname, ohe_sparse)
    print('OHE shape:', ohe_sparse.shape)
    # MT only
    oheMT_sparse = sp.csr_matrix(ohe[:, int(len(cols)/2):])
    print('OHE(MT) shape:', oheMT_sparse.shape)
    sp.save_npz(ohe_dir + output_fname.replace('.npz','MT.npz'), oheMT_sparse)


if __name__ == "__main__":
    main()