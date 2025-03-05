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
    data_folder = address_dict['PON-Sol2']
    input_fname = 'ponsol2_all_rev.csv'
    output_fname = 'ohe_rev.npz'
    ohe_subfolder = subfolders['ohe']
    df = pd.read_csv(data_folder + input_fname)
    sequence_base_list = df['sequence_base'].tolist()
    sequence_list = df['sequence'].tolist()
    seq_base_ohe = one_hot_encode_sequences(sequence_base_list)
    seq_ohe = one_hot_encode_sequences(sequence_list)
    ohe_3D = np.concatenate((seq_base_ohe,seq_ohe), axis=1)
    ohe = ohe_3D.reshape(ohe_3D.shape[0], ohe_3D.shape[1]*ohe_3D.shape[2])
    # W = seq_ohe.shape[1]
    # cols = [f'{i}_{aa}_{WT_or_MT}' for WT_or_MT in ['WT','MT'] for i in range(1,W+1) for aa in aaList]
    # ohe_df = pd.DataFrame(ohe, columns=cols, index=df['name'].tolist())

    # Convert to sparse format
    ohe_sparse = sp.csr_matrix(ohe)
    # Save to file
    sp.save_npz(data_folder + ohe_subfolder + output_fname, ohe_sparse)

if __name__ == "__main__":
    main()