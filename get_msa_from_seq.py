#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import os
import numpy as np
from collections import Counter
import pandas as pd
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt
from variables import address_dict, subfolders

def run_msa(seq_fname, msa_fname, method, data_folder, sequence_subfolder, msa_subfolder):
    import subprocess
    in_file = f'{data_folder}{sequence_subfolder}{seq_fname}'
    out_file = f'{data_folder}{msa_subfolder}{msa_fname}'

    if method=='mafft':
        mafft  = './msa/mafft-mac/mafft.bat'
        mafft_command = f'{mafft} {in_file} > {out_file}'
        print(mafft_command)
        subprocess.call(mafft_command, shell=True)

    # not working
    elif method=='clustalw2':
        from Bio.Align.Applications import ClustalwCommandline
        clustalw2_cline = ClustalwCommandline("clustalw2", infile=in_file)
        print(clustalw2_cline)

    elif method == 'clustalo':
        import subprocess
        seq_fasta_fpath = os.path.abspath(f'{data_folder}{sequence_subfolder}{seq_fname}')
        msa_fpath = os.path.abspath(f'{data_folder}{msa_subfolder}{msa_fname}')
        clustalo_fpath = os.path.abspath(f'./{msa_subfolder}/clustalo/clustalo.exe')
        cmd = f'{clustalo_fpath} -i {seq_fasta_fpath} -o {msa_fpath} --force'
        msa_file = open(msa_fpath, 'w')
        subprocess.run(cmd, stdout=msa_file, encoding="utf8")
        msa_file.close()

    # not working
    elif method == 'muscle':
        muscle = './msa/muscle_macOS'
        subprocess.run([muscle, "-align", in_file, "-out", out_file])

def get_consensus_at_position(msa_array, position=0):
    """Calculate the consensus residue at a given position in an MSA."""
    # Extract the column at the specified position
    column = list(msa_array[:,position])
    column = [x for x in column if x!='-']
    # Count frequency of residues
    residue_counts = Counter(column)
    # Identify the most frequent residue
    consensus_residue = max(residue_counts, key=residue_counts.get)
    return consensus_residue

def get_consensus_sequence(msa_array):
    """Computes the consensus sequence for the entire MSA."""
    consensus = [get_consensus_at_position(msa_array, i) for i in range(msa_array.shape[1])]
    return consensus

def get_consensus_scores(msa_array):
    """Compute the consensus score for each position in an MSA."""
    num_sequences = msa_array.shape[0]
    sequence_length = msa_array.shape[1]
    consensus_scores = []
    for i in range(sequence_length):
        # Extract column at position i
        column = list(msa_array[:,i])
        # Count frequencies of each residue
        residue_counts = Counter(column)
        # Identify the most frequent residue
        most_common_residue, max_count = residue_counts.most_common(1)[0]
        # Compute consensus score (frequency of the most common residue)
        consensus_score = max_count / num_sequences
        consensus_scores.append(consensus_score)
    return consensus_scores

def plot_msa_seaborn(msa_fpath, color_scheme='Taylor', plot_msa_pos_range=None,
                     wrap_length=300, xtick_interval=25, ytick_interval=100,
                     show_seq_names=False, label_residues=None, filter_by_refseq=None, savefig=None):
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    aa_to_cmap_color_mapping = {
        'Clustal': {
            '-': 0,  # Gaps
            'A': 1, 'V': 1, 'L': 1, 'I': 1, 'M': 1, 'F': 1, 'W': 1, 'C': 1,  # Hydrophobic (Green)
            'K': 2, 'R': 2,  # Positive charge (Blue)
            'D': 3, 'E': 3,  # Negative charge (Red)
            'S': 4, 'T': 4, 'N': 4, 'Q': 4,  # Polar (Cyan)
            'Y': 5, 'H': 5,  # Aromatic (Magenta)
            'G': 6,  # Glycine (Orange)
            'P': 7,  # Proline (Yellow)
            'X': 8
        },
        'Taylor': {
            '-': 0,  # Gaps (Light Gray)
            'A': 1, 'V': 1, 'L': 1, 'I': 1, 'M': 1,  # Hydrophobic (Green)
            'F': 2, 'Y': 2, 'W': 2,  # Aromatic (Blue)
            'K': 3, 'R': 3, 'H': 3,  # Positive Charge (Red)
            'D': 4, 'E': 4,  # Negative Charge (Magenta)
            'S': 5, 'T': 5, 'N': 5, 'Q': 5,  # Polar Uncharged (Cyan)
            'C': 6,  # Cysteine (Yellow)
            'G': 7,  # Glycine (Orange)
            'P': 8,  # Proline (Brown)
            'X': 9  # Ambiguous/Unknown (Black)
        }
    }
    palette = {
        'Clustal': [
            "#d3d3d3",  # Gaps (Light Gray)
            "#32CD32",  # Hydrophobic (Green)
            "#0000FF",  # Positive (Blue)
            "#FF0000",  # Negative (Red)
            "#00FFFF",  # Polar (Cyan)
            "#FF00FF",  # Aromatic (Magenta)
            "#FFA500",  # Glycine (Orange)
            "#FFFF00",  # Proline (Yellow)
            "#000000"  # Ambiguous (Black)
        ],
        'Taylor': [
            "#D3D3D3",  # Gaps (Light Gray)
            "#33FF00",  # Hydrophobic (Green)
            "#0099FF",  # Aromatic (Blue)
            "#FF0000",  # Positive Charge (Red)
            "#CC00FF",  # Negative Charge (Magenta)
            "#00FFFF",  # Polar Uncharged (Cyan)
            "#FFFF00",  # Cysteine (Yellow)
            "#FF9900",  # Glycine (Orange)
            "#996633",  # Proline (Brown)
            "#000000"  # Ambiguous (Black)
        ]
    }

    # Load MSA using Biopython
    alignment = AlignIO.read(msa_fpath, "fasta")
    # get sequence names
    seq_names = [record.id for record in alignment]
    # Convert MSA to a NumPy array
    msa_array = np.array([list(record.seq) for record in alignment])

    # filter MSA positions by reference sequence
    if filter_by_refseq is not None:
        ref_seq = list(msa_array[filter_by_refseq, :])
        idx_notgap = [i for i, x in enumerate(ref_seq) if x != '-']
        print(len(idx_notgap), idx_notgap)
        msa_array = msa_array[:, idx_notgap]
        ref_seq = np.array(ref_seq)[idx_notgap]

    # label residues to annotate (i.e. consensus sequence or reference sequence)
    consensus_seq = np.array(get_consensus_sequence(msa_array))
    if label_residues is None or label_residues=='consensus':
        annotate_seq = consensus_seq
    elif label_residues == 'ref' and filter_by_refseq is not None:
        annotate_seq = ref_seq
        diff_seq = [consensus_aa if consensus_aa!=ref_aa else '' for consensus_aa, ref_aa in zip(consensus_seq, ref_seq)]

    # get filtered range of positions to plot
    if plot_msa_pos_range is None:
        start_pos_offset = 0
    else:
        [start_pos, end_pos] = plot_msa_pos_range
        if end_pos is None:
            end_pos = msa_array.shape[1]
        msa_array = msa_array[:,start_pos:end_pos]
        annotate_seq = annotate_seq[start_pos:end_pos]
        start_pos_offset = start_pos

    # get MSA dimensions
    msa_len = msa_array.shape[1]
    num_sequences = msa_array.shape[0]
    num_rows = int(np.ceil(msa_len / wrap_length))

    # Map residues to color codes
    colors = aa_to_cmap_color_mapping[color_scheme]
    msa_numeric = np.vectorize(colors.get)(msa_array)

    # Create custom color map
    cmap = ListedColormap(palette[color_scheme])

    # Plot heatmap of MSA using Seaborn
    fig, ax = plt.subplots(num_rows, 1, figsize=(25,20))

    # plot MSA row by row
    for row_idx in range(num_rows):
        if num_rows==1:
            ax_row = ax
            start_pos = 0
            end_pos = msa_len
            row_len = msa_len
            msa_numeric_row = msa_numeric
            annotate_seq_row = annotate_seq
        else:
            ax_row = ax[row_idx]
            start_pos = row_idx*wrap_length
            end_pos = min((row_idx+1)*wrap_length, msa_len)
            row_len = end_pos-start_pos
            msa_numeric_row = msa_numeric[:,start_pos:end_pos]
            annotate_seq_row = annotate_seq[start_pos:end_pos]

        # plot msa segment
        sns.heatmap(msa_numeric_row, ax=ax_row, cmap=cmap, cbar=False, xticklabels=False, yticklabels=False)
        # add vertical lines
        for x in range(row_len):
            ax_row.axvline(x=x, linewidth=0.2, color='k')
        # add tick labels
        ax_row.set_xticks(np.arange(0, row_len, xtick_interval))
        ax_row.set_xticklabels(np.arange(start_pos, end_pos, xtick_interval)+start_pos_offset, fontsize=10)
        # annotate sequence
        if label_residues is not None:
            for res_idx, res in enumerate(annotate_seq_row):
                ax_row.annotate(res, (res_idx,0), fontsize=8)

        if show_seq_names:
            ax_row.set_yticks(np.arange(0, num_sequences, ytick_interval)+0.5)
            ax_row.set_yticklabels(seq_names, fontsize=10)
        else:
            ax_row.set_yticks(np.arange(0, num_sequences, ytick_interval))
            ax_row.set_yticklabels(np.arange(0, num_sequences, ytick_interval), fontsize=10)

    # add labels and title
    plt.suptitle("MSA Visualization", fontsize=18)
    plt.xlabel("Position", fontsize=14)
    plt.savefig(savefig, bbox_inches='tight')
    plt.show()

def visualize_msa(msa_fpath, how='seaborn', color_scheme='Taylor', plot_msa_pos_range=None,
                  wrap_length=300, xtick_interval=25, ytick_interval=100,
                  show_seq_names=False, label_residues=None, filter_by_refseq=None, savefig=None):
    # get figure save name
    if savefig is None:
        savefig = msa_fpath.replace('.fasta','.png')
    else:
        savefig = f'{savefig}.png'

    # visualize MSA using PyMSAViz
    if how=='pymsaviz':
        from pymsaviz import MsaViz
        mv = MsaViz(msa_fpath, color_scheme=color_scheme, wrap_length=wrap_length, show_grid=True, show_seq_char=False, show_consensus=True, show_count=True)
        mv.set_plot_params(ticks_interval=50, x_unit_size=0.04, show_consensus_char=False)
        fig = mv.plotfig()
        mv.savefig(savefig)
        plt.show()

    # visualize MSA using Seaborn
    elif how=='seaborn':
        plot_msa_seaborn(msa_fpath, color_scheme, plot_msa_pos_range,
                         wrap_length, xtick_interval, ytick_interval,
                         show_seq_names, label_residues, filter_by_refseq, savefig)


def main():
    data_folder = address_dict['ECOHARVEST']
    sequence_subfolder = subfolders['sequences']
    msa_subfolder = subfolders['msa']
    seq_fname = 'CALB_phmmer_uniprot_trembl_incE=1e-03_E=1e-03.fasta' # 'lipases_initialMSA.fasta' # 'CALB_patent_hits.fasta' # 'UPO.fasta' # 'UPO_filtered.fasta' #
    msa_method = 'mafft' # 'clustalo' #
    msa_fname = seq_fname[:seq_fname.find('.fa')] + f'_{msa_method}' + '.fasta'
    get_msa = False
    plot_msa = 'seaborn' # 'pymsaviz' #
    plot_msa_pos_range = [150,250]
    wrap_length = 400 # 600
    filter_by_refseq = 0 # None
    label_residues = 'ref' # 'consensus' #
    xtick_interval = 10
    ytick_interval = 100
    if ytick_interval==1:
        show_seq_names = True
    else:
        show_seq_names = False
    savefig = None

    # get MSA
    if get_msa:
        run_msa(seq_fname, msa_fname, method=msa_method, data_folder=data_folder, sequence_subfolder=sequence_subfolder, msa_subfolder=msa_subfolder, savefig=savefig)

    # visualize MSA
    if plot_msa is not None:
        msa_fpath = f'{data_folder}{msa_subfolder}{msa_fname}'
        visualize_msa(msa_fpath, how=plot_msa, color_scheme='Taylor', plot_msa_pos_range=plot_msa_pos_range,
                      wrap_length=wrap_length, xtick_interval=xtick_interval, ytick_interval=ytick_interval,
                      show_seq_names=show_seq_names, label_residues=label_residues, filter_by_refseq=filter_by_refseq, savefig=savefig)

if __name__ == "__main__":
    main()