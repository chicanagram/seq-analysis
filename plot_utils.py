import numpy as np
import matplotlib.pyplot as plt
from collections import Counter

def get_consensus_at_position(msa_array, position=0):
    """Calculate the consensus residue at a given position in an MSA."""
    # Extract the column at the specified position
    column = list(msa_array[:,position])
    column = [x for x in column if x!='-']
    # Count frequency of residues
    residue_counts = Counter(column)
    # Identify the most frequent residue
    if len(residue_counts)>0:
        consensus_residue = max(residue_counts, key=residue_counts.get)
    else:
        consensus_residue = ''
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
                     show_seq_names=False, label_residues=None, show_all_sequences=False, fontsize=8, filter_by_refseq=None, savefig=None):
    import seaborn as sns
    from matplotlib.colors import ListedColormap, BoundaryNorm
    from Bio import AlignIO

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
    palettes = {
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
    # Create custom color map
    colors = aa_to_cmap_color_mapping[color_scheme]
    palette = palettes[color_scheme]
    cmap = ListedColormap(palette)
    norm = BoundaryNorm(np.arange(len(palette) + 1) - 0.5, len(palette))

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
    if label_residues == 'ref' and filter_by_refseq is not None:
        annotate_seq = ref_seq
        diff_seq = [consensus_aa if consensus_aa!=ref_aa else '' for consensus_aa, ref_aa in zip(consensus_seq, ref_seq)]
    else:
        annotate_seq = consensus_seq
        diff_seq = ['']*len(annotate_seq)

    # get filtered range of positions to plot
    if plot_msa_pos_range is None:
        start_pos_offset = 0
    else:
        [start_res, end_res] = plot_msa_pos_range
        if end_res is None:
            end_res = msa_array.shape[1]+1
        msa_array = msa_array[:,start_res-1:end_res-1]
        annotate_seq = annotate_seq[start_res-1:end_res-1]
        diff_seq = diff_seq[start_res-1:end_res-1]
        start_pos_offset = start_res-1

    # get MSA dimensions
    msa_len = msa_array.shape[1]
    num_sequences = msa_array.shape[0]
    num_rows = int(np.ceil(msa_len / wrap_length))

    # Map residues to color codes
    msa_numeric = np.vectorize(colors.get)(msa_array)
    pos_counter_dict = {seq_num:0 for seq_num in range(len(seq_names))}

    # Plot heatmap of MSA using Seaborn
    fig, ax = plt.subplots(num_rows, 1, figsize=(25,20))

    # plot MSA row by row
    for row_idx in range(num_rows):
        if num_rows==1:
            ax_row = ax
            start_pos = 0
            end_pos = msa_len
            row_len = msa_len
        else:
            ax_row = ax[row_idx]
            start_pos = row_idx*wrap_length
            end_pos = min((row_idx+1)*wrap_length, msa_len)
            row_len = end_pos-start_pos
        msa_array_row = msa_array[:,start_pos:end_pos]
        msa_numeric_row = msa_numeric[:,start_pos:end_pos]
        annotate_seq_row = annotate_seq[start_pos:end_pos]
        diff_seq_row = diff_seq[start_pos:end_pos]
        print('diff_seq_row:', diff_seq_row)

        # plot msa segment
        sns.heatmap(msa_numeric_row, ax=ax_row, cmap=cmap, norm=norm, cbar=False, xticklabels=False, yticklabels=False)

        # add vertical lines
        for x in range(row_len):
            ax_row.axvline(x=x, linewidth=0.2, color='k')
        # add tick labels
        ax_row.set_xticks(np.arange(0, row_len, xtick_interval))
        ax_row.set_xticklabels(np.arange(start_pos+1, end_pos+1, xtick_interval)+start_pos_offset, fontsize=fontsize)
        # annotate sequence
        if label_residues is not None:
            for res_idx, (res, diff_res) in enumerate(zip(annotate_seq_row, diff_seq_row)):
                if res!='':
                    ax_row.annotate(res, (res_idx, 0), fontsize=fontsize, c=palette[colors[res]])
                if diff_res!='':
                    ax_row.annotate(diff_res, (res_idx, -num_sequences/80), fontsize=fontsize, annotation_clip=False, c=palette[colors[diff_res]])

        if show_seq_names:
            ax_row.set_yticks(np.arange(0, num_sequences, ytick_interval)+0.5)
            ax_row.set_yticklabels(seq_names, fontsize=10)
            if show_all_sequences:
                for seq_num in range(num_sequences):
                    for res_idx, res in enumerate(msa_array_row[seq_num]):
                        if res != '-':
                            pos_counter_dict[seq_num] += 1
                            pos = pos_counter_dict[seq_num]
                            if pos%10==0:
                                ax_row.annotate(str(pos)+'\n'+res, (res_idx, seq_num+0.5), fontsize=fontsize*0.75, c='k')
                            else:
                                ax_row.annotate(''+'\n'+res, (res_idx, seq_num + 0.5), fontsize=fontsize * 0.75, c='k')
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
                  show_seq_names=False, label_residues=None, show_all_sequences=False, fontsize=8, filter_by_refseq=None, savefig=None):
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
                         show_seq_names, label_residues, show_all_sequences, fontsize, filter_by_refseq, savefig)

def symlog(data):
    idx_pos = np.where(data > 0)
    idx_neg = np.where(data < 0)
    data_symlog = np.zeros((data.shape[0], data.shape[1]))
    data_symlog[:] = np.nan
    data_symlog[idx_pos] = np.log(data[idx_pos])
    data_symlog[idx_neg] = -np.log(-data[idx_neg])
    return data_symlog

def annotate_heatmap(array_2D, ax, ndecimals=2, fontsize=8):
    for (j,i),label in np.ndenumerate(array_2D):
        if ndecimals==0 and ~np.isnan(label):
            label = int(label)
        else:
            label = round(label,ndecimals)
        ax.text(i,j,label,ha='center',va='center', color='0.8', fontsize=fontsize, fontweight='bold')

def heatmap(array, c='viridis', ax=None, cbar_kw={}, cbarlabel="", datamin=None, datamax=None, logscale_cmap=False,
            annotate=None, row_labels=None, col_labels=None, show_gridlines=True, fontsize=8):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    if not ax:
        ax = plt.gca()
    cmap = getattr(plt.cm, c)

    # get array size and xy labels
    data = array.astype(float)
    ny, nx = data.shape

    # get row and column labels
    if row_labels is None:
        row_labels = list(np.arange(ny) + 1)
    if col_labels is None:
        col_labels = list(np.arange(nx) + 1)

    # get locations of nan values and negative values, replace values so these don't trigger an error
    naninds = np.where(np.isnan(data) == True)
    infinds = np.where(np.isinf(data) == True)
    if len(infinds[0]) > 0:
        data[infinds] = np.nan
    if len(naninds[0]) > 0:
        data[naninds] = np.nanmean(data)
    if len(infinds[0]) > 0:
        data[infinds] = np.nanmean(data)
    data_cmap = data.copy()

    # get min and max values
    if datamin is None:
        datamin = np.nanmin(data_cmap)
    if datamax is None:
        datamax = np.nanmax(data_cmap)

    # get colormap to plot
    if logscale_cmap:  # plot on logscale
        data_cmap = symlog(data_cmap)
        datamin, datamax = np.min(data_cmap), np.max(data_cmap)

    # get cmap gradations
    dataint = (datamax - datamin) / 100
    norm = plt.Normalize(datamin, datamax + dataint)
    # convert data array into colormap
    colormap = cmap(norm(data_cmap))

    # Set the positions of nan values in colormap to 'lime'
    colormap[naninds[0], naninds[1], :3] = 0, 1, 0
    colormap[infinds[0], infinds[1], :3] = 1, 1, 1

    # plot colormap
    im = ax.imshow(colormap, interpolation='nearest')

    # Create colorbar
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.07)
    cbar = ax.figure.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)

    if logscale_cmap == True:
        cbar_labels = cbar.ax.get_yticks()
        cbar.set_ticks(cbar_labels)
        cbar_labels_unlog = list(np.round(np.exp(np.array(cbar_labels)), 2))
        cbar.set_ticklabels(cbar_labels_unlog)

    # Turn off gridlines if required
    ax.tick_params(axis='both', which='both', length=0, gridOn=show_gridlines)

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontsize=fontsize, ha="right")
    ax.set_yticklabels(row_labels, fontsize=fontsize)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90,
             rotation_mode="anchor")

    # Annotate
    if annotate is not None:
        if isinstance(annotate, int):
            ndecimals = annotate
        else:
            ndecimals = 3
        annotate_heatmap(array, ax, ndecimals=ndecimals, fontsize=fontsize)

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    # set xticks
    ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar, ax

def plot_variant_heatmap(arr, seq, N_res_per_heatmap_row, aaList, seq_name=None, savefig=None, figtitle=None):
    import matplotlib.pyplot as plt
    # Visualize the heatmaps
    seq_len = len(seq)
    residue_num = list(np.arange(1, seq_len + 1))
    num_heatmaps = int(np.ceil(seq_len / N_res_per_heatmap_row))
    heatmap_min = np.min(arr)
    heatmap_max = np.max(arr)
    fig, ax = plt.subplots(num_heatmaps, 1, figsize=(N_res_per_heatmap_row / len(aaList) * 4, num_heatmaps * 4))
    for k in range(num_heatmaps):
        if num_heatmaps == 1:
            ax_k = ax
        else:
            ax_k = ax[k]
        residue_num_k = residue_num[k * N_res_per_heatmap_row:min((k + 1) * N_res_per_heatmap_row, seq_len)]
        start_idx = k * N_res_per_heatmap_row
        end_idx = min((k + 1) * N_res_per_heatmap_row, seq_len)
        heatmap_k = arr[:, start_idx:end_idx]
        wt_idxs_k = np.array([[aaList.index(wt_aa),res_idx] for res_idx, wt_aa in enumerate(seq[start_idx:end_idx])])
        im = ax_k.imshow(heatmap_k, cmap="viridis", aspect="auto", vmin=heatmap_min, vmax=heatmap_max)
        ax_k.scatter(wt_idxs_k[:,1], wt_idxs_k[:,0], c='r', s=4)
        ax_k.set_yticks(range(20), aaList)
        ax_k.set_xticks(range(len(residue_num_k)), residue_num_k, fontsize=7, rotation=45)
    fig.colorbar(im, orientation='vertical')
    if figtitle is not None:
        if seq_name is not None:
            figtitle = seq_name + ': ' + figtitle
        plt.suptitle(figtitle, y=0.93, fontsize=16)
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight')
    plt.show()
