from variables import address_dict, subfolders
from utils import run_msa, filter_and_save_sequences, run_clipkit, replace_edge_gaps

def run_pycdhit(data_folder, seq_fname, thres):
    import os
    import subprocess
    cwd = os.getcwd()
    os.chdir('../py-cdhit/')

    # Run the script inside the Conda environment
    command = f'conda run -n pycdhit python run_pycdhit.py --dir {data_folder} --f {seq_fname} --thres {thres}'
    subprocess.run(command, shell=True, executable="/bin/bash")
    os.chdir(cwd)


def main():
    data_folder = address_dict['ECOHARVEST']
    msa_dir = data_folder + subfolders['msa']
    seq_dir = data_folder + subfolders['sequences']
    seq_fname = 'CALA_phmmer_uniprot_trembl_incE1e-03_E1e-03.fasta' # 'lipases_initialMSA.fasta' # 'CALB_patent_hits.fasta' # 'UPO.fasta' # 'UPO_filtered.fasta' #
    msa_method = 'mafft' # 'clustalo' #
    msa_suffix = '_filt' # ''
    msa_fname = seq_fname[:seq_fname.find('.fa')] + f'_{msa_method}' + msa_suffix + '.fasta'

    # steps to execute
    filter_msa = None # ['CALA','<', 5, 0.5, msa_method] #
    cluster_msa = None # 0.99 #
    trim_msa = 'kpic-gappy' # 'smart-gap' # None #
    filter_by_refseq = 0  # None #
    remove_terminalX = False

    # plot parameters
    plot_msa = 'seaborn' # 'pymsaviz' # None #
    plot_msa_pos_range = [200,236] # None #
    wrap_length = 400 # 600
    label_residues = 'ref' # 'consensus' #
    ytick_interval = 100
    if ytick_interval==1:
        show_seq_names = True
    else:
        show_seq_names = False
    savefig = None
    if plot_msa_pos_range is None:
        fontsize = 8
        xtick_interval = 20
    else:
        fontsize = 20
        xtick_interval = 5

    # get input msa fname
    msa_input = msa_fname

    # FILTER MSA
    if filter_msa is not None:
        print('Filtering MSA by sequence identity...')
        seq_filt_fname = filter_and_save_sequences(msa_input, reference_id=filter_msa[0], inequality=filter_msa[1], id_threshold=filter_msa[2], len_threshold=filter_msa[3],
                                                   output_filtered_fasta=msa_input.replace('.fasta', '_filt.fasta'), msa_dir=msa_dir, seq_dir=seq_dir)
        # regenerate MSA
        msa_output = run_msa(seq_filt_fname, msa_input.replace('.fasta', '_filt.fasta'), filter_msa[4], seq_dir, msa_dir)
        msa_input = msa_output

    # CLUSTER MSA
    if cluster_msa is not None:
        run_pycdhit(seq_dir, seq_fname, cluster_msa)

    # TRIM MSA
    if trim_msa is not None:
        print('Trimming MSA...')
        msa_output = run_clipkit(msa_input, msa_input.replace('.fasta', '_trimmed.fasta'), msa_dir=msa_dir, keep_ref_seq_positions=filter_by_refseq, mode=trim_msa)
        msa_input = msa_output

    # REPLACE EDGE GAPS IN MSA
    if remove_terminalX:
        print('Removing terminal gaps in MSA...')
        msa_output = replace_edge_gaps(msa_input, msa_input.replace('.fasta', '_termX.fasta'), msa_dir=msa_dir)

if __name__ == "__main__":
    main()