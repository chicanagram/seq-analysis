import os
import multiprocessing
from blastp import run_blastp
from phmmer import run_phmmer
from jackhmmer import run_jackhmmer
from parse_seqsearch_output import parse_hmmer_output, parse_blastp_output
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta, run_msa
from variables import address_dict, subfolders
from convert_msa_format import convert_msa_format

def run_seqsearch(search_type, query_fasta, db_fasta, output_file, e_thres, incE_thres, max_target_seqs, num_cpu, query_dir, db_dir, output_dir):
    if search_type=='blastp':
        run_blastp(query_fasta, db_fasta, output_file, e_thres, max_target_seqs, query_dir, db_dir, output_dir)
    elif search_type=='phmmer':
        run_phmmer(query_fasta, db_fasta, output_file, e_thres, incE_thres, max_target_seqs, num_cpu, query_dir, db_dir, output_dir)
    elif search_type=='jackhmmer':
        run_jackhmmer(query_fasta, db_fasta, output_file, e_thres, incE_thres, max_target_seqs, num_cpu, query_dir, db_dir, output_dir)

def run_search_and_align_pipeline(
        query_seq,
        query_seq_name,
        db_name,
        e_thres,
        incE_thres,
        max_target_seqs,
        num_cpu,
        search_type,
        msa_method,
        msa_fmt,
        seqsearch_dir,
        sequences_dir,
        msa_dir,
        db_dir,
        i=None,
        num_queries=None,
        label_with_index=True
):
    query_seq_name_actual = query_seq_name
    if label_with_index:
        query_seq_name = str(i)

    if incE_thres is None or search_type == 'blastp':
        seqsearch_raw = query_seq_name + '_' + search_type + '_' + db_name + '_E' + "{:.0e}".format(e_thres)
    else:
        seqsearch_raw = query_seq_name + '_' + search_type + '_' + db_name + '_incE' + "{:.0e}".format(
            incE_thres) + '_E' + "{:.0e}".format(e_thres)

    print(seqsearch_raw)
    seqsearch_csv = seqsearch_raw + '.csv'
    seqsearch_fasta = seqsearch_raw + '.fasta'
    seqsearch_msa = f'{seqsearch_raw}_{msa_method}.fasta'
    incomplete_processing = []
    try:
        # temporarily generate fasta file for sequence
        query_fasta = query_seq_name + '_temp.fasta'
        write_sequence_to_fasta([query_seq], [query_seq_name], query_fasta.replace('.fasta', ''), sequences_dir)

        # run homolog search
        print(f'Running {search_type} homolog search for {i + 1}/{num_queries}: {query_seq_name_actual} ({len(query_seq)})')
        run_seqsearch(search_type, query_fasta, db_name, seqsearch_raw, e_thres, incE_thres, max_target_seqs, num_cpu,
                      sequences_dir, db_dir, seqsearch_dir)

        # parse homolog search output to get fasta input for MSA generation
        print(f'Parsing results from homolog search for {i + 1}/{num_queries}: {query_seq_name_actual} ({len(query_seq)})')
        if search_type in ['phmmer', 'jackhmmer']:
            parse_hmmer_output(seqsearch_raw, seqsearch_csv, seqsearch_fasta.replace('.fasta', ''), query_fasta,
                               sequences_dir, seqsearch_dir)
        elif search_type in ['blastp']:
            parse_blastp_output(seqsearch_raw, seqsearch_csv, seqsearch_fasta.replace('.fasta', ''), query_fasta,
                                sequences_dir, seqsearch_dir)
        # delete raw seqsearch file
        os.remove(seqsearch_dir+seqsearch_raw)
        print('Removed raw sequence search results file after parsing to FASTA and CSV:', seqsearch_raw)

        # MSA generation
        print(f'Performing {msa_method} alignment for {i + 1}/{num_queries}: {query_seq_name_actual} ({len(query_seq)})')
        run_msa(seqsearch_fasta, seqsearch_msa, method=msa_method, seq_dir=sequences_dir, msa_dir=msa_dir)
        if msa_fmt != 'fasta':
            convert_msa_format(seqsearch_msa, 'fas', msa_fmt, msa_dir, from_perl_or_python='perl', options=None)

        # delete temporarily generated FASTA file for query sequence
        os.remove(sequences_dir + query_fasta)
        print(f'Removed temporary query sequence FASTA file for {i}: {query_seq_name_actual}')

        # log
        logfile_complete = sequences_dir + 'run_search_and_align_COMPLETE.txt'
        if os.path.exists(logfile_complete):
            mode = 'a'
        else:
            mode = 'w'
        with open(logfile_complete, mode) as f:
            f.write(f'{i},{query_seq_name_actual}\n')
    except Exception as e:
        print(f'Did not finish processing {i}: {query_seq_name_actual}:', e)
        incomplete_processing.append((i, query_seq_name_actual))
        # log
        logfile_incomplete = sequences_dir + 'run_search_and_align_INCOMPLETE.txt'
        if os.path.exists(logfile_incomplete):
            mode = 'a'
        else:
            mode = 'w'
        with open(logfile_incomplete, mode) as f:
            f.write(f'{i},{query_seq_name_actual}\n')

    return incomplete_processing


if __name__=='__main__':
    search_type = 'blastp'
    msa_method = 'mafft'
    msa_fmt = 'a3m'
    data_folder = address_dict['PON-Sol2'] # address_dict['ECOHARVEST'] #
    seqsearch_dir = data_folder + subfolders['seqsearch']
    sequences_dir = data_folder + subfolders['sequences']
    msa_dir = data_folder + subfolders['msa']
    db_dir = address_dict['databases']
    input_seqs_fasta = 'ponsol2_WT.fasta' # 'CALA_UML_KLIP.fasta'
    db_name = 'uniprot_trembl'  #  'uniprot_sprot' # 'Pfam-A'
    e_thres = 1e-4  # 1e-3 #
    incE_thres = None # 0.001
    max_target_seqs = 750
    num_cpu = None
    num_processes = 3

    # get all sequences to build MSAs for
    query_seqs, query_seq_names, _ = fetch_sequences_from_fasta(sequences_dir+input_seqs_fasta)
    num_queries = len(query_seq_names)
    args_list = []
    i_offset = 0
    for i, (query_seq, query_seq_name) in enumerate(zip(query_seqs[:30],query_seq_names[:30])):
        args = (
            query_seq,
            query_seq_name,
            db_name,
            e_thres,
            incE_thres,
            max_target_seqs,
            num_cpu,
            search_type,
            msa_method,
            msa_fmt,
            seqsearch_dir,
            sequences_dir,
            msa_dir,
            db_dir,
            i+i_offset,
            num_queries
        )
        args_list.append(args)

    # Create multiprocessing pool
    with multiprocessing.Pool(processes=num_processes) as pool:
        res = pool.starmap(run_search_and_align_pipeline, args_list)
    incomplete_processing = [entry[0] for entry in res if len(entry)>0]

    print('INCOMPLETE:', len(incomplete_processing))
    for entry in incomplete_processing:
        print(entry)