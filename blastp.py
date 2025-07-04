import subprocess
import time
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta

def run_blastp(query_fasta, db_name, output_name, e_thres=1e-5, max_target_seqs=None, query_dir='', db_dir='', df_subfolder='', output_dir='', outfmt=0, other_params={}, run_multiprocessing=False):
    start_time = time.time()
    append_to_cmd = []
    for param_name, param_val in other_params.items():
        append_to_cmd.append('-'+param_name)
        append_to_cmd.append(str(param_val))
    cmd = [
        "blastp",
        "-query", query_dir + query_fasta,
        "-db", f'{db_dir}{df_subfolder}/{db_name}.fasta',
        "-out", output_dir + output_name + '.out',
        "-evalue", str(e_thres),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", str(outfmt)
    ] + append_to_cmd

    print(' '.join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    end_time = time.time()
    time_taken_min = (end_time - start_time) / 60
    if result.returncode == 0:
        print(f"BLASTP finished successfully. Time taken: {round(time_taken_min, 3)} min")
    else:
        print("BLASTP failed:")
        print(result.stderr)

def run_biopython_qblast(input, format_type='Text', hitlist_size=100, output_dir=None, output_name=None, db_name='nr', run_multiprocessing=False):
    from Bio import Blast
    from Bio.Seq import Seq

    if input.find('.fasta')>-1:
        seqs, seq_names, _ = fetch_sequences_from_fasta(input)

        for i, (seq, seq_name) in enumerate(zip(seqs, seq_names)):
            # if i<5: continue
            seq = Seq(seq.replace('-',''))
            print(f'Running BLASTp search on {i}, {seq_name}, {seq}...')
            start_time = time.time()
            result_stream = Blast.qblast("blastp", database=db_name, sequence=seq, ncbi_gi=True, format_type=format_type, hitlist_size=hitlist_size, alignments=hitlist_size, descriptions=hitlist_size)
            if output_name is None:
                out_fname = seq_name + f'_blastp_{db_name}.out'
            else:
                out_fname = output_name[i] + f'_{db_name}.out'
            output_fpath = output_dir + out_fname
            with open(output_fpath, "wb") as out_stream:
                out_stream.write(result_stream.read())
            result_stream.close()
            end_time = time.time()
            time_taken_min = (end_time - start_time) / 60
            print(f"BLASTP finished successfully on {i}, {seq_name}. Time taken: {round(time_taken_min, 3)} min")

def run_biopython_psiblast():
    from Bio.Blast.Applications import NcbipsiblastCommandline
    cline = NcbipsiblastCommandline(help=True)
    NcbipsiblastCommandline(cmd='psiblast', help=True)

if __name__ == '__main__':
    # Define input files and output file
    data_folder = address_dict['ECOHARVEST']
    search_type = 'blastp'
    data_subfolder = 'lipases' # 'sidestream_cocktail' #
    query_dir = data_folder + subfolders['sequences'] + data_subfolder + '/'
    output_dir = data_folder + subfolders['seqsearch'] + data_subfolder + '/'
    db_dir = address_dict['databases']
    query_fasta = 'RML_pdb.fasta' #  'queryseqs_enzcocktail.fasta' #  'CALA.fasta'
    max_target_seqs = 500
    run_local_blast = False
    run_multiprocessing = False

    # run blastp locally
    if run_local_blast:
        db_name = 'uniprot_trembl'
        outfmt = 0
        e_thres = 1e-5
        other_params = {'taxids': 4751}  # {}
        output_name = query_fasta.split('.')[0] + '_' + search_type + '_' + db_name + '_E' + "{:.0e}".format(e_thres)
        run_blastp(query_fasta, db_name, output_name, e_thres, max_target_seqs, query_dir, db_dir, db_name, output_dir, outfmt, other_params, run_multiprocessing=run_multiprocessing)
    # run qblast with biopython
    else:
        db_name = 'pataa'
        outfmt = 'Text'
        run_biopython_qblast(query_dir+query_fasta, format_type=outfmt, hitlist_size=max_target_seqs, output_dir=output_dir, output_name=None, db_name=db_name, run_multiprocessing=run_multiprocessing)