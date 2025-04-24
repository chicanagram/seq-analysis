import subprocess
import time
from variables import address_dict, subfolders
import os

def run_blastp(query_fasta, db_name, output_file, e_thres=1e-5, max_target_seqs=None, query_dir='', db_dir='', output_dir='', outfmt=0):
    start_time = time.time()
    cmd = [
        "blastp",
        "-query", query_dir + query_fasta,
        "-db", f'{db_dir}{db_name}/{db_name}.fasta',
        "-out", output_dir + output_file,
        "-evalue", str(e_thres),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", str(outfmt)
    ]
    print(' '.join(cmd))
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    end_time = time.time()
    time_taken_min = (end_time - start_time) / 60
    if result.returncode == 0:
        print(f"BLASTP finished successfully. Time taken: {round(time_taken_min, 3)} min")
    else:
        print("BLASTP failed:")
        print(result.stderr)

def run_biopython_blastp(protein_sequence, blast_fname, data_folder, blast_subfolder=subfolders['sequences'], expect=10.0, hitlist_size=50):
    """
    Run a blast search given an input sequence
    """
    from Bio.Blast import NCBIWWW
    print('Starting blast search...')
    result_handle = NCBIWWW.qblast("blastp", "nr", protein_sequence, expect=expect, hitlist_size=hitlist_size)
    print(f'Completed blast search.')
    # save result
    if blast_fname is not None:
        with open(f"{data_folder}{blast_subfolder}{blast_fname}.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()
        print()
    return result_handle

def run_biopython_psiblast():
    from Bio.Blast.Applications import NcbipsiblastCommandline
    cline = NcbipsiblastCommandline(help=True)
    NcbipsiblastCommandline(cmd='psiblast', help=True)

if __name__ == '__main__':
    # Define input files and output file
    search_type = 'blastp'
    query_dir = address_dict['ECOHARVEST'] + subfolders['sequences']
    output_dir = address_dict['ECOHARVEST'] + subfolders['seqsearch']
    db_dir = address_dict['databases']
    query_fasta = 'CALA.fasta'
    db_name = 'uniprot_trembl' # 'uniprot_sprot' #
    e_thres = 1e-5
    max_target_seqs = 1000
    output_file = query_fasta.split('.')[0] + '_' + search_type + '_' + db_name + '_E' + "{:.0e}".format(e_thres)
    outfmt = 0

    # run blastp
    run_blastp(query_fasta, db_name, output_file, e_thres, max_target_seqs, query_dir, db_dir, output_dir, outfmt)