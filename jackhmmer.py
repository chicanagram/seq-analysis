import subprocess
import time
from variables import address_dict, subfolders

def run_jackhmmer(query_fasta, db_name, output_file, e_thres=1e-5, incE_thres=1e-5, max_target_seqs=None, num_cpu=None, query_dir='', db_dir='', output_dir='', outfmt=0):
    # Define PHMMER command
    cmd = [
        "jackhmmer",
        "-o", output_dir + output_file,
        "-E", str(e_thres),
        "--incE", str(incE_thres),
        "--cpu", str(num_cpu),
              query_dir + query_fasta,
              f'{db_dir}{db_name}/{db_name}.fasta',
    ]
    print(' '.join(cmd))

    # Run PHMMER
    try:
        start = time.time()
        subprocess.run(cmd, check=True, encoding="utf-8")
        end = time.time()
        print(
            f"Jackhmmer search completed ({round((end - start) / 60, 3)} min). Results saved in {output_dir + output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error running Jackhmmer: {e}")
    except FileNotFoundError:
        print("Error: Jackhmmer is not installed or not in the system PATH.")


if __name__ == '__main__':
    # Define input files and output file
    search_type = 'jackhmmer'
    query_dir =  address_dict['ECOHARVEST'] + subfolders['sequences']
    output_dir = address_dict['ECOHARVEST'] + subfolders['seqsearch']
    db_dir = address_dict['databases']
    query_fasta = 'CALB.fasta'
    db_name = 'uniprot_trembl' # 'uniprot_sprot' # 'Pfam-A.full' # 'refProteomes'
    e_thres = 0.001 #
    incE_thres = 0.00001 # 0.1 #
    max_target_seqs = None
    num_cpu = 8
    output_file = query_fasta.split('.')[0] + '_' + search_type + '_' + db_name + '_incE'+"{:.0e}".format(incE_thres) + '_E'+"{:.0e}".format(e_thres)
    # run jackhmmer
    run_jackhmmer(query_fasta, db_name, output_file, e_thres, incE_thres, max_target_seqs, num_cpu, query_dir, db_dir, output_dir)

