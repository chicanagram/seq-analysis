import subprocess
import time
from variables import address_dict, subfolders

# Define input files and output file
search_type = 'jackhmmer'
query_dir =  address_dict['ECOHARVEST'] + subfolders['sequences']
output_dir = address_dict['ECOHARVEST'] + subfolders['seqsearch']
db_dir = address_dict['databases']
query_fasta = 'CALB.fasta'
db_fasta = 'uniprot_trembl.fasta' # 'uniprot_sprot.fasta' # 'Pfam-A.full' # 'refProteomes.fasta'
e_thres = 0.001 #
incE_thres = 0.00001 # 0.1 #
num_cpu = 8
output_file = query_fasta.split('.')[0] + '_' + search_type + '_' + db_fasta.split('.')[0] + '_incE'+"{:.0e}".format(incE_thres) + '_E'+"{:.0e}".format(e_thres)

# Define PHMMER command
jackhmmer_cmd = [
    "jackhmmer",
    "-o", output_dir + output_file+'.out',
    # "-A", output_dir + output_file+'.ali',
    # "--tblout", output_dir + output_file+'.tblout',
    # "--domtblout", output_dir + output_file+'.domtblout',
    "-E", str(e_thres),
    "--incE", str(incE_thres),
    "--cpu", str(num_cpu),
    query_dir + query_fasta,
    db_dir + db_fasta,
]
print(' '.join(jackhmmer_cmd))

# Run PHMMER
try:
    start = time.time()
    subprocess.run(jackhmmer_cmd, check=True, encoding="utf-8")
    end = time.time()
    print(f"Jackhmmer search completed ({round((end-start)/60,3)} min). Results saved in {output_dir + output_file}")
except subprocess.CalledProcessError as e:
    print(f"Error running Jackhmmer: {e}")
except FileNotFoundError:
    print("Error: Jackhmmer is not installed or not in the system PATH.")
