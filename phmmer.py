import subprocess
import time
from variables import address_dict, subfolders

# Define input files and output file
search_type = 'phmmer'
query_dir =  address_dict['ECOHARVEST'] + subfolders['sequences']
output_dir = address_dict['ECOHARVEST'] + subfolders['seqsearch']
db_dir = address_dict['databases']
query_fasta = 'CALA.fasta'
db_fasta = 'uniprot_trembl.fasta.gz' #  'Pfam-A.full' # 'uniprot_sprot.fasta' # 'Pfam-A.full.gz' # 'refProteomes.fasta'
e_thres = 0.001 # 0.01 #
incE_thres = 0.001 # 0.1 #
num_cpu = 4
output_file = query_fasta.split('.')[0] + '_' + search_type + '_' + db_fasta.split('.')[0] + '_incE='+"{:.0e}".format(incE_thres) + '_E='+"{:.0e}".format(e_thres)

# Define PHMMER command
phmmer_cmd = [
    "phmmer",
    "-o", output_dir + output_file+'.out',
    # "-A", output_file+'.ali',
    # "--tblout", output_file+'.tblout',
    "-E", str(e_thres),
    "--incE", str(incE_thres),
    "--cpu", str(num_cpu),
    query_dir + query_fasta,
    db_dir + db_fasta,
]
print(' '.join(phmmer_cmd))

# Run PHMMER
try:
    start = time.time()
    subprocess.run(phmmer_cmd, check=True, encoding="utf-8")
    end = time.time()
    print(f"PHMMER search completed ({round((end-start)/60,3)} min). Results saved in {output_file}")
except subprocess.CalledProcessError as e:
    print(f"Error running PHMMER: {e}")
except FileNotFoundError:
    print("Error: PHMMER is not installed or not in the system PATH.")
