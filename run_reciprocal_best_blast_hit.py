
#!/usr/bin/env python3

'''
Script for finding orthologs using reciprocal BLAST hits.

Sript usage is as follows:
    ./find_orthologs.py -i1 <Input file 1> -i2 <Input file 2> -o <Output file name> –t <Sequence type – n/p>

where "n" specifies a nucleotide sequence and "p" specifies a protein sequence.
'''
import os
import shutil
import argparse
import subprocess
from zipfile import ZipFile
from Bio import Entrez
from variables import address_dict, subfolders
from blastp import run_blastp
from parse_seqsearch_output import parse_blastp_output

class ReciprocaBestHit:
    def __init__(
            self,
            data_folder,
            data_subfolder,
            sequences_subfolder='sequences/',
            seqsearch_subfolder='seqsearch/',
            db_subfolder='by_taxon/'
    ):
        self.data_folder = data_folder
        self.data_subfolder = data_subfolder
        self.sequences_dir = data_folder + sequences_subfolder + data_subfolder + '/'
        self.seqsearch_dir = data_folder + seqsearch_subfolder + data_subfolder + '/'
        self.db_subfolder = db_subfolder
        self.cwd = os.getcwd()

    def get_NCBI_accession(self, query_str):

        cmd = 'esearch -db gene -query "xylanase AND Aspergillus niger[Organism]" | efetch -format docsum'
        cmd = 'esearch -db protein -query "P55329" | elink -target gene | efetch -format docsum'
        subprocess.run(cmd)
        return seq_record

    def get_orthologs(self, accession, output_fpath):
        cmd = f'datasets download gene accession {accession} --ortholog'
        subprocess.run(cmd)
        return None

    def get_taxon_data(self, species, output_fpath):
        """once we have the taxid, we can fetch the record"""
        # run datasets download
        output_dir = os.path.dirname(output_fpath)
        os.chdir(output_dir)
        cmd = f'datasets download gene taxon "{species}" --include protein'
        subprocess.run(cmd)
        os.chdir(self.cwd)

        # unzip downloaded data
        ncbi_dataset_zip = output_dir+'/'+'ncbi_dataset.zip'
        ncbi_dataset_fpath = output_dir+'/'+'ncbi_dataset'
        with ZipFile(ncbi_dataset_zip, 'r') as zip_object:
            zip_object.extractall(ncbi_dataset_fpath)

        # extract fasta file of proteins from that taxon
        source_fasta_fpath = f'{ncbi_dataset_fpath}/ncbi_dataset/data/protein.faa'
        shutil.copy(source_fasta_fpath, output_fpath)

        # delete downloaded ncbi directory (zipped and unzipped)
        shutil.rmtree(ncbi_dataset_fpath)
        os.remove(ncbi_dataset_zip)

    def blast_local_taxon_db(self, query_fasta, output_file, db_name, e_thres=1e-5, max_target_seqs=100):
        """
        TO EDIT
        """
        run_blastp(query_fasta, db_name, output_file, e_thres, max_target_seqs, self.sequences_dir, self.seqsearch_dir, self.db_subfolder, output_dir, outfmt, other_params)

        "-db", f'{db_dir}{db_name}/{db_name}.fasta',

    def get_reciprocal_hits(self, file_one,file_two, input_sequence_type, tmp_db_dir, output_file):

        a = []                  #List a stores the first index headers (db values)
        a1 = []                 #List a1 stores the second index headers (query headers)
        b = []                  #List b stores the first index headers (db values)
        b1 = []
        output_list_1 = []
        output_list_2 = []      #Inversing the list in order to find common things bw output_list_1 and output_list_2
        output_list = []        #Contains only reciprocal blast hits

        if input_sequence_type == "n":  #For nucleotide(n)

            file_1_db = "makeblastdb -in " + str(file_one) + f" -dbtype prot -out {tmp_db_dir}/db1"   #Creating db for input1
            file_1_db_subprocess = subprocess.check_output(file_1_db.split())

            result_1 = f"blastn -db {tmp_db_dir}/db1 -query " + str(file_two) + f" -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out {tmp_db_dir}/output_file_1" #Output for query = input2 and db = input1
            result_1_subprocess = subprocess.check_output(result_1.split())

            file_2_db = "makeblastdb -in " + str(file_two) + f" -dbtype prot -out {tmp_db_dir}/db2"   #Creating db for input2
            file_2_db_subprocess = subprocess.check_output(file_2_db.split())

            result_2 = f"blastn -db {tmp_db_dir}/db2 -query " + str(file_one) + f" -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out {tmp_db_dir}/output_file_2" #Output for query = input1 and db = input2
            result_2_subprocess = subprocess.check_output(result_2.split())


            with open("tmp/output_file_1") as fh1:
                for line in fh1.readlines():
                    if line.startswith("lcl"):
                        a.append(line.split()[0])           #List a stores the first index headers (db values)
                    if line.startswith("lcl"):
                        a1.append(line.split()[1])          #List a1 stores the second index headers (query headers)

            for i,j in zip(a,a1):
                output_list_1.append(i+"\t"+j+"\n")

            with open("tmp/output_file_2") as fh2:
                for line in fh2.readlines():
                    if line.startswith("lcl"):
                        b.append(line.split()[0])           #List b stores the first index headers (db values)
                    if line.startswith("lcl"):
                        b1.append(line.split()[1])          #List b1 stores the second index headers (query headers)

            for i,j in zip(b,b1):
                output_list_2.append(j+"\t"+i+"\n")         #Inversing the list in order to find common things bw output_list_1 and output_list_2

            for x in output_list_1:                         #Checking common hits beween output_list_1 and output_list_2
                if x in output_list_2:
                    output_list.append(x)                   #Contains only reciprocal blast hits

        elif input_sequence_type == "p":  #For protein(p)

            file_1_db = "makeblastdb -in " + str(file_one) + " -dbtype prot -out tmp/db1"   #Creating db for input1
            file_1_db_subprocess = subprocess.check_output(file_1_db.split())

            result_1 = "blastp -db tmp/db1 -query " + str(file_two) + " -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out "+str(output_file) #Output for query = input2 and db = input1
            result_1_subprocess = subprocess.check_output(result_1.split())

            file_2_db = "makeblastdb -in " + str(file_two) + " -dbtype prot -out tmp/db2"   #Creating db for input2
            file_2_db_subprocess = subprocess.check_output(file_2_db.split())

            result_2 = "blastp -db tmp/db2 -query " + str(file_one) + " -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out "+str(output_file) #Output for query = input1 and db = input2
            result_2_subprocess = subprocess.check_output(result_2.split())

            with open("tmp/output_file_1") as fh1:
                for line in fh1.readlines():
                    if line.startswith("lcl"):
                        a.append(line.split()[0])
                    if line.startswith("lcl"):
                        a1.append(line.split()[1])

            for i,j in zip(a,a1):
                output_list_1.append(i+"\t"+j+"\n")

            with open("tmp/output_file_2") as fh2:
                for line in fh2.readlines():
                    if line.startswith("lcl"):
                        b.append(line.split()[0])
                    if line.startswith("lcl"):
                        b1.append(line.split()[1])

            for i,j in zip(b,b1):
                output_list_2.append(j+"\t"+i+"\n")

            for x in output_list_1:         #Checking common hits beween output_list_1 and output_list_2
                if x in output_list_2:
                    output_list.append(x)   #Contains only reciprocal blast hits

        return(output_list)

def main():
    # Argparse code

    parser = argparse.ArgumentParser()
    parser.add_argument("-i1", help = "Takes in input sequence 1.")
    parser.add_argument("-i2", help = "Takes in input sequence 2.")
    parser.add_argument("-o", help = "Your output file.")
    parser.add_argument("-t", help = "n for nucleotide sequence, p for protein sequence.")
    args = parser.parse_args()

    file_one = args.i1             #Populating the variables
    file_two = args.i2
    output_file = args.o
    input_sequence_type = args.t

    # Subprocess code

    '''
    output_list is a list of reciprocal BLAST hits. Each element is a tab
    separated pair of gene names. Eg:
    ["lcl|AM421808.1_cds_CAM09336.1_10	lcl|AE002098.2_cds_NMB0033_33", "lcl|AM421808.1_cds_CAM09337.1_11	lcl|AE002098.2_cds_NMB0034_34", "lcl|AM421808.1_cds_CAM09338.1_12	lcl|AE002098.2_cds_NMB0035_35", "lcl|AM421808.1_cds_CAM09339.1_13	lcl|AE002098.2_cds_NMB0036_36", ...]
    '''
    data_folder = address_dict['ECOHARVEST']
    tmp_db_dir = data_folder + subfolders['seqsearch'] + 'tmp'
    make_dir = f'mkdir {tmp_db_dir}'
    subprocess.call(make_dir.split())

    output_list = get_reciprocal_hits(file_one, file_two, input_sequence_type, tmp_db_dir, output_file)
    with open(output_file, 'w') as output_fh:
        for ortholog_pair in output_list:
            output_fh.write(ortholog_pair)

    remove_dir = f"rm -rf {tmp_db_dir}"
    subprocess.call(remove_dir.split())


if __name__ == "__main__":

    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'sidestream_cocktail'
    seqsearch_dir = data_folder + subfolders['seqsearch'] + data_subfolder + '/'
    sequences_dir = data_folder + subfolders['sequences'] + data_subfolder + '/'
    species_name = 'Aspergillus niger'
    format_type = 'XML'
    hitlist_size = 5
    output_fpath = seqsearch_dir + f'{species_name}_qblast.xml'
    # taxid  = get_tax_id(species_name)
    # print(taxid)

    # get sequence record
    sequence = 'MKVTAAFAGLLVTAFAAPVPEPVLVSRSAGINYVQNYNGNLGDFTYDESAGTFSMYWEDGVSSDFVVGLGWTTGSSKAITYSAEYSASGSSSYLAVYGWVNYPQAEYYIVEDYGDYNPCSSATSLGTVYSDGSTYQVCTDTRTNEPSITGTSTFTQYFSVRESTRTSGTVTVANHFNFWAQHGFGNSDFNYQVMAVEAWSGAGSASVTISS'
    seq_record = get_seq_id(sequence)
    print(seq_record.id)



    # # get taxonomic ID
    # taxid = 5061
    # out = get_tax_data(taxid, seqsearch_dir+f'{species_name}_{taxid}.fasta')

    # main()
