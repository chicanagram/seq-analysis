from variables import address_dict, subfolders


def fetch_sequence_by_accession(accession, fasta_fpath=None, email_id="charmaine_chia@bii.a-star.edu.sg"):
    from Bio import Entrez, SeqIO
    Entrez.email = email_id  # required by NCBI
    handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()

    if fasta_fpath is not None:
        SeqIO.write(record, fasta_fpath, "fasta")
        print(f"Sequence saved as {fasta_fpath}")


if __name__=='__main__':
    data_folder = address_dict['influenza-resistance']
    data_subfolder = 'IC50benchmark'
    accession_fname = 'IC50benchmark_Zanamivir_ACC_deduped.txt'
    accession_fpath = data_folder + subfolders['sequences'] + data_subfolder + '/' + accession_fname
    missed_accessions = []
    # get accession ids
    with open(accession_fpath, 'r') as f:
        accession_ids = f.read()
    accession_ids = accession_ids.split(', ')
    print(len(accession_ids), accession_ids)

    for accession in accession_ids:
        try:
            fasta_fpath = data_folder + subfolders['sequences'] + data_subfolder + '/' + accession + '.fasta'
            fetch_sequence_by_accession(accession, fasta_fpath)
        except:
            missed_accessions.append(accession)
            print('Unable to fetch the sequence for ', accession)
    print('Missed accessions:', len(missed_accessions), missed_accessions)
