from Bio import Entrez, SeqIO
from search_sequence_by_accession import fetch_sequence_by_accession

Entrez.email = "charmaine_chia@bii.a-star.edu.sg"

def search_protein(strain_name, gene_name=None, geography=None, collection_date=None, max_records=10):
    query_parts = [f'"{strain_name}"[All Fields]']
    if gene_name:
        query_parts.append(f'{gene_name}[Protein name]')
    if geography:
        query_parts.append(f'{geography}[All Fields]')
    if collection_date:
        query_parts.append(f'{collection_date}[All Fields]')
    query = " AND ".join(query_parts)

    handle = Entrez.esearch(db="protein", term=query, retmax=max_records)
    record = Entrez.read(handle)
    handle.close()
    return record['IdList']

def fetch_protein_sequences(ids, filename=None):
    handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
    records = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    if filename:
        SeqIO.write(records, filename, "fasta")
        print(f"Saved {len(records)} protein sequences to {filename}")
    return records

# Example usage
strain = 'A/Schweiz/1/2023' # "A/Tasmania/503/2020"
geo = None # "Australia"
gene = None # "neuraminidase"
collection_year = None # "2020"

ids = search_protein(strain, gene_name=gene, geography=geo, collection_date=collection_year)

print(f"\nProtein IDs found for {strain} ({gene}, {geo}, collected {collection_year}):\n", ids)

if ids:
    proteins = fetch_protein_sequences(ids, filename="Tasmania_503_2020_NA_Australia.fasta")
    for protein in proteins:
        print(protein.id, protein.description)
        print(len(protein.seq), protein.seq)
else:
    print("No matching protein sequences found.")