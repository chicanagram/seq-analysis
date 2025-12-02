#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
from variables import address_dict, subfolders
data_folder = address_dict['examples']
database_folder = address_dict['databases']
sequence_subfolder = subfolders['sequences']
msa_subfolder = subfolders['msa']
pdb_subfolder = subfolders['pdb']

def read_pfam(file_handle,query_number):

    from collections import OrderedDict

    line = file_handle.readline()
    query_find=False #the indicate whether we have found the query

    seqs = {}
    ids = OrderedDict()
    gs = {}
    gr = {}
    gf = {}
    gc = {}
    passed_end_alignment = False

    #assume the query only appears once in the document
    while line and not query_find: #find the line where the query number would be there
        if line[:7] == "#=GF AC":
            print(line)
            AC = line[7:].strip().split(".",1)[0]
            if AC == query_number:
                query_find = True
                break
        try:
            line=file_handle.readline()
        except UnicodeDecodeError:
            print(f'Reached end of database. {query_number} not found.')

    # when the chunk the alignment is found
    while not passed_end_alignment and query_find and line:
        if line[:2] == "//": #the end of the alignment and leave the loop
            passed_end_alignment = True
            break
        elif line == "": # blank line, ignore
            pass
        elif line[0] != "#":
            # Sequence
            # Format: "<seqname> <sequence>"
            assert not passed_end_alignment
            parts = [x.strip() for x in line.split(" ", 1)]
            if len(parts) != 2:
                # This might be someone attempting to store a zero length sequence?
                raise ValueError(
                    "Could not split line into identifier "
                    "and sequence:\n" + line)
            seq_id, seq = parts
            if seq_id not in ids:
                ids[seq_id] = True
                seqs.setdefault(seq_id, '')
                seqs[seq_id] += seq.replace(".", "-")
        elif len(line) >= 5:
            # Comment line or meta-data
            if line[:5] == '#=GC ':
            # Generic per-Column annotation, exactly 1 char per column
            # Format: "#=GC <feature> <exactly 1 char per column>"
                feature, text = line[5:].strip().split(None, 2)
                if feature not in gc:
                    gc[feature] = ""
                gc[feature] += text.strip()  # append to any previous entry
                    # Might be interleaved blocks, so can't check length yet
            elif line[:5] == '#=GS ':
            # Generic per-Sequence annotation, free text
            # Format: "#=GS <seqname> <feature> <free text>"
                seq_id, feature, text = line[5:].strip().split(None, 2)
                # if seq_id not in ids:
                #    ids.append(seq_id)
                if seq_id not in gs:
                    gs[seq_id] = {}
                if feature not in gs[seq_id]:
                    gs[seq_id][feature] = [text]
                else:
                    gs[seq_id][feature].append(text)
            elif line[:5] == "#=GR ":
            # Generic per-Sequence AND per-Column markup
            # Format: "#=GR <seqname> <feature> <exactly 1 char per column>"
                seq_id, feature, text = line[5:].strip().split(None, 2)
                # if seq_id not in ids:
                #    ids.append(seq_id)
                if seq_id not in gr:
                    gr[seq_id] = {}
                if feature not in gr[seq_id]:
                    gr[seq_id][feature] = ""
                gr[seq_id][feature] += text.strip()  # append to any previous entry
        line=file_handle.readline()

    assert len(seqs) <= len(ids)

    return ids,seqs
def pfam_access_functions(query):

    import os
    import pickle
    from prody import searchPfam, fetchPfamMSA, parsePfamPDBs
    from Bio import AlignIO
    fn = query['query_fn']
    val = query['query_val']

    # NOT WORKING
    if fn=='searchPfam':
        # query format: uniprot or PDB id e.g. :uniprot:'PIWI_ARCFU', :pdb:'3luc'
        res = searchPfam(val)

    # NOT WORKING
    elif fn == 'fetchPfamMSA':
        res = fetchPfamMSA(val, format='stockholm', folder=data_folder + msa_subfolder) # formatting is not working
        in_file = f'{data_folder}{msa_subfolder}{val}_full.sth'
        out_file = f'{data_folder}{msa_subfolder}{val}_alignment.fasta'
        _ = AlignIO.convert(in_file, 'stockholm', out_file, 'fasta')
        os.remove(in_file)

    # NOT WORKING
    elif fn == 'parsePfamPDBs':
        cwd = os.getcwd()
        os.chdir(f'{data_folder}{pdb_subfolder}')
        res = parsePfamPDBs(val)
        os.chdir(cwd)

    elif fn == 'parsePfamFull':
        import gzip
        pfam_fpath = f'{database_folder}Pfam-A.full.gz'
        pfam = gzip.open(pfam_fpath, "rt")
        ids, seqs = read_pfam(pfam, val)
        res = {}
        for seq_id in ids:
            seq = seqs[seq_id]
            res.update({seq_id: seq})
            print(seq_id + "\t" + seq)
        # save alignment as MSA in fasta format

    return res


def main():
    query_dict = {
        # 'query_fn': 'searchPfam', 'query_val': '3luc', # 'PIWI_ARCFU',
        # 'query_fn': 'fetchPfamMSA', 'query_val': 'PF02171', #'PF00754' #
        'query_fn': 'parsePfamPDBs', 'query_val': 'PIWI_ARCFU',
        # 'query_fn': 'parsePfamFull', 'query_val': 'PF00754' #'PF16121' ## 'PF02171' # ## 'PF00041'
    }

    import os
    print(os.getcwd())
    res = pfam_access_functions(query_dict)
    for k, v in res.items():
        print(k, v)

if __name__=='__main__':
    main()