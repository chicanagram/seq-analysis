#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  7 07:58:44 2024

@author: charmainechia
"""
import pyhmmer
import pandas as pd
pd.set_option('display.max_columns', None)
import re
import matplotlib.pyplot as plt
from variables import address_dict, subfolders
data_folder = address_dict['examples']
database_folder = address_dict['databases']
hmm_subfolder = subfolders['hmm']

def parse_database(database_name='swissprot.fasta', database_folder=database_folder, convert_to_dataframe=False):
    # build index from downloaded sequences
    taxonomy_by_name = {}
    organism_by_name = {}
    description_by_name = {}
    if convert_to_dataframe:
        df = []
    # iterate through database, parse taxonomy / organism / description
    with pyhmmer.easel.SequenceFile(f'{database_folder}{database_name}', digital=True) as seq_file:
        database = seq_file.read_block()
    for seq in database:
        match = re.search(b"(.*) OS=(.*) OX=(\d+)", seq.description)
        taxonomy = int(match.group(3))
        organism = match.group(2)
        description = match.group(1)
        taxonomy_by_name[seq.name] = taxonomy
        organism_by_name[seq.name] = organism
        description_by_name[seq.name] = description
        if convert_to_dataframe:
            df.append({'name': seq.name.decode('utf-8'), 'description':description.decode('utf-8'), 'taxonomy':taxonomy, 'organism':organism.decode('utf-8')})
    print('# of sequences in database:', len(database))
    if convert_to_dataframe:
        df = pd.DataFrame(df)
    else:
        df = None
    return database, taxonomy_by_name, organism_by_name, description_by_name, df
def find_database_match(df, match_by_dict, print_df_match=False):
    for k, v in match_by_dict.items():
        if v is not None:
            df = df[df[k].str.contains(v, case=False)]
    seq_names_matched = df.name.tolist()
    print(*enumerate(seq_names_matched))
    if print_df_match:
        print(df)
    return df, seq_names_matched

def get_hmm_cutoff(hits, hmm):
    from math import floor, ceil
    # noise cutoff
    noise_score = max(hit.score for hit in hits if not hit.included)
    noise_domscore = max(hit.best_domain.score for hit in hits if not hit.included)
    hmm.cutoffs.noise = ceil(noise_score), ceil(noise_domscore)
    # gathering cutoff
    gathering_score = min(hit.score for hit in hits if hit.included)
    gathering_domscore = min(hit.best_domain.score for hit in hits if hit.included)
    hmm.cutoffs.gathering = floor(gathering_score), floor(gathering_domscore)
    # trusted cutoff
    diff = [
        hits[i-1].score - hits[i].score
        for i in range(1, len(hits))
        if hits[i-1].included
    ]
    argmax = max(range(len(diff)), key=diff.__getitem__)
    hit = hits[argmax]
    hmm.cutoffs.trusted = floor(hit.score), floor(hit.best_domain.score)
    # record results
    cutoff_dict = {
        'noise': hmm.cutoffs.noise,
        'gathering': hmm.cutoffs.gathering,
        'trusted': hmm.cutoffs.trusted
    }
    return cutoff_dict

def run_simple_search(ref_seq, database, plot_score_distribution=True, filter_hits=100, print_filtered_hits=False):
    """
    Perform a simple search to find close homologs without losing specificity
    """
    # find close homologs without losing specificity
    alphabet = pyhmmer.easel.Alphabet.amino()
    pli = pyhmmer.plan7.Pipeline(alphabet)
    hits = pli.search_seq(ref_seq, database)
    print('# of hits:', len(hits))
    # visualize score distribution to determine inclusion threshold
    if plot_score_distribution:
        plt.plot([hit.score for hit in hits], "o-")
        plt.xlabel("Hit rank")
        plt.ylabel("Score (bits)")
        plt.show(block=False)
        plt.pause(0.001)
    # filter hits
    if filter_hits is not None:
        hits_filtered = []
        for hit in hits:
            if hit.score >= filter_hits:
                hits_filtered.append(hit)
                if print_filtered_hits:
                    print(
                        f"{hit.name.decode():<20}"
                        f"\t{hit.score:5.1f}"
                        f"\t{description_by_name[hit.name].decode():34}"
                        f"\t{organism_by_name[hit.name].decode()}"
                    )
        print('# of hits (filtered):', len(hits_filtered))
    else:
        hits_filtered = hits
    return hits_filtered
def run_iterative_hmm_search(ref_seq, database, bitscore_cutoff=250, filter_by_taxonomy=None, plot_score_distribution=True, print_included_hits=True, print_excluded_hits=False, save_selected_seq_to_fasta=False, save_hmm_fname=None):
    # initialize search pipeline
    alphabet = pyhmmer.easel.Alphabet.amino()
    pli = pyhmmer.plan7.Pipeline(alphabet, incT=bitscore_cutoff, incdomT=bitscore_cutoff)

    if filter_by_taxonomy is not None:
        # get taxonomy database
        import taxopy
        taxdb = taxopy.TaxDb(keep_files=False)
        # definite taxonomy filter function
        def taxonomy_filter_function(top_hits):
            for hit in top_hits:
                taxid = taxonomy_by_name[hit.name]
                taxon = taxopy.Taxon(taxid, taxdb)
                if hit.included and taxon.rank_name_dictionary.get('class') != filter_by_taxonomy: # e.g. filter_by_taxonomy='Halobacteria'
                    hit.dropped = True
        # initialize search class
        search = pli.iterate_seq(ref_seq, database, select_hits=taxonomy_filter_function)

    else:
        search = pli.iterate_seq(ref_seq, database)

    # iteratively search database
    iterations = []
    while not search.converged:
        iteration = next(search)
        iterations.append(iteration)
        print(
            f"Iteration={iteration.iteration} Hits={len(iteration.hits):2} Included={len(iteration.hits.included):2} Converged={search.converged}")

    # visualize score distribution on each iteration
    if plot_score_distribution:
        for iteration in iterations:
            plt.plot([hit.score for hit in iteration.hits], "o-", label=f"Iteration {iteration.iteration}")
        plt.legend()
        plt.xlabel("Hit rank")
        plt.ylabel("Score (bits)")
        plt.show(block=False)
        plt.pause(0.001)

    # get final hmm
    hits = iterations[-1].hits
    hmm = iterations[-1].hmm
    # get cutoffs
    cutoffs_dict = get_hmm_cutoff(hits, hmm)
    print('Cutoffs:', *cutoffs_dict.items())

    # print included hits (final iteration)
    print('# of hits (final iteration):', len(iterations[-1].hits.included))
    if print_included_hits:
        print('HITS:')
        for hit in iterations[-1].hits:
            if hit.included:
                print(
                    f"{hit.name.decode():<20}"
                    f"\t{hit.score:5.1f}"
                    f"\t{description_by_name[hit.name].decode():34}"
                    f"\t{organism_by_name[hit.name][:60].decode()}"
                )
    # print hits excluded (final iteration)
    if print_excluded_hits:
        print()
        print('EXCLUDED HITS:')
        for hit in iterations[-1].hits:
            if not hit.included and hit.score > 50.0:
                print(
                    f"{hit.name.decode():<20}"
                    f"\t{hit.score:5.1f}"
                    f"\t{description_by_name[hit.name].decode():34}"
                    f"\t{organism_by_name[hit.name][:60].decode()}"
                )

    # save the HMM
    if save_hmm_fname is not None:
        hmm.name = save_hmm_fname.encode('utf-8')
        hmm.description = None
        hmm.accession = None
        with open(f"{data_folder}{hmm_subfolder}{save_hmm_fname}.hmm", "wb") as f:
            hmm.write(f)

    return iterations, hmm


# parse database
database, taxonomy_by_name, organism_by_name, description_by_name, df = parse_database(database_name='swissprot.fasta', convert_to_dataframe=True, database_folder=database_folder)

# find ref_seq by matching key words in name, description, organism, taxonomy
match_by_dict = {
    'name': 'B0R2U4', #None, #
    'description': None, # 'peroxidase',
    'organism': None,
    'taxonomy': None,
}
_, seq_names_matched = find_database_match(df, match_by_dict, print_df_match=False)

# get target ref sequence object in database
if len(seq_names_matched)>0:
    ref_seq_name_str = seq_names_matched[0]
    ref_seq_name = ref_seq_name_str.encode('utf-8')
    ref_seq = next(seq for seq in database if ref_seq_name in seq.name)
else: print('No matches found.')

# # perform a simple search to find close homologs without losing specificity --> use this to determine threshold
# hits = run_simple_search(ref_seq, database, plot_score_distribution=True, filter_hits=100, print_filtered_hits=True)

# perform iterative HMM search
iterations, hmm = run_iterative_hmm_search(
    ref_seq, database,
    bitscore_cutoff=100, filter_by_taxonomy='Halobacteria',
    plot_score_distribution=True, print_included_hits=True, print_excluded_hits=True,
    save_hmm_fname=match_by_dict['name'], save_selected_seq_to_fasta=False
)

# EXTRACT & SAVE RESULTS

# show plots
plt.show()

def main():
    return None

if __name__ == "__main__":
    main()