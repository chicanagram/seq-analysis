import pandas as pd
pd.set_option('display.max_columns', None)
from variables import address_dict, subfolders

if __name__=='__main__':
    # Input and output file paths
    data_folder = address_dict['ECOHARVEST']
    data_subfolder = 'sidestream_cocktail'
    seqsearch_dir = data_folder + subfolders['seqsearch'] + data_subfolder + '/'
    sequences_dir = data_folder + subfolders['sequences'] + data_subfolder + '/'
    query_seq_fname = 'queryseqs_enzcocktail.fasta' # 'XYNA_ASPNG.fasta' # 'MmCAR.fasta'
    csv_fpath = seqsearch_dir + data_subfolder + '.csv'


    # define enztype to seqname mapping
    enztype_to_seqname = {
        ('xylanase', 'Aspergillus niger'): ['sp|P55329|XYNA_ASPNG', ],
        ('endoglucanase', 'Trichoderma reesei'): ['sp|O14405|LP9A_HYPJE', 'sp|A0A024SH20|GUN2_HYPJR', 'sp|A0A024SNB7|GUN1_HYPJR', 'tr|A0A024RXC1|A0A024RXC1_HYPJR', 'tr|A0A024SNE0|A0A024SNE0_HYPJR'],
        ('exoglucanase', 'Trichoderma reesei'): ['sp|A0A024RXP8|GUX1_HYPJR', 'sp|A0A024SH76|GUX2_HYPJR'],
        ('endoglucanase', 'Trichoderma harzianum'): ['tr|A0A0F9WYH5|A0A0F9WYH5_TRIHA', 'sp|P53626|E13B_TRIHA', 'tr|A0A0F9X8Y8|A0A0F9X8Y8_TRIHA', 'tr|A0A0F9ZZ66|A0A0F9ZZ66_TRIHA', 'tr|A0A0F9XWP3|A0A0F9XWP3_TRIHA'],
        ('exoglucanase', 'Trichoderma harzianum'): ['sp|Q9P8P3|GUX1_TRIHA'],
    }
    # reverse the mapping
    seqname_to_enzytype = {}
    for enztype, seqname_list in enztype_to_seqname.items():
        for seqname in seqname_list:
            seqname_to_enzytype[seqname] = enztype
            print(seqname, enztype)

    # load csv
    df = pd.read_csv(csv_fpath, index_col=0)
    # drop duplicate target sequences
    print(len(df))
    df = df.drop_duplicates(subset=['Target Organism', 'Target Seq'], ignore_index=True)
    print(len(df))

    # add query enzyme class and query organism labels
    if 'Query Enz Class' not in df and 'Query Organism' not in df:
        seq_name_list = df['Query Seq Name'].tolist()
        query_enz_class_list = []
        query_organism_list = []
        for seq_name in seq_name_list:
            query_enz_class_list.append(seqname_to_enzytype[seq_name][0])
            query_organism_list.append(seqname_to_enzytype[seq_name][1])
        df.insert(2, 'Query Enz Class', query_enz_class_list)
        df.insert(3, 'Query Organism', query_organism_list)
        df.insert(4, 'Query_EnzClass_Organism', [f'{query_enz_class}_{query_organism}' for query_enz_class, query_organism in zip(query_enz_class_list, query_organism_list)])
    # sort data
    df = df.sort_values(by=['Target Organism', 'Query_EnzClass_Organism', 'Identities'], ascending=[True,True,False], ignore_index=True)

    # get data stats
    df_stats = []
    stats_to_obtain = ['Target Organism', 'Query Enz Class', 'Seq Name', 'Query_EnzClass_Organism', 'Query Seq Name']
    target_organism_list = list(set(df['Target Organism'].tolist()))
    target_organism_list.sort()
    print('# of unique target organisms:', len(target_organism_list))

    for target_organism in target_organism_list:
        df_targetOrganism = df[df['Target Organism'] == target_organism]
        df_stats_targetOrganism = {'Target Organism': target_organism}
        for statscol in stats_to_obtain[1:]:
            lst_deduped_statscol = list(set(df_targetOrganism[statscol].tolist()))
            df_stats_targetOrganism.update({
                'Num Hits '+statscol: len(lst_deduped_statscol),
                'Hits '+statscol: ', '.join(lst_deduped_statscol)
            })
        df_stats.append(df_stats_targetOrganism)

    df_stats = pd.DataFrame(df_stats)
    df_stats = df_stats.sort_values(by=['Num Hits Query Enz Class', 'Num Hits Seq Name', 'Num Hits Query_EnzClass_Organism', 'Target Organism'], ascending=[False, False, False, True], ignore_index=True)
    df_stats.to_csv(csv_fpath.replace('.csv','_stats.csv'))
    print(df_stats)

    # sort data
    target_organism_sorted = df_stats['Target Organism'].tolist()
    for i, target_organism in enumerate(target_organism_sorted):
        df_targetOrganism = df[df['Target Organism']==target_organism]
        if i==0:
            df_sorted = df_targetOrganism.copy()
        else:
            df_sorted = pd.concat([df_sorted, df_targetOrganism], ignore_index=True)
    # save final sorted data
    df_sorted.to_csv(csv_fpath.replace('.csv', '_sorted.csv'))
    print(df)



