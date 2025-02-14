
import argparse
import yaml
from collections import defaultdict
import os
import sys
#from tqdm import tqdm
import copy
import time
import itertools
import numpy as np
#from scipy import sparse
import pandas as pd
# add this file's directory to the path so these imports work from anywhere
sys.path.insert(0,os.path.dirname(__file__))
#from FastSinkSource.run_eval_algs import setup_runners
from FastSinkSource.src.plot import plot_utils
from FastSinkSource.src.utils import config_utils
from FastSinkSource.src.algorithms import runner
from setup_datasets import parse_mapping_file


def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        #config_map = yaml.load(conf, Loader=yaml.FullLoader)
        config_map = yaml.load(conf)
    # TODO check to make sure the inputs are correct in config_map

    return config_map, kwargs


def setup_opts():
    ## Parse command line args.
    parser = argparse.ArgumentParser(description="Script to pull together predictions from multiple datasets and/or algorithms, and sort them by med rank")

    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="fss_inputs/config_files/stringv11/400.yaml",
                       help="Configuration file used when running FastSinkSource")
    group.add_argument('--alg', '-A', dest='algs', type=str, action="append",
                       help="Algorithms for which to get results. Must be in the config file. " +
                       "If not specified, will get the list of algs with 'should_run' set to True in the config file")
    group.add_argument('--id-mapping-file', type=str, default="datasets/mappings/human/uniprot-reviewed-status.tab.gz",
                       help="Table downloaded from UniProt to map to gene names. Expected columns: 'Entry' and 'Gene names'")
    group.add_argument('--num-pred-to-write', '-W', type=int, default=100,
            help="Number of predictions to keep. " +
            "If 0, none will be written. If -1, all will be written. Default=100")
    group.add_argument('--factor-pred-to-write', '-N', type=float, 
            help="Keep the predictions <factor>*num_pos. " +
            "For example, if the factor is 2, a term with 5 annotations would get the nodes with the top 10 prediction scores written to file.")
    group.add_argument('--out-pref', type=str, default="outputs/stats/pred-table-",
                       help="Output prefix for writing the table of predictions. Default=outputs/stats/pred-table-")
    group.add_argument('--round', type=int, default=3,
                       help="Round to the given number of decimal places in the output file")
    #group.add_argument('--download-only', action='store_true', default=False,
    #                   help="Stop once files are downloaded and mapped to UniProt IDs.")

#    # additional parameters
#    group = parser.add_argument_group('Additional options')
#    group.add_argument('--forcealg', action="store_true", default=False,
#            help="Force re-running algorithms if the output files already exist")
#    group.add_argument('--forcenet', action="store_true", default=False,
#            help="Force re-building network matrix from scratch")
#    group.add_argument('--verbose', action="store_true", default=False,
#            help="Print additional info about running times and such")

    return parser


def main(config_map, **kwargs):
    """
    *config_map*: everything in the config file
    *kwargs*: all of the options passed into the script
    """
    
    input_settings = config_map['input_settings']
    input_dir = input_settings['input_dir']
    alg_settings = config_map['algs']
    output_settings = config_map['output_settings']
    if config_map.get('eval_settings'):
        kwargs.update(config_map['eval_settings'])
    postfix = kwargs.get("postfix")
    postfix = "" if postfix is None else postfix
    algs = plot_utils.get_algs_to_run(alg_settings, **kwargs)
    del kwargs['algs']  # remove algs from kwargs
    uniprot_to_gene = None
    # also add the protein name
    uniprot_to_prot_names = None
    if kwargs.get('id_mapping_file'):
        df = pd.read_csv(kwargs['id_mapping_file'], sep='\t', header=0) 
        ## keep only the first gene for each UniProt ID
        uniprot_to_gene = {p: genes.split(' ')[0] for p, genes in zip(df['Entry'], df['Gene names'].astype(str))}
        if 'Protein names' in df.columns:
            uniprot_to_prot_names = dict(zip(df['Entry'], df['Protein names'].astype(str)))

    # for each dataset, get the prediction scores from each method
    # store a single df for each alg
    df_all = pd.DataFrame()
    for dataset in input_settings['datasets']:
        # load the positive examples and remove them from the predictions
        pos_neg_file = "%s/%s" % (input_dir, dataset['pos_neg_file'])
        df = pd.read_csv(pos_neg_file, sep='\t')
        orig_pos = df[df['2020-03-sarscov2-human-ppi'] == 1]['prots']

        dataset_name = config_utils.get_dataset_name(dataset) 
        alg_pred_files = config_utils.get_dataset_alg_prediction_files(
            output_settings['output_dir'], dataset, alg_settings, algs, **kwargs)
        for alg, pred_file in alg_pred_files.items():
            if not os.path.isfile(pred_file):
                print("Warning: %s not found. skipping" % (pred_file))
                continue
            print("reading: %s" % (pred_file))
            df = pd.read_csv(pred_file, sep='\t')
            # remove the original positives
            df = df[~df['prot'].isin(orig_pos)]
            df.reset_index(inplace=True, drop=True)

            df = df[['prot', 'score']]
            # reset the index again to store the current rank as a column
            df.reset_index(inplace=True)
            df.columns = ['Rank', 'Prot', 'Score']
            df['Rank'] += 1

            # now set the index as the uniprot ID 
            df.set_index('Prot', inplace=True)

            # add the dataset-specific settings to the column headers
            #df.columns = ["%s-%s" % (curr_name, col) for col in ["Rank", "Score"]]
            # UPDATE: add levels to the column index
            tuples = [(dataset_name, alg, col) for col in df.columns]
            index = pd.MultiIndex.from_tuples(tuples)
            df.columns = index
            print(df.head())

            # now add the columns to the master dataframe
            df_all = pd.concat([df_all, df], axis=1)
            #for col in df.columns:
            #    alg_dfs[alg][col] = df[col]

    print(df_all.head())

    # sort the dataframe by the average rank
    # first compute the ranks, then add them to the dataframe later
    print("Sorting genes by average rank")
    # for each dataset, replace NA with the highest rank in that network
    df_ranks = df_all[[col for col in df_all.columns if 'Rank' in col]]
    df_ranks.fillna(df_ranks.max(axis=0), inplace=True)
    # not exactly sure what the "bottom" option does from pandas
    #           .rank(method='average', na_option='bottom')

    #nan_rank_val = len(orig_pos)*2
    #print("\tempty ranks getting a value of %d" % (nan_rank_val))
    #df_ranks = df_all[[col for col in df_all.columns if 'Rank' in col]].fillna(nan_rank_val)
    #print(df_ranks.head())
    med_rank = df_ranks.median(axis=1).sort_values()
    med_rank = med_rank.round(1)
    #print(med_rank)
    #print(len(med_rank))

    # figure out how many predictions to keep
    num_pred_to_write = kwargs.get('num_pred_to_write', 100) 
    if kwargs.get('factor_pred_to_write') is not None:
        num_pred_to_write = kwargs['factor_pred_to_write']*len(orig_pos)

    if num_pred_to_write == -1:
        print("Keeping all predictions from each method")
    else:
        print('Keeping %d predictions (with ties) from each alg' % (num_pred_to_write))
        # for each dataset, keep the top k predictions, with ties
        for dataset in input_settings['datasets']:
            name = dataset['plot_exp_name']
            curr_algs = df_all[name].columns.levels[0]
            for alg in curr_algs:
                sub_df = df_all[(name, alg)]
                sub_df.sort_values(by='Rank', inplace=True)

                topk = int(min(len(sub_df), num_pred_to_write))
                score_topk = sub_df.iloc[topk-1]['Score']
                sub_df = sub_df[sub_df['Score'] >= score_topk]
                # TODO for some reason this doesn't work if that are nan or inf in the column.
                sub_df['Score'] = sub_df['Score'].round(kwargs.get('round', 3))
                # use object to keep the value as an integer in the output
                sub_df['Rank'] = sub_df['Rank'].astype(int).astype(object)
                df_all[[(name, alg, col) for col in ['Rank', 'Score']]] = sub_df

    # now remove rows that are NA for all values
    df_all.dropna(how='all', inplace=True)

    # and the average rank
    df_all.insert(0, 'MedianRank', med_rank)
    df_all.sort_values(by='MedianRank', inplace=True)
    # add the protein names
    if uniprot_to_prot_names is not None:
        df_all.insert(0, 'ProteinNames', df_all.index.map(uniprot_to_prot_names))
    # add the gene names
    if uniprot_to_gene is not None:
        df_all.insert(0, 'GeneName', df_all.index.map(uniprot_to_gene))

    print(df_all.head())

    # make sure the output directory exists
    os.makedirs(os.path.dirname(kwargs['out_pref']), exist_ok=True)
    #for alg, df in alg_dfs.items():
    out_file = "%s%s.tsv" % (kwargs['out_pref'], '-'.join(algs))
    print("\nWriting %s" % (out_file))
    df_all.to_csv(out_file, sep='\t')


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    main(config_map, **kwargs)

