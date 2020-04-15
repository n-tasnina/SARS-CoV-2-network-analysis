import scipy.stats as stats
import pandas as pd
import argparse
import yaml
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from textwrap import wrap
from src.setup_datasets import parse_gmt_file
sys.path.insert(0, os.path.dirname(__file__))


def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)
    return config_map, kwargs


def setup_opts():
    # Parse command line args.
    parser = argparse.ArgumentParser(description="Script for Fisher's exact test")
    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="config.yaml",
                       help="Configuration file for this script.")
    # NURE: allow multiple values to this argument since we want to trace the p-value as we increase k.
    # We can discuss on Slack. The code to plot the p-values as in my 2011 paper on HIV dependency factors can also be part of this script.
    # group.add_argument('--k', type=int, help="top k predictions to consider")

    group.add_argument('--ks', type=list, default=[100, 2000, 100], help="value of k will be from list[0] to list[1] increasing by amount of list[2]")

    return parser


def parse_geneset(config_map):
    # this will filter out the positive Krogan gene sets
    # return* a mapping of gene_set_name to uniprotids of genes belonging to the protein set

    uniprotids, prot_descriptions = parse_gmt_file(config_map['reference_gene_set']['uniprotids'])
    for key in prot_descriptions:
        if "Krogan" in key:
            uniprotids.pop(key)
    return uniprotids


# def find_total_protein_in_network(config_map):
#
#     protein_space = {}
#     network_info_dir = config_map['network_info_dir']
#
#     for dirpath, dirs, files in os.walk(network_info_dir):
#
#         for filename in files:
#
#             fname = os.path.join(dirpath, filename)
#
#             if 'mapping-stats' in fname:
#
#                 network_info = pd.read_csv(fname, sep='\t')
#                 network_info.rename(columns={'Unnamed: 0': 'network_type'}, inplace=True)
#                 dirpath_split = dirpath.split('/')
#
#                 for i in range(len(network_info['network_type'])):
#
#                     network_source_name = dirpath_split[dirpath_split.__len__() - 1] + '_' + network_info.values[i][0]
#                     protein_space[network_source_name] = network_info.values[i][3]
#
#     return protein_space

def calc_Fisher_exact_test_score(geneset_protein_list, predicted_protein_list, k):

    # Calculate Fisher's exact coefficient
    # returns:
    #       1. oddratio
    #       2. pvalue
    #       3. 4 values from contingency table

    total_number_of_protein = len(predicted_protein_list)

    # get the list of top-k predictions. Check if there are at least k predictions.
    if k < len(predicted_protein_list):
        predicted_top_k_proteins = predicted_protein_list[:k]
    else:
        predicted_top_k_proteins = predicted_protein_list
        k = len(predicted_top_k_proteins)

    # NURE: Add a comment explaining why you are doing this.
    reduced_geneset_protein_list = list(set(geneset_protein_list).intersection(set(predicted_protein_list)))

    proteins_common_in_reduced_geneset_and_top_k_prediction = list(set(reduced_geneset_protein_list).intersection(set(predicted_top_k_proteins)))
    number_of_proteins_common_in_reduced_geneset_and_top_k_prediction = len(proteins_common_in_reduced_geneset_and_top_k_prediction)
    a = number_of_proteins_common_in_reduced_geneset_and_top_k_prediction
    b = len(reduced_geneset_protein_list) - a
    c = len(predicted_top_k_proteins) - a
    d = total_number_of_protein - a - b - c

    oddsratio, pvalue = stats.fisher_exact([[a, b], [c, d]],'greater')

    # fraction_of_predicted_postive_by_random_model = (a+c)/total_number_of_protein
    # if(a+b>0):
    #     fraction_of_predicted_postive_by_this_model = a/(a+b)
    # else:
    #     fraction_of_predicted_postive_by_this_model = 0
    #
    # difference_with_random_prediction = fraction_of_predicted_postive_by_this_model - fraction_of_predicted_postive_by_random_model
    return oddsratio, pvalue,a, b, c, d


def get_list_of_positive_proteins(config_map):
    # returns: the list of Krogan proteins that have been used as positive examples in prediction
    positive_proteins = pd.read_csv(config_map['positive_protein_file'], '\t')
    return list(positive_proteins['prots'])


def plot(rank, pvalue, geneset_name, algo_dirpath):

    plt.show()
    threshold = 0.05
    pvalue = np.asarray(pvalue)
    rank = np.asarray(rank)
    significant = np.ma.masked_where(pvalue <= threshold, pvalue)
    non_significant = np.ma.masked_where(pvalue > threshold, pvalue)

    fig, ax = plt.subplots()
    ax.plot(rank, significant, '--o', rank, non_significant,'--o')
    # plt.plot(rank, pvalue,'--o')
    plt.xlabel('rank')
    plt.ylabel('pvalue')
    plt.title("\n".join(wrap('overlap between ' + geneset_name + ' and proteins predicted by ' + os.path.basename(algo_dirpath))))
    plt.xticks(np.arange(0, max(rank) + 1, 200))


    geneset_name = geneset_name.replace("/", "_")
    # print(geneset_name)
    plt.savefig(algo_dirpath +"/Plots/"+geneset_name, format='png')

    # plt.show()
    # plt.savefig(algo_dirpath + "/" + '52: SARS coronavirus hypothetical protein sars9b from Virus-Host PPI P-HIPSTer 2020', format='png')
    # plt.savefig( algo_dirpath + "/" + '57: SARS coronavirus nsp7-pp1a/pp1ab (gene: orf1ab) from Virus-Host PPI P-HIPSTer 2020', ormat='png')


def main(config_map, **kwargs):

    # k = top k predictions to consider to find overlap
    # k = kwargs.get('k')
    k_list = []

    if (kwargs.get('ks')):
        ks = kwargs.get('ks')
        for k in range(ks[0], ks[1] + 1, ks[2]):
            k_list.append(k)

    geneset_uniprotids_dict = parse_geneset(config_map)

    # this is the Krogan protein list
    positive_protein_list = get_list_of_positive_proteins(config_map)

    predicted_prot_dir = config_map['predicted_prot_dir']

    # oddratio, pvalue, difference_with_random_prediction, fraction_of_predicted_postive_by_random_model, fraction_of_predicted_postive_by_this_model, a, b, c, d = \
    #     calc_Fisher_exact_test_score(['Q6UX04','8IWA5','O60885','O00203'], ['Q6UX04','8IWA5','O60885','O00203','P25440','Q86VM9','ABFSS',"JKJ",'ajisiwdhid'],3)
    # print(difference_with_random_prediction,fraction_of_predicted_postive_by_random_model,fraction_of_predicted_postive_by_this_model,pvalue)

    for dirpath, dirs, files in os.walk(predicted_prot_dir):

        for filename in files:

            fname = os.path.join(dirpath, filename)

            # this code assume that pred-scores contain score for all the proteins in the network.
            # have to consider removing this dependency if suggested
            if 'pred-scores' in fname and 'stats' not in fname:

                # filtering out the positive proteins from predicted list of proteins
                predicted_prot_info = pd.read_csv(fname, sep='\t')
                predicted_protein_list = predicted_prot_info[~predicted_prot_info['prot'].isin(positive_protein_list)]['prot']
                print(len(predicted_protein_list))

                Fishers_exact_test_output_file = dirpath + '/' + os.path.basename(dirpath)+'_Fishers_exact_test_score.csv'
                os.mkdir(dirpath + '/Plots')

                oddratio_list = []
                pvalue_list = []
                a_list=[]
                b_list=[]
                c_list=[]
                d_list=[]

                reference_prot_set_list = []
                rank = []

                for gene_set_name in geneset_uniprotids_dict:

                    for k in k_list:

                        oddratio, pvalue, a, b, c, d = \
                            calc_Fisher_exact_test_score(geneset_uniprotids_dict[gene_set_name], predicted_protein_list, k)

                        reference_prot_set_list.append(gene_set_name)
                        rank.append(k)
                        a_list.append(a)
                        b_list.append(b)
                        c_list.append(c)
                        d_list.append(d)
                        oddratio_list.append(oddratio)
                        pvalue_list.append(pvalue)

                    plot(rank[-len(k_list):], pvalue_list[-len(k_list):], gene_set_name, dirpath)

                print('writing to' + Fishers_exact_test_output_file)
                scores = pd.DataFrame({ 'protein_set_name': pd.Series(reference_prot_set_list),
                                        'rank':pd.Series(rank),
                                        'a': pd.Series(a_list),
                                        'b': pd.Series(b_list),
                                        'c': pd.Series(c_list),
                                        'd': pd.Series(d_list),
                                        # 'fraction_by_random_predictor': pd.Series(fraction_by_random_predictor_list),
                                        # 'fraction_by_this_predictor': pd.Series(fraction_by_this_predictor_list),
                                        # 'difference_with_random_prediction': pd.Series(difference_with_random_prediction_list),
                                        'pvalue': pd.Series(pvalue_list),
                                        'oddratio': pd.Series(oddratio_list)})
                scores.to_csv(Fishers_exact_test_output_file, '\t')


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    main(config_map, **kwargs)
