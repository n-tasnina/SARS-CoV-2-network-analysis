
import argparse
import yaml
import itertools
from collections import defaultdict
import os
import sys
from tqdm import tqdm, trange
import time
import numpy as np
from scipy import sparse as sp
# packages in this repo
# add this file's directory to the path so these imports work from anywhere
sys.path.insert(0,os.path.dirname(__file__))
#from src import setup_sparse_networks as setup
import src.setup_sparse_networks as setup
import src.algorithms.alg_utils as alg_utils
import src.algorithms.runner as runner
import src.evaluate.eval_utils as eval_utils
import src.evaluate.cross_validation as cross_validation
import src.evaluate.eval_leave_one_species_out as eval_loso
import src.utils.file_utils as utils
import src.utils.string_utils as string_utils
import src.utils.config_utils as config_utils


def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        #config_map = yaml.load(conf, Loader=yaml.FullLoader)
        config_map = yaml.load(conf)
    # TODO check to make sure the inputs are correct in config_map

    # for each argument not set in opts, remove it from kwargs.
    # That way the default will be used instead of None
    for key in list(kwargs.keys()):
        if kwargs[key] is None:
            del kwargs[key]

    return config_map, kwargs


def setup_opts():
    ## Parse command line args.
    parser = argparse.ArgumentParser()  #description='')

    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="config-files/config.yaml",
            help="Configuration file")
    group.add_argument('--term', '-T', type=str, action="append",
            help="Specify the terms to use. Can use this option multiple times")

    # evaluation parameters
    group = parser.add_argument_group("Evaluation options (can be set under 'eval_settings' in config file)")
    group.add_argument('--only-eval', action="store_true",
            help="Perform evaluation only (i.e., skip prediction mode)")
    group.add_argument('--cross-validation-folds', '-C', type=int,
            help="Perform cross validation using the specified # of folds. Usually 5")
    group.add_argument('--num-reps', type=int, 
            help="Number of times to repeat the CV process. Default=1")
    group.add_argument('--cv-seed', type=int, 
            help="Seed to use for the random number generator when splitting the annotations into folds. " + \
            "If --num-reps > 1, the seed will be incremented by 1 each time. Should only be used for testing purposes")
    group.add_argument('--sample-neg-examples-factor', type=float, 
            help="If specified, sample negative examples randomly without replacement from the protein universe equal to <sample_neg_examples_factor> * # positives")
    group.add_argument('--loso', action="store_true",
            help="Leave-One-Species-Out evaluation. For each species, leave out all of its annotations " +
            "and evaluate how well they can be recovered from the annotations of the other species. " +
            "Must specify the 'taxon_file' in the config file")
    group.add_argument('--taxon', '-S', dest="taxons", type=str, action='append',
            help="Specify the species taxonomy ID for which to evaluate. Multiple may be specified. Otherwise, all species will be used")
    group.add_argument('--write-prec-rec', action="store_true",
            help="Also write a file containing the precision and recall for every positive example. " + \
            "If a single term is given, only the prec-rec file, with the term in its name, will be written.")
    group.add_argument('--early-prec', '-E', type=str, action="append", default=["k1"],
            help="Report the precision at the specified recall value (between 0 and 1). " + \
            "If prefixed with 'k', for a given term, the precision at (k * # ann) # of nodes is given. Default: k1")

    # additional parameters
    group = parser.add_argument_group('Additional options')
    group.add_argument('--num-pred-to-write', '-W', type=int,
            help="Number of predictions to write to the file for predictions mode (meaning if --only-eval isn't specified). " + \
            "If 0, none will be written. If -1, all will be written. Default=10")
    group.add_argument('--factor-pred-to-write', '-N', type=float, 
            help="Write the predictions <factor>*num_pos for each term to file. " +
            "For example, if the factor is 2, a term with 5 annotations would get the nodes with the top 10 prediction scores written to file.")
    group.add_argument('--postfix', type=str, 
            help="String to add to the end of the output file name(s)")
    group.add_argument('--forcealg', action="store_true",
            help="Force re-running algorithms if the output files already exist")
    group.add_argument('--forcenet', action="store_true",
            help="Force re-building network matrix from scratch")
    group.add_argument('--verbose', action="store_true",
            help="Print additional info about running times and such")

    return parser


def run(config_map, **kwargs):
    input_settings, input_dir, output_dir, alg_settings, kwargs \
        = config_utils.setup_config_variables(config_map, **kwargs)

    for dataset in input_settings['datasets']:
        # add options specified for this dataset to kwargs  
        # youngs_neg: for a term t, a gene g cannot be a negative for t if g shares an annotation with any gene annotated to t 
        #kwargs['youngs_neg'] = dataset.get('youngs_neg') 

        net_obj, ann_obj, eval_ann_obj = setup_dataset(dataset, input_dir, alg_settings, **kwargs) 
        # if there are no annotations, then skip this dataset
        if len(ann_obj.terms) == 0:
            print("No terms found. Skipping this dataset")
            continue
        # the outputs will follow this structure:
        # outputs/<net_version>/<exp_name>/<alg_name>/output_files
        out_dir = "%s/%s/%s/" % (output_dir, dataset['net_version'], dataset['exp_name'])
        alg_runners = setup_runners(alg_settings, net_obj, ann_obj, out_dir, **kwargs)

        # first run prediction mode since it is the fastest
        if kwargs.get('only_eval',False) is False:
            # run algorithms in "prediction" mode 
            run_algs(alg_runners, **kwargs) 
            # if specified, write the SWSN combined network to a file
            save_net = dataset['net_settings'].get('save_net', None) if 'net_settings' in dataset else None
            if net_obj.weight_swsn is True and save_net is not None:
                out_file = "%s/%s/%s" % (input_dir, dataset['net_version'], save_net)
                # the SWSN network is part of the runner object. Need to organize that better
                net_obj.save_net(out_file)

            # if a pos_neg_file_eval was passed in (e.g., for temporal holdout validation),
            # use it to evaluate the predictions
            if eval_ann_obj is not None:
                exp_type = "eval"
                for run_obj in alg_runners:
                    out_file = "%s/%s-%s%s.txt" % (
                        run_obj.out_dir, exp_type, run_obj.params_str, kwargs.get("postfix", ""))
                    utils.checkDir(os.path.dirname(out_file))
                    eval_utils.evaluate_ground_truth(
                        run_obj, eval_ann_obj, out_file, **kwargs)

        if kwargs.get('cross_validation_folds') is not None:
            # run cross validation
            cross_validation.run_cv_all_terms(alg_runners, ann_obj, folds=kwargs['cross_validation_folds'], **kwargs)

        if kwargs.get('loso') is True:
            # add the taxon file paths for this dataset to kwargs
            for arg in ['taxon_file', 'only_taxon_file']:
                kwargs[arg] = "%s/%s" % (input_dir, dataset[arg]) 
            # now run the leave-one-species-out eval
            eval_loso.eval_loso(alg_runners, ann_obj, eval_ann_obj=eval_ann_obj, **kwargs)


def setup_net(input_dir, dataset, **kwargs):
    # load the network matrix and protein IDs
    net_files = None
    if 'net_files' in dataset:
        net_files = ["%s/%s/%s" % (input_dir, dataset['net_version'], net_file) for net_file in dataset['net_files']]
    unweighted = dataset['net_settings'].get('unweighted', False) if 'net_settings' in dataset else False
    # if multiple networks are passed in, then set multi_net to True automatically
    if (net_files is not None and len(net_files) > 1) or 'string_net_files' in dataset:
        if dataset.get('multi_net') is False:
            print("WARNING: multiple networks were passed in. Setting 'multi_net' to True")
        dataset['multi_net'] = True

    # parse and store the networks 
    if dataset.get('multi_net') is True: 
        # if multiple file names are passed in, then map each one of them
        if net_files is not None or 'string_net_files' in dataset:
            string_net_files = ["%s/%s/%s" % (input_dir, dataset['net_version'], string_net_file) for string_net_file in dataset['string_net_files']]
            string_nets = None 
            if 'string_nets' in dataset['net_settings']:
                string_nets = string_utils.convert_string_naming_scheme(dataset['net_settings']['string_nets'])
            # they all need to have the same rows and columns, which is handled by this function
            # this function also creates the multi net file if it doesn't exist
            string_cutoff = dataset['net_settings'].get('string_cutoff', 150) 
            out_pref = "%s/sparse-nets/c%d-" % (os.path.dirname(string_net_files[0]), string_cutoff)
            utils.checkDir(os.path.dirname(out_pref))
            sparse_nets, net_names, prots = setup.create_sparse_net_file(
                    out_pref, net_files=net_files, string_net_files=string_net_files, 
                    string_nets=string_nets,
                    string_cutoff=string_cutoff,
                    forcenet=kwargs.get('forcenet',False))
        else:
            # if a .mat file with multiple sparse matrix networks inside of it is passed in, read that here
            net_names_file = "%s/%s/%s" % (input_dir, dataset['net_version'], dataset['net_settings']['net_names_file'])
            node_ids_file  = "%s/%s/%s" % (input_dir, dataset['net_version'], dataset['net_settings']['node_ids_file'])
            sparse_nets, net_names, prots = alg_utils.read_multi_net_file(net_file, net_names_file, node_ids_file)

        weight_method = dataset['net_settings']['weight_method'].lower()
        net_obj = setup.Sparse_Networks(
            sparse_nets, prots, net_names=net_names, weight_method=weight_method,
            unweighted=unweighted, verbose=kwargs.get('verbose',False)
        )
    else:
        if net_files is None:
            print("ERROR: no net files specified in the config file. Must provide either 'net_files', or 'string_net_files'")
            sys.exit()
        W, prots = alg_utils.setup_sparse_network(net_files[0], forced=kwargs.get('forcenet',False))
        net_obj = setup.Sparse_Networks(
            W, prots, unweighted=unweighted, verbose=kwargs.get('verbose',False))
    return net_obj


def load_annotations(prots, dataset, input_dir, **kwargs):
    # limit the terms to whatever is specified either in the only_functions_file,
    # or the --term command-line option
    only_functions_file = None
    # TODO add a 'test' option or something to be able to limit to a single term
    if 'only_functions_file' in dataset and dataset['only_functions_file'] != '':
        only_functions_file = "%s/%s" % (input_dir, dataset['only_functions_file'])
    selected_terms = alg_utils.select_terms(only_functions_file=only_functions_file) 

    # write/load the processed annotation matrix in a pos_neg_file version of the network folder
    # TODO this is hacky, but works for now
    pos_neg_str = dataset['pos_neg_file'] \
                       .replace('/','-').replace('pos-neg-','') \
                       .replace('.txt','').replace('.tsv','').replace('.gz','')
    pos_neg_str += '-Yneg' if kwargs.get('youngs_neg') else ''
    net_files = dataset['net_files'] if 'net_files' in dataset else dataset['string_net_files']
    sparse_ann_file = "%s/%s/sparse-anns/%s-%s.npz" % (
        input_dir, dataset['net_version'], net_files[0].split('.')[0], pos_neg_str)
    # now build the annotation matrix
    pos_neg_file = "%s/%s" % (input_dir, dataset['pos_neg_file'])
    ann_obj = setup.create_sparse_ann_and_align_to_net(
            pos_neg_file, sparse_ann_file, prots, **kwargs)

    if selected_terms is not None:
        ann_obj.limit_to_terms(selected_terms)
    else:
        selected_terms = ann_obj.terms

    eval_ann_obj = None
    # also check if an evaluation pos_neg_file was given
    if dataset.get('pos_neg_file_eval', '') != '':
        pos_neg_str = dataset['pos_neg_file_eval'] \
                        .replace('/','-').replace('pos-neg-','') \
                        .replace('.txt','').replace('.tsv','').replace('.gz','')
        pos_neg_str += '-Yneg' if kwargs.get('youngs_neg') else ''
        sparse_ann_file = "%s/%s/sparse-anns/%s.npz" % (input_dir, dataset['net_version'], pos_neg_str)
        pos_neg_file_eval = "%s/%s" % (input_dir, dataset['pos_neg_file_eval'])
        eval_ann_obj = setup.create_sparse_ann_and_align_to_net(
                pos_neg_file_eval, sparse_ann_file, prots, **kwargs)
        # also limit the terms in the eval_ann_obj to those from the pos_neg_file
        eval_ann_obj.limit_to_terms(selected_terms)
    ann_obj.selected_terms = selected_terms
    return selected_terms, ann_obj, eval_ann_obj


def get_algs_to_run(alg_settings):
    # these are the algs to run
    algs = []
    for alg in alg_settings:
        if alg_settings[alg].get('should_run', [True])[0] is True:
            algs.append(alg)
    return algs


def setup_dataset(dataset, input_dir, alg_settings, **kwargs):
    if kwargs.get('verbose'):
        utils.print_memory_usage()
    # setup the network matrix first
    net_obj = setup_net(input_dir, dataset, **kwargs)
    if kwargs.get('verbose'):
        utils.print_memory_usage()
    #ann_obj = setup_annotations(input_dir, dataset, **kwargs)
    selected_terms, ann_obj, eval_ann_obj = load_annotations(net_obj.nodes, dataset, input_dir, **kwargs)
    if kwargs.get('verbose'):
        utils.print_memory_usage()

    #algs = get_algs_to_run(alg_settings)
    return net_obj, ann_obj, eval_ann_obj


def setup_runners(alg_settings, net_obj, ann_obj, out_dir, **kwargs):
    # these are the algs with 'should_run' set to [True]
    algs = get_algs_to_run(alg_settings)
    print("Algs to run: %s" % (', '.join(algs)))

    alg_runners = []
    for alg in algs:
        params = alg_settings[alg]
        combos = [dict(zip(params, val))
            for val in itertools.product(
                *(params[param] for param in params))]
        for combo in combos:
            run_obj = runner.Runner(alg, net_obj, ann_obj, out_dir, combo, **kwargs)
            alg_runners.append(run_obj) 

    print("\t%d total runners" % (len(alg_runners)))
    return alg_runners


def run_algs(alg_runners, **kwargs):
    """
    Runs all of the specified algorithms with the given network and annotations.
    Each runner should return the GO term prediction scores for each node in a sparse matrix.
    """
    # if no negative examples are given and an algorithm needs negative examples,
    # sample negative examples 
    neg_factor = kwargs.get('sample_neg_examples_factor')
    num_reps = kwargs.get('num_reps', 1)
    # first check to see if the algorithms have already been run
    # and if the results should be overwritten
    for run_obj in alg_runners:
        out_file = run_obj.out_pref + ".txt"
        run_obj.out_file = out_file
    if kwargs.get('forcealg') is True or kwargs.get('num_pred_to_write') == 0:
        runners_to_run = alg_runners
    else:
        runners_to_run = []
        for run_obj in alg_runners:
            if os.path.isfile(run_obj.out_file):
                print("%s already exists. Use --forcealg to overwite" % (run_obj.out_file))
            else:
                runners_to_run.append(run_obj)

    params_results = {}

    print("Generating inputs and running the algorithms")
    # now setup the inputs for the runners
    for run_obj in runners_to_run:
        # TODO make a better way of indicating an alg needs negative examples than just having 'plus' in the name
        if 'plus' not in run_obj.name and neg_factor is not None:
            orig_ann_obj = run_obj.ann_obj
            # sample negative examples, repeat the method the given number of times,
            # and then average the resulting scores
            avg_term_scores = sp.lil_matrix(orig_ann_obj.ann_matrix.shape, dtype=np.float)
            print("Sampling negative examples and averaging the scores over %d runs of %s" % (num_reps, run_obj.name))
            for rep in trange(1,num_reps+1):
                new_ann_matrix = eval_utils.sample_neg_examples(
                    orig_ann_obj, sample_neg_examples_factor=neg_factor)
                run_obj.ann_obj.ann_matrix = new_ann_matrix
                run_obj.setupInputs()
                run_obj.run()
                avg_term_scores += run_obj.term_scores 
            # finally, divide by the number of repetitions to get the average
            avg_term_scores = avg_term_scores / num_reps 
            run_obj.term_scores = avg_term_scores
        else:
            # run the method like normal
            run_obj.setupInputs()
            # TODO storing all of the runners scores simultaneously could be costly (too much RAM).
            run_obj.run()
            params_results.update(run_obj.params_results)

    # parse the outputs. Only needed for the algs that write output files
    for run_obj in runners_to_run:
        run_obj.setupOutputs()

        # write to file if specified
        num_pred_to_write = kwargs.get('num_pred_to_write',10)
        if kwargs.get('factor_pred_to_write') is not None:
            # make a dictionary with the # ann*factor for each term
            num_pred_to_write = {} 
            for i in range(run_obj.ann_matrix.shape[0]):
                y = run_obj.ann_matrix[i,:]
                positives = (y > 0).nonzero()[1]
                num_pred_to_write[run_obj.terms[i]] = len(positives) * kwargs['factor_pred_to_write']
        if num_pred_to_write != 0:
            utils.checkDir(os.path.dirname(run_obj.out_file)) 
            alg_utils.write_output(run_obj.term_scores, run_obj.ann_obj.terms, run_obj.ann_obj.prots,
                         run_obj.out_file, num_pred_to_write=num_pred_to_write)

    #eval_loso.write_stats_file(runners_to_run, params_results)
    print(params_results)
    print("Finished")


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    run(config_map, **kwargs)
