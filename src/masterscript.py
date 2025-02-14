
import argparse
import yaml
from collections import defaultdict
import os
import sys
#from tqdm import tqdm
import copy
import time
#import numpy as np
#from scipy import sparse
import pandas as pd
#import subprocess
# packages in this repo
# add this file's directory to the path so these imports work from anywhere
sys.path.insert(0,os.path.dirname(__file__))
from setup_datasets import setup_dataset_files, run_command


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
    parser = argparse.ArgumentParser(description="Script to download and parse input files, and (TODO) run the FastSinkSource pipeline using them.")

    # general parameters
    group = parser.add_argument_group('Main Options')
    group.add_argument('--config', type=str, default="config-files/master-config.yaml",
                       help="Configuration file for this script.")
    group.add_argument('--download-only', action='store_true', default=False,
                       help="Stop once files are downloaded and mapped to UniProt IDs.")
    group.add_argument('--force-download', action='store_true', default=False,
                       help="Force re-downloading and parsing of the input files")

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
    dataset_settings = config_map['dataset_settings']

    datasets_dir = dataset_settings['datasets_dir']

    # Download and parse the ID mapping files 
    # Download, parse, and map (to uniprot) the network files 
    setup_dataset_files(datasets_dir, dataset_settings['datasets_to_download'], dataset_settings.get('mappings'), **kwargs)

    if kwargs.get('download_only'):
        return

    # Now setup the config file to run the FastSinkSource pipeline using the specified networks
    # For now I will assume some of the essential structure is already setup (i.e., committed to the repo) 
    # TODO setup the pos-neg file(s) automatically (the file containing the positive and negative examples for which to use when running FSS)
    fss_settings = config_map['fastsinksource_pipeline_settings']
    config_files = setup_fss_config(datasets_dir, fss_settings) 

    # Now run FSS on each of the config files
    # TODO allow for parallelization of multiple networks / algorithms
    print("\nRunning the FastSinkSource pipeline on %d config files" % (len(config_files)))
    for config_file in config_files:
        log_file = config_file.replace('.yaml', '.log')
        command = "python -u src/FastSinkSource/run_eval_algs.py "  + \
                  " --config %s " % (config_file) + \
                  " >> %s 2>&1 " % (log_file)
        run_command(command) 


def setup_fss_config(datasets_dir, fss_settings):
    """
    Setup the config file(s) to run the FastSinkSource pipeline

    *returns*: a list of config files that are ready to be used to run FSS
    """
    fss_dir = fss_settings['input_dir']
    config_dir = "%s/config_files" % (fss_dir)
    # start setting up the config map settings that will remain the same for all datasets
    config_map = {
        'input_settings': {'input_dir': fss_dir,
            'datasets': [],},
        'output_settings': {'output_dir': fss_settings['output_dir']},
        'algs': fss_settings['algs'],
    }
    eval_str = ""
    if fss_settings.get('eval_settings'):
        eval_s = fss_settings.get('eval_settings')
        config_map['eval_settings'] = eval_s
        # append a string of the evaluation settings specified to the yaml file so you can run multiple in parallel.
        # TODO also add a postfix parameter
        eval_str = "%s%s%s%s" % (
            "-cv%s" % eval_s['cross_validation_folds'] if 'cross_validation_folds' in eval_s else "",
            "-nf%s" % eval_s['sample_neg_examples_factor'] if 'sample_neg_examples_factor' in eval_s else "",
            "-nr%s" % eval_s['num_reps'] if 'num_reps' in eval_s else "",
            "-seed%s" % eval_s['cv_seed'] if 'cv_seed' in eval_s else "",
            )

    base_dataset_settings = {
        'exp_name': fss_settings['exp_name'],
        'pos_neg_file': fss_settings['pos_neg_file']}

    config_files = []
    for net in fss_settings['networks_to_run']:
        net_config_map = copy.deepcopy(config_map)
        name = net['name']
        # if specified, use the given net_version
        net_version = net.get('net_version')
        if net_version is not None:
            name = net_version
        net_settings = net.get('net_settings')
        if net.get('network_collection'):
            base_dir = "%s/networks/%s" % (datasets_dir, name)
            base_net_files = "%s/net-files.txt" % (base_dir)
            net_files = pd.read_csv(base_net_files, header=None, squeeze=True)
            net_files = ["%s/%s" % (base_dir, f) for f in net_files]

            for net_file in net_files:
                curr_settings = setup_net_settings(
                    fss_dir, copy.deepcopy(base_dataset_settings), name,
                    net_file, net_settings=net_settings)
                net_config_map['input_settings']['datasets'].append(curr_settings) 
        else:
            # this is a single file
            net_file = "%s/networks/%s/%s" % (datasets_dir, net['name'], net['file_name'])
            curr_settings = setup_net_settings(
                fss_dir, copy.deepcopy(base_dataset_settings), name,
                net_file, string_networks=net.get('string_networks'),
                net_settings=net_settings)
            net_config_map['input_settings']['datasets'].append(curr_settings) 
        # now write the config file 
        config_file = "%s/%s%s.yaml" % (config_dir, name, eval_str)
        write_yaml_file(config_file, net_config_map)
        config_files.append(config_file)
    return config_files


def write_yaml_file(yaml_file, config_map):
    # make sure the directory exists first
    os.makedirs(os.path.dirname(yaml_file), exist_ok=True)
    print("Writing %s" % (yaml_file))
    with open(yaml_file, 'w') as out:
        yaml.dump(config_map, out, default_flow_style=False)


def setup_net_settings(
        input_dir, dataset_settings, net_version,
        net_file, string_networks=False, net_settings=None):
    """
    Add the network file and all other network settings to this config file
    Will create a symbolic link from the original dataset 
    """
    file_name = os.path.basename(net_file)
    new_net_file = "%s/networks/%s/%s" % (input_dir, net_version, file_name)
    if not os.path.isfile(new_net_file):
        print("\tadding symlink from %s to %s" % (net_file, new_net_file))
        os.makedirs(os.path.dirname(new_net_file), exist_ok=True)
        os.symlink(os.path.abspath(net_file), new_net_file)
    dataset_settings['net_version'] = "networks/"+net_version
    # The FSS pipeline uses this to distinguish the file with all individual string channels
    if string_networks is True:
        dataset_settings['string_net_files'] = [file_name]
    else:
        dataset_settings['net_files'] = [file_name]
    dataset_settings['exp_name'] += "-%s" % (file_name.split('.')[0])
    if net_settings is not None:
        dataset_settings['net_settings'] = net_settings
    return dataset_settings


if __name__ == "__main__":
    config_map, kwargs = parse_args()
    main(config_map, **kwargs)
