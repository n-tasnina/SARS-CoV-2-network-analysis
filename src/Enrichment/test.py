import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd
import argparse
import yaml
import sys
import os

def nCr(n, r):
    return (fact(n) / (fact(r)
                       * fact(n - r)))


# Returns factorial of n
def fact(n):
    res = 1

    for i in range(2, n + 1):
        res = res * i

    return res


def find_probability(a,b,c,d):
    x = nCr(a+c, a)*nCr(b+d, b)/ nCr(a+b+c+d, a+b)
    print(x)
    return x

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


def parse_args():
    parser = setup_opts()
    args = parser.parse_args()
    kwargs = vars(args)
    with open(args.config, 'r') as conf:
        config_map = yaml.load(conf, Loader=yaml.FullLoader)
    return config_map, kwargs



def main(config_map, **kwargs):

    predicted_prot_dir = config_map['predicted_prot_dir']

    for dirpath, dirs, files in os.walk(predicted_prot_dir):

        network_name = dirpath.split('/')[3]
        print(network_name)

        # print(dirpath)
        for filename in files:

            fname = os.path.join(dirpath, filename)
            print(fname)



if __name__ == "__main__":
    config_map, kwargs = parse_args()
    main(config_map, **kwargs)
