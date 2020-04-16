import matplotlib.pyplot as plt
import numpy as np
from textwrap import wrap
import pandas as pd
import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def plot(rank, pvalue, geneset_name, algo_dirpath):

    plt.show()
    threshold = 0.05
    pvalue = np.asarray(pvalue)
    rank = np.asarray(rank)
    significant = np.ma.masked_where(pvalue <= threshold, pvalue)
    non_significant = np.ma.masked_where(pvalue > threshold, pvalue)

    fig, ax = plt.subplots()
    ax.plot(rank, significant, '--o', rank, non_significant,'--o')
    plt.xlabel('rank')
    plt.ylabel('pvalue')
    plt.title("\n".join(wrap('overlap between ' + geneset_name + ' and proteins predicted by ' + os.path.basename(algo_dirpath))))
    plt.xticks(np.arange(0, max(rank) + 1, 200))

    geneset_name = geneset_name.replace("/", "_")
    plt.savefig(algo_dirpath +"/Plots/"+geneset_name, format='png')

