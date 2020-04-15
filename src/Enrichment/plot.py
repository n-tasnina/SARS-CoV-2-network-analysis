import matplotlib.pyplot as plt
import pandas as pd

def plot(rank, pvalue, geneset_name,algorithm):
    plt.plot(rank,pvalue)
    plt.xlabel('rank')
    plt.ylabel('pvalue')
    plt.title('overlap between '+geneset_name+ ' and proteins predicted by '+algorithm)
    plt.show()


