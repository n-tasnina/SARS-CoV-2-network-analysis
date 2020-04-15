import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats

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

def main():

    # total = 9000
    # a_b = 100
    # c_d = total-a_b
    #
    # for a_c in range(100, 500, 100):
    #     b_d = total - a_c
    #     for a in range (0,5,1):
    #         b = a_b - a
    #         c= a_c - a
    #         d = b_d - b
    #         print('for rank = ', a_c, ' and overlap = ', a )
    #         oddratio, pvalue1 = stats.fisher_exact([[a, b], [c, d]], 'greater')
    #         print(pvalue1)
    #
    #         oddratio, pvalue2 = stats.fisher_exact([[a, b], [c, d]],'less')
    #         print(pvalue2)
    #
    #         prob = find_probability(a,b,c,d)
    #
    #         print(pvalue1+pvalue2-prob)
    a = '57: SARS coronavirus nsp7-pp1a/pp1ab (gene: orf1ab) from Virus-Host PPI P-HIPSTer 2020'
    a = a.replace('/','_')
    print(a, type(a))





if __name__ == "__main__":
    main()

