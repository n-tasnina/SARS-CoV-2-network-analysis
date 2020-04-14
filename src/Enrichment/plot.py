import matplotlib.pyplot as plt
import pandas as pd

def plot(file,geneset_name):
    df = pd.read_csv(file,'\t')
