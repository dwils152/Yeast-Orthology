import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():

    hits_df = pd.read_csv(sys.argv[1], sep="\t", header=None)
    hits_df.columns = ['ORF', 'Query', 'E-value']

    # Plot the data
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.kdeplot(data=hits_df, x='E-value', color='#3C5488FF')
    plt.axvline(x=0.1, color='r', linestyle='--', linewidth=0.5)
    
    # Add axis labels and a title
    ax.set_xlabel('E-value')
    ax.set_ylabel('Frequency')

    # Save the plot
    plt.savefig('e-vals.png', dpi=600)

if __name__ == "__main__":
    main()