import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def main():

    hits_df = pd.read_csv(sys.argv[1], sep="\t", header=None)
    hits_df.columns = ['ORF', 'Query', 'E-value']

    # Plot the data
    fig, ax = plt.subplots(figsize=(7, 12))
    sns.histplot(data=hits_df, x='E-value', ax=ax, color='#3C5488FF', bins=10, kde=True)
    
    # Add axis labels and a title
    ax.set_xlabel('E-value')
    ax.set_ylabel('Frequency')
    #ax.tick_params(axis='y', which='major', labelsize=5)
    #xfmt = ScalarFormatter(useMathText=True)
    #xfmt.set_powerlimits((-1,4))

    # Save the plot
    plt.savefig('e-vals.png', dpi=600)

if __name__ == "__main__":
    main()