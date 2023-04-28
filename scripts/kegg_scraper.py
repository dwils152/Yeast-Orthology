import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import seaborn as sns

def main():

    """Reads in the KEGG html file and parses the table to get the number of genes for each organism."""

    # Read in the html file
    kegg_table = pd.read_html('../data/genes_statistics.html')
    kegg_df = kegg_table[1]
    fungal_df = kegg_df[kegg_df['Grp']['Grp'].str.contains('E.Fun')]
    fun_genes = fungal_df.loc[:, [('Org', 'Org'), ('Protein', 'KO'), ('Protein', 'All')]]
    fun_genes = fun_genes.sort_values(by=[('Protein', 'All')], ascending=False)
    fun_genes.columns = fun_genes.rename(columns={('Org', 'Org'): 'Organism', ('Protein', 'KO'): 'KO Annotated Genes', ('Protein', 'All'): 'Genes'}).columns.droplevel(0)
    
    # Plot the data
    fig, ax = plt.subplots(figsize=(7, 12))
    sns.barplot(data=fun_genes, x='All', y='Org', ax=ax, color='#3C5488FF', legend=True, label='All Genes')
    sns.barplot(data=fun_genes, x='KO', y='Org', ax=ax, color='#8491B4FF', legend=True, label='KO Annotated Genes')
    
    
    # Add axis labels and a title
    ax.set_xlabel('Genes')
    ax.set_ylabel('Species')
    ax.set_title('KEGG: Fungal Genes')  
    ax.tick_params(axis='y', which='major', labelsize=5)
    xfmt = ScalarFormatter(useMathText=True)
    xfmt.set_powerlimits((-1,4))

    # Show the plot
    plt.savefig('fungal_genes.png', dpi=600)





if __name__ == "__main__":
    main()