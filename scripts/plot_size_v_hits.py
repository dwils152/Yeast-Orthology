import pandas as pd
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def main():

    genome_sizes = pd.read_csv(sys.argv[1], sep='\t', header=None, names=["Genome", "Size"])
    hits = pd.read_csv(sys.argv[2], sep=' ', header=None, names=["Hits", "Genome"])
    hits["Genome"] = hits["Genome"].str.replace('.orfs.hmmsearch.tsv', '.fna')
    merged = pd.merge(hits, genome_sizes, on="Genome")
    
    fig, ax = plt.subplots(figsize=(7, 7))
    sns.jointplot(data=merged, x="Hits", y="Size",
                  kind="reg", truncate=False,
                  color="#3C5488FF", height=7)

    ax.set_xlabel("Genome Size (bp)", fontsize=35)
    ax.set_ylabel("Hmmer Hits", fontsize=35)
    plt.savefig("size_v_hits.png", dpi=600)
    

if __name__ == "__main__":
    main()