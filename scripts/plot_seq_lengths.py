from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import sys

def main():
    # Define the two fasta files
    fasta1 = sys.argv[1]
    fasta2 = sys.argv[2]

    # Define the function to compute the sequence lengths from a fasta file
    def get_seq_lengths(file):
        lengths = []
        for record in SeqIO.parse(file, "fasta"):
            lengths.append(len(record.seq))
        return lengths

    # Get the sequence lengths for each file
    lengths1 = get_seq_lengths(fasta1)
    lengths2 = get_seq_lengths(fasta2)

    # Create a two-panel figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

    # Plot the length distributions for each file using Seaborn
    sns.histplot(lengths1, bins=50, ax=ax1, color="#6A6599FF")
    ax1.set_xlabel("Sequence Length")
    ax1.set_ylabel("Frequency")
    ax1.set_title("All Hits")

    sns.histplot(lengths2, bins=50, ax=ax2, color="#79AF97FF")
    ax2.set_xlabel("Sequence Length")
    ax2.set_ylabel("Frequency")
    ax2.set_title("All Hits w. Annotations")

    # Show the plot
    plt.tight_layout()
    plt.savefig(sys.argv[3], dpi=600)

if __name__ == "__main__":
    main()