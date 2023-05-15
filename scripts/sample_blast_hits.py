import sys
import pandas as pd

def main():
    # Get the input and output filenames from the command line
    yeast_tree = sys.argv[1]
    blast_hits = sys.argv[2]


    with open(yeast_tree, 'r') as f:
        data = [line.strip().split('\t') for line in f]

    df = pd.DataFrame(data, columns=['ID', 'Classification'])
    df[['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']] = df['Classification'].str.split(';', expand=True)
    df.drop('Classification', axis=1, inplace=True)
    df['ID'] = df['ID'].astype(int)

    # Read the blast hits file into a DataFrame
    blast_df = pd.read_csv(blast_hits, sep='\t', header=None, names=['qseqid', 'sseqid', 'sseq', 'staxid'])

    merged_df = pd.merge(df, blast_df, left_on='ID', right_on='staxid')

    #randomly sample a maximum of 5 hits per family
    sampled_df = merged_df.groupby('Family').apply(lambda x: x.sample(min(len(x), 5)))
    sampled_df.to_csv(sys.argv[3], columns=['qseqid', 'sseqid', 'sseq', 'staxid'], index=False, sep='\t', header=False)

if __name__ == '__main__':
    main()
