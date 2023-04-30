import pandas as pd
import numpy as np

def main():

    genome_metadata = pd.read_csv('assembly_summary.txt', sep='\t', header=1)
    accesion_id = pd.read_csv('ncbi_dataset.tsv', sep='\t', header=0)

    # Merge the two dataframes
    merged_df = pd.merge(genome_metadata, accesion_id, left_on='# assembly_accession', right_on='Assembly Accession')

    # Add the file to the end of the path and write to a list
    with open('ftp_links.txt', 'w') as f:
        ftp_paths = merged_df['ftp_path'].tolist()
        for path in ftp_paths:
            file_name = path.split('/')[-1]
            f.write(f'{path}/{file_name}_genomic.fna.gz\n')

if __name__ == "__main__":
    main()