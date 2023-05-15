import sys
from Bio import SeqIO


def main():
    
    # Define the input fasta file and tsv file names
    fasta_file = sys.argv[1]
    tsv_file = sys.argv[2]

    # Read the tsv file and store the data in a dictionary
    tsv_data = {}
    with open(tsv_file) as tsv:
        for line in tsv:
            fields = line.strip().split('\t')
            # If the current key already exists in the dictionary, add the new value to the existing list of values
            if fields[0] in tsv_data:
                tsv_data[fields[0]].append(fields[1:])
            else:
                tsv_data[fields[0]] = [fields[1:]]

    # Create a set of unique headers to match to the input fasta file
    matching_headers = set()
    for key in tsv_data:
        matching_headers.add(key)

    # Create a list of sequence records for matching sequences
    matching_seqs = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        # If the current sequence header matches a key in the set of unique headers, add the sequence record to the matching_seqs list
        if record.id.split()[0] in matching_headers:
            matching_seqs.append(record)
            # Remove the current header from the set of unique headers to prevent redundancy
            matching_headers.remove(record.id.split()[0])

    # Write the matching sequences to a new fasta file
    SeqIO.write(matching_seqs, sys.argv[3], 'fasta')
    
if __name__ == "__main__":
    main()
