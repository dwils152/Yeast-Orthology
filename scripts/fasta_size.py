import sys
import os

def read_fasta_file(filename):
    with open(filename, 'r') as fin:
        sequences = []
        current_sequence = ''
        for line in fin:
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_sequence:
            sequences.append(current_sequence)
    return sequences

def count_total_length(sequences):
    return sum(len(seq) for seq in sequences)

def main():

    fasta = sys.argv[1]
    basename = os.path.basename(fasta)
    length = count_total_length((read_fasta_file(fasta)))
    print(f'{basename}\t{length}')

if __name__ == "__main__":
    main()