import pandas as pd
import numpy as np
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():

    fasta_file = sys.argv[1]
    hmmer_hits = sys.argv[2]

    # Read in fasta file
    records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    #print(records)

    hmmer_dict = {}
    with open(hmmer_hits, 'r') as fin:
        for line in fin:
            coords = line.split('\t')[0]
            identifier, start, end = coords.split(':')
            hmmer_dict[identifier] = (int(start), int(end))

    # Extract sequences
    hit_records = list()
    for identifier, start_end in hmmer_dict.items():
        fasta_key = identifier.split('_')[1]
        start, end = start_end
        start, end = int(start), int(end)
        hit_record = records[fasta_key]
        hit_seq = hit_record.seq
        if start > end:
            hit_seq = hit_seq[end+1:start+1].reverse_complement() #makes start always inclusive regardless of strand
        else:
            hit_seq = hit_seq[start:end] 
        seq_record = SeqRecord(hit_seq, id=identifier)
        hit_records.append(seq_record)

    # Write out fasta file
    SeqIO.write(hit_records, sys.argv[3], "fasta")

if __name__ == '__main__':
    main()