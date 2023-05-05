import sys
import numpy as np
from Bio import SeqIO
import itertools
import random
from collections import OrderedDict

def main():
    
    records = SeqIO.to_dict(SeqIO.parse(sys.argv[1], 'fasta'))
    #sort the records by the keys
    records = OrderedDict(sorted(records.items()))

    #create an empty numpy array that how many rows as there are sequences in the fasta file and as many columns as there are kmers
    kmer_array = np.zeros((len(records), 4**3))

    non_standard_base_dict = {
        'R': ['A', 'G'],
        'Y': ['C', 'T'],
        'S': ['G', 'C'],
        'W': ['A', 'T'],
        'K': ['G', 'T'],
        'M': ['A', 'C'],
        'B': ['C', 'G', 'T'],
        'D': ['A', 'G', 'T'],
        'H': ['A', 'C', 'T'],
        'V': ['A', 'C', 'G'],
        'N': ['A', 'C', 'G', 'T']
    }

    for record in records.values():

        record.seq = record.seq.upper()
        #print(record.seq)
        kmer_dict = get_kmer_dict(3)
        for i in range(len(record.seq)-3):
            kmer = str(record.seq[i:i+3])
            #exclude the KMER if any of the bases are not A, T, C, or G
            if kmer not in kmer_dict.keys():
                new_kmer = ''
                for nt in kmer:
                    if nt not in 'ATCG':
                        new_kmer += random.choice(non_standard_base_dict[nt])
                    else:
                        new_kmer += nt
                kmer = new_kmer
            kmer_dict[kmer] += 1

        #add the values from the kmer_dict to the numpy array
        kmer_array[i] = list(kmer_dict.values())
    
    #write the numpy array to a file
    np.save(sys.argv[2], kmer_array)


def get_kmer_dict(k):
    kmer_dict = {}
    for kmer in itertools.product('ATCG', repeat=k): #add Ns? # Why does a sequence have S ._.
        kmer_dict[''.join(kmer)] = 0
    #sort the dictionary by the keys
    kmer_dict = OrderedDict(sorted(kmer_dict.items()))
    return kmer_dict


if __name__ == '__main__':
    main()