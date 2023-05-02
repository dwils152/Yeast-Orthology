import pandas as pd
import numpy as np
import sys

def main():
    try:
        hits = pd.read_csv(sys.argv[1], delim_whitespace=True,
                            header=None, usecols=[0,2,4], comment='#')
        hits.columns = ['ORF', 'Target', 'E-value']
        hits.to_csv(sys.argv[1]+".tsv", sep='\t', header=False, index=False)

        bed = hits['ORF'].str.split('_', expand=True)[1]
        bed = bed.str.split(':', expand=True)
        bed.columns = ['gi', 'start', 'end']
        bed.to_csv(sys.argv[1]+".bed", sep='\t', header=False, index=False)
    except:
        ...

if __name__ == '__main__':
    main()