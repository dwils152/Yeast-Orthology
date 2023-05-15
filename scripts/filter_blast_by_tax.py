import pandas as pd
import sys

def main():

    blast = pd.read_csv(sys.argv[1], sep = "\t", header = None,
                        names=['qseqid', 'sseqid' ,'sseq', 'staxid'])
    
    yeast_taxs = pd.read_csv(sys.argv[2], sep = "\t", header = None, names=['taxid', 'name'])

    blast_yeast = blast[blast['staxid'].isin(yeast_taxs['taxid'])]

    blast_yeast.to_csv(sys.argv[3], sep = "\t", index = False, header = False)

if __name__ == "__main__":
    main()