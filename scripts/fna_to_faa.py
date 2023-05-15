from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys

def main(input_file, output_file):
    # Read the input nucleic acid FASTA file
    nucleic_acid_records = SeqIO.parse(input_file, 'fasta')

    # Translate the sequences to amino acids
    amino_acid_records = []
    for record in nucleic_acid_records:
        aa_seq = record.seq.translate()
        aa_record = SeqRecord(aa_seq, id=record.id, description=record.description)
        amino_acid_records.append(aa_record)

    # Write the translated sequences to an output amino acid FASTA file
    SeqIO.write(amino_acid_records, output_file, 'fasta')


if __name__ == '__main__':

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    main(input_file, output_file)
