from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_gb_file(file_path: str) -> SeqRecord:
    """
    Function to return the objects in the genbank file

    :param file_path:
    :return str:
    """
    gb_obj = SeqIO.read(file_path, 'gb')
    return gb_obj


def read_genome_file(file_path: str) -> str:
    """
    Function to read a genome fasta file

    :param str file_path: file path of the genome fasta file
    :return str:
    """
    fasta_sequences = SeqIO.parse(open(file_path), 'fasta')
    genome_sequence = [seq for seq in fasta_sequences][0]
    return genome_sequence