from Bio import SeqIO
import gzip
import sys
from Bio.Seq import Seq
import os


def mk_or_mv_dir(directorypath):
    try:
        os.makedirs(directorypath)
        return directorypath
    except FileExistsError:
        return directorypath


def list_files(directory):
    filenames = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            filenames.append(os.path.join(root, file))
    return filenames


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def load_data(infile):
    """
    Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)

    Parameters:
        infile (str): Path to the input file.

    Returns:
        dict or bool: Dictionary containing sequence information if successful, False otherwise.
    """
    try:
        print("Reading data from input file...", file=sys.stderr)
        if infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(infile, "rt", encoding="latin-1")
            if (
                infile.endswith(".fasta.gz")
                or infile.endswith(".fa.gz")
                or infile.endswith(".fas.gz")
                or infile.endswith(".fna.gz")
            ):
                # Read file as fasta if it is fasta
                records = SeqIO.parse(y, "fasta")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as value
                    seq_dict.update({record.id: str(record.seq)})
                y.close()
                return seq_dict
        # Read file directly as fasta if it is not zipped fastq
        elif (
            infile.endswith(".fasta")
            or infile.endswith(".fa")
            or infile.endswith(".fas")
            or infile.endswith(".fna")
        ):
            with open(infile, "r") as y:
                records = SeqIO.parse(y, "fasta")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as value
                    seq_dict.update({record.id: str(record.seq)})
                y.close()
                return seq_dict
        else:
            raise ValueError("File is the wrong format")
        print("Done", file=sys.stderr)
    except FileNotFoundError:
        print(
            f"File {infile} does not exist, proceeding with the analysis...",
            file=sys.stderr,
        )
        return False
