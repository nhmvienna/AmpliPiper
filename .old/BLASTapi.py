from Bio.Blast import NCBIWWW
from argparse import ArgumentParser
from Bio import SeqIO
import gzip
from Bio.Blast import NCBIXML
import os

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to input file with sequences", required=True)
argparse.add_argument(
    "-o", "--outfile", help="Path to output file with sequences", required=True)


args = argparse.parse_args()


inf = args.infile


def load_data(infile):
    """Load data from infile if it is in fasta format (after having unzipped it, if it is zipped)"""
    if infile.endswith(".gz"):  # If file is gzipped, unzip it
        y = gzip.open(infile, "rt", encoding="latin-1")
        # Read file as fasta if it is fasta
        if infile.endswith(".fasta.gz") or infile.endswith(".fna.gz") or infile.endswith(".fas.gz") or infile.endswith(".fa.gz"):
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
        else:
            y.close()
            raise ValueError("File is the wrong format")
    # Read file directly as fasta if it is a not zipped fasta: handle also more uncommon extensions :-)
    elif infile.endswith(".fasta") or infile.endswith(".fna") or infile.endswith(".fas") or infile.endswith(".fa"):
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fasta")
            sequences = {}
            for record in records:
                sequences.update({str(record.id): str(record.seq)})
            y.close()
            return sequences
    else:
        raise ValueError("File is the wrong format")


if __name__ == "__main__":
    print("Loading data...")

    # Load query sequences from the input file
    query_sequences = load_data(inf)
    print("Loaded")

    # Prepare the output file for BLAST results
    result_file_path = os.path.join(args.outfile)
    result_file = open(result_file_path, "w")

    # Write the header to the output file
    result_file.write(
        "QUERY_ID,HIT_ID,HIT_LENGTH,PERC_IDENTITY,GAPS,QUERY_COVERAGE,EVAL\n")
    result_file.close()

    # Iterate through each query sequence
    for query_sequence in list(query_sequences.keys()):

        # Perform a BLAST search against the NCBI nr database
        print(f"Starting BLAST for sequence {query_sequence}")
        result_handle = NCBIWWW.qblast(
            "blastn", "nt", query_sequences[query_sequence])

        # Parse and print the result
        blast_record = NCBIXML.read(result_handle)

        # Get the length of the query sequence
        qlen = blast_record.query_length

        # Open the result file in append mode
        result_file = open(result_file_path, "a")
        print(f"Writing BLAST for sequence {query_sequence}")

        # Process each alignment and its associated HSPs
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                title = alignment.title.replace(",", " -")

                # Calculate percentage identity, query coverage, and round the values
                perc_identity = 100 * \
                    round(int(hsp.identities) / int(hsp.align_length), 4)
                query_coverage = 100 * \
                    round((int(hsp.query_end) -
                          int(hsp.query_start) + 1) / int(qlen), 4)

                # Write the result to the output file
                result_file.write(
                    f"{query_sequence},{title},{alignment.length},{perc_identity},{hsp.gaps}/{int(hsp.align_length)},{query_coverage},{hsp.expect}\n")

        result_handle.close()
        result_file.close()
        print("Finished")


# https://drive.google.com/file/d/1-ysk9prcplU0vnciteZSwQ22YKcS4wea/view?usp=sharing
