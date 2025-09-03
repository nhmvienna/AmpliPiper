import edlib
from math import ceil
from module import FileSystem
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import os
from argparse import ArgumentParser
from datetime import datetime
import gzip
import sys

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument(
    "-p", "--primers", help="Path to the csv file containing primers", required=True)
argparse.add_argument(
    "-o", "--output", help="Path to the output folder", required=True)
argparse.add_argument(
    "-rp", "--reads_percentage", help="Providing this command, the demultiplexer will output the top nr%% of reads", type=float, required=False, default=1.0)
argparse.add_argument(
    "-th", "--kthreshold", help="0.5 means that there should not be any differences between primer and target sequence. The higher this value the more liberal. Should not exceed 0.6", type=float, required=False, default=0.5)


args = argparse.parse_args()

inf = args.infile
prim = args.primers
TH = args.kthreshold
readper = args.reads_percentage


def load_data(infile):
    """Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)"""
    try:
        print("Reading data from input file...", file=sys.stderr)
        if infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(infile, "rt", encoding="latin-1")
            if infile.endswith(".fastq.gz"):  # Read file as fastq if it is fastq
                records = SeqIO.parse(y, "fastq")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
                    seq_dict.update(
                        {record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                y.close()
                return seq_dict
        # Read file directly as fastq if it is a not zipped fastq
        elif infile.endswith(".fastq"):
            with open(infile, "r") as y:
                records = SeqIO.parse(y, "fastq")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
                    seq_dict.update(
                        {record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                y.close()
                return seq_dict
        else:
            raise ValueError("File is the wrong format")
        print("Done", file=sys.stderr)
    except FileNotFoundError:
        print(
            f"File {infile} does not exist, proceeding with the analysis...", file=sys.stderr)
        return


def read_primer_table(csv):
    """Read primers table, provided as a csv file whose separatore MUST be comma and whose fields MUST be named and ordered as follows: ID,FWD,REV 
    ID: contains the demultiplexing unit to which the primers are referred
    FWD and REV: idicate respectively the forward and reverse primer sequences"""
    print("Reading primers table...", file=sys.stderr)
    if os.path.splitext(csv)[1] == ".csv":
        try:
            primers = pd.read_csv(csv)  # read primers table thanks to pandas
            fp = primers["FWD"]  # Forward sequences
            rp = primers["REV"]  # Reverse sequences
            ids = primers["ID"]  # Demultiplexing unit ID
            fwd_dict = {}
            rev_dict = {}
            for i in range(len(ids)):
                # Assign to each demultiplexing unit ID its reverse and forward primers
                fwd_dict.update({ids[i]: fp[i]})
                rev_dict.update({ids[i]: rp[i]})
            print("Done", file=sys.stderr)
            return fwd_dict, rev_dict
        except KeyError or ValueError:
            raise KeyError(
                "Fields of the csv do not comply with the requirements")
    else:
        raise ValueError("Input file is the wrong format")


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def demultiplex(infile, csv):
    """Demultiplexing function"""
    sequences_dict = load_data(infile)  # Load data from infile
    seq_keys = list(sequences_dict.keys())  # Extract headers
    sequences = [seq[0]
                 for seq in list(sequences_dict.values())]  # Extract sequences
    correspondences = [("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("M", "A"), ("M", "C"), ("K", "G"), ("K", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), (
        "H", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]  # Set correspondences between degenerated bases for edlib to know what to retain as equivalent while aligning: option is suggested by the PyPi page of the edlib package
    fwd_dict, rev_dict = read_primer_table(csv)  # Read primers table
    demultiplexed = {}
    print("Total number of Sequences:"+str(len(sequences)))
    for i in list(fwd_dict.keys()):
        # Set a dictionary with the demultiplexing units IDs as key: the list passed as value of the dictionary will contain the reads assigned to every demultiplexed unit
        demultiplexed.update({i: {}})
    for i in list(demultiplexed.keys()):
        for j in range(len(sequences)):  # Loop through the sequences loaded with load_data
            # Set the maximum numbers of bps at the beginning of the sequences that will be aligned with the forward primer
            nf = len(fwd_dict[i])
            # align the forward primer with the first nf bp of the sequence, if the edit distance is bigger than k, the alignment is aborted and the functions returns -1 as editDistance
            aln1 = edlib.align(sequences[j][:nf], fwd_dict[i], task="path",
                               mode="NW", additionalEqualities=correspondences)

            # if aln1["editDistance"] <= nf*TH:
            #     nice = edlib.getNiceAlignment(
            #         aln1, sequences[j][:nf], fwd_dict[i])
            #     print(nice)

            # If the alignement is significant (alignment process not aborted at aln1)
            if aln1["editDistance"] <= nf*TH:
                # Set the maximum numbers of bps at the end of the sequences that will be aligned with the reverse primer
                nr = len(rev_dict[i])
                # the reverse primer is reverse complemented, the alignment goes as before
                aln2 = edlib.align(sequences[j][-nr:], reverse_complement(rev_dict[i]),
                                   task="path", mode="NW", additionalEqualities=correspondences)
                # If aln1 and aln2 were significant (alignment processes not aborted), the sequence is assigned to the demultiplexing unit
                if aln2["editDistance"] <= nr*TH:
                    # Update dictionary with header, sequence and base quality
                    demultiplexed[i].update(
                        {seq_keys[j]: [sequences[j], sequences_dict[seq_keys[j]][1]]})
                else:
                    continue
            else:  # Check for reverse strand
                nf = len(rev_dict[i])
                aln3 = edlib.align(sequences[j][:nf], rev_dict[i], task="path",
                                   mode="NW", additionalEqualities=correspondences)
                if aln3["editDistance"] <= nf*TH:
                    nr = len(fwd_dict[i])
                    # the reverse primer is reverse complemented, the alignment goes as before
                    aln4 = edlib.align(sequences[j][-nr:], reverse_complement(fwd_dict[i]), k=nr*TH,
                                       task="path", mode="NW", additionalEqualities=correspondences)
                    if aln4["editDistance"] <= nr*TH:
                        # Update dictionary with header, sequence and base quality
                        demultiplexed[i].update(
                            {seq_keys[j]: [sequences[j], sequences_dict[seq_keys[j]][1]]})
                    else:
                        continue
                else:
                    continue
        print("%d reads assigned to %s (%g of the total)" % (len(demultiplexed[i]), i, round(
            len(demultiplexed[i])/len(sequences), 4)), file=sys.stderr)
    # Return the dictionary with the demultiplexing units as keys and lists of their assigned reads as values
    return demultiplexed


def write_demult_files(output, inf, demultiplexed):
    print("Writing temporary demultiplexed files...", file=sys.stderr)
    FileSystem.makedir_orchange(output)
    demultiplexed_files = []
    for key in list(demultiplexed.keys()):
        with open(os.path.join(output, key+"_tmp.fastq"), "w") as f:
            c = 0
            demultiplexed_files.append(os.path.join(output, key+"_tmp.fastq"))
            for k in list(demultiplexed[key].keys()):
                # Write demultiplexed file as fasta file, with fasta header and sequence
                # Write the header as from the original file
                f.write("@"+k+"\n")
                f.write(demultiplexed[key][k][0]+"\n")  # Write the sequence
                f.write("+\n")
                # Write the base quaality string
                f.write(demultiplexed[key][k][1]+"\n")
        f.close()
    print("Done", file=sys.stderr)
    return demultiplexed_files


def select_perc_reads(demultiplexed_files, percentage: float, output):
    """Select the top percentage of reads in terms of quality"""
    print("Selecting reads and writing definitive demultiplexed files...", file=sys.stderr)
    for infile in demultiplexed_files:
        with open(infile, "r") as y:
            records = SeqIO.parse(y, "fastq")
            seq_dict = {}  # Create a dictionary to store everything from the file
            bq_dict = {}  # Create a dictionary to store everything from the file
            for record in records:
                # ignore sequences <300bp:
                if len(record.seq) < 300:
                    continue
                # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
                seq_dict.update(
                    {record.id: [str(record.seq), record.format("fastq").split("\n")[3]]})
                # Update dictionary with header as key, average base quality as value
                bq_dict.update({record.id: round(sum(list(record.letter_annotations["phred_quality"]))/len(
                    list(record.letter_annotations["phred_quality"])), 2)})
        y.close()
        bqs = [value for value in list(bq_dict.values())]
        bqs.sort()
        bqs.reverse()
        # now added the possibility to consider also a fixed number
        if percentage > 1 and percentage <= len(bqs):
            nr = ceil(percentage)
        elif percentage < 1:
            nr = ceil(percentage*len(bqs))
        else:
            nr = len(bqs)
        print("Selecting %d reads from %s" % (nr, infile), file=sys.stderr)
        bqs = bqs[:nr]
        final_list = []
        with open(os.path.join(output, FileSystem.get_base_dir(infile)[1].replace("_tmp", "")+".fastq"), "w") as f:
            for key in list(seq_dict.keys()):
                if len(final_list) < nr and bq_dict[key] in bqs:
                    # Write the header as from the original file
                    f.write("@"+key+"\n")
                    f.write(seq_dict[key][0]+"\n")  # Write the sequence
                    f.write("+\n")
                    # Write the base quaality string
                    f.write(seq_dict[key][1]+"\n")
                    final_list.append(1)
                else:
                    pass
        f.close()
        os.remove(infile)
    print("Done", file=sys.stderr)


if __name__ == "__main__":
    # Call the functions
    start = datetime.now()
    d = demultiplex(inf, prim)
    demultiplexed_files = write_demult_files(args.output, inf, d)
    select_perc_reads(demultiplexed_files, readper, args.output)
    end = datetime.now()
    print('Duration of demultiplexing process: {}'.format(
        end - start), file=sys.stderr)
