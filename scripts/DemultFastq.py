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
from collections import defaultdict as d
from collections import Counter

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
argparse.add_argument(
    "-sr", "--sizerange", help="+/- size ", type=float, required=False, default=0.5)
argparse.add_argument(
    "-mr", "--minreads", help="minimum number of reads", type=float, required=False, default=0.5)


args = argparse.parse_args()
inf = args.infile
prim = args.primers
TH = args.kthreshold
readper = args.reads_percentage
sizerange = args.sizerange
minreads = args.minreads


def load_data(infile):
    """Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)"""
    try:
        print("Reading data from input file...", file=sys.stderr)
        if infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(infile, "rt", encoding="latin-1")
        elif infile.endswith(".fastq"):
            y = open(infile, "r")
        else:
            raise ValueError("File is the wrong format")
        records = SeqIO.parse(y, "fastq")
        seq_dict = []  # Create a dictionary to store everything from the file
        for record in records:
            # Update dictionary with header as key, sequence as [0] element of the value list, and base quality string as [1] element of the value list
            seq_dict.append([record.id,
                             str(record.seq),
                            record.format("fastq").split("\n")[3],
                            round(sum(list(record.letter_annotations["phred_quality"]))/len(
                                list(record.letter_annotations["phred_quality"])), 2)])
        y.close()
        print("Done", file=sys.stderr)
        return seq_dict

    except FileNotFoundError:
        print(
            f"File {infile} does not exist, proceeding with the analysis...", file=sys.stderr)
        return


def read_primer_table(csv):
    """Read primers table, provided as a csv file whose separatore MUST be comma and whose fields MUST be named and ordered as follows: ID,FWD,REV
    ID: contains the demultiplexing unit to which the primers are referred
    FWD and REV: idicate respectively the forward and reverse primer sequences"""
    print("Reading primers table...", file=sys.stderr)
    try:
        primers = pd.read_csv(csv)  # read primers table thanks to pandas
        fp = primers["FWD"]  # Forward sequences
        rp = primers["REV"]  # Reverse sequences
        ids = primers["ID"]  # Demultiplexing unit ID
        size = primers["SIZE"]  # Demultiplexing unit ID
        primer_dict = d(list)
        for i in range(len(ids)):
            # Assign to each demultiplexing unit ID its reverse and forward primers
            primer_dict[ids[i]] = [fp[i].upper(), rp[i].upper(), size[i]]
        print("Done", file=sys.stderr)
        return primer_dict
    except KeyError or ValueError:
        raise KeyError(
            "Fields of the csv do not comply with the requirements")


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def do_alignment(F, R, SEQ, TH):
    ''' test fwd and rev alignment'''
    correspondences = [("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("M", "A"), ("M", "C"), ("K", "G"), ("K", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), (
        "H", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]  # Set correspondences between degenerated bases for edlib to know what to retain as equivalent while aligning: option is suggested by the PyPi page of the edlib package
    aln1 = edlib.align(
        F, SEQ,
        task="path",
        mode="HW",
        additionalEqualities=correspondences)
    # If the alignement is significant (alignment process not aborted at aln1) also try reverse
    if aln1["editDistance"] > len(F)*TH:
        return False, "NA"
    aln2 = edlib.align(
        R, SEQ,
        task="path",
        mode="HW",
        additionalEqualities=correspondences)
    if aln2["editDistance"] <= len(R)*TH:
        return True, (aln1["locations"][0], aln2["locations"][0])

    else:
        # else return false
        return False, "NA"


def demultiplex(infile, csv):
    """Demultiplexing function"""
    sequences_dict = load_data(infile)  # Load data from infile
    # seq_keys = list(sequences_dict.keys())  # Extract headers
    # sequences = [seq[0]
    #              for seq in list(sequences_dict.values())]  # Extract sequences
    primer_dict = read_primer_table(csv)  # Read primers table
    demultiplexed = d(list)
    print("Total number of Sequences:"+str(len(sequences_dict)))
    for primername, primers in primer_dict.items():
        # Set a dictionary with the demultiplexing units IDs as key: the list passed as value of the dictionary will contain the reads assigned to every demultiplexed unit
        for SeqDATA in sequences_dict:  # Loop through the sequences loaded with load_data
            # Set the maximum numbers of bps at the beginning of the sequences that will be aligned with the forward primer
            FWD, REV, SIZE = primers
            # exclude sequences outside sizerange
            if len(SeqDATA[1]) < SIZE-int(sizerange) or len(SeqDATA[1]) > SIZE+int(sizerange):
                continue
            T1, D1 = do_alignment(FWD, reverse_complement(REV), SeqDATA[1], TH)
            T2, D2 = do_alignment(REV, reverse_complement(FWD), SeqDATA[1], TH)

            if T1:
                ALL = sorted([i for sub in D1 for i in sub])
                LENGTH = max(ALL)-min(ALL)
                # print(primername, D1, LENGTH, SIZE, int(sizerange), SeqDATA)
                if LENGTH > SIZE-int(sizerange):
                    demultiplexed[primername].append((
                        SeqDATA, (min(ALL), max(ALL))))
            elif T2:
                ALL = sorted([i for sub in D2 for i in sub])
                LENGTH = max(ALL)-min(ALL)
                if LENGTH > SIZE-int(sizerange):
                    demultiplexed[primername].append((
                        SeqDATA, (min(ALL), max(ALL))))
                # print(SeqDATA.extend((min(ALL), max(ALL))))
    return demultiplexed


def sort_filter(demultiplexed, ReadCount, Output):
    for k, v in demultiplexed.items():

        BQdict = d(list)
        for items, RANGE in v:
            BQdict[items[0]] = items[3]
        ToPrint = [x for x, y in Counter(
            dict(BQdict)).most_common(int(ReadCount))]
        if len(ToPrint) < int(minreads):
            continue
        print(k, len(ToPrint))
        OUT = open(Output+"/"+k+".fastq", "wt")
        for items, RANGE in v:
            ID, SEQ, BQ, AvBQ = items
            if ID not in ToPrint:
                continue
            OUT.write("@"+ID+"\n"+SEQ[RANGE[0]:RANGE[1]] +
                      "\n+\n"+BQ[RANGE[0]:RANGE[1]]+"\n")
        OUT.close()


if __name__ == "__main__":
    # Call the functions
    start = datetime.now()
    demult = demultiplex(inf, prim)
    # demultiplexed_files = write_demult_files(args.output, inf, d)
    sort_filter(demult, readper, args.output)
    end = datetime.now()
    print('Duration of demultiplexing process: {}'.format(
        end - start), file=sys.stderr)
