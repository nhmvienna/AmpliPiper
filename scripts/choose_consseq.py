from argparse import ArgumentParser
import os
import sys
from Bio import SeqIO
from module import FileSystem

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to the input files directory", required=True)

args = argparse.parse_args()

inf = args.infile


def has_two_lines(infile):
    ##Checks if there are two lines because, if there are more, it means there is more than one consensus sequence
    """Checks if the file has more than two lines: if so, returns True, else returns False"""
    with open(infile) as f:
        lines = f.readlines()
    f.close()
    if len(lines) == 2:
        return True
    return False


def read_infile(infile):
    """Reads the provided file, only if it is in fasta format and has more than one sequence in it, and returns a dictionary with the sequences as values and the headers as keys. Otherwise, returns False (the file contains only one seq) or raises ValueError (the file is the wrong format)"""
    print("Reading input file...", file=sys.stderr)
    if os.path.splitext(infile)[1] == ".fasta" and not has_two_lines(infile): ##Checks if there is more than one consensus sequence
        consseqs = {}
        with open(infile, "r") as fastafile:
            for record in SeqIO.parse(fastafile, "fasta"):
                consseqs.update({str(record.description): str(record.seq)}) ##Update the dictionary with the fasta header of the seq and the seq itself
        fastafile.close()
        print("Done", file=sys.stderr)
        return consseqs
    return False


def choose_the_most_supported(infile):
    """Creates and writes a file named selectedconsensus.fasta (in the same directory of the provided input file) with the most supported consensus sequence (if there are more than 1), or that is a copy of the input file (if there is only a consensus sequence)"""
    consseqs = read_infile(infile)
    if consseqs != False:
        support = []
        keys = list(consseqs.keys())
        for key in keys:
            supporting_reads = int(key.split("(")[1].split(")")[0]) ##Split the headers to get the read depth: >consensus_S7_0_0_(1492) becomes int(1492)##
            support.append(supporting_reads)
        seq = consseqs[keys[support.index(max(support))]] ##return the consensus sequence that correspond to the maximum read depth
        basedir = FileSystem.get_base_dir(infile)[0]
        flname = os.path.join(basedir, "selectedconsensus.fasta") ##Save the selected consensus sequence in a file named selectedconsensus.fasta
        with open(flname, "w") as f:
            f.write(">"+keys[support.index(max(support))].split("_")[1]+"\n")  ##Split the headers to get the name of the consensus sequence (>consensus_S7_0_0_(1492) becomes S7) and use it as header in the new fasta
            f.write(seq+"\n")
        f.close()
    elif consseqs == False and os.path.splitext(infile)[1] == ".fasta": ##If the original file contains only one consensus sequence, copy it in selectedconsensus.fasta
        with open(infile, "r") as f:
            lines = f.readlines()
        f.close()
        basedir = FileSystem.get_base_dir(infile)[0]
        flname = os.path.join(basedir, "selectedconsensus.fasta")
        with open(flname, "w") as fp:
            fp.write(">"+lines[0].split("_")[1]+"\n")
            fp.write(lines[1])
        fp.close()
    else:
        raise ValueError("Input file is the wrong format")


if __name__ == "__main__":
    try:
        choose_the_most_supported(inf)
    except FileNotFoundError:
        print(f"File {inf} dos not exist, proceeding with the analysis...", file=sys.stderr)
