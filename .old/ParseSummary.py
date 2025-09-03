from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import glob
import os

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)




def parse_AS(IN, OUT):
    # (1) ampliconsorter
    AS = f"{IN}/consensus_seqs/*/*/results.txt"
    AS_out = open(OUT+"_AS.csv", "wt")
    AS_out.write(
        "ID,Locus,HaplotypeID,TotalReads,ReadCount,Frequency,HaplotypeLength\n")
    for i in glob.glob(AS):
        ID = i.split("/")
        IND = ID[-3]
        LOCUS = ID[-2]
        FILE = load_data(i)
        RD, PC, HID = [], [], []
        for l in FILE:
            if l.startswith("-->"):
                HID.append(l.split("-> ")[1].split(".fasta contains")[0])
                RD.append(l.split("contains ")[1].split(" sequences")[0])
                PC.append(
                    str(round(float(l.split("sequences (")[1].split("% of total")[0])/100, ndigits=3)))
            if l.startswith("- used_reads"):
                URline = l.rstrip().split("= ")
                if len(URline) > 1:
                    UR = URline[1]
                else:
                    UR = "0"
        # length of sequences
        FASTA = f"{IN}/consensus_seqs/{IND}/{LOCUS}/{LOCUS}_consensussequences.fasta"
        if os.path.exists(FASTA):
            FILE = load_data(FASTA)
            SEQL = []
            for l in FILE:
                if l.startswith(">"):
                    continue
                SEQL.append(str(len(l.rstrip())))
        else:
            SEQL = []
        for i in range(len(RD)):
            AS_out.write(",".join([IND,
                                   LOCUS,
                                   HID[i],
                                   UR,
                                   RD[i],
                                   PC[i],
                                   SEQL[i]]) +
                         "\n")


if __name__ == "__main__":
    parse_AS(options.IN, options.OUT)
