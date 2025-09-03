import edlib
from Bio.Seq import Seq
from argparse import ArgumentParser
import sys
from collections import defaultdict as d

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument(
    "-s", "--sanger", help="Path to the csv file containing primers", required=True)
argparse.add_argument(
    "-p", "--primer", help="Primer as comma-sep list", required=True)
argparse.add_argument(
    "-d", "--divergence", help="Divergence to ref", required=True)
argparse.add_argument(
    "-f", "--frequency", help="Frequency", required=True)
argparse.add_argument(
    "-r", "--replicate", help="Repolicate", required=True)
argparse.add_argument(
    "-c", "--SC", help="SimilarityConsensus", required=True)


args = argparse.parse_args()
cons = args.infile
sanger = args.sanger
FWD, REV = args.primer.split(",")


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def do_alignment(Short, Long):
    ''' test fwd and rev alignment'''
    aln = edlib.align(
        Short, Long,
        task="path",
        mode="HW")

    raln = edlib.align(
        Short, reverse_complement(Long),
        task="path",
        mode="HW")
    # print(raln, aln)
    # If the alignement is significant (alignment process not aborted at aln1) also try reverse
    if aln["editDistance"] > raln["editDistance"]:
        return raln["editDistance"], raln["editDistance"]/len(Short)
    else:
        return aln["editDistance"], aln["editDistance"]/len(Short)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


def readFASTA(x):
    DATA = d(str)
    for l in x:
        if l.startswith(">"):
            ID = l[1:].rstrip()
            continue
        DATA[ID] += l.rstrip()
    return DATA


Cons = readFASTA(load_data(cons))
Sanger = readFASTA(load_data(sanger))
# remove primers:
for k, v in Sanger.items():
    Sanger = v[len(FWD):-len(REV)]

TEST = d()
for ID, ConsSeq, in Cons.items():
    SangerSeq = Sanger
    if len(ConsSeq) >= len(SangerSeq):
        Full, rel = do_alignment(SangerSeq, ConsSeq)
    else:
        Full, rel = do_alignment(ConsSeq, SangerSeq)
    TEST[Full] = rel
    # print(SangerSeq)
    # print(ConsSeq)
print(args.divergence, args.replicate, args.SC, args.frequency, len(Cons), min(
    TEST.keys()), TEST[min(TEST.keys())], sep=",")
