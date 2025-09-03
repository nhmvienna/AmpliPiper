import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--counts", dest="CO", help="count of output reads")

(options, args) = parser.parse_args()
parser.add_option_group(group)


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


SEQ = []
COUNTS = options.CO.split(",")

for l in load_data(options.IN):
    if l.startswith(">"):
        continue
    SEQ.append(l.rstrip())

for i in range(len(COUNTS)):
    COUNT = int(COUNTS[i])
    C = 0
    for j in range(COUNT):
        print(">"+str(i)+"_"+str(j)+"\n"+SEQ[i])
