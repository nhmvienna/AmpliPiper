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
parser.add_option("--gblocks", dest="GB", help="Output file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

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


def meanstdv(x):
    ''' calculate mean, stdev and standard error : x must be a list of numbers'''
    from math import sqrt
    n, mean, std, se = len(x), 0, 0, 0
    if len(x) == 0:
        return "na", "na", "na"
    for a in x:
        mean = mean + a
    mean = mean / float(n)
    if len(x) > 1:
        for a in x:
            std = std + (a - mean)**2
        std = sqrt(std / float(n-1))
        se = std/sqrt(n)
    else:
        std = 0
        se = 0
    return mean, std, se


GB = d(int)
for l in load_data(options.GB):
    if l.startswith(">"):
        ID = l.rstrip()[1:]
        continue
    GB[ID] += l.replace(" ", "").count("-")

M, SD = meanstdv(GB.values())[:2]

EXCL = []
for k, v in GB.items():
    if v > M+SD:
        EXCL.append(k)

print(EXCL)

for l in load_data(options.IN):
    if l.startswith(">"):
        ID = l.rstrip()[1:]
    if ID in EXCL:
        continue
    print(l.rstrip())
