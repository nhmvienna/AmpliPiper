from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
from module import load_data

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--Alignment", dest="Al", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)





Rep = []

for l in load_data(options.IN):
    if l.startswith(">"):
        Rep.append(l.rstrip())
C = 0
for l in load_data(options.Al):
    if l.startswith(">"):
        print(Rep[C])
        C += 1
        continue
    print(l.rstrip())
