from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import glob

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

(options, args) = parser.parse_args()
parser.add_option_group(group)


# '/media/inter/mkapun/projects/HAPLOTYPES/Susi_test5/results/haplotypes/*/*_aln.fasta'


LIST = d(lambda: d(str))
Primers = d(str)

S = 1
Lengths = d(int)
for i in glob.glob(options.IN):
    FILE = load_data(i)
    C = 0
    for l in FILE:
        if l.startswith(">"):
            C += 1
            ID = l.rstrip()
            continue
        if not ID.endswith("-1"):
            continue
        Primers[ID]
        if C < 2:
            Lengths[S] += len(l.rstrip())
        LIST[i][ID] += l.rstrip()
    S += 1

Start = 1
out2 = open(options.OUT+".part", "wt")
for k, v in sorted(Lengths.items()):
    out2.write("DNA, part"+str(k)+" = "+str(Start)+"-"+str(v+Start-1)+"\n")
    Start = v+Start

Global = d(str)

for k, v in LIST.items():
    LEN = len(list(v.values())[0])
    for PR in Primers.keys():
        if PR not in LIST[k]:
            Global[PR] += "-"*LEN
        else:
            Global[PR] += LIST[k][PR]
out1 = open(options.OUT+".fasta", "wt")
for k, v in Global.items():
    if v.count("-")/len(v) > 0.3:
        continue
    out1.write(k+"\n"+v+"\n")
