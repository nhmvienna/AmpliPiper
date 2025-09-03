from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--primername", dest="PI", help="Input file")
parser.add_option("--names", dest="names", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

NAME = d(list)
for l in load_data(options.names):
    if l.startswith("SAMPLE"):
        continue
    a = l.rstrip().rstrip(",").split(",")
    ID = "_".join(a[0].split("_")[:-1])
    # ID = a[0]
    if len(a) == 1:
        NAME[ID].append("NA")
        continue
    HASH = d(lambda: d(lambda: d(str)))
    for sample in a[1:]:
        Spec = sample.split(" (")[0].replace(" ", "_")
        Spec = Spec.split("_")[0][:3]+"."+Spec.split("_")[1][:3]
        Sim = float(sample.split(" (")[1].split("%)")[0])
        count = int(sample.split("count = ")[1])
        HASH[ID][Sim][count] = [Spec, str(round(Sim/100, 2))]
    Sim = max(HASH[ID].keys())
    count = max(HASH[ID][Sim].keys())
    NAME[ID].append(HASH[ID][Sim][count])

NAME2 = d(str)
# now test if
if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
    for k, v in NAME.items():
        for i in range(len(v)):
            NAME2[k+"_"+str(i+1)] = "_".join(v[i])
else:
    for k, v in NAME.items():
        Spec, Sim = list(zip(*v))
        if len(list(set(Spec))) > 1:
            NAME2[k] = "NA"
        else:
            NAME2[k] = "_".join(v[0])


for l in load_data(options.IN):
    if l.startswith(">"):
        if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
            ID = l[1:].rstrip()
        else:
            ID = "_".join(l[1:].rstrip().split("_")[:-1])
        if ID in NAME2:
            print(l.rstrip()+"-"+NAME2[ID])
        else:
            print(l.rstrip())
        continue
    print(l.rstrip())
