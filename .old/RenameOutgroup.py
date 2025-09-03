from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import os
import sys

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--haplotypes", dest="HA", help="Input file")
parser.add_option("--primername", dest="PI", help="Input file")
parser.add_option("--samplenames", dest="SA", help="Input file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

SA = options.SA
if SA == "no":
    print("no")
    sys.exit()

SAh = []
for l in load_data(options.HA):
    if l.startswith(">"):
        for i in SA.split(","):
            if i in l:
                SAh.append(l.rstrip()[1:])
                continue
SAh = ",".join(SAh)

if os.path.isfile(options.IN):
    NAME = d(list)
    for l in load_data(options.IN):
        if l.startswith("SAMPLE"):
            continue
        a = l.rstrip().rstrip(",").split(",")
        ID = "-".join(a[0].split("-")[:-1])
        # ID = a[0]
        if len(a) == 1:
            NAME[ID].append("NA")
            continue
        HASH = d(lambda: d(lambda: d(str)))
        for sample in a[1:]:
            Spec = sample.split(" (")[0].replace(" ", "_")
            Spec = Spec.split("_")[0][:3]+"."+Spec.split("_")[1]
            Sim = float(sample.split(" (")[1].split("%)")[0])
            count = int(sample.split("count = ")[1])
            HASH[ID][Sim][count] = [Spec, str(round(Sim/100, 2))]
        Sim = max(HASH[ID].keys())
        count = max(HASH[ID][Sim].keys())
        NAME[ID].append(HASH[ID][Sim][count])
    # print(NAME)
    NAME2 = d(lambda: d(str))
    # now test if
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME.items():
            for i in range(len(v)):
                NAME2[k][str(i+1)] = "_".join(v[i])
    else:
        for k, v in NAME.items():
            Spec, Sim = list(zip(*v))
            if len(list(set(Spec))) > 1:
                NAME2[k]["N"] = "NA"
            else:
                NAME2[k]["N"] = "_".join(v[0])

    tree = load_data(options.IN).readline()
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME2.items():
            for I, v1 in v.items():
                SAh = SAh.replace(k+"-"+I, k+"-"+v1+"-"+I)
    else:
        for k, v in NAME2.items():
            SAh = SAh.replace(k, k+"-"+v["N"])
    print(SAh)
else:
    print(SAh)
