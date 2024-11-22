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
    if "Fq_" in a[0].split("_")[-1]:
        ID = "_".join(a[0].split("_")[:-2])
    else:
        ID = "_".join(a[0].split("_")[:-1])
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
            tree = tree.replace(k+"_"+I, k+"-"+v1+"_"+I)
else:
    for k, v in NAME2.items():
        tree = tree.replace(k, k+"-"+v["N"])

out = open(options.IN, "w")
out.write(tree+"\n")
