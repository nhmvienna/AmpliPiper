from module import load_data
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import sys
import sys

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, "< put description here >")

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--primername", dest="PI", help="Input file")
parser.add_option("--names", dest="names", help="Output file")
parser.add_option("--outgroup", dest="OG", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)

OG = options.OG.split(",")
OGnew = []

TEST = 0
NAME = d(list)
for l in load_data(options.names):

    # skip header
    if l.startswith("SAMPLE"):
        continue

    # split by comma
    a = l.rstrip().rstrip(",").split(",")

    # First test if Filename contains "Freq_" indicating --freqthreshold 0 in this case do not split by different Consensus sequences and do not append Species names to loci other than the barcode loci
    if "Freq_" in a[0]:
        TEST = 1
        ID = a[0]
    else:
        ID = "-".join(a[0].split("-")[:-1])
        # ID = a[0]

    # If no Species was identified, use ["NA", "0"]
    if len(a) == 1:
        # Only append if this ID hasn't been processed yet
        if ID not in NAME or len(NAME[ID]) == 0:
            NAME[ID].append(["NA", "0"])
        continue

    # If at least one Species identified, parse Species table and only keep best Hit
    HASH = d(lambda: d(lambda: d(str)))
    for sample in a[1:]:
        # skip if sample is empty
        if sample == "":
            continue

        # skip if sample is "N/A" or "NA"
        if "N/A" in sample or "NA" in sample:
            continue

        # replace spaces with "_"
        Spec = sample.split(" (")[0].replace(" ", "_")

        # only keep three letters from genus name and connect with dot to species names
        Spec = Spec.split("_")[0][:3]+"."+Spec.split("_")[1]

        # keep similarity values
        Sim = float(sample.split(" (")[1].split("%)")[0])

        # keep count values
        count = int(sample.split("count=")[1])

        # fill hash and recode similarity in percent
        HASH[ID][Sim][count] = [Spec, str(round(Sim/100, 2))]

    # obtain highest counts and similarities
    if len(HASH[ID]) == 0:
        # Only append if this ID hasn't been processed yet
        if ID not in NAME or len(NAME[ID]) == 0:
            NAME[ID].append(["NA", "0"])
        continue
    else:
        Sim = max(HASH[ID].keys())
        count = max(HASH[ID][Sim].keys())

        # for each ID, retain Spec with highest counts/similarity
        # Only append if this ID hasn't been processed yet or if this is a better hit
        if ID not in NAME or len(NAME[ID]) == 0:
            NAME[ID].append(HASH[ID][Sim][count])
        else:
            # Check if current hit is better than existing one
            existing_sim = float(NAME[ID][0][1]) if len(
                NAME[ID][0]) > 1 and NAME[ID][0][1] != "0" else 0
            current_sim = Sim / 100
            if current_sim > existing_sim:
                NAME[ID] = [HASH[ID][Sim][count]]  # Replace with better hit

    # rename outgroup samples
    for outgroup in OG:
        if outgroup not in a[0]:
            continue
        ID = "-".join(a[0].split("-")[:-1])
        EXT = a[0].split("-")[-1]
        SPEC = HASH[ID][Sim][count]
        OGnew.append(ID+"-"+"_".join(SPEC)+"-"+EXT)


# if --freqthreshold==0, only append names to COX1, etc.
if TEST == 1:
    tree = load_data(options.IN).readline()
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        # Keep track of replacements to avoid double-replacement
        replacements_made = set()
        for k, v in NAME.items():
            # split into ID and number of consensus
            ID = "-".join(k.split("-")[:-1])
            EXT = k.split("-")[-1]
            old_id = ID + "-" + EXT

            # Handle different formats of v[0]
            if len(v) > 0:
                if isinstance(v[0], list) and len(v[0]) > 0:
                    species_name = v[0][0] if v[0][0] != "NA" else "unidentified"
                elif isinstance(v[0], str):
                    species_name = v[0] if v[0] != "NA" else "unidentified"
                else:
                    species_name = "unidentified"
            else:
                species_name = "unidentified"

            new_id = ID + "-" + species_name + "-" + EXT

            # Only replace if we haven't already replaced this ID
            if old_id not in replacements_made and old_id in tree:
                tree = tree.replace(old_id, new_id)
                replacements_made.add(old_id)

else:
    # make more specific dictionary for each consensus per locus separately
    NAME2 = d(lambda: d(str))
    # if barcoding locus, keep species ID specific for each consensus number
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        for k, v in NAME.items():
            for i in range(len(v)):
                NAME2[k][str(i+1)] = "_".join(v[i])
    else:
        # if there are more than one consensus sequences for the barcoding locus, skip adding species name as a whole
        for k, v in NAME.items():
            if len(v) == 0:
                NAME2[k]["N"] = "NA"
                continue

            # Handle case where v contains strings instead of lists
            if isinstance(v[0], str):
                NAME2[k]["N"] = v[0] if v[0] != "NA" else "NA"
                continue

            # Extract species and similarity info
            try:
                Spec, Sim = list(zip(*v))
                if len(list(set(Spec))) > 1:
                    NAME2[k]["N"] = "NA"
                else:
                    NAME2[k]["N"] = "_".join(v[0]) if isinstance(
                        v[0], list) else str(v[0])
            except (ValueError, TypeError):
                # Handle malformed data
                NAME2[k]["N"] = "NA"
    # OK, now read the tree file
    tree = load_data(options.IN).readline()

    # If barcoding locus, attach species names for each consensus
    if options.PI in ["COX1", "ITS", "MATK_RBCL"]:
        # Keep track of replacements to avoid double-replacement
        replacements_made = set()
        for k, v in NAME2.items():
            for I, v1 in v.items():
                old_id = k + "-" + I
                new_id = k + "-" + v1 + "-" + I

                # Only replace if we haven't already replaced this ID
                if old_id not in replacements_made and old_id in tree:
                    tree = tree.replace(old_id, new_id)
                    replacements_made.add(old_id)
    else:
        # for all other, attach species name if only single consensus of diagnostic barcoding locus with species names
        replacements_made = set()
        for k, v in NAME2.items():
            old_id = k
            new_id = k + "-" + v["N"]

            # Only replace if we haven't already replaced this ID
            if old_id not in replacements_made and old_id in tree:
                tree = tree.replace(old_id, new_id)
                replacements_made.add(old_id)

# print new tree
out = open(options.IN, "w")
out.write(tree+"\n")

# print new IDs in outgrouplist
print(",".join(OGnew))

# # Add debugging output to identify potential issues
# print(f"Processing completed. Summary of NAME dictionary:", file=sys.stderr)
# for k, v in NAME.items():
#     if len(v) > 1:
#         print(
#             f"Warning: ID {k} has {len(v)} species assignments: {v}", file=sys.stderr)
#     elif len(v) == 1 and isinstance(v[0], list) and len(v[0]) > 0:
#         print(f"ID {k}: {v[0][0]}", file=sys.stderr)
