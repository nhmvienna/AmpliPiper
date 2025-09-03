from optparse import OptionParser, OptionGroup
import numpy as np
from scipy.stats import multinomial
from collections import defaultdict as d
import os
import gzip
import sys

# Author: Martin Kapun (modified for robustness and additional functionality)

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file"
parser = OptionParser(usage=usage)
group = OptionGroup(parser, 'Description')
group.add_option("--input", dest="IN",
                 help="Input CSV file containing genomic data.")
group.add_option("--path", dest="PA", help="Path to consensus sequence files.")
group.add_option("--FreqTH", dest="FT",
                 help="Minimum frequency threshold for filtering haplotypes.", default=0.0, type=float)
group.add_option("--output", dest="OUT", help="Output directory for results.")
parser.add_option_group(group)

(options, args) = parser.parse_args()

#########################################################   FUNCTIONS   #########################################################################


def load_data(x):
    """
    Load data from a file or STDIN. Supports gzipped files.
    Args:
        x (str): Path to the file or '-' for STDIN.
    Returns:
        file object: File handle for reading.
    """
    if x == "-":
        return sys.stdin
    elif x.endswith(".gz"):
        return gzip.open(x, "rt", encoding="latin-1")
    else:
        return open(x, "r", encoding="latin-1")


def likelihood(observed, expected):
    """
    Calculate multinomial likelihood for observed counts given expected proportions.
    Args:
        observed (list): Observed allele counts.
        expected (list): Expected allele proportions.
    Returns:
        float: Multinomial likelihood.
    """
    return multinomial.pmf(observed, n=np.sum(observed), p=expected)


def estimatePloidy(OBS, output, Locus, ID):
    """
    Estimate ploidy based on observed haplotype counts.
    Args:
        OBS (list): Observed haplotype counts.
        output (file): File handle for writing likelihood data.
        Locus (str): Locus identifier.
        ID (str): Sample identifier.
    Returns:
        defaultdict: Ploidy likelihood results.
    """
    PLO = d(lambda: d(list))
    OBS = [int(x) for x in OBS]
    COUNTl = len(OBS)
    if COUNTl not in TEST:
        return PLO  # Return empty if ploidy model is undefined
    for k, v in TEST[COUNTl].items():
        for EXP in v:
            EXP /= np.sum(EXP)
            LH = likelihood(sorted(OBS), EXP)
            PLO[LH]["samples"] = sorted(OBS)
            PLO[LH]["PLOIDY"].append([str(COUNTl), str(k)])
            output.write(",".join(
                [ID, Locus, "/".join([str(x) for x in OBS]), str(COUNTl), str(k), str(LH)]) + "\n")
    return PLO


#########################################################   TEST MODELS   #########################################################################
TEST = {
    2: {1: [np.array([0.0, 1.0])],
        2: [np.array([0.5, 0.5])],
        3: [np.array([1/3, 2/3])],
        4: [np.array([0.5, 0.5]), np.array([0.25, 0.75])]},
    3: {1: [np.array([0.0, 0.0, 1.0])],
        2: [np.array([0.0, 0.5, 0.5])],
        3: [np.array([1/3, 1/3, 1/3])],
        4: [np.array([0.25, 0.25, 0.5])]},
    4: {1: [np.array([0.0, 0.0, 0.0, 1.0])],
        2: [np.array([0.0, 0.0, 0.5, 0.5])],
        3: [np.array([0.0, 1/3, 1/3, 1/3])],
        4: [np.array([0.25, 0.25, 0.25, 0.25])]}
}

# TEST = {
#     2: {2: [np.array([0.5, 0.5])],
#         3: [np.array([1/3, 2/3])],
#         4: [np.array([0.5, 0.5]), np.array([0.25, 0.75])],
#         5: [np.array([1/5, 4/5]), np.array([2/5, 3/5])],
#         6: [np.array([1/6, 5/6]), np.array([2/6, 4/6]), np.array([3/6, 3/6])]},
#     3: {3: [np.array([1/3, 1/3, 1/3])],
#         4: [np.array([0.25, 0.25, 0.5])],
#         5: [np.array([1/5, 1/5, 3/5]), np.array([1/5, 2/5, 2/5])],
#         6: [np.array([1/6, 1/6, 4/6]), np.array([1/6, 2/6, 3/6])]},
#     4: {4: [np.array([0.25, 0.25, 0.25, 0.25])],
#         5: [np.array([1/5, 1/5, 1/5, 2/5])],
#         6: [np.array([1/6, 1/6, 1/6, 3/6]), np.array([1/6, 1/6, 2/6, 2/6])]},
#     5: {5: [np.array([1/5, 1/5, 1/5, 1/5, 1/5])],
#         6: [np.array([1/6, 1/6, 1/6, 1/6, 2/6])]},
#     6: {6: [np.array([1/6, 1/6, 1/6, 1/6, 1/6, 1/6])]}
# }

#########################################################   MAIN SCRIPT   #########################################################################
FreqTH = float(options.FT)
output = open(options.IN.split(".csv")[0] + ".likelihoods", "wt")
output.write("ID,Locus,Reads,Alleles,Ploidy,Likelihood\n")

KEEP = d(lambda: d(list))
PLOIDY = d(lambda: d(list))
FREQ = d(lambda: d(lambda: d(int)))
for l in load_data(options.IN):
    if l.startswith("ID,"):
        continue  # Skip header

    ID, Locus, HaplotypesCount, TotalReads, ReadCount, FrequencyOfTotal, HaplotypesLengths = l.rstrip().split(",")
    if ReadCount == "NA":
        continue

    Reads = ReadCount.split("/")
    READID = dict(zip(Reads, range(len(Reads))))

    GrandTotal = sum(map(int, Reads))
    Keep = []
    NewTotal = 0
    for C in range(len(Reads)):
        if int(Reads[C]) / GrandTotal >= FreqTH:
            Keep.append(C)
            NewTotal += int(Reads[C])

    # Handle single allele case
    if len(Keep) == 1:
        KEEP[Locus][ID].append(Keep[0])
        # Use "NA" to indicate no ploidy testing
        PLOIDY[Locus][ID].append([1, 1])
        FREQ[Locus][ID][0] = 1
        continue

    # Skip ploidy testing if more than 4 haplotypes pass the threshold
    if len(Keep) > 4:
        for C in Keep:
            KEEP[Locus][ID].append(C)
            FREQ[Locus][ID][C] = int(Reads[C]) / NewTotal
        PLOIDY[Locus][ID].append([len(Keep), "NA"])
        continue

    # Estimate ploidy otherwise
    NEW = estimatePloidy([Reads[x] for x in Keep], output, Locus, ID)
    BEST = max(NEW.keys(), default=None)

    if BEST is not None:
        for sample in list(set(NEW[BEST]["samples"])):
            KEEP[Locus][ID].append(READID[str(sample)])
            FREQ[Locus][ID][READID[str(sample)]] = int(
                Reads[READID[str(sample)]]) / NewTotal
        for ploidy in NEW[BEST]["PLOIDY"]:
            PLOIDY[Locus][ID].append(ploidy)

# Write ploidy results
Out = open(options.IN + ".ploidy", "wt")
for l in load_data(options.IN):
    if l.startswith("ID,"):
        Out.write(l.rstrip() + ",NoAlleles,ExpectedPloidy\n")
        continue

    ID, Locus, HaplotypesCount, TotalReads, ReadCount, FrequencyOfTotal, HaplotypesLengths = l.rstrip().split(",")
    if ReadCount == "NA":
        Out.write(l.rstrip() + ",NA,NA\n")
        continue

    if ID in PLOIDY[Locus]:
        NoA, EP = list(zip(*PLOIDY[Locus][ID]))
        Out.write(l.rstrip() + "," + "/".join(list(set([str(x) for x in NoA]))) + "," +
                  "/".join([str(x) for x in EP]) + "\n")
    else:
        Out.write(l.rstrip() + ",NA,NA\n")

Out.close()

# Write filtered FASTA files
for Locus, v in KEEP.items():
    directory = os.path.join(options.OUT, Locus)
    os.makedirs(directory, exist_ok=True)
    Out = open(os.path.join(directory, Locus + ".fasta"), "wt")
    for ID, v1 in v.items():
        FILE = os.path.join(options.PA, ID, Locus, Locus +
                            "_consensussequences.fasta")
        FASTA = d(str)
        C = ""
        for l in load_data(FILE):
            if l.startswith(">"):
                if C == "":
                    C = 0
                else:
                    C += 1
                continue
            FASTA[C] += l
        C = 1
        for key in v1:
            # print(key, FREQ[Locus][ID][key])
            if FreqTH == 0:
                Out.write(">" + ID + "-Freq_" +
                          str(round(FREQ[Locus][ID][key], 2)) + "-" + str(C)+"\n")
            else:
                Out.write(">" + ID + "-" + str(C) + "\n")
            Out.write(FASTA[key])
            C += 1
    Out.close()
