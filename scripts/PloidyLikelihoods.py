import numpy as np
from scipy.stats import multinomial
import sys
from collections import defaultdict as d


def likelihood(observed, expected):
    # Multinomial likelihood: multinomial.pmf(observed counts | expected proportions)
    # print(observed, np.sum(observed), expected)
    return multinomial.pmf(observed, n=np.sum(observed), p=expected)


def estimatePloidy(OBS):
    PLOIDY = d(lambda: d(list))
    COUNTl = len(OBS)
    for k, v in TEST[COUNTl].items():
        for EXP in v:
            EXP /= np.sum(EXP)
            LH = likelihood(sorted(OBS), EXP)
            PLOIDY[LH]["samples"] = sorted(OBS)
            PLOIDY[LH]["PLOIDY"] = [COUNTl, k]

    while (len(OBS) >= 3):
        Total = sum(OBS)
        OBS = sorted(OBS)[1:]
        Freq = [x/sum(OBS) for x in OBS]
        OBS = [round(x*Total, 0) for x in Freq]
        print(OBS)
        COUNTl = len(OBS)

        for k, v in TEST[COUNTl].items():
            for EXP in v:
                EXP /= np.sum(EXP)
                LH = likelihood(sorted(OBS), EXP)
                PLOIDY[LH]["samples"].append(sorted(OBS))
                PLOIDY[LH]["PLOIDY"].append([COUNTl, k])
    return PLOIDY


TEST = {2:
        {1: [np.array([0.0, 1.0])],
         2: [np.array([1/2, 1/2])],
         3: [np.array([1/3, 2/3])],
         4: [np.array([1/2, 1/2]), np.array([1/4, 3/4])]},
        3:
        {1: [np.array([0.0, 0.0, 1.0])],
         2: [np.array([0.0, 1/2, 1/2])],
            3: [np.array([1/3, 1/3, 1/3])],
            4: [np.array([1/4, 1/4, 1/2])]},
        4: {1: [np.array([0.0, 0.0, 0.0, 1.0])],
            2: [np.array([0.0, 0.0, 1/2, 1/2])],
            3: [np.array([0, 1/3, 1/3, 1/3])],
            4: [np.array([1/4, 1/4, 1/4, 1/4])]}}


# Observed haplotype frequencies
OBS = [1203, 1232, 2]  # e.g., from read counts

if len(OBS) > 4:
    OBS = sorted(OBS)[:4]

NEW = estimatePloidy(OBS)
print(NEW.keys())
BEST = max(NEW.keys())
print(NEW[BEST])


# dictionary to store the consensus IDs to keep for each locus and sample
KEEP = d(lambda: d(list))
PLOIDY = d(lambda: d(list))
for l in load_data(options.IN):
    # skip header
    if l.startswith("ID,"):
        continue
    ID, Locus, HaplotypesCount, TotalReads, ReadCount, FrequencyOfTotal, HaplotypesLengths = l.rstrip().split(",")
    # skip if NA
    if ReadCount == "NA":
        continue
    # keep single haplotype if homozygous
    Reads = ReadCount.split("/")
    if len(Reads) == 1:
        KEEP[Locus][ID].append(0)
        continue
    # calculate total reads across all consensus sequences
    GrandTotal = 0
    for C in Reads:
        GrandTotal += int(C)
    # calculate frequencies and remove consensus sequences with rel. Freq <10%
    Freq = d(float)
    for C in range(len(Reads)):
        RF = int(Reads[C])/GrandTotal
        print(C, RF, int(Reads[C]))
        if RF < 0.1:
            continue
        Freq[C] = int(Reads[C])
    # if only one haplotype left, keep and proceed
    if len(Freq) == 1:
        KEEP[Locus][ID].append(0)
        continue

    Total = sum(Freq.values())
    print(Total)
    Ploidies = 0
