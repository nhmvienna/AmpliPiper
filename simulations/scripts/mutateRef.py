import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import random

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--primers", dest="PR", help="fwd and rev comma-sep")
parser.add_option("--HowManyAs", dest="HMA", help="fwd and rev comma-sep")
parser.add_option("--divergence", dest="DIV",
                  help="list of differences, comma-sep")
parser.add_option("--output", dest="OUT", help="Output file")

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


def mutate_sequence(sequence, mutation_rate, seed=None):
    """
    Introduces random mutations into a given sequence with an option to set a seed for reproducibility.

    Parameters:
    - sequence (str): The original sequence (e.g., DNA sequence).
    - mutation_rate (float): The proportion of the sequence to mutate (between 0 and 1).
    - seed (int, optional): Seed for the random number generator. Default is None.

    Returns:
    - str: The mutated sequence.
    """

    if not 0 <= mutation_rate <= 1:
        raise ValueError("Mutation rate must be between 0 and 1.")

    # Set the seed for reproducibility if provided
    if seed is not None:
        random.seed(seed)

    # Possible bases or characters for mutation
    bases = list(set(sequence))

    # Convert sequence to a list to allow mutations
    seq_list = list(sequence)
    seq_len = len(seq_list)

    # Calculate the number of mutations based on the mutation rate
    num_mutations = int(seq_len * mutation_rate)

    # Randomly select positions to mutate
    mutation_positions = random.sample(range(seq_len), num_mutations)

    for pos in mutation_positions:
        original_base = seq_list[pos]
        # Select a random base that is different from the original
        new_base = random.choice([b for b in bases if b != original_base])
        seq_list[pos] = new_base

    return ''.join(seq_list)


# Example usage
# sequence = "ACGTACGTACGT"
# mutation_rate = 0.2  # 20% of the sequence will mutate
# seed = 42  # Setting seed for reproducibility

# mutated_sequence = mutate_sequence(sequence, mutation_rate, seed=seed)
# print("Original sequence:", sequence)
# print("Mutated sequence:", mutated_sequence)

Fwd, Rev = options.PR.split(",")
As = "".join(["A"*int(options.HMA)])

SEQ = ""

for l in load_data(options.IN):
    if l.startswith(">"):
        ID = l.rstrip()[1:]
        continue
    SEQ += l.rstrip()

SEQ2 = SEQ[len(Fwd):-len(Rev)]
FWD = SEQ[:len(Fwd)]
REV = SEQ[-len(Rev):]

for i in options.DIV.split(","):
    SEQ3 = mutate_sequence(SEQ2, float(i), 42)
    out1 = open(options.OUT+"_short_"+i+".fa", "wt")
    out2 = open(options.OUT+"_long_"+i+".fa", "wt")
    ToPrintShort = FWD+SEQ3+REV
    out1.write(">"+ID+"\n"+ToPrintShort+"\n")
    ToPrintLong = As+FWD+SEQ3+REV+As
    out2.write(">"+ID+"\n"+ToPrintLong+"\n")
