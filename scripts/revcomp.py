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

(options, args) = parser.parse_args()
parser.add_option_group(group)

code = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C', '-': ''}

DATA = d(str)
try:
    with open(options.IN, 'r') as in_file:
        for line in in_file:
            if line != '\n':  # skip empty lines
                line = line.strip()  # Remove new line character (I'm working on windows)
                if line.startswith('>'):
                    head = line.rstrip()
                else:
                    DATA[head] += line.rstrip()

    for k, v in DATA.items():
        print(k)
        print(''.join(code[nuc.upper()] for nuc in v[::-1]))
except FileNotFoundError:
    print(f"File {options.IN} does not exist, proceeding with the analysis...", file=sys.stderr)
