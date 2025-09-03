import sys
import csv
import os
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--path", dest="IN", help="Input path")
parser.add_option("--primer", dest="PR", help="primers file")
parser.add_option("--samples", dest="SA", help="samples file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--id-column", dest="ID_COL", help="Name of ID column (default: auto-detect)", default=None)
parser.add_option("--delimiter", dest="DELIM", help="CSV delimiter (default: comma)", default=",")
parser.add_option("--verbose", dest="VERBOSE", help="Verbose output", action="store_true", default=False)

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or uncompressed file or from STDIN'''
    import gzip
    if x == "-":
        return sys.stdin
    elif x.endswith(".gz"):
        return gzip.open(x, "rt", encoding="utf-8-sig")
    else:
        return open(x, "r", encoding="utf-8-sig")


def get_csv_reader(file_handle, delimiter=","):
    '''Create a CSV reader with proper dialect detection'''
    return csv.reader(file_handle, delimiter=delimiter)


def detect_id_column(csv_file, delimiter=",", id_column_name=None):
    '''Detect the ID column in a CSV file, handling BOM and various formats'''
    try:
        with load_data(csv_file) as f:
            reader = get_csv_reader(f, delimiter)
            header = next(reader)
            
            # Clean up headers (strip whitespace and potential BOM artifacts)
            clean_header = [col.strip().strip('\ufeff').strip('ï»¿') for col in header]
            
            if options.VERBOSE:
                print(f"Raw header: {header}", file=sys.stderr)
                print(f"Clean header: {clean_header}", file=sys.stderr)
            
            # If specific column name provided, find it
            if id_column_name:
                for i, col in enumerate(clean_header):
                    if col.upper() == id_column_name.upper():
                        return i
                raise ValueError(f"Column '{id_column_name}' not found in {clean_header}")
            
            # Auto-detect ID column (look for common ID column names)
            id_candidates = ['ID', 'NAME', 'SAMPLE', 'PRIMER', 'LOCUS']
            for candidate in id_candidates:
                for i, col in enumerate(clean_header):
                    if col.upper() == candidate:
                        return i
            
            # Default to first column
            return 0
            
    except Exception as e:
        if options.VERBOSE:
            print(f"Error detecting ID column in {csv_file}: {e}", file=sys.stderr)
        return 0


def getIDS(csv_file, delimiter=",", id_column_name=None):
    ''' get IDS from CSV file with flexible column detection '''
    ids = []
    
    try:
        id_col_index = detect_id_column(csv_file, delimiter, id_column_name)
        
        with load_data(csv_file) as f:
            reader = get_csv_reader(f, delimiter)
            header = next(reader)  # Skip header
            
            for row_num, row in enumerate(reader, start=2):
                if len(row) > id_col_index:
                    id_value = row[id_col_index].strip()
                    if id_value:  # Only add non-empty IDs
                        ids.append(id_value)
                elif options.VERBOSE:
                    print(f"Warning: Row {row_num} has insufficient columns: {row}", file=sys.stderr)
                    
    except Exception as e:
        print(f"Error reading CSV file {csv_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    if options.VERBOSE:
        print(f"Found {len(ids)} IDs in {csv_file}: {ids[:5]}{'...' if len(ids) > 5 else ''}", file=sys.stderr)
    
    return ids


TEST = {2:
        {2: [(0.5, 0.5)],
         3: [(0.3333333333333333, 0.6666666666666666)],
         4: [(0.5, 0.5), (0.25, 0.75)]},
        3:
        {3: [
            (0.3333333333333333, 0.3333333333333333, 0.3333333333333333)],
            4: [(0.25, 0.25, 0.5)]},
        4: {4: [(0.25, 0.25, 0.25, 0.25)]}}

# Perform chi-square goodness-of-fit test
# chi2_stat, p_val = chisquare(observed, f_exp=expected)

# Validate required arguments
if not options.IN:
    print("Error: --path is required", file=sys.stderr)
    sys.exit(1)
if not options.PR:
    print("Error: --primer is required", file=sys.stderr)
    sys.exit(1)
if not options.SA:
    print("Error: --samples is required", file=sys.stderr)
    sys.exit(1)

# Check if files exist
for file_path, name in [(options.PR, "primers"), (options.SA, "samples")]:
    if not os.path.exists(file_path):
        print(f"Error: {name} file not found: {file_path}", file=sys.stderr)
        sys.exit(1)

if options.VERBOSE:
    print(f"Reading primers from: {options.PR}", file=sys.stderr)
    print(f"Reading samples from: {options.SA}", file=sys.stderr)
    print(f"Using delimiter: '{options.DELIM}'", file=sys.stderr)

print("ID,Locus,HaplotypesCount,TotalReads,ReadCount,FrequencyOfTotal,HaplotypesLengths")
for LOCUS in getIDS(options.PR, options.DELIM, options.ID_COL):
    for IND in getIDS(options.SA, options.DELIM, options.ID_COL):
        if not os.path.exists(f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/results.txt"):
            RD, PC, UR, SEQL = [], "NA", "NA", "NA"
            print(",".join([IND,
                            LOCUS,
                            "NA",
                            "NA",
                            "NA",
                            "NA",
                            "NA"]))
        else:
            FILE = load_data(
                f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/results.txt")
            RD, PC = [], []
            for l in FILE:
                if l.startswith("-->"):
                    RD.append(l.split("contains ")[1].split(" sequences")[0])
                    PC.append(
                        str(round(float(l.split("sequences (")[1].split("% of total")[0])/100, ndigits=3)))
                if l.startswith("- used_reads"):
                    URline = l.rstrip().split("= ")
                    if len(URline) > 1:
                        UR = URline[1]
                    else:
                        UR = "0"
            if os.path.exists(f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/{LOCUS}_consensussequences.fasta"):
                FILE = load_data(
                    f"{options.IN}/results/consensus_seqs/{IND}/{LOCUS}/{LOCUS}_consensussequences.fasta")
                SEQL = []
                for l in FILE:
                    if l.startswith(">"):
                        continue
                    SEQL.append(str(len(l.rstrip())))
            else:
                SEQL = ["NA"]
            HAP = len(RD)
            if RD == []:
                RD = ["NA"]
            if PC == []:
                PC = ["NA"]

            print(",".join([IND,
                            LOCUS,
                            str(HAP),
                            UR,
                            "/".join(RD),
                            "/".join(PC),
                            "/".join(SEQL)]))
