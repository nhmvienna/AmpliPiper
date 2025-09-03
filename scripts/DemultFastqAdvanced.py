#!/usr/bin/env python3

"""
DemultFastqAdvanced.py - Advanced demultiplexing with internal primer detection and multi-amplicon splitting

This enhanced version of DemultFastq.py can:
1. Detect primer pairs anywhere in the read (not just at the ends)
2. Identify multiple amplicons within a single read
3. Split reads containing multiple amplicons into separate sequences

Author: AmpliPiper development team
Date: August 2025
"""

import edlib
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from argparse import ArgumentParser
from datetime import datetime
import gzip
import sys
from collections import defaultdict as d
from collections import Counter

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to the file with raw reads", required=True)
argparse.add_argument(
    "-p", "--primers", help="Path to the csv file containing primers", required=True)
argparse.add_argument(
    "-o", "--output", help="Path to the output folder", required=True)
argparse.add_argument(
    "-rp", "--reads_percentage", help="Providing this command, the demultiplexer will output the top nr%% of reads", type=float, required=False, default=1.0)
argparse.add_argument(
    "-th", "--kthreshold", help="0.5 means that there should not be any differences between primer and target sequence. The higher this value the more liberal. Should not exceed 0.6", type=float, required=False, default=0.5)
argparse.add_argument(
    "-sr", "--sizerange", help="+/- size ", type=float, required=False, default=0.5)
argparse.add_argument(
    "-mr", "--minreads", help="minimum number of reads", type=float, required=False, default=0.5)
argparse.add_argument(
    "--disable-split-multi", help="Disable multi-amplicon splitting (for legacy compatibility)", action="store_true", default=False)
argparse.add_argument(
    "--min-amplicon-length", help="Minimum length for split amplicons", type=int, default=100)
argparse.add_argument(
    "--threads", help="Number of threads to use for processing", type=int, default=1)

args = argparse.parse_args()
inf = args.infile
prim = args.primers
TH = args.kthreshold
readper = args.reads_percentage
sizerange = args.sizerange
minreads = args.minreads
split_multi = not args.disable_split_multi
min_amplicon_length = args.min_amplicon_length
threads = args.threads


def load_data(infile):
    """Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)"""
    try:
        print("Reading data from input file...", file=sys.stderr)
        if infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(infile, "rt", encoding="latin-1")
        elif infile.endswith(".fastq"):
            y = open(infile, "r")
        else:
            raise ValueError("File is the wrong format")

        records = SeqIO.parse(y, "fastq")
        seq_dict = []  # Create a list to store everything from the file

        # Process records with progress bar (no pre-counting for speed)
        print("Loading sequences...", file=sys.stderr)
        for record in records:
            # Extract quality scores more efficiently
            quality_scores = record.letter_annotations["phred_quality"]
            avg_quality = sum(quality_scores) / len(quality_scores)

            # Store sequence data
            seq_dict.append([
                record.id,
                str(record.seq),
                record.format("fastq").split("\n")[3],
                round(avg_quality, 2)
            ])

        y.close()
        print(f"Loaded {len(seq_dict)} sequences", file=sys.stderr)
        return seq_dict

    except FileNotFoundError:
        print(
            f"File {infile} does not exist, proceeding with the analysis...", file=sys.stderr)
        return


def read_primer_table(csv):
    """Read primers table, provided as a csv file whose separatore MUST be comma and whose fields MUST be named and ordered as follows: ID,FWD,REV
    ID: contains the demultiplexing unit to which the primers are referred
    FWD and REV: idicate respectively the forward and reverse primer sequences"""
    print("Reading primers table...", file=sys.stderr)
    try:
        primers = pd.read_csv(csv)  # read primers table thanks to pandas
        fp = primers["FWD"]  # Forward sequences
        rp = primers["REV"]  # Reverse sequences
        ids = primers["ID"]  # Demultiplexing unit ID
        size = primers["SIZE"]  # Demultiplexing unit ID
        primer_dict = d(list)
        for i in range(len(ids)):
            # Assign to each demultiplexing unit ID its reverse and forward primers
            primer_dict[ids[i]] = [fp[i].upper(), rp[i].upper(), size[i]]
        print("Done", file=sys.stderr)
        return primer_dict
    except KeyError or ValueError:
        raise KeyError(
            "Fields of the csv do not comply with the requirements")


def reverse_complement(seq: str):
    """Returns the reverse complementary of a DNA sequence (also degenerate)"""
    dna = Seq(seq)
    return str(dna.reverse_complement())


def do_alignment_fast(F, R, SEQ, TH):
    """Fast version of the original do_alignment function"""
    correspondences = [("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("M", "A"), ("M", "C"), ("K", "G"), ("K", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("B", "C"), ("B", "G"), ("B", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), (
        "H", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]

    # Test forward primer alignment
    aln1 = edlib.align(
        F, SEQ,
        task="path",
        mode="HW",
        additionalEqualities=correspondences)

    # If the alignment is not significant, return false
    if aln1["editDistance"] == -1 or aln1["editDistance"] > len(F)*TH:
        return False, "NA"

    # Test reverse primer alignment
    aln2 = edlib.align(
        R, SEQ,
        task="path",
        mode="HW",
        additionalEqualities=correspondences)

    # Check if reverse primer also matches within threshold
    if aln2["editDistance"] != -1 and aln2["editDistance"] <= len(R)*TH:
        # Both primers match - return the locations
        if aln1["locations"] and aln2["locations"]:
            return True, (aln1["locations"][0], aln2["locations"][0])

    # If we get here, only forward primer matched or reverse primer failed
    return False, "NA"


def find_all_primer_positions_fast(primer, sequence, TH):
    """Fast primer matching using edlib's HW mode for speed"""
    max_mismatches = int(len(primer) * TH)

    # IUPAC correspondences for degenerate bases (same as original script)
    correspondences = [("R", "A"), ("R", "G"), ("Y", "C"), ("Y", "T"), ("M", "A"), ("M", "C"), ("K", "G"), ("K", "T"), ("S", "G"), ("S", "C"), ("W", "A"), ("W", "T"), ("B", "C"), (
        "B", "G"), ("B", "T"), ("D", "A"), ("D", "G"), ("D", "T"), ("H", "A"), ("H", "C"), ("H", "T"), ("V", "A"), ("V", "C"), ("V", "G"), ("N", "A"), ("N", "C"), ("N", "G"), ("N", "T")]

    # Use edlib's HW (semi-global) mode with IUPAC support
    aln = edlib.align(
        primer, sequence,
        task="locations",
        mode="HW",  # Semi-global alignment - much faster
        k=max_mismatches,  # Limit search space
        additionalEqualities=correspondences  # Handle IUPAC codes efficiently
    )

    positions = []
    if aln["editDistance"] != -1 and aln["editDistance"] <= max_mismatches:
        for start, end in aln["locations"]:
            positions.append((start, end + 1, aln["editDistance"]))

    return positions


def find_amplicons_in_read_fast(sequence, primer_dict, TH, expected_size, size_tolerance):
    """Fast amplicon finding with optimized search for multi-amplicon detection"""
    amplicons = []

    for primer_id, (fwd_primer, rev_primer, size) in primer_dict.items():
        # Pre-compute reverse complements
        fwd_rc = reverse_complement(fwd_primer)
        rev_rc = reverse_complement(rev_primer)

        # Find primer positions using fast method
        fwd_positions = find_all_primer_positions_fast(
            fwd_primer, sequence, TH)
        fwd_rc_positions = find_all_primer_positions_fast(fwd_rc, sequence, TH)
        rev_positions = find_all_primer_positions_fast(
            rev_primer, sequence, TH)
        rev_rc_positions = find_all_primer_positions_fast(rev_rc, sequence, TH)

        # Early exit if no primers found
        if not (fwd_positions or fwd_rc_positions or rev_positions or rev_rc_positions):
            continue

        # Check combinations more efficiently
        all_combinations = [
            (fwd_positions, rev_rc_positions, "fwd+rev_rc"),
            (fwd_rc_positions, rev_positions, "fwd_rc+rev"),
        ]

        for fwd_pos_list, rev_pos_list, orientation in all_combinations:
            if not fwd_pos_list or not rev_pos_list:
                continue

            for fwd_start, fwd_end, fwd_dist in fwd_pos_list:
                for rev_start, rev_end, rev_dist in rev_pos_list:
                    # Quick boundary check
                    if fwd_end <= rev_start:  # Forward before reverse
                        amplicon_start = fwd_start
                        amplicon_end = rev_end
                        insert_start = fwd_end
                        insert_end = rev_start
                    elif rev_end <= fwd_start:  # Reverse before forward
                        amplicon_start = rev_start
                        amplicon_end = fwd_end
                        insert_start = rev_end
                        insert_end = fwd_start
                    else:
                        continue  # Overlapping primers, skip

                    amplicon_length = amplicon_end - amplicon_start
                    insert_length = insert_end - insert_start

                    # Calculate size tolerance for this specific primer pair
                    if size_tolerance < 1.0:
                        # If < 1.0, treat as proportion of expected size
                        tolerance = int(size * size_tolerance)
                    else:
                        # If >= 1.0, treat as absolute bp value
                        tolerance = int(size_tolerance)

                    # Quick size check
                    if (amplicon_length >= size - tolerance and
                        amplicon_length <= size + tolerance and
                            insert_length >= min_amplicon_length):

                        amplicons.append({
                            'primer_id': primer_id,
                            'start': amplicon_start,
                            'end': amplicon_end,
                            'insert_start': insert_start,
                            'insert_end': insert_end,
                            'length': amplicon_length,
                            'insert_length': insert_length,
                            'orientation': orientation,
                            'fwd_distance': fwd_dist,
                            'rev_distance': rev_dist,
                            'expected_size': size
                        })

    # Sort by position and remove overlaps efficiently
    if amplicons:
        amplicons.sort(key=lambda x: x['start'])
        non_overlapping = [amplicons[0]]  # Start with first amplicon

        for amp in amplicons[1:]:
            # Check only against last non-overlapping amplicon
            if amp['start'] >= non_overlapping[-1]['end']:
                non_overlapping.append(amp)

        return non_overlapping

    return []


def demultiplex_multi_amplicon(sequences_dict, primer_dict, TH, sizerange, split_multi, min_amplicon_length):
    """Advanced demultiplexing with multi-amplicon detection and splitting"""
    demultiplexed = d(list)
    multi_amplicon_stats = d(int)
    multi_amplicon_details = []  # Store detailed information about multi-amplicon reads

    for seq_data in sequences_dict:
        seq_id, sequence, quality, avg_quality = seq_data

        # Find all amplicons in this read
        # Calculate size tolerance properly - use the sizerange parameter directly
        # If sizerange >= 1.0, treat as absolute bp; if < 1.0, treat as proportion
        all_amplicons = find_amplicons_in_read_fast(
            sequence, primer_dict, TH, 0, sizerange)

        if not all_amplicons:
            continue

        # Track multi-amplicon reads
        if len(all_amplicons) > 1:
            multi_amplicon_stats[len(all_amplicons)] += 1

            # Store detailed information about this multi-amplicon read
            loci_order = [amp['primer_id'] for amp in all_amplicons]
            positions = [(amp['start'], amp['end']) for amp in all_amplicons]
            lengths = [amp['length'] for amp in all_amplicons]
            insert_lengths = [amp['insert_length'] for amp in all_amplicons]

            multi_amplicon_details.append({
                'read_id': seq_id,
                'num_amplicons': len(all_amplicons),
                'loci_order': ';'.join(loci_order),
                'positions': ';'.join([f"{start}-{end}" for start, end in positions]),
                'amplicon_lengths': ';'.join(map(str, lengths)),
                'insert_lengths': ';'.join(map(str, insert_lengths)),
                'read_length': len(sequence),
                'avg_quality': avg_quality
            })

        if split_multi and len(all_amplicons) > 1:
            # Split multi-amplicon reads
            for i, amplicon in enumerate(all_amplicons):
                # Create new sequence ID for split amplicon
                new_seq_id = f"{seq_id}_amplicon_{i+1}_{amplicon['primer_id']}"

                # Extract amplicon sequence and quality
                amp_sequence = sequence[amplicon['start']:amplicon['end']]
                amp_quality = quality[amplicon['start']:amplicon['end']]

                # Calculate new average quality
                if amp_quality:
                    amp_avg_quality = sum(
                        ord(c) - 33 for c in amp_quality) / len(amp_quality)
                else:
                    amp_avg_quality = avg_quality

                # Create new sequence data
                new_seq_data = [new_seq_id, amp_sequence,
                                amp_quality, amp_avg_quality]

                # Add to demultiplexed results
                demultiplexed[amplicon['primer_id']].append((
                    new_seq_data, (0, len(amp_sequence))
                ))
        else:
            # Use the best amplicon (first one after sorting)
            best_amplicon = all_amplicons[0]
            demultiplexed[best_amplicon['primer_id']].append((
                seq_data, (best_amplicon['start'], best_amplicon['end'])
            ))

    # Print multi-amplicon statistics
    if multi_amplicon_stats:
        print("\n=== Multi-amplicon Read Statistics ===", file=sys.stderr)
        total_multi = sum(multi_amplicon_stats.values())
        print(
            f"Total reads with multiple amplicons: {total_multi}", file=sys.stderr)
        for num_amplicons, count in sorted(multi_amplicon_stats.items()):
            print(
                f"  Reads with {num_amplicons} amplicons: {count}", file=sys.stderr)

    return demultiplexed, multi_amplicon_details


def demultiplex_like_original(sequences_dict, primer_dict, TH, sizerange):
    """Demultiplexing using the original algorithm but with optimizations"""
    demultiplexed = d(list)
    matched_sequences = 0

    # Use the same loop structure as original
    for primername, primers in primer_dict.items():
        for SeqDATA in sequences_dict:
            FWD, REV, SIZE = primers

            # exclude sequences outside sizerange (same as original)
            # Convert sizerange properly - it's a proportion, so multiply by SIZE
            size_tolerance = int(
                SIZE * sizerange) if sizerange < 1.0 else int(sizerange)

            if len(SeqDATA[1]) < SIZE - size_tolerance or len(SeqDATA[1]) > SIZE + size_tolerance:
                continue

            # Use the same alignment calls as original
            T1, D1 = do_alignment_fast(
                FWD, reverse_complement(REV), SeqDATA[1], TH)
            T2, D2 = do_alignment_fast(
                REV, reverse_complement(FWD), SeqDATA[1], TH)

            if T1:
                ALL = sorted([i for sub in D1 for i in sub])
                LENGTH = max(ALL)-min(ALL)
                if LENGTH > SIZE - size_tolerance:
                    demultiplexed[primername].append((
                        SeqDATA, (min(ALL), max(ALL))))
                    matched_sequences += 1
            elif T2:
                ALL = sorted([i for sub in D2 for i in sub])
                LENGTH = max(ALL)-min(ALL)
                if LENGTH > SIZE - size_tolerance:
                    demultiplexed[primername].append((
                        SeqDATA, (min(ALL), max(ALL))))
                    matched_sequences += 1

    print(
        f"Matched {matched_sequences} sequence-primer combinations", file=sys.stderr)
    return demultiplexed


def do_alignment(F, R, SEQ, TH):
    """Legacy function for compatibility - now uses fast alignment"""
    return do_alignment_fast(F, R, SEQ, TH)


def demultiplex_advanced(infile, csv):
    """Advanced demultiplexing function - chooses algorithm based on split_multi flag"""
    sequences_dict = load_data(infile)  # Load data from infile
    primer_dict = read_primer_table(csv)  # Read primers table

    print("Total number of Sequences:" +
          str(len(sequences_dict)), file=sys.stderr)
    print("Number of primer pairs:" + str(len(primer_dict)), file=sys.stderr)

    multi_amplicon_details = []  # Initialize empty list for legacy compatibility

    if not split_multi:
        print(
            f"Using legacy algorithm (multi-amplicon splitting disabled)", file=sys.stderr)

        # Use the original algorithm for compatibility
        demultiplexed = demultiplex_like_original(
            sequences_dict, primer_dict, TH, sizerange)
    else:
        print(f"Using advanced multi-amplicon algorithm with splitting enabled", file=sys.stderr)
        print(f"  Min amplicon length: {min_amplicon_length}", file=sys.stderr)

        # Use advanced algorithm with multi-amplicon detection
        demultiplexed, multi_amplicon_details = demultiplex_multi_amplicon(
            sequences_dict, primer_dict, TH, sizerange, split_multi, min_amplicon_length)

    # Print results summary
    total_matches = sum(len(v) for v in demultiplexed.values())
    print(f"Total matches found: {total_matches}", file=sys.stderr)
    for primer_name, reads in demultiplexed.items():
        if reads:
            print(
                f"Primer pair {primer_name}: {len(reads)} reads", file=sys.stderr)

    return demultiplexed, multi_amplicon_details


def demultiplex(infile, csv):
    """Main demultiplexing function - now uses original algorithm"""
    return demultiplex_advanced(infile, csv)


def sort_filter(demultiplexed, ReadCount, Output):
    """Sort and filter demultiplexed reads - fixed to match original behavior"""
    print("Processing and writing output files...", file=sys.stderr)

    # Pre-filter valid primers (same as original)
    valid_primers = [(k, v) for k, v in demultiplexed.items()
                     if len(v) >= int(minreads)]

    if not valid_primers:
        print("No primer pairs met minimum read threshold", file=sys.stderr)
        return

    for k, v in valid_primers:
        # Use the same logic as original script
        BQdict = d(list)
        for items, RANGE in v:
            BQdict[items[0]] = items[3]

        # Calculate number of reads to output
        # If ReadCount is <= 1.0, treat as percentage; otherwise as absolute number
        if ReadCount <= 1.0:
            num_reads_to_output = max(1, int(len(v) * ReadCount))
        else:
            num_reads_to_output = int(ReadCount)

        ToPrint = [x for x, y in Counter(
            dict(BQdict)).most_common(num_reads_to_output)]

        print(f"{k}: {len(ToPrint)} reads", file=sys.stderr)

        # Write output file (same as original)
        output_path = f"{Output}/{k}.fastq"
        with open(output_path, "wt") as OUT:
            for items, RANGE in v:
                ID, SEQ, BQ, AvBQ = items
                if ID in ToPrint:
                    OUT.write("@"+ID+"\n"+SEQ[RANGE[0]:RANGE[1]] +
                              "\n+\n"+BQ[RANGE[0]:RANGE[1]]+"\n")


def write_multi_amplicon_details(multi_amplicon_details, output_dir):
    """Write detailed information about multi-amplicon reads to a CSV file"""
    if not multi_amplicon_details:
        print("No multi-amplicon reads found to write.", file=sys.stderr)
        return

    output_file = f"{output_dir}/multi_amplicon_reads.csv"

    try:
        with open(output_file, 'w') as f:
            # Write header
            f.write(
                "Read_ID,Num_Amplicons,Loci_Order,Positions,Amplicon_Lengths,Insert_Lengths,Read_Length,Avg_Quality\n")

            # Write data for each multi-amplicon read
            for read_info in multi_amplicon_details:
                f.write(
                    f"{read_info['read_id']},{read_info['num_amplicons']},{read_info['loci_order']},")
                f.write(
                    f"{read_info['positions']},{read_info['amplicon_lengths']},{read_info['insert_lengths']},")
                f.write(
                    f"{read_info['read_length']},{read_info['avg_quality']:.2f}\n")

        print(
            f"Multi-amplicon details written to: {output_file}", file=sys.stderr)
        print(
            f"Total multi-amplicon reads: {len(multi_amplicon_details)}", file=sys.stderr)

        # Print summary statistics
        loci_combinations = {}
        for read_info in multi_amplicon_details:
            loci_order = read_info['loci_order']
            loci_combinations[loci_order] = loci_combinations.get(
                loci_order, 0) + 1

        print("Multi-amplicon loci combinations:", file=sys.stderr)
        for combination, count in sorted(loci_combinations.items(), key=lambda x: x[1], reverse=True):
            print(f"  {combination}: {count} reads", file=sys.stderr)

    except Exception as e:
        print(f"Error writing multi-amplicon details: {e}", file=sys.stderr)


if __name__ == "__main__":
    # Call the functions
    start = datetime.now()

    if not split_multi:
        print(f"Starting demultiplexing with legacy algorithm (multi-amplicon splitting disabled)...", file=sys.stderr)
    else:
        print(f"Starting advanced demultiplexing with multi-amplicon detection and splitting...", file=sys.stderr)

    print(f"Parameters:", file=sys.stderr)
    print(f"  K-threshold: {TH}", file=sys.stderr)
    print(f"  Size range: ±{int(sizerange)}", file=sys.stderr)
    print(f"  Min reads: {int(minreads)}", file=sys.stderr)
    print(f"  Multi-amplicon splitting: {split_multi}", file=sys.stderr)
    print(f"  Min amplicon length: {min_amplicon_length}", file=sys.stderr)

    demult, multi_amplicon_details = demultiplex(inf, prim)
    sort_filter(demult, readper, args.output)

    # Write multi-amplicon details to file
    write_multi_amplicon_details(multi_amplicon_details, args.output)

    end = datetime.now()
    print('Duration of demultiplexing process: {}'.format(
        end - start), file=sys.stderr)

    # Summary statistics
    total_amplicons = sum(len(v) for v in demult.values())
    print(f"\n=== Final Summary ===", file=sys.stderr)
    print(f"Total amplicons found: {total_amplicons}", file=sys.stderr)
    print(f"Primer pairs with hits:", file=sys.stderr)
    for primer_id, amplicons in demult.items():
        if amplicons:
            print(f"  {primer_id}: {len(amplicons)} amplicons", file=sys.stderr)
