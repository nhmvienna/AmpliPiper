from pygenomeviz import GenomeViz
from argparse import ArgumentParser
import pandas as pd
from Bio import SeqIO
import gzip
import sys
import os

# Create an ArgumentParser to handle command-line arguments
argparse = ArgumentParser()
argparse.add_argument(
    "-v", "--vcf_file", help="Path to the vcf file with SNPs variants", required=True
)
argparse.add_argument(
    "-r",
    "--reference",
    help="Path to the fasta file containing the reference sequence",
    required=True,
)
argparse.add_argument("-o", "--output", help="Path to the output file", required=True)

# Parse the command-line arguments
args = argparse.parse_args()

vcf = args.vcf_file
ref = args.reference
out = args.output


def load_data(infile):
    """
    Load data from infile if it is in fastq format (after having unzipped it, if it is zipped)

    Parameters:
        infile (str): Path to the input file.

    Returns:
        dict or bool: Dictionary containing sequence information if successful, False otherwise.
    """
    try:
        print("Reading data from input file...", file=sys.stderr)
        if infile.endswith(".gz"):  # If file is gzipped, unzip it
            y = gzip.open(infile, "rt", encoding="latin-1")
            if (
                infile.endswith(".fasta.gz")
                or infile.endswith(".fa.gz")
                or infile.endswith(".fas.gz")
                or infile.endswith(".fna.gz")
            ):
                # Read file as fasta if it is fasta
                records = SeqIO.parse(y, "fasta")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as value
                    seq_dict.update({record.id: str(record.seq)})
                y.close()
                return seq_dict
        # Read file directly as fasta if it is not zipped fastq
        elif (
            infile.endswith(".fasta")
            or infile.endswith(".fa")
            or infile.endswith(".fas")
            or infile.endswith(".fna")
        ):
            with open(infile, "r") as y:
                records = SeqIO.parse(y, "fasta")
                seq_dict = {}  # Create a dictionary to store everything from the file
                for record in records:
                    # Update dictionary with header as key, sequence as value
                    seq_dict.update({record.id: str(record.seq)})
                y.close()
                return seq_dict
        else:
            raise ValueError("File is the wrong format")
        print("Done", file=sys.stderr)
    except FileNotFoundError:
        print(
            f"File {infile} does not exist, proceeding with the analysis...",
            file=sys.stderr,
        )
        return False


def read_vcf(vcf_file: str):
    """
    Read data from the VCF file and extract relevant information.

    Parameters:
        vcf_file (str): Path to the VCF file.

    Returns:
        tuple or bool: Tuple containing variant dictionary and sample name if successful, False otherwise.
    """
    try:
        # Check if the VCF file is not gzipped
        if not vcf_file.endswith(".gz"):
            # Open the VCF file for reading and writing
            v = open(vcf_file, "r+")
            # Read all lines from the VCF file
            lines = v.readlines()
            # Separate lines with and without '#' (header lines)
            without_ash = [line for line in lines if not line.startswith("#")]
            with_ash = [line for line in lines if line.startswith("#")]
            # Extract the last header line and remove '#'
            headerline = with_ash[len(with_ash) - 1].replace("#", "")
            # Insert the modified header line at the beginning of the file
            without_ash.insert(0, headerline)
            # Move the file pointer to the beginning and truncate the file
            v.seek(0)
            v.truncate()
            # Write the modified lines back to the file
            for line in without_ash:
                v.write(line)
            # Close the file
            v.close()
            # Read the VCF data using pandas
            vcf = pd.read_csv(vcf_file, delimiter="\t")
            # Extract relevant information from the VCF data
            name = list(set(list(vcf["CHROM"])))[0]
            pos = list(vcf["POS"])
            refs = list(vcf["REF"])
            alts = list(vcf["ALT"])
            if len(pos) == 1:
                phases = [
                    (
                        sample.split(":")[0].split("/")[0],
                        sample.split(":")[0].split("/")[1],
                    )
                    for sample in list(vcf["SAMPLE"])
                ]
            else:
                phases = [
                    (
                        sample.split(":")[0].replace("/", "|").split("|")[0],
                        sample.split(":")[0].replace("/", "|").split("|")[1],
                    )
                    for sample in list(vcf["SAMPLE"])
                ]
            # Create a dictionary containing variant information
            vars_dict = {
                pos[i]: [(refs[i], phases[i][0]), (alts[i], phases[i][1])]
                for i in range(len(pos))
            }
            # Return the variant dictionary and sample name
            return vars_dict, vcf_file.split("/")[-4] + "_" + name
        else:
            # If the VCF file is gzipped, open it with gzip
            y = gzip.open(vcf_file, "rt", encoding="latin-1")
            # Read all lines from the gzipped VCF file
            lines = y.readlines()
            # Create a temporary VCF file without gzipping
            temporary_vcf = os.path.join(
                os.path.dirname(vcf_file),
                os.path.basename(vcf_file.replace(".vcf.gz", ".tmp.vcf")),
            )
            tmpvcf = open(temporary_vcf, "w")
            # Write the lines to the temporary VCF file
            for line in lines:
                tmpvcf.write(line)
            # Close the temporary VCF file
            tmpvcf.close()
            # Recursively call read_vcf on the temporary VCF file
            vars_dict, name = read_vcf(temporary_vcf)
            # Remove the temporary VCF file
            os.remove(temporary_vcf)
            # Return the variant dictionary and sample name
            return vars_dict, name
    except FileNotFoundError:
        # If the file is not found, return False
        return False, False


def is_a_phase(vars_dict: dict):
    """
    Identify phases in the variant dictionary.

    Parameters:
        vars_dict (dict): Dictionary containing variant information.

    Returns:
        list: List of unique phases.
    """
    phases = []
    # Iterate over each key in the variant dictionary
    for i in list(vars_dict.keys()):
        # Iterate over each item in the value (list of tuples) for the key
        for j in vars_dict[i]:
            # Iterate over each element in the tuple
            for k in j:
                try:
                    # Try to convert the element to an integer
                    h = int(k)
                    # If successful, consider it as a phase and add it to the phases list
                    phases.append(k)
                except Exception as e:
                    # If conversion fails, continue to the next element
                    continue
    # Return the list of unique phases
    return list(set(phases))


def render_vcf(inf, vcf, output):
    """
    Render VCF data onto a GenomeViz plot.

    Parameters:
        inf (str): Path to the reference sequence file.
        vcf (str): Path to the VCF file.
        output (str): Path to the output file.
    """
    # Read variant data from the VCF file
    var_dicts, name = read_vcf(vcf)

    # Check if variant data exists
    if var_dicts is not False:
        # Load sequence data from the reference file
        seqdict = load_data(inf)

        # Check if sequence data exists
        if seqdict is not False:
            # Determine the length of the sequence track
            track_length = len(list(seqdict.values())[0])

            # Create a GenomeViz object
            gv = GenomeViz(tick_style="axis")

            # Identify phases from the variant data
            phases = is_a_phase(var_dicts)

            # Iterate over phases to create feature tracks
            for phase in phases:
                track = gv.add_feature_track(name=name + "_" + phase, size=track_length)

                # Iterate over variant positions to add features
                for key in list(var_dicts.keys()):
                    label = []
                    # Iterate over variants at the current position
                    for j in var_dicts[key]:
                        # Check if the variant belongs to the current phase
                        if j[1] == phase:
                            label.append(j[0])
                    # Add feature to the track for the current position
                    track.add_feature(
                        start=int(key),
                        end=int(key) + 1,
                        strand=1,
                        plotstyle="box",
                        label=str(key) + ": " + ",".join(label),
                        labelrotation=0,
                        labelha="center",
                        facecolor="blue",
                        edgecolor="black",
                        linewidth=1,
                    )

            # Generate the GenomeViz plot
            fig = gv.plotfig()

            # Save the plot to the output file
            fig.savefig(output)
        else:
            # Print a message if the reference file does not exist
            print(
                f"Reference file {inf} does not exist, proceeding with the analysis..."
            )
    else:
        # Print a message if the VCF file does not exist
        print(f"Reference file {vcf} does not exist, proceeding with the analysis...")


if __name__ == "__main__":
    # Execute the main rendering function
    render_vcf(ref, vcf, out)
