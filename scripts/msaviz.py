from pymsaviz import MsaViz
from argparse import ArgumentParser

# Create an ArgumentParser to handle command-line arguments
argparse = ArgumentParser()
argparse.add_argument(
    "-m",
    "--msa_file",
    help="Path to the MSA file with aligned sequences in fasta format",
    required=True,
)

argparse.add_argument(
    "-o", "--output", help="Path to the output file", required=True)

# Parse the command-line arguments
args = argparse.parse_args()

msa = args.msa_file
out = args.output


def generate_msa_visualization(msa_file, output):
    mf = open(msa_file)
    lines = mf.readlines()
    mf.close()
    newmsa_file = msa_file + "_capitalized"
    nms = open(newmsa_file, "w")
    for line in lines:
        nms.write(line.upper())
    nms.close()
    mv = MsaViz(
        newmsa_file,
        wrap_length=200,
        show_count=True,
        show_consensus=True,
        consensus_color="#348ABD",
    )
    mv.savefig(output)


if __name__ == "__main__":
    generate_msa_visualization(msa, out)
