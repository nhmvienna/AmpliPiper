from argparse import ArgumentParser
import numpy as np
import pandas as pd
from Bio.Seq import Seq
import edlib
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for rendering without a display


# Argument parsing
argparse = ArgumentParser()
argparse.add_argument(
    "-p", "--primers", help="Path to the csv file containing primers", required=True)
argparse.add_argument(
    "-o", "--output", help="Path to the output csv file", required=True)
args = argparse.parse_args()


def reverse_complement(seq: str):
    """Return the reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def read_primer_table(csv_file):
    """Read primers from a CSV file."""
    primers = pd.read_csv(csv_file)
    return {row["ID"]: [row["FWD"], row["REV"]] for _, row in primers.iterrows()}


def calculate_edit_distance(seq1, seq2):
    """Calculate normalized edit distance between two sequences."""
    aln = edlib.align(seq1, seq2, task="path", mode="NW")
    return aln["editDistance"], f'{seq1}\n{seq2}'


def create_distance_matrix(primers):
    """Create a pairwise distance matrix for all primers."""
    primer_ids = list(primers.keys())
    all_primers = {f"{id}_fwd": primers[id][0] for id in primer_ids}
    all_primers.update({f"{id}_rev": primers[id][1] for id in primer_ids})
    all_primers.update({f"rc_{id}": reverse_complement(seq)
                       for id, seq in all_primers.items()})

    primer_names = list(all_primers.keys())
    matrix_size = len(primer_names)
    matrix = np.zeros((matrix_size, matrix_size))
    lens = [len(seq) for seq in all_primers.values()]

    for i in range(matrix_size):
        for j in range(matrix_size):
            if i != j:
                dist, _ = calculate_edit_distance(
                    all_primers[primer_names[i]], all_primers[primer_names[j]])
                matrix[i, j] = dist / max(lens)

    return matrix, primer_names


def plot_heatmap(matrix, labels, output_file):
    """Plot and save a heatmap from a distance matrix (lower triangle along the other diagonal with flipped x-axis)."""
    # Flip the matrix to orient the lower triangle along the other diagonal
    flipped_matrix = np.fliplr(np.flipud(matrix))  # Flip rows and columns

    # Retain only the lower triangle
    mask_upper_triangle = np.triu_indices_from(flipped_matrix, k=1)
    flipped_matrix[mask_upper_triangle] = np.nan  # Set upper triangle to NaN

    # Custom color palette
    custom_colors = ["#7ac1dc", "#ff8e3e"]
    cmap = LinearSegmentedColormap.from_list("custom_palette", custom_colors)

    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(flipped_matrix, cmap=cmap, aspect="auto")

    # Add labels and flip the x-axis
    ax.set_xticks(np.arange(len(labels)), labels=labels)
    ax.set_yticks(np.arange(len(labels)), labels=labels)  # Keep y-axis labels
    ax.invert_xaxis()  # Flip the x-axis
    # Rotate x-axis labels by 90 degrees
    plt.setp(ax.get_xticklabels(), rotation=90)
    cbar = ax.figure.colorbar(im, ax=ax, label="Normalized Edit Distance")

    plt.tight_layout()
    plt.savefig(output_file, dpi=200)

    return flipped_matrix


# Main script
if __name__ == "__main__":
    primers = read_primer_table(args.primers)
    matrix, labels = create_distance_matrix(primers)

    # Save the distance matrix as CSV
    pd.DataFrame(matrix, index=labels, columns=labels).to_csv(args.output)

    # Save the heatmap as PNG
    flipped_matrix = plot_heatmap(
        matrix, labels, args.output.replace(".csv", ".png"))

    # Calculate and print the minimum distance and corresponding sequences
    lower_triangle = flipped_matrix[np.tril_indices_from(flipped_matrix, k=-1)]
    min_distance = np.nanmin(lower_triangle)  # Ignore NaN values
    # Get the indices of the min distance
    min_index = np.where(flipped_matrix == min_distance)

    # Extract sequence names
    seq1 = labels[min_index[0][0]]
    seq2 = labels[min_index[1][0]]

    print(f"The minimum normalized edit distance is: {min_distance:.4f}")
    print(f"The sequences with the minimum distance are: {seq1} and {seq2}")
