from ete3 import Tree
import subprocess as sp
from argparse import ArgumentParser
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# Create an ArgumentParser to handle command-line arguments
argparse = ArgumentParser()
argparse.add_argument(
    "-t",
    "--trees",
    help="Path to the folder containing tree files",
    required=True,
)

argparse.add_argument(
    "-a",
    "--astral",
    help="Path to the ASTRAL tree file",
    required=True,
)

argparse.add_argument(
    "-od",
    "--outputdir",
    help="Path of the output directory (without further instructions, program will output a csv named treecompare.csv and a png named treecompare.png in the outputdir)",
    required=True,
)

argparse.add_argument(
    "-on",
    "--outputname",
    help="Customize the name of the output file in the output directory (default: treecompare)",
    required=False,
    default="treecompare",
)

# Parse the command-line arguments
args = argparse.parse_args()

treef = args.trees
astral = args.astral
outd = args.outputdir
outn = args.outputname


def heatmap(
    data, row_labels, col_labels, ax=None, cbar_kw=None, cbarlabel="", **kwargs
):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="black", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=("black", "white"),
    threshold=None,
    **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) < threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def find_all_trees(treefold: str, astraltree: str):
    """
    This function takes a directory path (`treefold`) and searches for all files with the extension ".treefile"
    in the directory and its subdirectories. It returns a list of the full paths of all found files.

    Args:
        treefold (str): The directory path to search for ".treefile" files.
        astraltree (str): ASTRAL tree file path.

    Returns:
        list: A list containing the full paths of all ".treefile" files found in the specified directory (and its subdirectories) and the path to the ASTRAL tree
    """
    filenames = []
    # Enumerate all the files in the directory with the walk function
    for root, dirs, files in os.walk(treefold):
        for file in files:
            # Add files if they contain a phylogenetic tree
            if file.endswith(".treefile"):
                filenames.append(os.path.join(root, file))
    filenames.append(astraltree)
    return filenames


def calculate_all_distances(treef: str, astraltree: str, outdir: str, outname: str):
    """
    Calculate Robinson-Foulds distances among gene trees and save results to a CSV file and a heatmap.

    Args:
        treef (list): List of treefile paths.
        astraltree (str): ASTRAL tree file path.
        outdir (str): Path to the output directory.
        outname (str): Name of the output files.
    Returns:
        None
    """
    # Find all treefiles in the specified directory
    trees = find_all_trees(treef, astraltree)

    # Create a temporary distance file for each comparison
    output = os.path.join(outdir, outname + ".csv")

    # Initialize the output CSV file
    out = open(output, "w")
    out.write("TREE1,TREE2,DIST\n")
    out.close()

    # Initialize lists to store distances, labels, and compared couples
    compared_couples = []

    # Initialize an empty distance matrix
    distance_matrix = np.zeros((len(trees), len(trees)))
    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix[i])):
            distance_matrix[i][j] = np.nan
    # Loop through all pairs of treefiles
    for i in range(len(trees)):
        for j in range(len(trees)):
            if i != j:
                # Check if the pair has not been compared before
                if [trees[i], trees[j]] not in compared_couples:
                    compared_couples.append([trees[i], trees[j]])
                    compared_couples.append([trees[j], trees[i]])

                    # Run Robinson-Foulds command and retrieve the distance
                    tree1 = Tree(trees[i])
                    tree2 = Tree(trees[j])
                    (
                        rf,
                        max_rf,
                        common_leaves,
                        parts_t1,
                        parts_t2,
                        discard_t1,
                        discard_t2,
                    ) = tree1.robinson_foulds(tree2, unrooted_trees=True)
                    dist = rf

                    # Append the distance to the distance matrix
                    distance_matrix[i, j] = dist

                    # Append the distance and labels to the output CSV file
                    out = open(output, "a")
                    out.write(
                        f"{trees[i].split('/')[-2]},{trees[j].split('/')[-2]},{dist}\n"
                    )

                else:
                    continue
            else:
                continue
    fig, ax = plt.subplots()
    trees1 = [t.split("/")[-2] for t in trees]
    im, cbar = heatmap(
        distance_matrix,
        trees1,
        trees1,
        ax=ax,
        cmap="coolwarm",
        cbarlabel="Robinson-Foulds dist",
    )
    texts = annotate_heatmap(im, valfmt="{x:.1f}")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, outname + ".png"))


if __name__ == "__main__":
    calculate_all_distances(treef, astral, outd, outn)
