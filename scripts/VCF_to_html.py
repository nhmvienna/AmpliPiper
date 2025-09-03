import subprocess
import os
import gzip
from argparse import ArgumentParser
import sys

# Create an ArgumentParser to handle command-line arguments
argparse = ArgumentParser()
argparse.add_argument(
    "-rf",
    "--results_folder",
    help="Path to the folder where the results are stored (the same passed with --output argument in the main pipeline)",
    required=True,
)

argparse.add_argument(
    "-py",
    "--python_script",
    help="Path to the python script to convert VCF to PNG",
    required=True,
)

argparse.add_argument(
    "-sf",
    "--stylesheet_file",
    help="Path to the CSS stylesheet file",
    required=True,
)

argparse.add_argument(
    "-hf",
    "--html_file",
    help="Path to the output html file",
    required=True,
)


# Parse the command-line arguments
args = argparse.parse_args()

resultsdir = args.results_folder
pyscript = args.python_script
html = args.html_file
css = args.stylesheet_file


def generate_images_recursively(resultsfolder, script):
    """
    Generate images recursively for variant-called genes using a specified script.

    Parameters:
    - resultsfolder (str): The base folder containing results.
    - script (str): The path to the Python script for image generation.

    Returns:
    - generated_imgs (list): A list of paths to the generated images.
    """

    # Construct the base path for consensus sequences
    basepath = os.path.join(resultsfolder, "results/consensus_seqs")

    # Get a list of sample directories
    samples = [
        os.path.join(basepath, sample)
        for sample in os.listdir(basepath)
        if os.path.isdir(os.path.join(basepath, sample))
    ]

    # Initialize a list to store paths of generated images
    generated_imgs = []

    # Iterate over each sample
    for sample in samples:
        # Get paths for variant-called genes
        variantcalled_genes = [
            os.path.join(sample, gene + f"/{gene}/variant_calls.snps.phased.vcf.gz")
            for gene in os.listdir(sample)
            if (
                os.path.isdir(os.path.join(sample, gene))
                and os.path.exists(
                    os.path.join(
                        sample, gene + f"/{gene}/variant_calls.snps.phased.vcf.gz"
                    )
                )
            )
        ]

        # Get paths for reference genes
        reference_genes = [
            os.path.join(sample, gene + f"/selectedconsensus.fasta")
            for gene in os.listdir(sample)
            if (
                os.path.isdir(os.path.join(sample, gene))
                and os.path.exists(
                    os.path.join(sample, gene + f"/selectedconsensus.fasta")
                )
                and os.path.exists(
                    os.path.join(
                        sample, gene + f"/{gene}/variant_calls.snps.phased.vcf.gz"
                    )
                )
            )
        ]

        # Iterate over each gene
        if len(variantcalled_genes) > 0 and len(variantcalled_genes) == len(
            reference_genes
        ):
            print("Equal!")
        else:
            print("Different!")
        for i in range(len(variantcalled_genes)):
            # Read variant-called genes file
            f = gzip.open(variantcalled_genes[i], mode="rt", encoding="latin-1")
            lines = [line for line in f.readlines() if not line.startswith("#")]
            f.close()

            # Check if there are variants in the file
            if len(lines) > 0:
                # Run the image generation script using subprocess
                subprocess.run(
                    f"python3 {script} -v {variantcalled_genes[i]} -r {reference_genes[i]} -o {os.path.join(os.path.dirname(reference_genes[i]), sample+'_'+os.path.basename(os.path.dirname(reference_genes[i]))+'.png')}",
                    shell=True,
                )

                # Print the command for reference
                print(
                    f"python3 {script} -v {variantcalled_genes[i]} -r {reference_genes[i]} -o {os.path.join(os.path.dirname(reference_genes[i]), sample+'_'+os.path.basename(os.path.dirname(reference_genes[i]))+'.png')}"
                )

                # Append the generated image path to the list
                generated_imgs.append(
                    os.path.join(
                        os.path.dirname(reference_genes[i]),
                        sample
                        + "_"
                        + os.path.basename(os.path.dirname(reference_genes[i]))
                        + ".png",
                    )
                )
            else:
                continue

        # Print a message indicating completion for the current sample
        print(f"Finished generating images for {sample}", file=sys.stderr)

    # Print a message indicating completion for all samples
    print("Finished generating all images", file=sys.stderr)

    # Return the list of generated image paths
    return generated_imgs


def display_images_recursively(generated_imgs, htmlpath, stylesheetpath) -> None:
    """
    Create an HTML file to display generated images with links for each image.

    Parameters:
    - generated_imgs (list): A list of paths to the generated images.
    - htmlpath (str): The path to save the HTML file.
    - stylesheetpath (str): The path to the CSS stylesheet.

    Returns:
    - None
    """

    # HTML header with metadata and stylesheet link
    htmlheader = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="{stylesheetpath}">
    <title>Variant analysis visualization</title>
    """

    # JavaScript script for displaying images
    htmlheader1 = """\n<script>
		function displayImage(imagePath) {
      		document.getElementById('imageContainer').src = imagePath;
    	}
  	</script>
</head>
<body>
<aside>
    <nav>
    <a href="master.html"><p> Analysis Overview (Go Back).</p></a>
    <center><h1>Variant analysis visualization</h1></center>
    <br>
    <br>
    """

    # HTML closer with image container
    closer = """</nav>
</aside>
<main>
  <img id="imageContainer" src="" alt="Selected Image">
</main>
</body>
"""

    # Write the HTML file
    html = open(htmlpath, "w")
    html.write(htmlheader)
    html.write(htmlheader1)

    # Write links for each generated image
    for img in generated_imgs:
        html.write(
            f"<a href=\"#\" onclick=\"displayImage('{img}')\" >{os.path.basename(img).split('.')[0]}</a>\n"
        )

    # Write the HTML closer
    html.write(closer)
    html.close()

    # Print a completion message
    print("Everything's done!", file=sys.stderr)


if __name__ == "__main__":
    generated_images = generate_images_recursively(resultsdir, pyscript)
    display_images_recursively(generated_images, html, css)
