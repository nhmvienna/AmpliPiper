import subprocess
import os
from argparse import ArgumentParser
import sys

# Create an ArgumentParser to handle command-line arguments
argparse = ArgumentParser()
argparse.add_argument(
    "-hap",
    "--haplotypes_folder",
    help="Path to the 'haplotypes' folder within the results folder",
    required=True,
)

argparse.add_argument(
    "-py",
    "--python_script",
    help="Full command for the python script to convert VCF to PNG",
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

resultsdir = args.haplotypes_folder
pyscript = args.python_script
html = args.html_file
css = args.stylesheet_file


def list_files(directory):
    filenames = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            filenames.append(os.path.join(root, file))
    return filenames


def generate_images_recursively(resultsfolder, script):
    """
    Generate images recursively for variant-called genes using a specified script.

    Parameters:
    - resultsfolder (str): The base folder containing results.
    - script (str): The path to the Python script for image generation.

    Returns:
    - generated_imgs (list): A list of paths to the generated images.
    """

    # List directory and retrieve MSA files
    filenames = list_files(resultsfolder)
    samples = [f for f in filenames if f.endswith("_aln.fasta")]
    # Initialize a list to store paths of generated images
    generated_imgs = []

    # Iterate over each sample
    for sample in samples:
        sampleout = sample + "_visualized.png"
        subprocess.run(
            f"{script} -m {sample} -o {sampleout}",
            shell=True,
        )
        # Append the generated image path to the list
        generated_imgs.append(sampleout)

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
    with open(htmlpath, 'w') as html_file:
        html_file.write(htmlheader)
        html_file.write(htmlheader1)
        # Write links for each generated image
        for img in generated_imgs:
            html_file.write(
                f"<a href=\"#\" onclick=\"displayImage('{img}')\" >{os.path.basename(img).split('.')[0]}</a>\n"
            )
        # Write the HTML closer
        html_file.write(closer)


# Print a completion message
print("Everything's done!", file=sys.stderr)


if __name__ == "__main__":
    generated_images = generate_images_recursively(resultsdir, pyscript)
    display_images_recursively(generated_images, html, css)
