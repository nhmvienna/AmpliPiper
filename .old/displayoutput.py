import cairosvg
from argparse import ArgumentParser
import webbrowser
import csv
import os


argparse = ArgumentParser()
argparse.add_argument("-r", "--resultsfolder",
                      help="Path to the results directory", required=True)
argparse.add_argument("-out", "--output",
                      help="Path to the outputfolder.", required=True)
argparse.add_argument("-loci", nargs='*', default=[],
                      help='Loci variables to be appended to the array')
argparse.add_argument("-samples", nargs='*', default=[],
                      help='Samples variables to be appended to the array')


args = argparse.parse_args()
resultsdir = args.resultsfolder
outdir = args.output
loci = args.loci
samples = args.samples
loci.append('Astral')


def get_svg_files_by_directory(input_path):
    for root, dirs, files in os.walk(input_path):
        svg_files = [file for file in files if file.lower().endswith('.svg')]
        if svg_files:
            for svg_file in svg_files:
                svg_path = os.path.join(root, svg_file)
                png_path = os.path.join(
                    root, svg_file.rsplit('.', 1)[0] + ".png")
                if os.stat(svg_path).st_size != 0:
                    print("Converting SVG files to PNG files for results display.")
                    cairosvg.svg2png(url=svg_path, write_to=png_path)


# def get_analysis_folders(base_path):
#    analysis_folders = [folder for folder in os.listdir(
#        base_path) if os.path.isdir(os.path.join(base_path, folder))]
#    return analysis_folders


def get_png_files(base_path):
    png_paths = []
    for root, _, files in os.walk(base_path):
        for file in files:
            if file.lower().endswith('.png'):
                png_paths.append(os.path.join(root, file))
    # print(pdf_paths)
    return png_paths


def list_directories(path):
    directories = []
    for entry in os.listdir(path):
        if entry != 'html':
            full_path = os.path.join(path, entry)
            if os.path.isdir(full_path):
                directories.append(entry)
    return directories


directories_in_dir = list_directories(resultsdir)

headnb = """<!DOCTYPE html>
<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<style>
body {
  font-family: "Lato", sans-serif;
}

/* Fixed sidenav, full height */
.sidenav {
  height: 100%;
  width: 500px;
  position: fixed;
  z-index: 1;
  top: 0;
  left: 0;
  background-color: #111;
  overflow-x: hidden;
  padding-top: 20px;
}

/* Style the sidenav links and the dropdown button */
.sidenav a, .dropdown-btn {
  padding: 6px 8px 6px 16px;
  text-decoration: none;
  font-size: 20px;
  color: #818181;
  display: block;
  border: none;
  background: none;
  width: 100%;
  text-align: left;
  cursor: pointer;
  outline: none;
}

/* On mouse-over */
.sidenav a:hover, .dropdown-btn:hover {
  color: #f1f1f1;
}

/* Main content */
.main {
  margin-left: 520px; /* Same as the width of the sidenav */
  font-size: 20px; /* Increased text to enable scrolling */
  padding: 0px 10px;
  max-width: 100%; /* Limit the maximum width of the container */
  min-height: 200px; /* Set a minimum height for the container */
  overflow: auto; /* Enable scrolling if content exceeds container size */
}

/* Images within the main content */
.main img {
  max-width: 100%; /* Make sure images don't exceed their container width */
  max-height: calc(100vh - 40px); /* Limit image height to fit within the viewport minus padding */
  min-height: 100px; /* Set a minimum height for the images */
  height: auto; /* Maintain aspect ratio */
}

/* Add an active class to the active dropdown button */
.active {
  background-color: green;
  color: white;
}

/* Dropdown container (hidden by default). Optional: add a lighter background color and some left padding to change the design of the dropdown content */
.dropdown-container {
  display: none;
  background-color: #262626;
  padding-left: 1px;
}

/* Optional: Style the caret down icon */
.fa-caret-down {
  float: right;
  padding-right: 8px;
}

/* Some media queries for responsiveness */
@media screen and (max-height: 450px) {
  .sidenav {padding-top: 15px;}
  .sidenav a {font-size: 18px;}
}
</style>
</head>
<body>
<div class="sidenav">

"""

htmlbutton = """<li>
                <button class="dropdown-btn">{} 
                    <i class="fa fa-caret-down"></i>
                </button>
                <div class="dropdown-container">
                    <ul>
                        {}
                    </ul>
                </div>
            </li>"""

onclicklink = """<a href="#" onclick="displayImage('{}', '{}')">{}</a>\n"""
showtablellink = """<li>
                            <button class="dropdown-btn" onclick="toggleTable()">{}</button>
                        </li>"""


mastercontent_closer = """
    </div>
<script>
function displayImage(imagePath, explanatoryText) {
    // Display the image and explanatory text
    var imageContainer = document.getElementById("imageContainer");
    var explanatoryTextContainer = document.getElementById("explanatoryText");
    imageContainer.src = imagePath;
    explanatoryTextContainer.textContent = explanatoryText;
    toggleTable2();
}

/* Loop through all dropdown buttons to toggle between hiding and showing its dropdown content - This allows the user to have multiple dropdowns without any conflict */
var dropdown = document.getElementsByClassName("dropdown-btn");
var i;

for (i = 0; i < dropdown.length; i++) {
  dropdown[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var dropdownContent = this.nextElementSibling;
    if (dropdownContent.style.display === "block") {
      dropdownContent.style.display = "none";
    } else {
      dropdownContent.style.display = "block";
    }
  });
}

function toggleTable() {
    var tableContainer = document.getElementById("tableContainer");
    var imageContainer = document.getElementById("imageContainer");
    var explanatoryText = document.getElementById("explanatoryText");

    if (tableContainer.style.display === "none") {
        // Show the table and hide other elements
        tableContainer.style.display = "block";
        imageContainer.style.display = "none";
        explanatoryText.style.display = "none";
    } else {
        // Hide the table and show other elements
        tableContainer.style.display = "none";
        imageContainer.style.display = "block";
        explanatoryText.style.display = "block";
    }
}

function toggleTable2() {
    var tableContainer = document.getElementById("tableContainer");
    var imageContainer = document.getElementById("imageContainer");
    var explanatoryText = document.getElementById("explanatoryText");
        // Hide the table and show other elements
        tableContainer.style.display = "none";
        imageContainer.style.display = "block";
        explanatoryText.style.display = "block";
}

</script>

</body>
</html> """


mastercontent_main = """
        </ul>
    </div>
    
    
    <div class="main">
        <p id="explanatoryText"></p>
        <img id="imageContainer" src="" alt="">
<div class="table-container" id="tableContainer" style="display: none;">"
<p id="explanatoryText">{}</p>
{}
</div>"""


def csv_to_html_table(csv_file):
    html = "<table border='1'>\n"
    with open(csv_file, 'r', newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        for row in csvreader:
            html += "<tr>\n"
            for cell in row:
                html += f"<td>{cell}</td>\n"
            html += "</tr>\n"
    html += "</table>"
    return html


def writehtml(dir, directories, mastercontent_main):
    mastercontent = ""
    master = "master_desaster.html"
    mastercontent += headnb
    for analysis in directories:
        if analysis != "SpeciesID":
            path = os.path.join(dir, analysis)
            links = ""  # Initialize links for each analysis
            for file in os.listdir(path):
                if file.lower().endswith('.html') or file.lower().endswith('.png'):
                    filepath = os.path.join(dir, analysis, file).replace(
                        (resultsdir+"/"), "")
                    # rel_path = os.path.join(analysis, file)
                    links += onclicklink.format(filepath, analysis, analysis)
            link_inner = ""
            link_inner_inner = ""
            # print(analysis)
            for elem in list_directories(path):
                if elem in loci:
                    # print(elem)
                    new_path = os.path.join(path, elem)
                    x = os.listdir(new_path)
                    numfile = sum(1 for file in x if file.endswith(
                        '.png') or file.endswith('.html'))
                    if numfile < 2:
                        for file in os.listdir(new_path):
                            if file.lower().endswith('.html') or file.lower().endswith('.png'):
                                filepath = os.path.join(new_path, file).replace(
                                    (resultsdir+"/"), "")
                                link_inner_inner += onclicklink.format(
                                    filepath, f'Results for {analysis} analysis at locus {elem}', elem)
                        link_inner = link_inner_inner
                    else:
                        link_n = ""
                        for file in os.listdir(new_path):
                            if file.lower().endswith('.html') or file.lower().endswith('.png'):
                                filepath = os.path.join(new_path, file).replace(
                                    (resultsdir+"/"), "")
                                link_n += onclicklink.format(
                                    filepath, f'Results for {analysis} analysis at locus: {elem}', file)
                        link_inner_inner = links + \
                            htmlbutton.format(elem, link_n)
                        link_inner += link_inner_inner
                elif elem in samples:
                    htmlbutton_o = ""
                    link_sample = ""
                    subdir = os.path.join(path, elem)
                    link_lo = ""
                    for subelem in list_directories(subdir):
                        locus = subelem
                        # print(locus)
                        sample = elem
                        if subelem in loci:
                            subelem_path = os.path.join(subdir, subelem)
                            link_locus = ""
                            for file in os.listdir(subelem_path):
                                if file.lower().endswith('consensus.fasta'):
                                    filepath = os.path.join(subelem_path, file)
                                    with open(filepath, 'r') as f:
                                        file_content = f.read().replace('\n', '\\n')
                                    link_locus = onclicklink.format(
                                        '', file_content, locus)
                                    link_sample += link_locus
                            buttonsample = htmlbutton.format(
                                sample, link_sample)
                    link_inner += buttonsample
                else:
                    for dirpath, dirs, files in os.walk(path):
                        for file in files:
                            if file.lower().endswith('.html') or file.lower().endswith('.png'):
                                filepath = os.path.join(
                                    dirpath, file).replace(resultsdir, "")
                                link_inner += f'<li><a href="{filepath}">{file}</a></li>'
                                mastercontent += "\n"
        else:
            path = os.path.join(dir, analysis)
            for dirpath, dirs, files in os.walk(path):
                for file in files:
                    if file.lower().endswith('final.csv'):
                        filepath = os.path.join(dirpath, file)
                        method = os.listdir(path)[0]
                        # print("FOUND csv.file", filepath)
                        print(dirpath, dirs, files)
                        x = csv_to_html_table(filepath)
                        csvlink = showtablellink.format(file)
                        link_inner = htmlbutton.format(method, csvlink)
                        mastercontent_main = mastercontent_main.format(
                            f"Summarized results for the Species Identification approach with {method}. ", x)
        link_inner = links+link_inner
        RES1 = htmlbutton.format(analysis, link_inner)
        mastercontent += RES1
        mastercontent += "\n"
    with open(os.path.join(outdir, "results.html"), 'w') as outmaster:
        outmaster.write(mastercontent)
        outmaster.write(mastercontent_main)
        outmaster.write(mastercontent_closer)


if __name__ == "__main__":
    get_svg_files_by_directory(resultsdir)
    writehtml(resultsdir, directories_in_dir, mastercontent_main)
    print("DONE!")
    print("***********HTML file generated successfully.*****************")
    webbrowser.open(f'{outdir}/results.html')
