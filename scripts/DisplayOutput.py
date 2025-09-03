# new display output
import cairosvg
from argparse import ArgumentParser
import webbrowser
import csv
import os
import shutil


argparse = ArgumentParser()
argparse.add_argument("-p", "--parameters", default="",
                      help="Parameters used to run the pipeline", required=True)
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
# loci.append('astraltree')
parameters = args.parameters


def get_png_files(base_path):
    png_paths = []
    for root, _, files in os.walk(base_path):
        for file in files:
            if file.lower().endswith('.png'):
                png_paths.append(os.path.join(root, file))
    # print(pdf_paths)
    return png_paths


def get_svg_files_by_directory(input_path):
    for root, dirs, files in os.walk(input_path):
        svg_files = [
            file for file in files if file.lower().endswith('updated.svg')]
        if svg_files:
            for svg_file in svg_files:
                svg_path = os.path.join(root, svg_file)
                png_path = os.path.join(
                    root, svg_file.rsplit('.', 1)[0] + ".png")
                if os.stat(svg_path).st_size != 0:
                    print("Converting SVG files to PNG files for results display.")
                    cairosvg.svg2png(url=svg_path, write_to=png_path)


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


def list_directories(path):
    directories = []
    for entry in os.listdir(path):
        if entry != 'html':
            full_path = os.path.join(path, entry)
            if os.path.isdir(full_path):
                directories.append(entry)
    return directories


onclicklink = """<a href="#" onclick="displayImage('{}', '{}')"{}</a>\n"""
showlink = """<a href="#" onclick="showHaplotype('{}', '{}','consensus_{}')">{}</a>\n"""

showtablellink = """<li>
                            <button class="dropdown-btn" onclick="toggleTable()">{}</button>
                        </li>"""
showGridLink = """<a href="#" onclick="showGrid('{}')">{}</a>"""

showPicLinkButton = """<div><li><button class="dropdown-btn" onclick="displayImage('{}', '{}')" id="" > {} </button></li></div>"""


gridcontainer = """ <div class="grid-container" id="{}" style="display: none;">
<div class="click-zoom"><p style="padding-left: 1ch;">{}</p><label><input type="checkbox"><img src="{}"</label></div>"""

gridcontainer2 = """ <div class="grid-container" id="{}" style="display: none;">
<div class="click-zoom2"><p style="padding-left: 1ch;">{}</p><label><input type="checkbox"><img src="{}"</label>
</div>"""

consensuscontainer_id = """<div class="consensusContainer" id="{}" style="display: none;">
<button class="copy-btn" onclick="copyAllToClipboard()">Copy All</button>"""

# conensus_button="""<button class="copy-btn" onclick="copyToClipboard()">Copy {} </button>"""
conensus_button = """<button class="copy-btn" id="{}" value="{}" onclick="copyToClipboard('{}')" >COPY {} </button>"""


headnb = """<!DOCTYPE html>
<html>
<head>
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Nunito:wght@400;700&display=swap">
<style>

table {
          float: left; /* Aligns the table to the left */
          border-collapse: collapse; /* Ensures borders are collapsed */
          margin-right: 20px; /* Adds some space between the table and other content */
      }
      th, td {
          color: #333;
          border: 1px solid black; /* Adds a solid border around table cells */
          padding: 15px; /* Increases space between cell content and border */
      }

  body {
    font-family: "Lato", sans-serif;
  }
  
  /* Fixed sidenav, full height */
  .sidenav {
    height: 100%;
    width: 300px;
    position: fixed;
    z-index: 1;
    top: 0;
    left: 0;
    background-color: #C3E6DA;
    overflow-x: hidden;
    padding-top: 20px;
  }
  
  /* Style the sidenav links and the dropdown button */
  .sidenav a, .dropdown-btn {
      padding: 6px 8px 6px 16px;
      text-decoration: none;
      font-size: 20px;
      color: #21614A;
      display: block;
      border: none;
      background: none;
      width: 100%;
      text-align: left;
      cursor: pointer;
      outline: none;
    }

    /* Style for the second type of dropdown button */
    .dropdown-btn-secondary {
      padding: 6px 8px 6px 16px;
      text-decoration: none;
      font-size: 20px;
      color: #ffffff;
      display: block;
      background-color: #FF8E3E; /* Orange */
      border: none;
      width: 100%;
      text-align: left;
      cursor: pointer;
      outline: none;
    }

  
  /* On mouse-over */
  .sidenav a:hover, .dropdown-btn:hover {
    color: #FF8E3E;
  }
        
  body {
    font-family: "Lato", sans-serif;
  }
  
  /* Main content */
  .main {
    margin-left: 360px; /* Same as the width of the sidenav */
    font-size: 20px; /* Increased text to enable scrolling */
    /*padding: 0px 10px;
    max-width: 100%; /* Limit the maximum width of the container */
    overflow: scroll; /* Enable scrolling if content exceeds container size */
    height: 50%;
    padding-bottom: 50px;
  }
  
  /* Images within the main content */
  .main img {
    max-width: 100%; /* Make sure images don't exceed their container width */
    height: auto; /* Maintain aspect ratio */
    margin-bottom: 10px;
    overflow: scroll;
  }
  
  
  /* Fancy enhancements */
  .main h1 {
    font-size: 56px;
    color: #333;
    text-align: center;
    margin-bottom: 30px;
  }
  
  .main p {
    font-size: 24px;
    line-height: 1.6;
    color: #21614A;
    text-align: justify;
    margin-bottom: 20px;
  }

  .landingpagecontent {
    font-size: 24px;
    line-height: 1.6;
    color: #21614A;
    text-align: justify;
    margin-bottom: 20px;
  }
  
  
  /* Add an active class to the active dropdown button */
  .active {
    background-color: #C3E6DA;
    color: white;
  }
  
  /* Dropdown container (hidden by default). Optional: add a lighter background color and some left padding to change the design of the dropdown content */
  .dropdown-container {
    display: none;
    background-color: #C3E6DA;
    padding-left: 1px;
  }
  
  /* Optional: Style the caret down icon */
  .fa-caret-down {
    float: right;
    padding-right: 8px;
  }
  

  .logo {
    padding: 20px; /* Adjust padding as needed */
    text-align: left;
  }

  .logo img {
    width: 180px; /* Adjust image width */
    height: auto; /* Maintain aspect ratio */
  }


  /* This is the HOME/Amplipiper Button that reloads the page when clicked */

  .btn {
    background-color: #FF8E3E;
    border: none;
    color: white;
    padding: 12px 16px;
    font-size: 19px;
    cursor: pointer;
    font-family: 'Nunito extra bold', sans-serif;
  }
  
  /* Darker background on mouse-over */
  .btn:hover {
    background-color: #FF8E3E;
  }

  .small {
      font-size: 14px;
    }
    .medium {
      font-size: 18px;
    }
    .large {
      font-size: 40px;
    }
  
   .main-container {
    margin-top: 10px;
    width: 100%;
    height: auto;
    background: white;
    overflow: visible;
  }

  .grid-container {
    margin-left: 20px;
    width: 95%;
    height: auto;
    background: white;
    overflow: scroll;
  }

  .maindiv {
    width: 100%;
    height: 100%;
    background: white;
    transform-origin:0% 0%;
    margin-bottom: 10px;
  }
  
  .click-zoom {
  margin-left: 0;
  /* Specify the dimensions of the container */
  width: 90%; /* Adjust as needed */
  height: 90%; /* Adjust as needed */
  /* Ensure that the overflow content is scrollable */
  overflow: visible; /* or overflow: scroll; */
  /* Ensure the container stays on top of other elements */
  position: relative;
}


/* Hide the checkbox */
.click-zoom input[type=checkbox] {
  display: none;
  width: 100%;
}

/* Style for the image */
.click-zoom img {
  transform-origin: 0 0;
  transition: transform 0.25s ease;
  cursor: zoom-in;
  /* Ensure the image stays on top of other elements when zoomed */
  position: relative;
  z-index: 1;
  margin-bottom: 200px; 
}

/* Apply zoom effect when checkbox is checked */
.click-zoom input[type=checkbox]:checked ~ img {
  transform: scale(1.8); /* Adjust zoom level as needed */
  cursor: zoom-out;
  /* Increase z-index further to ensure the zoomed image stays on top */
  z-index: 2;
}
  
/* zoooom2 for trees*/
.click-zoom2 {
  width: 30%; /* Adjust as needed */
  height: auto; /* Adjust as needed */
  /* Ensure that the overflow content is scrollable */
  /* Ensure the container stays on top of other elements */
  position: relative;
}

.click-zoom2 input[type=checkbox] {
  display: none;
  width: 90%;
}

/* Style for the image */
.click-zoom2 img {
  transform-origin: 0 0;
  transition: transform 0.25s ease;
  cursor: zoom-in;
  /* Ensure the image stays on top of other elements when zoomed */
  position: relative;
  z-index: 1;
  width:100%;
  margin-bottom: 200px; 
}

/* Apply zoom effect when checkbox is checked */
.click-zoom2 input[type=checkbox]:checked ~ img {
  transform: scale(1.8); /* Adjust zoom level as needed */
  cursor: zoom-out;
  /* Increase z-index further to ensure the zoomed image stays on top */
  z-index: 2;
}
  
  </style>
</head>
<body>
<div class="sidenav">
<div class="logo">
    <img src="../.logos/tettris.png" alt="Logo" onclick="showExplanation()">
</div>
<li>
  <button class="dropdown-btn-secondary" onclick=""  >Summary
      <i class="fa fa-caret-down"></i>
  </button>
  <div class="dropdown-container">
      <ul>
          <a href="#" onclick="showPloidy()"> Table</a>
          <a href="#" onclick="showGrid('SummaryHaplotypes') "> Reconstructed Haplotypes </a>
          <a href="#" onclick="showGrid('SummaryPloidy') "> Modelled Ploidy </a>
          <a href="#" onclick="showGrid('SummaryPrimers') "> Primers </a>
      </ul>
  </div>
</li> 
"""

# dir="/media/inter/ssteindl/Amplicon/latest_HAPLOTYPES/HAPLOTYPES/lepidoptera_test_output/results"
directories_in_dir = ["summary", "tree",
                      "haplotypes", "SpeciesID", "SpeciesDelim"]

# for dirc in directories_in_dir:
#    full_path = os.path.join(dir, dirc)
#    #print(full_path)
#    get_png_files(full_path)

htmlbutton_empty = """<li>
                <button class="dropdown-btn" {} >{} 
                    <i class="fa fa-caret-down"></i>
                </button>
                <div class="dropdown-container">
                    <ul>
                        {}
                    </ul>
                </div>
            </li> """

# allparameters=f"""quality = {parameters.split(",")[0]},
# similarconsensus= {parameters.split(",")[1]}, nreads= {parameters.split(",")[2]}, sizerange= {parameters.split(",")[3]},
# minreads= {parameters.split(",")[4]}, threads= {parameters.split(",")[5]}, kthres= {parameters.split(",")[6]},
# force= {parameters.split(",")[7]}, blast= {parameters.split(",")[8]}, partition= {parameters.split(",")[9]},
# """


allparameters = f"""
<table border="1" style="text-align:left; padding: 15px;">
  <tr>
    <th>Parameters  </th>
    <th>Values  </th>
  </tr>
  <tr>
    <td>Quality threshold  </td>
    <td>{parameters.split(",")[0]}</td>
  </tr>
  <tr>
    <td>Similar consensus threshold  </td>
    <td>{parameters.split(",")[1]}</td>
  </tr>
  <tr>
    <td>Number of reads</td>
    <td>{parameters.split(",")[2]}</td>
  </tr>
  <tr>
    <td>Size range</td>
    <td>{parameters.split(",")[3]}</td>
  </tr>
  <tr>
    <td>Minimum reads required</td>
    <td>{parameters.split(",")[4]}</td>
  </tr>
  <tr>
    <td>Number of threads</td>
    <td>{parameters.split(",")[5]}</td>
  </tr>
  <tr>
    <td>K-threshold</td>
    <td>{parameters.split(",")[6]}</td>
  </tr>
  <tr>
    <td>Force Flag</td>
    <td>{parameters.split(",")[7]}</td>
  </tr>
  <tr>
    <td>BLAST usage</td>
    <td>{parameters.split(",")[8]}</td>
  </tr>
  <tr>
    <td>Partition strategy</td>
    <td>{parameters.split(",")[9]}</td>
  </tr>
  <tr>
    <td>Outgroup Definition</td>
    <td>{parameters.split(",")[10]}</td>
  </tr>
</table> 
"""

mastercontent_main = f"""
        </ul>
        <div class="logo">
          <img src="../.logos/nhm.svg.png" alt="Logo">
        </div>
    </div>
    
    <div class="main">
    <div class="landingpagecontent"
      <p id="LandingPage"> <span class="large"><br> Your analysis has finished.</span><br>
      <span class=""> Click on the tabs on the side to show the corresponding analysis results.</span><br>
      <br>
      <span class="large"><br> The following parameters were used for the analysis:</span><br></p>
      {allparameters}  </div>
        <p id="explanatoryText"> </p>
        <div class="main-container">
          <div  class="maindiv">
              <img id="imageContainer" src="" alt="">
          </div>
        </div>

<div class="table-container" id="tableContainer" style="display: none;">
<p id="explanatoryText">{{}}</p>
{{}}
</div>
"""

ploidyContainerEmpty = """
<div class="ploidy-container" id="ploidyContainer" style="display: none;">
<p id="PloidyContainer">{}</p>
{}
</div>
"""


mastercontent_closer = """
</div>

<script>
  function displayImage(imagePath, explanatoryText) {
      // Display the image and explanatory text
      hideLandingPageContent();
      var imageContainer = document.getElementById("imageContainer");
      var explanatoryTextContainer = document.getElementById("explanatoryText");
      imageContainer.src = imagePath;
      var mainContainer = document.querySelector(".main img");
      explanatoryTextContainer.textContent = explanatoryText;
      toggleTable2();
      hideLandingPageContent();
      hidePloidy();

      var gridContainers = document.querySelectorAll('.grid-container');
      // Hide all grid containers
      gridContainers.forEach(function(container) {
      container.style.display = 'none';
    });

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


  var dropdown2 = document.getElementsByClassName("dropdown-btn-secondary");
  var i;
  for (i = 0; i < dropdown2.length; i++) {
    dropdown2[i].addEventListener("click", function() {
      this.classList.toggle("active");
      var dropdownContent2 = this.nextElementSibling;
      if (dropdownContent2.style.display === "block") {
        dropdownContent2.style.display = "none";
      } else {
        dropdownContent2.style.display = "block";
      }
    });
  }


  function showTable(){
      var tableContainer = document.getElementById("tableContainer");
      tableContainer.style.display = "block";
  }
  


  function showTable(){
      var tableContainer = document.getElementById("tableContainer");
      tableContainer.style.display = "block";
  }
  
  function toggleTable() {
      var tableContainer = document.getElementById("tableContainer");
      var ploidyContainer = document.getElementById("ploidyContainer");
      var imageContainer = document.getElementById("imageContainer");
      var explanatoryText = document.getElementById("explanatoryText");
      var ZoomButton = document.getElementById("zoomheader");
      var gridContainers = document.querySelectorAll('.grid-container');
      var elements = document.querySelectorAll('.LandingPage');
        // Hide all grid containers
        gridContainers.forEach(function(container) {
        container.style.display = 'none'; 
        });
  
      if (tableContainer.style.display === "none") {
          // Show the table and hide other elements
          tableContainer.style.display = "block";
          imageContainer.style.display = "none";
          ploidyContainer.style.display = "none";
          explanatoryText.style.display = "none";
          LandingPage.style.display = "none";
          ZoomButton.style.display = "none";
      } else {
          // Hide the table and show other elements
          tableContainer.style.display = "none";
          ploidyContainer.style.display = "none";
          imageContainer.style.display = "block";
          explanatoryText.style.display = "block";
          LandingPage.style.display = "none";
          hideLandingPageContent();
      }
      
      gridContainers.forEach(function(container) {
        container.style.display = 'none';
        hideLandingPageContent();
        hidePloidy();
    });
  }


  function showPloidy() {
    var ploidyContainer = document.getElementById("ploidyContainer");
    ploidyContainer.style.display = "block";
    var imageContainer = document.getElementById("imageContainer");
    imageContainer.style.display = "none";
    var tableContainer = document.getElementById("tableContainer");
    tableContainer.style.display = "none";
    var gridContainers = document.getElementsByClassName("grid-container");
    for (var i = 0; i < gridContainers.length; i++) {
      gridContainers[i].style.display = "none";
    }
    var explanatoryText = document.getElementById("explanatoryText");
    explanatoryText.style.display = "none";
    var LandingPage = document.getElementById('LandingPage');
    LandingPage.style.display = "none";
  }

  function showExplanation() {
        var gridContainers = document.querySelectorAll('.grid-container');
        var imageContainer = document.getElementById("imageContainer");
        var tableContainer = document.getElementById("tableContainer");
        var explanatoryText = document.getElementById("explanatoryText");
        var elements = document.querySelectorAll('.LandingPage');
          tableContainer.style.display = "none";
          explanatoryText.style.display = "none";
          imageContainer.style.display = "none";
          LandingPage.style.display = "none";
        hideConsensus();
        hidePloidy();

        // Hide all grid containers
        gridContainers.forEach(function(container) {
            container.style.display = 'none';
        });
        var explanation = document.getElementById("LandingPage");
        if (explanation.style.display === "none" || explanation.style.display === "") {
            explanation.style.display = "block";
        } else {
            explanation.style.display = "none";
        }
    }

  function hideLandingPageContent() {
          var elements = document.querySelectorAll('.LandingPage');
          elements.forEach(function(element) {
              element.style.display = 'none';
          });
      }

  function hideConsensus() {
          var consbuttons = document.querySelectorAll('.consensusContainer');
          consbuttons.forEach(function(element) {
              element.style.display = 'none';
          });
      }

  function hidePloidy() {
    var ploidyContainer = document.getElementById("ploidyContainer");
    ploidyContainer.style.display = "none";
  }
  
  function toggleTable2() {
      var tableContainer = document.getElementById("tableContainer");
      var imageContainer = document.getElementById("imageContainer");
      var explanatoryText = document.getElementById("explanatoryText");
          // Hide the table and show other elements
          tableContainer.style.display = "none";
          imageContainer.style.display = "block";
          explanatoryText.style.display = "block";
          LandingPage.style.display = "none";
          hideLandingPageContent();
          hideConsensus();
          hidePloidy();
      var gridContainers = document.querySelectorAll('.grid-container');
        // Hide all grid containers
        gridContainers.forEach(function(container) {
        container.style.display = 'none'; 
        });
  }
  
  


function showGrid(Id) {
    var gridContainers = document.querySelectorAll('.grid-container');
    var imageContainer = document.getElementById("imageContainer");
    var tableContainer = document.getElementById("tableContainer");
    var explanatoryText = document.getElementById("explanatoryText");
    var elements = document.querySelectorAll('.LandingPage');
      tableContainer.style.display = "none";
      explanatoryText.style.display = "none";
      LandingPage.style.display = "none";
    hideConsensus();
    hidePloidy();
    
    // Hide all grid containers
    gridContainers.forEach(function(container) {
        container.style.display = 'none';
    });

    // Show the grid container with the specified ID
    var gridContainer = document.getElementById(Id);
    if (gridContainer) {
        gridContainer.style.display = 'block';
        imageContainer.style.display = "none";
    }

  }


  function showButton(elem) {
    var consContainer = document.getElementById(elem);
    if (consContainer.style.display === "none") {
    consContainer.style.display = 'block';
    hidePloidy();
  }  }


const checkbox = document.getElementById('checkbox');
let zoomLevel = 1; // Initial zoom level

checkbox.addEventListener('click', () => {
  if (zoomLevel === 1) {
    // First click, zoom in once
    zoomLevel = 1.8; // Adjust zoom level as needed
  } else {
    // Second click, zoom in further
    zoomLevel *= 2.8; // Zoom in by multiplying the current zoom level
  }
  // Apply the zoom level to the image
  document.querySelector('.click-zoom img').style.transform = `scale(${zoomLevel})`;
});

</script>
</body>
</html> """


def extrafiles(path, analysis):
    for file in os.listdir(path):
        if file.lower().endswith('comparison.html') or file.lower().endswith('comparison.png'):
            treelink = os.path.join(analysis, file)
            overview = """ onclick="displayImage('""" + treelink + \
                """', 'Distance of different trees displayed by heatmap.')" """
            return overview
    return None

# Function to filter out non-PNG files


def ignore_non_png_files(dir, files):
    return [file for file in files if not (file.endswith('.png') or file.endswith('.pdf') or file.endswith('summary.csv') or os.path.isdir(os.path.join(dir, file)))]


def writehtml(dir, directories, mastercontent_main):
    mastercontent = ""
    mastercontent += headnb
    RES3 = showPicLinkButton.format(
        "haplotypes/distance_matrices.png", "Genetic Distances", "Genetic Distances")
    mastercontent += RES3
    for analysis in directories:
        path = os.path.join(dir, analysis)
        if os.path.exists(path) and analysis != "SpeciesDelim":
            shutil.copytree(path, os.path.join(
                outdir, analysis), ignore=ignore_non_png_files)
        else:
            print("Directory", analysis, "not found.")
        link_inner = ""
        links = ""
        if analysis == "summary":
            sfile = os.path.join(dir, analysis, "summary.csv")
            ploidy = csv_to_html_table(sfile)
            ploidyContainer = ploidyContainerEmpty.format("", ploidy)
            mastercontent_main += ploidyContainer
            for analysisname, filename in {'SummaryHaplotypes': 'summary/summary.csv.png', 'SummaryPrimers': 'summary/primers/primers_dist.png', 'SummaryPloidy': 'summary/summary.csv_ExpectedPloidy.png'}.items():
                gridc = gridcontainer.format(analysisname, "", filename)
                mastercontent_main += gridc
                mastercontent_main += """</div>"""
            # gridc=gridcontainer.format(analyis_el, filepath)
        if analysis == "tree":
            # addinfo = extrafiles(path, analysis)
            htmlbutton = htmlbutton_empty.format(
                "", "Phylogenetic Trees", "{}")
            allinks = ""
            # print(htmlbutton)
            # link_inner = "" #old position of "link_inner"
            link_inner_inner = ""
            # astraltreepath = os.path.join(dir, "astraltree")
            for elem in list_directories(path):
                # print("ELEM IN LOCI! PATH IS:", path, elem)
                new_path = os.path.join(path, elem)
                x = os.listdir(new_path)
                analyis_el = analysis+elem
                for file in os.listdir(new_path):
                    if file.lower().endswith('.html') or file.lower().endswith('.png'):
                        filepath = os.path.join(
                            new_path, file).replace((dir+"/"), "")
                # link_inner_inner += onclicklink.format(filepath, f'Results for {analysis} analysis at locus {elem}', elem)
                # showgridTree here
                        # piclink = showGridLink.format(analyis_el,elem)
                        # allinks += piclink
                        gridc = gridcontainer2.format(
                            analyis_el, elem, filepath)
                mastercontent_main += gridc
                link_inner_inner += showGridLink.format(analyis_el, elem)
                mastercontent_main += """</div>"""
            # for file in os.listdir(astraltreepath):
            #     if file.lower().endswith('.html') or file.lower().endswith('.png'):
            #         filepath = os.path.join(
            #             astraltreepath, file).replace((dir+"/"), "")
            #         gridc = gridcontainer2.format(
            #             'astraltree', "ASTRAL", filepath)
            #         mastercontent_main += gridc
            #         mastercontent_main += """</div>"""
            #         link_inner_inner += showGridLink.format(
            #             'astraltree', 'ASTRAL Tree')
            #     else:
            #         continue
    # print(link_inner_inner)
            RES1 = htmlbutton.format(link_inner_inner)
            mastercontent += RES1
        if (analysis == "SpeciesDelim" and os.path.exists(path)) or analysis == "haplotypes":
            if analysis == "SpeciesDelim":
                get_svg_files_by_directory(path)
                shutil.copytree(path, os.path.join(
                    outdir, analysis), ignore=ignore_non_png_files)
                variable = "Species Delimitation"
            else:
                variable = "Consensus Alignment"
            allinks = ""
            htmlbutton = htmlbutton_empty.format("", variable, "{}")
            for elem in list_directories(path):
                new_path = os.path.join(path, elem)
                analyis_el = analysis+elem
                # piclink = showGridLink.format(elem,elem)
                # gridc = gridcontainer.format(elem)
                # print(elem)
                # allinks += piclink
                for file in os.listdir(new_path):
                    if file.lower().endswith('updated.png') or file.lower().endswith('visualized.png'):
                        piclink = showGridLink.format(analyis_el, elem)
                        allinks += piclink
                        filepath = os.path.join(
                            new_path, file).replace((dir+"/"), "")
                        # print(filepath)
                        # insert heading here
                        gridc = gridcontainer.format(
                            analyis_el, elem, filepath)
                        # gridc = gridcontainer.format(analyis_el, filepath)
                    else:
                        continue
                mastercontent_main += gridc
                # mastercontent_main += """ <p style="padding-left: 1ch;">"""+elem+"""</p></div>"""
                mastercontent_main += """</div>"""
            RES2 = htmlbutton.format(allinks)
            mastercontent += RES2
            # link_n += onclicklink.format(elem)
            # mastercontent_main += gridc
        if analysis == "SpeciesID" and os.path.exists(path):
            variable = "Species Identification"
            for dirpath, dirs, files in os.walk(path):
                for file in files:
                    if file.lower().endswith('final.csv'):
                        filepath = os.path.join(dirpath, file)
                        method = os.listdir(path)[0]
                        # print("FOUND csv.file", filepath)
                        x = csv_to_html_table(filepath)
                        csvlink = showtablellink.format(variable)
                        # link_inner = htmlbutton_empty.format(analysis, method, csvlink)
                        mastercontent_main = mastercontent_main.format(
                            f"Summarized results for the Species Identification approach with {method}. ", x)
            # link_inner = links+link_inner
            # RES3 = htmlbutton_empty.format("", analysis, csvlink)
            mastercontent += csvlink
        if analysis == "SpeciesID" and not os.path.exists(path):
            print("SpeciesID Folder does not exist")
            # RES1 = htmlbutton.format(analysis, link_inner)
            # mastercontent += RES1
    # pairwise_distances=os.path.join(dir, "haplotypes/distance_matrices.png")
    # insert=onclicklink.format(pairwise_distances, "", "")
    # RES3 = htmlbutton_empty.format(insert,"GeneticDistances", "")
    # htmlbutton = htmlbutton_empty.format(addinfo, analysis, "{}")
    # RES3 = showPicLinkButton.format("haplotypes/distance_matrices.png","Genetic Distances", "Genetic Distances")
    # mastercontent += RES3
    with open(os.path.join(outdir, "results.html"), 'w') as outmaster:
        outmaster.write(mastercontent)
        outmaster.write(mastercontent_main)
        outmaster.write(mastercontent_closer)


if __name__ == "__main__":
    writehtml(resultsdir, directories_in_dir, mastercontent_main)
    print("DONE!")
    print("***********HTML file generated successfully.*****************")
    webbrowser.open(f'{outdir}/results.html')
