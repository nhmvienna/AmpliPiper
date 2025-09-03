import os
import cairosvg
from argparse import ArgumentParser
import webbrowser

argparse = ArgumentParser()
argparse.add_argument("-r", "--resultsfolder",
                      help="Path to the results directory", required=True)
argparse.add_argument("-st", "--stylesheetfolder",
                      help="Path to the directory where stylesheet is stored.", required=True)
argparse.add_argument("-out", "--output",
                      help="Path to the outputfolder.", required=True)


args = argparse.parse_args()
resultsdir = args.resultsfolder
stylesheet = args.stylesheetfolder
outdir = args.output

# define functions for reading results directory


def get_svg_files_by_directory(input_path):
    for root, dirs, files in os.walk(input_path):
        svg_files = [file for file in files if file.lower().endswith('.svg')]
        if svg_files:
            for svg_file in svg_files:
                svg_path = os.path.join(root, svg_file)
                png_path = os.path.join(
                    root, svg_file.rsplit('.', 1)[0] + ".png")
                print("Converting SVG files to PNG files for results display.")
                cairosvg.svg2png(url=svg_path, write_to=png_path)


def get_analysis_folders(base_path):
    analysis_folders = [folder for folder in os.listdir(
        base_path) if os.path.isdir(os.path.join(base_path, folder))]
    return analysis_folders


def get_png_files(base_path):
    png_paths = []
    for root, _, files in os.walk(base_path):
        for file in files:
            if file.lower().endswith('.png'):
                png_paths.append(os.path.join(root, file))
    # print(pdf_paths)
    return png_paths


# functions for html
hideAllPages = """{
validPageIds.forEach(function (validId) {
        var page = document.getElementById(validId);
        if (page) {
          page.style.display = 'none';
        }
      });
    }"""

showPageId = """{
hideAllPages();
      // Show the selected page
      var targetPage = document.getElementById(pageId);
      if (targetPage) {
        targetPage.style.display = 'block';
      } else {
        // Display an error or fallback behavior if pageId is not found
        console.error("Invalid pageId:", pageId);
        // Optionally, show a default or home page
        showHomePage();
      }
    }"""

showHomePage = """{
hideAllPages();
}"""

showsubPage = """{
    event.preventDefault(); // Prevent the default behavior of the anchor tag
    hideAllPages();
    // Show the selected sub-page
    var targetPage = document.getElementById(pageId);
    if (targetPage) {
      targetPage.style.display = 'block';
    } else {
      // Display an error or fallback behavior if pageId is not found
      console.error("Invalid pageId:", pageId);
      // Optionally, show a default or home page
      showHomePage();
    }
  }"""

## main function to generate the htmls ###


def generate_html(analysis_folders, base_path, out):
    htmlHead = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Analysis for {analysis_folders}</title>
    <link href="static/vendor/bootstrap-main/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="static/css/style.css" rel="stylesheet">
    <script src="static/vendor/bootstrap-main/dist/js/bootstrap.bundle.min.js"></script>
</head>
<nav class="navbar bg-primary fixed-top">
  <div class="container-fluid">
    <a class="navbar-brand" href="#">
      <img src="static/assets/logo.svg" alt="Logo" width="30" height="24" class="d-inline-block align-text-top">
        <span class="text-white">Molester Analysis for {analysis_folders}</span>
    </a>
  </div>
</nav>
"""

    htmlBody = f"""
<nav class="navbar bg-primary fixed-top">
  <div class="container-fluid">
    <a class="navbar-brand" href="#">
      <img src="static/assets/logo.svg" alt="Logo" width="30" height="24" class="d-inline-block align-text-top">
        <span class="text-white">Molester Analysis for {analysis_folders}</span>
    </a>
  </div>
</nav>
<body>
<div class="row">
    <div class="col-4">
        <nav id="side-nav" class="h-100 flex-column align-items-stretch pe-4 border-end position-fixed">
            <nav class="nav nav-pills flex-column" role="tablist">
                <a class="nav-link" href="#item-1">Item 1</a>
                <nav class="nav nav-pills flex-column">
                    <a class="nav-link ms-3 my-1" href="#item-1-1">Item 1-1</a>
                    <a class="nav-link ms-3 my-1" href="#item-1-2">Item 1-2</a>
                </nav>
                <a class="nav-link" href="#item-2">Item 2</a>
                <a class="nav-link" href="#item-3">Item 3</a>
                <nav class="nav nav-pills flex-column">
                    <a class="nav-link ms-3 my-1" href="#item-3-1">Item 3-1</a>
                    <a class="nav-link ms-3 my-1" href="#item-3-2">Item 3-2</a>
                </nav>
            </nav>
        </nav>
    </div>
"""
    for analysis_folder in analysis_folders:
        html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{analysis_folder} Analysis</title>
  <link rel="stylesheet" href="{stylesheet}">
  <script>
    function displayImage(imagePath) {{
      document.getElementById('imageContainer').src = imagePath;
    }}
  </script>
</head>
<body>
  <aside>
    <nav>
    <a href="master.html"><p> Analysis Overview (Go Back).</p></a>
    <h1>{analysis_folder} analysis </h1>
"""
        searchpath = os.path.join(base_path, analysis_folder)
        directories = [d for d in os.listdir(
            searchpath) if os.path.isdir(os.path.join(searchpath, d))]
        for dir in directories:
            # print("looking for png files in",base_path, analysis_folder,dir)
            result = get_png_files(os.path.join(
                base_path, analysis_folder, dir))
            if result:
                if len(result) == 1:
                    for file in result:
                        html_content += f"""
<a href="#" onclick="displayImage('{file}')">Locus {dir}</a>
"""
                elif len(result) >= 1:
                    html_content += f"""
<div id="{dir}_{analysis_folder}">
<h2>{analysis_folder} - Locus {dir}</h2>
"""
                    for file in result:
                        root, extension = os.path.splitext(file)
                        parts = root.split('.')
                        last_extension = parts[-1]
                        html_content += f"""
<a href="#" onclick="displayImage('{file}')">{last_extension}</a>
"""
        html_content += """
</nav>
</aside>
<main>
  <img id="imageContainer" src="" alt="Selected Image">
</main>
</body>
</html>
"""
        with open(f'{out}/{analysis_folder}_analysis.html', 'w') as html_file:
            html_file.write(html_content)
            master_html_content += f"""
<li><a href="{analysis_folder}_analysis.html">{analysis_folder} Analysis</a></li>
"""
    master_html_content += """
</body>
</html>
"""
    with open(f'{out}/master.html', 'w') as master_html_file:
        master_html_file.write(master_html_content)


if __name__ == "__main__":
    get_svg_files_by_directory(resultsdir)
    analysis_folders = get_analysis_folders(resultsdir)
    generate_html(analysis_folders, resultsdir, outdir)
    print("DONE!")
    print("***********HTML file generated successfully.*****************")
    webbrowser.open(f'{outdir}/master.html')
