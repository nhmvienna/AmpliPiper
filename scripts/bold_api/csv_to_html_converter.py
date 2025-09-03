import pandas as pd
from statistics import mean
from math import isnan


def is_float_convertible(string):
    """
    Checks if a given string can be converted to a float.

    Args:
        string (str): The input string to be checked.

    Returns:
        bool: True if the string can be converted to a float, False otherwise.
    """
    try:
        # Attempt to convert the string to a float and check if it is NaN
        return isnan(float(string))
    except ValueError:
        # If the conversion to float raises a ValueError, the string is not convertible
        return False


def identified_species_name(taxonomic_dictionary: dict):
    """
    Extracts identified species names from a taxonomic dictionary.

    Args:
        taxonomic_dictionary (dict): A dictionary containing taxonomic information.

    Returns:
        list: A list of identified species names.
    """
    names = []
    for i in range(len(taxonomic_dictionary["Genera"])):
        # Check if Genera is not convertible to float or is empty/None
        if (
            is_float_convertible(taxonomic_dictionary["Genera"][i])
            or taxonomic_dictionary["Genera"][i] == " "
            or taxonomic_dictionary["Genera"][i] == ""
            or taxonomic_dictionary["Genera"][i] == None
        ):
            # Check if Families is not convertible to float or is empty/None
            if (
                is_float_convertible(taxonomic_dictionary["Families"][i])
                or taxonomic_dictionary["Families"][i] == " "
                or taxonomic_dictionary["Families"][i] == ""
                or taxonomic_dictionary["Families"][i] == None
            ):
                # Check if Orders is not convertible to float or is empty/None
                if (
                    is_float_convertible(taxonomic_dictionary["Orders"][i])
                    or taxonomic_dictionary["Orders"][i] == " "
                    or taxonomic_dictionary["Orders"][i] == ""
                    or taxonomic_dictionary["Orders"][i] == None
                ):
                    # Check if Classes is not convertible to float or is empty/None
                    if (
                        is_float_convertible(taxonomic_dictionary["Classes"][i])
                        or taxonomic_dictionary["Classes"][i] == " "
                        or taxonomic_dictionary["Classes"][i] == ""
                        or taxonomic_dictionary["Classes"][i] == None
                    ):
                        names.append(str(taxonomic_dictionary["Phyla"][i]))
                    else:
                        names.append(str(taxonomic_dictionary["Classes"][i]))
                else:
                    names.append(str(taxonomic_dictionary["Orders"][i]))
            else:
                names.append(str(taxonomic_dictionary["Families"][i]))
        else:
            # Check if Species is not convertible to float or is empty/None
            if (
                is_float_convertible(taxonomic_dictionary["Species"][i])
                or taxonomic_dictionary["Species"][i] == " "
                or taxonomic_dictionary["Species"][i] == ""
                or taxonomic_dictionary["Species"][i] == None
            ):
                names.append(f"{taxonomic_dictionary['Genera'][i]} spp.")
            else:
                # Check if Subspecies is not convertible to float or is empty/None
                if (
                    is_float_convertible(taxonomic_dictionary["Subspecies"][i])
                    or taxonomic_dictionary["Subspecies"][i] == " "
                    or taxonomic_dictionary["Subspecies"][i] == ""
                    or taxonomic_dictionary["Subspecies"][i] == None
                ):
                    names.append(
                        f"{taxonomic_dictionary['Genera'][i]} {taxonomic_dictionary['Species'][i]}"
                    )
                else:
                    names.append(
                        f"{taxonomic_dictionary['Genera'][i]} {taxonomic_dictionary['Species'][i]} {taxonomic_dictionary['Subspecies'][i]}"
                    )
    return names


def retrieve_important_info(filepath: str):
    """
    Reads taxonomic information from a TSV file and processes it to create a dictionary with species information for each sample.

    Args:
        filepath (str): The file path to the CSV file containing taxonomic information.

    Returns:
        dict: A dictionary containing species information for each sample.
    """
    try:
        # Create a temporary csv with headers to store data
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(filepath, sep=",")
        print(df)
        # Extract relevant columns from the DataFrame
        seqids = list(df["seqid"])
        genera = list(df["Genus"])
        species = list(df["Species"])
        phyla = list(df["Phylum"])
        classes = list(df["Class"])
        orders = list(df["Order"])
        families = list(df["Family"])
        subspecies = list(df["Subspecies"])
        similarities = list(df["Similarity"])

        # Create a taxonomic dictionary
        taxonomic_dictionary = {
            "Phyla": phyla,
            "Classes": classes,
            "Orders": orders,
            "Families": families,
            "Genera": genera,
            "Species": species,
            "Subspecies": subspecies,
        }
        print(taxonomic_dictionary)

        # Get unique sample IDs
        samples = list(set(seqids))

        # Identify species names for each entry in taxonomic dictionary
        complete_name_species = identified_species_name(taxonomic_dictionary)

        # Initialize a dictionary to store species information for each sample
        specsdict = {i: [] for i in samples}

        # Process taxonomic information for each sample
        for sample in samples:
            tmp = {}
            for j in range(len(seqids)):
                if sample == seqids[j]:
                    # Update the temporary dictionary with species information and similarities
                    if complete_name_species[j] not in list(tmp.keys()):
                        tmp.update(
                            {complete_name_species[j]: [[float(similarities[j])], 1]}
                        )
                    else:
                        tmp[complete_name_species[j]][0].append(float(similarities[j]))
                        tmp[complete_name_species[j]][1] += 1
                else:
                    continue

            # Create a list of species names with average similarities
            defin = [
                key
                + " "
                + "("
                + str(round(mean(tmp[key][0]), 2))
                + "%); count = "
                + str(tmp[key][1])
                for key in list(tmp.keys())
            ]

            # Update the dictionary with species information for the current sample
            specsdict[sample] = defin
        return specsdict
    except FileNotFoundError:
        return False


def write_csv(specsdict: dict, csvpath: str) -> None:
    """
    Writes species information from a dictionary to a CSV file.

    Args:
        specsdict (dict): A dictionary containing species information for each sample.
        csvpath (str): The file path to the CSV file to be created.
    Returns:
        None
    """
    # Get the maximum length of species information lists
    lens = [len(val) for val in specsdict.values()]
    maxlength = max(lens)

    # Create a string for the species columns in the CSV header
    specsstring = ",".join(f"Taxon{i} (sim%)" for i in range(1, maxlength + 1))

    # Open the CSV file in write mode
    with open(csvpath, "w") as csv:
        # Write the header to the CSV file
        csv.write(f"SAMPLE,{specsstring}\n")

        # Write species information for each sample to the CSV file
        for key in list(specsdict.keys()):
            samplelen = maxlength - len(specsdict[key])

            # If there are missing values, fill with "NA"
            if samplelen > 0:
                samplestr = ",".join(["NA" for i in range(samplelen)])
                csv.write(f"{key},{','.join(specsdict[key])},{samplestr}\n")
            else:
                # If no missing values, write the sample information directly
                csv.write(f"{key},{','.join(specsdict[key])}\n")


def html_header(htmlpath, stylesheetpath):
    """
    Update the HTML file at 'htmlpath' by adding a custom header and styling it with the specified stylesheet.

    Parameters:
    - htmlpath (str): Path to the HTML file to be updated.
    - stylesheetpath (str): Path to the CSS stylesheet to be linked in the HTML header.
    """

    # Open the HTML file in read mode
    html = open(htmlpath, "r")
    htmlines = html.readlines()
    html.close()

    # Define the custom HTML header
    htmlheader = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="stylesheet" href="{stylesheetpath}">
    <title>Species identification table</title>
</head>
<body>
<aside>
    <nav>
    <a href="master.html"><p> Analysis Overview (Go Back).</p></a>
    <center><h1>Species identification table</h1></center>
    <br>
    <br>
    </nav>
</aside>
    """

    # Open the HTML file in write mode to update it
    html = open(htmlpath, "w")

    # Write the custom HTML header
    html.write(htmlheader)

    # Process the existing content
    for i in htmlines:
        if (
            i.replace("\n", "") != '<table border="1" class="dataframe">'
            and i.replace("\n", "") != "</table>"
        ):
            html.write(i)
        elif i.replace("\n", "") == '<table border="1" class="dataframe">':
            html.write(f"<center>{i}")
        elif i.replace("\n", "") == "</table>":
            html.write(i.replace("\n", "") + "</center>" + "\n")

    # Add the closing tags to complete the HTML structure
    closer = """</body>
</html>"""
    html.write(closer)

    # Close the file
    html.close()


def convert_csv_to_html(csvpath, htmlpath):
    """
    Convert a CSV file at 'csvpath' to an HTML file at 'htmlpath' using pandas.

    Parameters:
    - csvpath (str): Path to the CSV file to be converted.
    - htmlpath (str): Path to the HTML file to be generated.
    """

    # Read the CSV file using pandas
    dataframe = pd.read_csv(csvpath)

    # Convert the DataFrame to an HTML table and write it to the specified HTML file
    dataframe.to_html(htmlpath)
