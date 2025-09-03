import time
import sys
import subprocess as sp
import os
from argparse import ArgumentParser
import pandas as pd
from statistics import mean
import xml.etree.ElementTree as ET
from csv_to_html_converter import *
from html_table_reader import *
from data_loader import *
from multiprocessing import Pool
import shutil

argparse = ArgumentParser()
argparse.add_argument(
    "-i", "--infile", help="Path to input file with sequences", required=True
)
argparse.add_argument(
    "-o", "--outdir", help="Path to output directory", required=True)
argparse.add_argument(
    "-p",
    "--primername",
    help="Name of the primer/locus to search [COX1, ITS or MATK_RBCL]",
    required=True,
)
argparse.add_argument(
    "-n",
    "--nhits",
    help="Number of maximum hits, between 0 and 100",
    required=False,
    default=10,
    type=int,
)
argparse.add_argument(
    "-c", "--cssfile", help="Path to the css file", required=True)


args = argparse.parse_args()


inf = args.infile
out = args.outdir
prim = args.primername
nh = args.nhits
cssfile = args.cssfile


if not os.path.exists(out):
    print("Outdir does not exist, creating it...")
    os.makedirs(out)

if nh > 100:
    print(
        "You set the number of top hits over the available threshold, setting them to 100..."
    )
    nh = 100

primertodb = {"ITS": "ITS", "COX1": "COX1_SPECIES_PUBLIC",
              "MATK_RBCL": "MATK_RBCL"}
primertotab = {"ITS": "fungiTabPane", "MATK_RBCL": "plantTabPane"}


def submit_to_api(sequence, header, outdir, primer, tophits, css):
    global primertotab, primertodb
    print(
        f"Submitted sequence: {header} at time {time.time()}", file=sys.stderr)
    if primer == "COX1":
        outputfile = os.path.join(outdir, header + "_identified.xml")

        command = sp.run(
            f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}',
            shell=True,
        )
        print(f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}'
              )
        while command.returncode != 0:
            time.sleep(10)
            command = sp.run(
                f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}',
                shell=True,
            )
        # Parse the XML file
        tree = ET.parse(outputfile)
        root = tree.getroot()
        # Find the "matches" element
        matches_element = root.find("match")

        # Check if the "matches" element is empty
        if type(matches_element) == type(None):
            print(
                f"Sequence {header} does not match with anything, proceeding with the analysis..."
            )
            return
        csvfile = outputfile.replace("identified.xml", "summarized_tmp.csv")
        c = open(csvfile, "w")
        c.write("Sample,Taxonomicidentification,Similarity\n")
        # Iterate over 'taxonomicidentification' and 'similarity' elements simultaneously
        counter = 0
        for elem1, elem2 in zip(
            root.iter("taxonomicidentification"), root.iter("similarity")
        ):
            c.write(
                f"{header.replace(',','-')},{elem1.text.replace(',','-')},{elem2.text.replace(',','-')}\n"
            )
            counter += 1
            if counter == tophits:
                break
        c.close()
        complete_name_species = list(pd.read_csv(
            csvfile)["Taxonomicidentification"])
        samples = list(set(list(pd.read_csv(csvfile)["Sample"])))
        unset_samples = list(pd.read_csv(csvfile)["Sample"])
        similarities = list(pd.read_csv(csvfile)["Similarity"])
        specsdict = {i: [] for i in samples}
        # Process taxonomic information for each sample
        for sample in samples:
            tmp = {}
            for j in range(len(complete_name_species)):
                if sample == unset_samples[j]:
                    # Update the temporary dictionary with species information and similarities
                    if complete_name_species[j] not in list(tmp.keys()):
                        tmp.update(
                            {complete_name_species[j]: [
                                [float(similarities[j])], 1]}
                        )
                    else:
                        tmp[complete_name_species[j]][0].append(
                            float(similarities[j]))
                        tmp[complete_name_species[j]][1] += 1
                else:
                    continue

            # Create a list of species names with average similarities
            defin = [
                key
                + " "
                + "("
                + str(round((mean(tmp[key][0]) * 100), 2))
                + "%); count = "
                + str(tmp[key][1])
                for key in list(tmp.keys())
            ]

            # Update the dictionary with species information for the current sample
            specsdict[sample] = defin
        write_csv(specsdict, outputfile.replace(
            "identified.xml", "summarized.csv"))
        convert_csv_to_html(
            outputfile.replace("identified.xml", "summarized.csv"),
            outputfile.replace("identified.xml", "summarized.html"),
        )
        html_header(outputfile.replace(
            "identified.xml", "summarized.html"), css)
    elif primer == "ITS" or primer == "MATK_RBCL":
        outputfile = os.path.join(outdir, header + "_identified.html")
        command = sp.run(
            f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={sequence}" > {outputfile}',
            shell=True,
        )
        while command.returncode != 0:
            time.sleep(10)
            command = sp.run(
                f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={sequence}" > {outputfile}',
                shell=True,
            )
        with open(outputfile, "r", encoding="utf-8") as html_file:
            html_content = html_file.read()
        html_file.close()
        if not "No match" in html_content:
            datatable_url = (
                "https://www.boldsystems.org"
                + html_content.split('<span style="text-decoration: none" result="')[
                    1
                ].split('">')[0]
            )
            outputfile = outputfile.replace(
                "_identified.html", "_datatable.html")
            command = sp.run(
                f'curl -s "{datatable_url}" > {outputfile}', shell=True)
            while command.returncode != 0:
                time.sleep(10)
                command = sp.run(
                    f'curl -s "{datatable_url}" > {outputfile}', shell=True)
            boldf = parse_html_for_topN_hits(
                outputfile,
                outputfile.replace(".html", ".csv"),
                nhits=tophits,
                seqid=header,
            )
            specsdict = retrieve_important_info(boldf)
            if specsdict != False:
                write_csv(
                    specsdict, boldf.replace(
                        "_reformatted.csv", "_summarized.csv")
                )
                convert_csv_to_html(
                    boldf.replace("_reformatted.csv", "_summarized.csv"),
                    htmlpath=boldf.replace(
                        "_reformatted.csv", "_summarized.html"),
                )
                html_header(
                    htmlpath=boldf.replace(
                        "_reformatted.csv", "_summarized.html"),
                    stylesheetpath=css,
                )
                print(
                    f"{boldf.replace('_reformatted.csv', '_summarized.html')} file created successfully"
                )
            else:
                print(
                    f"File {boldf} does not exist, proceeding with the analysis...")
        else:
            print(
                f"Sequence {header} does not match with anything, proceeding with the analysis..."
            )
            return
    else:
        raise ValueError(
            f"BOLD cannot be searched with {primer}: user COX1 (animals), ITS (fungi) or MATK_RBCL (plants)"
        )


def submit_chunk(chunk):
    with Pool() as pool:
        pool.starmap(
            submit_to_api,
            [
                (sequence, header, outdir, primer, tophits, css)
                for sequence, header, outdir, primer, tophits, css in chunk
            ],
        )


def submit_first(sequence, header, outdir, primer):
    global primertotab, primertodb

    if primer == "COX1":
        outputfile = os.path.join(outdir, header + "_identifiedfirst.xml")
        command = sp.run(
            f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}',
            shell=True,
        )
        print(
            f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}')
        C = 0
        while command.returncode != 0:
            if C > 100:
                break
            time.sleep(10)
            command = sp.run(
                f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={sequence}" | xmllint --format - > {outputfile}',
                shell=True,
            )
            C += 1
        # Parse the XML file
        sim = 0
        tree = ET.parse(outputfile)
        root = tree.getroot()
        matches_element = root.find("match")
        if type(matches_element) == type(None):
            print(
                f"Sequence {header} in fwd direction does not match with anything, trying rev..."
            )
            outputfile = os.path.join(outdir, header + "_identifiedfirst.xml")
            command = sp.run(
                f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={reverse_complement(sequence)}" | xmllint --format - > {outputfile}',
                shell=True,
            )
            print(
                f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={reverse_complement(sequence)}" | xmllint --format - > {outputfile}')
            while command.returncode != 0:
                time.sleep(10)
                command = sp.run(
                    f'curl -s "https://www.boldsystems.org/index.php/Ids_xml?db={primertodb[primer]}&sequence={reverse_complement(sequence)}" | xmllint --format - > {outputfile}',
                    shell=True,
                )
            # Parse the XML file
            sim = 0
            tree = ET.parse(outputfile)
            root = tree.getroot()
            matches_element = root.find("match")
            if type(matches_element) == type(None):
                print(
                    f"Sequence {header} does not match with anything both fwd and rev, proceeding with the analysis..."
                )
                return "nomatch", False
            counter = 0
            # Iterate over 'taxonomicidentification' and 'similarity' elements simultaneously
            for elem1, elem2 in zip(
                root.iter("taxonomicidentification"), root.iter("similarity")
            ):
                if counter == 0:
                    sim += float(elem2.text)
                    counter += 1
                if counter > 0:
                    break
            if sim >= 0.8:
                return "rev", True
            return "fwd", True
        counter = 0
        # Iterate over 'taxonomicidentification' and 'similarity' elements simultaneously
        for elem1, elem2 in zip(
            root.iter("taxonomicidentification"), root.iter("similarity")
        ):
            if counter == 0:
                sim += float(elem2.text)
                counter += 1
            if counter > 0:
                break
        if sim >= 0.8:
            return "fwd", True
        return "rev", True
    elif primer == "ITS" or primer == "MATK_RBCL":
        outputfile = os.path.join(outdir, header + "_identifiedfirst.html")
        command = sp.run(
            f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={sequence}" > {outputfile}',
            shell=True,
        )
        print(
            f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={sequence}" > {outputfile}')
        while command.returncode != 0:
            time.sleep(10)
            command = sp.run(
                f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={sequence}" > {outputfile}',
                shell=True,
            )
        with open(outputfile, "r", encoding="utf-8") as html_file:
            html_content = html_file.read()
        html_file.close()
        if not "No match" in html_content and not (html_content == "" or html_content == "\n"):
            datatable_url = (
                "https://www.boldsystems.org"
                + html_content.split('<span style="text-decoration: none" result="')[
                    1
                ].split('">')[0]
            )
            outputfile = outputfile.replace(
                "_identifiedfirst.html", "_datatablefirst.html"
            )
            command = sp.run(
                f'curl -s "{datatable_url}" > {outputfile}', shell=True)
            while command.returncode != 0:
                time.sleep(10)
                command = sp.run(
                    f'curl -s "{datatable_url}" > {outputfile}', shell=True)
            boldf = parse_html_for_topN_hits(
                outputfile, outputfile.replace(".html", ".csv"), nhits=1, seqid=header
            )
            df = pd.read_csv(boldf)
            sim = list(df["Similarity"])[0]
            if float(sim) >= 80:
                return "fwd", True
            return "rev", True
        else:
            print(
                f"Sequence {header} in fwd direction does not match with anything, trying rev..."
            )
            outputfile = os.path.join(outdir, header + "_identifiedfirst.html")
            command = sp.run(
                f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={reverse_complement(sequence)}" > {outputfile}',
                shell=True,
            )
            print(
                f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={reverse_complement(sequence)}" > {outputfile}')
            while command.returncode != 0:
                time.sleep(10)
                command = sp.run(
                    f'curl -s -X POST "https://www.boldsystems.org/index.php/IDS_BlastRequest" -d "tabtype={primertotab[primer]}" -d "searchdb={primertodb[primer]}" -d "sequence={reverse_complement(sequence)}" > {outputfile}',
                    shell=True,
                )
            with open(outputfile, "r", encoding="utf-8") as html_file:
                html_content = html_file.read()
            html_file.close()
            if not "No match" in html_content and not (html_content == "" or html_content == "\n"):
                datatable_url = (
                    "https://www.boldsystems.org"
                    + html_content.split('<span style="text-decoration: none" result="')[
                        1
                    ].split('">')[0]
                )
                outputfile = outputfile.replace(
                    "_identifiedfirst.html", "_datatablefirst.html"
                )
                command = sp.run(
                    f'curl -s "{datatable_url}" > {outputfile}', shell=True)
                while command.returncode != 0:
                    time.sleep(10)
                    sp.run(
                        f'curl -s "{datatable_url}" > {outputfile}', shell=True)
                boldf = parse_html_for_topN_hits(
                    outputfile, outputfile.replace(".html", ".csv"), nhits=1, seqid=header
                )
                df = pd.read_csv(boldf)
                sim = list(df["Similarity"])[0]
                if float(sim) >= 80:
                    return "rev", True
                return "fwd", True
            else:
                print(
                    f"Sequence {header} does not match with anything both fwd and rev, proceeding with the analysis..."
                )
                return "nomatch", False


if __name__ == "__main__":
    seqsdict = load_data(inf)
    sequences = list(seqsdict.values())
    headers = list(seqsdict.keys())

    direction, match = submit_first(
        sequences[0].upper(), headers[0], out, prim)
    c = 0

    while match == False:
        c += 1
        direction, match = submit_first(
            sequences[c].upper(), headers[c], out, prim)
        if match == True:
            break
        if c == len(sequences) - 1:
            print("No sequence matched with the database, quitting...")
            sys.exit(0)
    if direction == "fwd":
        data_chunks = [
            (sequences[i].upper(), headers[i], out, prim, nh, cssfile)
            for i in range(len(headers))
        ]

        for i in range(0, len(data_chunks), 10):
            chunk = data_chunks[i: i + 10]
            submit_chunk(chunk)
    else:
        data_chunks = [
            (
                reverse_complement(sequences[i]).upper(),
                headers[i],
                out,
                prim,
                nh,
                cssfile,
            )
            for i in range(len(headers))
        ]

        for i in range(0, len(data_chunks), 10):
            chunk = data_chunks[i: i + 10]
            submit_chunk(chunk)
    flnms = list_files(out)
    lines_to_last_csv = []
    for f in flnms:
        if (
            not f.endswith("_reformatted.csv")
            and not f.endswith("_summarized.csv")
            and not f.endswith("_summarized.html")
        ):
            os.remove(f)
        elif f.endswith("_reformatted.csv"):
            movedirpath = mk_or_mv_dir(os.path.join(out, "bold_outputs"))
            try:
                shutil.move(f, movedirpath)
            except Exception as e:
                print(e)
        elif f.endswith("_summarized.csv"):
            if len(lines_to_last_csv) == 0:
                c = open(f)
                lines = c.readlines()
                for line in lines:
                    lines_to_last_csv.append(line)
                c.close()
            else:
                c = open(f)
                lines = c.readlines()
                if len(lines[0].split(",")) > len(lines_to_last_csv[0].split(",")):
                    lines_to_last_csv[0] = lines[0]
                for i in range(1, len(lines)):
                    lines_to_last_csv.append(lines[i])
                c.close()
            movedirpath = mk_or_mv_dir(os.path.join(out, "summarized_outputs"))
            try:

                shutil.move(f, movedirpath)
            except Exception as e:
                print(e)
        elif f.endswith("_summarized.html"):
            movedirpath = mk_or_mv_dir(os.path.join(out, "htmls"))
            try:
                shutil.move(f, movedirpath)
            except Exception as e:
                print(e)
    last_csv = os.path.join(out, "summarized_outputs/final.csv")
    last_html = os.path.join(out, "htmls/final.html")
    lc = open(last_csv, "w")
    for line in lines_to_last_csv:
        lc.write(line)
    lc.close()
    # convert_csv_to_html(
    #        last_csv,
    #        last_html,
    #    )
    # html_header(last_html, cssfile)
