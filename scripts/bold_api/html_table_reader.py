from bs4 import BeautifulSoup
import csv


def parse_html_for_topN_hits(file_path, csv_filename, seqid, nhits=10):
    # Read HTML content from the file
    with open(file_path, "r", encoding="utf-8") as html_file:
        html_content = html_file.read()
    # Parse the HTML content with BeautifulSoup
    soup = BeautifulSoup(html_content, "html.parser")

    # Find the table and extract its content
    table = soup.find("table", class_="table resultTable noborder")

    # Create a CSV file for writing
    with open(csv_filename, "w", newline="", encoding="utf-8") as csvfile:
        csv_writer = csv.writer(csvfile)

        # Extract header row
        header_row = table.find("tr", class_="taxonHeader").find_all("h6")
        headers = [header.get_text(strip=True) for header in header_row]
        csv_writer.writerow(headers)

        # Extract data rows
        data_rows = table.find_all("tr", class_=lambda x: x and "Row" in x)
        for row in data_rows:
            columns = row.find_all(["td", "th"])
            data = [column.get_text(strip=True) for column in columns]
            csv_writer.writerow(data)
        cs = open(csv_filename, "r")
        lines = cs.readlines()
        cs.close()
        c = open(csv_filename.replace(".csv", "_reformatted.csv"), "w")
        for i in range(nhits + 1):
            if i == 0:
                c.write("seqid," + ",".join(lines[i].split(",")[:12]))
            else:
                c.write(str(seqid) + "," + ",".join(lines[i].split(",")[:12]) + "\n")
        c.close()

    print(
        f"CSV file '{csv_filename.replace('.csv','_reformatted.csv')}' created successfully."
    )
    return csv_filename.replace(".csv", "_reformatted.csv")
