import os
import re
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument("-f", "--filefolder",
                      help="Path to the svg files directory including subdirectories", required=True)
args = argparse.parse_args()
filefolder = args.filefolder
#scorepattern=r'monospace;" >[([0-9.])/([0-9.])] </text>'
scorepattern = r'monospace;" >\s*\[([0-9]+(?:\.[0-9]+)?)\/([0-9]+(?:\.[0-9]+)?)\]\s*</text>'

#x_pattern = r'x="([0-9.]+)"'

def update_svg_x_values2(file_path, output_path):
    with open(file_path, 'r') as file:
        svg_lines = file.readlines()
    updated_lines = []
    for line in svg_lines:
        # Remove lines containing 'asap score' or matching the scorepattern
        if 'asap score' in line or re.search(scorepattern, line):
            continue  # Skip this line
        # Change font size where applicable
        if 'style="color:black;font-size :12;font-family = monospace;"' in line:
            line = re.sub(r'font-size :12', 'font-size :14', line)
        if 'style="color:black;font-size :10;font-family = monospace;"' in line:
            line = re.sub(r'font-size :10', 'font-size :14', line)
        # Append the updated (or unchanged) line to the list
        updated_lines.append(line)
    with open(output_path, 'w') as file:
        file.writelines(updated_lines)
# Function to process all SVG files in the folder
def process_all_files_in_folder(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            if file_name.endswith('.groups.svg'):
                file_path = os.path.join(root, file_name)
                base_name, _ = os.path.splitext(file_name)
                output_file_name = base_name + '.updated.svg'
                output_path = os.path.join(root, output_file_name)
                update_svg_x_values2(file_path, output_path)
                print(f"Processed {file_path} -> {output_path}")

if __name__ == "__main__":
    print("Correcting SVG file parameters.")
    process_all_files_in_folder(filefolder)
