import os
import re
from argparse import ArgumentParser

argparse = ArgumentParser()
argparse.add_argument("-f", "--filefolder",
                      help="Path to the svg files directory including subdirectories", required=True)
args = argparse.parse_args()
filefolder = args.filefolder

def update_svg_x_values2(file_path, output_path):
    with open(file_path, 'r') as file:
        svg_lines = file.readlines()
    # Patterns for identifying x, x1, x2, and cx values
    x_pattern = r'x="([0-9.]+)"'
    x1_pattern = r'x1="([0-9.]+)"'
    x2_pattern = r'x2="([0-9.]+)"'
    cx_pattern = r'cx="([0-9.]+)"'
    y_pattern = r'y="([0-9.]+)"' 
    increment = 50
    additional_increment = 10

    def update_x(match, is_monospace, add_value=0):
        nonlocal increment
        value = float(match.group(1))
        if value > 10 and 'cx="' not in line:
            if is_monospace:
                increment += additional_increment  # For monospace
            updated_value = value + increment + add_value
            if updated_value > 2000:
                updated_value = 1700
            return f'x="{updated_value}"'
        return match.group(0)

    def update_x1_x2(match, add_value=250):
        value = float(match.group(1))
        updated_value = value + add_value
        if updated_value > 2000:
            updated_value = 1700
        return f'{match.group(0)[:3]}"{updated_value}"'

    def update_cx(match):
        value2 = float(match.group(1))
        updated_value = value2 + 250
        if updated_value > 2000:
            updated_value = 1700
        return f'cx="{updated_value}"'

    def update_y_values(line):
        y_pattern = r'(y|y1|y2)="([0-9.]+)"'
        return re.sub(y_pattern, lambda match: f'{match.group(1)}="{float(match.group(2)) + 10}"', line)

    updated_lines = []
    for line in svg_lines:
        is_monospace = 'family = monospace;' in line 
        line = re.sub(x_pattern, lambda match: update_x(match, is_monospace), line)

        is_line = 'line' in line
        if is_line:
            line = re.sub(x1_pattern, lambda match: update_x1_x2(match), line)
            line = re.sub(x2_pattern, lambda match: update_x1_x2(match), line)
        is_circle = 'circle' in line
        if is_circle:
            line = re.sub(cx_pattern, lambda match: update_cx(match), line)
        line = update_y_values(line)
        # Handling x1 == x2 and adjusting the x2 value when they are equal
        if 'x1="' in line and 'x2="' in line and 'y1="' in line and 'y2="' in line:
            x1_match = re.search(r'x1="([0-9.]+)"', line)
            x2_match = re.search(r'x2="([0-9.]+)"', line)
            y1_match = re.search(r'y1="([0-9.]+)"', line)
            y2_match = re.search(r'y2="([0-9.]+)"', line)
            if x1_match and x2_match and y1_match and y2_match:
                x1_value = float(x1_match.group(1))
                x2_value = float(x2_match.group(1))
                y1_value = float(y1_match.group(1))
                y2_value = float(y2_match.group(1))
                if x1_value == x2_value and x1_value > 1500 and y1_value == y2_value:
                    # Set x2 to be slightly greater than x1 (create the correct ending point for lines)
                    new_x2_value = x1_value + 100
                    line = re.sub(r'x2="([0-9.]+)"', f'x2="{new_x2_value}"', line)
        if 'color:black;font-size :10px' in line or 'Rank' in line:
            current_y_match = re.search(y_pattern, line)
            if current_y_match:
                target_y_value = current_y_match.group(1)  # Extract the y value
                new_y_value= float(target_y_value) - 12
                #print(line)
                updated_line = line.replace(f'y="{target_y_value}"', f'y="{new_y_value}"')
                #print(updated_line)
                line = re.sub(r'10px', '16px', updated_line)
                #print(line)
        #y_match = re.search(r'y="([0-9.]+)"', line)
        #if y_match:
            #previous_y_value = y_match.group(1)
        if 'style="color:black;font-size :12px;font-family = monospace;"' in line or 'asap' in line:
            current_y_match = re.search(y_pattern, line)
            ny=re.findall(y_pattern, line)
            #print("Old",line)
            if current_y_match and len(ny) <=1:
                current_y_value = current_y_match.group(1)
                new_y_value = float(current_y_value) + 15 
                updated_line = re.sub(y_pattern, f'y="{new_y_value}"', line)
                #print("Updated:",updated_line)
            else:
                #print(ny)
                #print(line)
                #print(type(ny))
                #replacement = 'SUBSTITUTED'  # Change this to your desired replacement
                second_match = ny[1]  # Get the second match
                replacement=float(second_match)+15
                # Replace the second match in the original line
                #print(second_match, replacement)
                updated_line = re.sub(f'y="{second_match}"', f'y="{replacement}"', line)  # Only replace the first occurrence of the second match
                #print("Updated line:", updated_line)
            line = re.sub(r'12px', '14px', updated_line)
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
