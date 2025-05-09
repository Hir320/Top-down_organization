import os
import chardet
import re
import pandas as pd

def find_file_recursively(root_directory, filename):
    # Walk the directory tree
    for dirpath, dirnames, filenames in os.walk(root_directory):
        if filename in filenames:
            return os.path.join(dirpath, filename)
    return None

def find_csv_in_same_directory(abf_file_path):
    if abf_file_path:
        directory = os.path.dirname(abf_file_path)
        pattern = re.compile(r'\d{4}\.\d{2}\.\d{2}\.csv')  # Regex pattern for YYYY.MM.DD.csv
        for file in os.listdir(directory):
            if pattern.match(file):
                return os.path.join(directory, file)
    return None

def find_jpgs_in_same_directory(abf_file_path):
    jpg_files = []
    if abf_file_path:
        directory = os.path.dirname(abf_file_path)
        pattern = re.compile(r'.*\.jpg$', re.IGNORECASE)  # Regex pattern for any file ending with .jpg (case insensitive)
        for file in os.listdir(directory):
            if pattern.match(file):
                jpg_files.append(os.path.join(directory, file))
    return jpg_files

def find_csv_files_recursively(root_directory, pattern=r'\d{4}\.\d{2}\.\d{2}\.csv'):
    """
    Find files with filenames matching the pattern in a directory recursively.

    Args:
    root_directory (str): The root directory to search.
    pattern (str): The regex pattern to match filenames.

    Returns:
    list: A list of matching file paths.
    """
    matching_files = []
    regex = re.compile(pattern)

    # Walk the directory tree
    for dirpath, dirnames, filenames in os.walk(root_directory):
        for filename in filenames:
            if regex.match(filename):
                matching_files.append(os.path.join(dirpath, filename))
    
    return matching_files

def find_line_with_abf_filename(filename, abf_filename):
    # Return None if abf_filename is None
    if abf_filename is None:
        return None

    with open(filename, 'rb') as file:
        raw_content = file.read()
    
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    file_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    abf_filename_pattern = re.compile(r'({})'.format(abf_filename))

    for line_num, line in enumerate(file_content):
        # Find the line containing the given abf filename
        if abf_filename_pattern.search(line):
            return line

    return None

def find_lines_with_abf_filenames(filename, abf_filenames):
    # Return an empty dictionary if abf_filenames is None
    if abf_filenames is None:
        return {}

    with open(filename, 'rb') as file:
        raw_content = file.read()
    
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    file_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    lines = {}

    for abf_filename in abf_filenames:
        abf_filename_pattern = re.compile(r'({})'.format(abf_filename))

        for line_num, line in enumerate(file_content):
            # Find the line containing the given abf filename
            if abf_filename_pattern.search(line):
                # Split the line at the abf filename and take the last part
                lines[abf_filename] = line.split(abf_filename)[-1].strip()
                break

    return lines
def find_all_lines_with_abf_filenames(filename):
    with open(filename, 'rb') as file:
        raw_content = file.read()
    
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    file_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    abf_filename_pattern = re.compile(r'([\w\-]+\.abf)')

    matching_lines = {}

    for line_num, line in enumerate(file_content):
        matches = abf_filename_pattern.findall(line)
        if matches and 'started' in line and 'completed' not in line and 'aborted' not in line:
            for match in matches:
                # Extracting only the relevant part after .abf filename
                relevant_info = line.split(match)[-1].strip()
                matching_lines[match] = f"{match}: started at {relevant_info}"
    # Remove "started at" from the strings in matching_lines
    for key, value in matching_lines.items():
        matching_lines[key] = value.replace("started at", "").strip()

    return matching_lines


def show_all_lines_with_abf_filenames(root_directory, log_filename):
    # Find the file recursively
    log_filepath = find_file_recursively(root_directory, log_filename)
    
    # Return if the log file is not found
    if log_filepath is None:
        print(f"The log file {log_filename} was not found in {root_directory}")
        return
    
    # Find the lines containing the abf filenames
    lines = find_all_lines_with_abf_filenames(log_filepath)
    return lines

def show_all_lines_with_two_stim(data_dict, regex_pattern):
    matched_files = {key: value for key, value in data_dict.items() if re.search(regex_pattern, value)}
    return matched_files

def extract_cell_id_from_log(filename, abf_filename):
    with open(filename, 'rb') as file:
        raw_content = file.read()
        
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    file_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    cell_id_pattern = re.compile(r'(## Cell\d+|## PVneuron-\d+)')
    abf_filename_pattern = re.compile(r'({})'.format(abf_filename))

    cell_id = None
    abf_line_number = None

    for line_num, line in enumerate(file_content):
        # Find the line containing the given abf filename
        if abf_filename_pattern.search(line):
            abf_line_number = line_num
            #print("ABF filename found at line: ", line_num)
            break

    if abf_line_number is not None:
        # Find the last heading of cells before the abf filename
        for line in reversed(file_content[:abf_line_number]):
            if cell_id_pattern.match(line):
                cell_id = line.strip()
                #print("Cell ID found: ", cell_id)
                break
    else:
        print("ABF filename not found in the log file.")

    return cell_id

def extract_slice_id_from_log(filename, abf_filename):
    with open(filename, 'rb') as file:
        raw_content = file.read()
        
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    file_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    slice_id_pattern = re.compile(r'(# Slice\d+)')
    abf_filename_pattern = re.compile(r'({})'.format(abf_filename))

    slice_id = None
    abf_line_number = None

    for line_num, line in enumerate(file_content):
        # Find the line containing the given abf filename
        if abf_filename_pattern.search(line):
            abf_line_number = line_num
            #print("ABF filename found at line: ", line_num)
            break

    if abf_line_number is not None:
        # Find the last heading of cells before the abf filename
        for line in reversed(file_content[:abf_line_number]):
            if slice_id_pattern.match(line):
                slice_id = line.strip()
                #print("Slice ID found: ", slice_id)
                break
    else:
        print("ABF filename not found in the log file.")

    return slice_id


def two_colour_check(filename, opsin_one, opsin_two):
    with open(filename, 'rb') as file:
        raw_content = file.read()
    # Detect file encoding
    detected_encoding = chardet.detect(raw_content)['encoding']
    
    # Read the content with the detected encoding
    text_content = raw_content.decode(detected_encoding, errors='replace').split('\n')

    # Flag to check if both "ChrimsonR" and "ChR2" are present in the text
    chrimsonr_present = opsin_one in text_content
    chr2_present = opsin_two in text_content
    both_present = chrimsonr_present and chr2_present

    return both_present

# Define a function to extract Optostimulus power and stimulus duration
def extract_stimulus_info(row):
    # Regular expression pattern to find "XXX% XX ms" after "R_XA" or "B_XA"
    pattern = r'(R|B)_(\d+\.?\d*) A (\d+)% (\d+\.?\d*) ms'
    matches = re.findall(pattern, row['log'])

    # Initialize variables
    r_power, r_duration, b_power, b_duration = None, None, None, None

    # Iterate through matches and extract information
    for match in matches:
        if match[0] == 'R':
            r_power = match[2]
            r_duration = match[3]
        elif match[0] == 'B':
            b_power = match[2]
            b_duration = match[3]

    return pd.Series([r_power, r_duration, b_power, b_duration])

def extract_stimulus_info_2(row):
    # Regular expression pattern
    pattern = r'(R|B)_\d+ A (\d+\.?\d*) (\d+\.?\d*) ms'
    matches = re.findall(pattern, row['log'])

    # Initialize variables
    r_power, r_duration, b_power, b_duration = None, None, None, None

    # Iterate through matches and extract information
    for match in matches:
        if match[0] == 'R':
            r_power = match[1]
            r_duration = match[2]
        elif match[0] == 'B':
            b_power = match[1]
            b_duration = match[2]

    return pd.Series([r_power, r_duration, b_power, b_duration])


def check_logfilename_contains(csv_data, specified_string):
    """
    Checks if the 'logfilename' column in the CSV data contains the specified string.
    
    :param csv_data: DataFrame containing the CSV data.
    :param specified_string: The string to check for in the 'logfilename' column.
    :return: DataFrame with a new boolean column indicating if the specified string is found.
    """
    # Adding a boolean column to indicate if 'logfilename' contains the specified string
    csv_data['contains_specified_string'] = csv_data['logfilename'].str.contains(specified_string, na=False)

    return csv_data
def check_all_false(csv_data, column_name):
    """
    Checks if all values in a specified column of the DataFrame are False.

    :param csv_data: DataFrame to check.
    :param column_name: The name of the column to check.
    :return: True if all values are False, False otherwise.
    """
    # Check if all values in the specified column are False
    return not csv_data[column_name].any()

def determine_br_order(row):
    log_text = row['log']
    b_index = log_text.find('B_')
    r_index = log_text.find('R_')

    if b_index == -1 and r_index == -1:
        return 'None'
    elif b_index == -1:
        return 'R'
    elif r_index == -1:
        return 'B'
    else:
        return 'B' if b_index < r_index else 'R'

    # Function to extract interval time
def extract_interval_time(log_entry):
    # Regular expression to find 'IT_XXX ms'
    match = re.search(r'IT_(\d+) ms', log_entry)
    if match:
        # Return the numeric part as an integer
        return int(match.group(1))
    else:
        
        # Return None or a default value if the pattern is not found
        return None
def check_substrings(log_entry):
    substances = []
    if 'PTX' in log_entry:
        substances.append('PTX')
    if 'NBQX' in log_entry:
        substances.append('NBQX')
    if 'AP-V' in log_entry:
        substances.append('AP-V')
    return ', '.join(substances) if substances else 'None'

def extract_rmp_data(row):
    # Check if the log contains 'RMP'
    rmp_match = re.search(r'RMP\s(-?\d+)\smV', row['log'])
    if rmp_match:
        # Extract the resting membrane potential value
        rmp_value = rmp_match.group(1)
        return pd.Series([row['logfilename'], row['SliceID'], row['CellID'], rmp_value])
    else:
        return pd.Series([None, None, None, None])
    
    