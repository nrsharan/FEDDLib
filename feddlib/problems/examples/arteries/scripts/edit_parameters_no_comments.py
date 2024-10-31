import xml.etree.ElementTree as ET
import argparse
import shutil
from datetime import datetime

def modify_parameters(file_path, param_list_numbers, modifications):
    # Load the XML file
    tree = ET.parse(file_path)
    root = tree.getroot()

    # Backup the original file with the current date and time
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file_path = f"{file_path}_{current_time}.bak"
    shutil.copy(file_path, backup_file_path)
    print(f"Backup created: {backup_file_path}")

    # Iterate through each ParameterList and modify the specified parameters
    for param_list in root.findall('.//ParameterList'):
        if param_list.get('name') in param_list_numbers:
            for param in param_list.findall('.//Parameter'):
                param_name = param.get('name')
                if param_name in modifications:
                    param.set('value', str(modifications[param_name]))
                    print(f"Modified {param_name} in ParameterList {param_list.get('name')} to {modifications[param_name]}")

    # Write the modified XML back to the file
    tree.write(file_path)
    print(f"Modifications saved to {file_path}")

if __name__ == "__main__":
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description='Modify XML parameters.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the XML file')
    args = parser.parse_args()

    # Define the ParameterList numbers and modifications
    param_list_numbers = ['1', '4', '5']
    modifications = {
        'Density': '1.1e0',  # Example modification
        'Alpha1': '15.0'     # Example modification
    }

    # Call the function with the provided arguments
    modify_parameters(args.file, param_list_numbers, modifications)
