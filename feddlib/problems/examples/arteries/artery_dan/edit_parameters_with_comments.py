import argparse
import shutil
from datetime import datetime

def modify_parameters(file_path, param_list_numbers, modifications):
    # Read the XML file as text
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Backup the original file with the current date and time
    current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_file_path = f"{file_path}_{current_time}.bak"
    shutil.copy(file_path, backup_file_path)
    print(f"Backup created: {backup_file_path}")

    # Modify the specified parameters directly in the text
    modified_lines = []
    in_target_param_list = False
    for line in lines:
        if '<ParameterList name="' in line:
            # Check if this ParameterList is one of the targets
            for number in param_list_numbers:
                if f'name="{number}"' in line:
                    in_target_param_list = True
                    break
            else:
                in_target_param_list = False

        if in_target_param_list:
            for param_name, new_value in modifications.items():
                if f'<Parameter name="{param_name}"' in line:
                    # Replace the value attribute
                    parts = line.split('value="')
                    if len(parts) > 1:
                        before_value, after = parts[1].split('"', 1)
                        line = parts[0] + f'value="{new_value}"' + after
                        print(f"Modified {param_name} to {new_value}")

        modified_lines.append(line)

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)
    print(f"Modifications saved to {file_path}")

if __name__ == "__main__":
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description='Modify XML parameters while preserving comments and their positions.')
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
