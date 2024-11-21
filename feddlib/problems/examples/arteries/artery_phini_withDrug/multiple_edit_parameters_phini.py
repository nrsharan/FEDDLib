import argparse
import shutil
from datetime import datetime

def modify_parameters(file_path, modification_sets):
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
    current_param_list = None
    for line in lines:
        if '<ParameterList name="' in line:
            # Check if this ParameterList is one of the targets
            current_param_list = line.split('name="')[1].split('"')[0]
            if any(current_param_list == number for mod_set in modification_sets for number in mod_set['param_list_numbers']):
                in_target_param_list = True
            else:
                in_target_param_list = False

        if in_target_param_list:
            for mod_set in modification_sets:
                if current_param_list in mod_set['param_list_numbers']:
                    modifications = mod_set['modifications']
                    for param_name, new_value in modifications.items():
                        if f'<Parameter name="{param_name}"' in line:
                            # Replace the value attribute
                            parts = line.split('value="')
                            if len(parts) > 1:
                                before_value, after = parts[1].split('"', 1)
                                line = parts[0] + f'value="{new_value}"' + after
                                print(f"Modified {param_name} to {new_value} in ParameterList {current_param_list}")

        modified_lines.append(line)

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(modified_lines)
    print(f"Modifications saved to {file_path}")

    # Verify the modifications
    verify_modifications(file_path, modification_sets)

def verify_modifications(file_path, modification_sets):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    verification_passed = True
    in_target_param_list = False
    current_param_list = None
    for line in lines:
        if '<ParameterList name="' in line:
            current_param_list = line.split('name="')[1].split('"')[0]
            if any(current_param_list == number for mod_set in modification_sets for number in mod_set['param_list_numbers']):
                in_target_param_list = True
            else:
                in_target_param_list = False

        if in_target_param_list:
            for mod_set in modification_sets:
                if current_param_list in mod_set['param_list_numbers']:
                    modifications = mod_set['modifications']
                    for param_name, expected_value in modifications.items():
                        if f'<Parameter name="{param_name}"' in line:
                            actual_value = line.split('value="')[1].split('"')[0]
                            if actual_value != expected_value:
                                print(f"Verification failed for {param_name}: expected {expected_value}, found {actual_value}")
                                verification_passed = False

    if verification_passed:
        print("All modifications verified successfully.")
    else:
        print("Some modifications did not match the expected values.")

if __name__ == "__main__":
    # Setup command line argument parsing
    parser = argparse.ArgumentParser(description='Modify XML parameters while preserving comments and their positions.')
    parser.add_argument('-f', '--file', type=str, required=True, help='Path to the XML file')
    args = parser.parse_args()

    # ParameterList numbers and domains:
    # 1 - Adventitia
    # 2 - Media
    # 3 - DegenMedia
    # 4 - Lipid
    # 5 - Calc1
    # 6 - Calc2
    # 7 - Calc3
    # 8 - FibCap

    # Define the modification sets
    modification_sets = [
        # {
        #     'param_list_numbers': ['1'], # Adventitia
        #     'modifications': {
        #         'MuA': '0.0',
        #         'Alpha2': '198.654',
        #         'Alpha1': '15.084',
        #         'Alpha4': '6.807',
        #         'KMin': '0.015'
        #     }
        # },
        {
            'param_list_numbers': ['4'], # Lipid
            'modifications': {
                # 'MuA': '0.0',
                'Alpha2': '200.0'
                # ,
                # 'Alpha1': '1.3'
            }
        },
        # {
        #     'param_list_numbers': ['5','6','7'], # Calc1, Calc2, Calc3
        #     'modifications': {
        #         # 'MuA': '0.0',
        #         'Alpha2': '151.73775'
        #     }
        # },
        # {
        #     'param_list_numbers': ['8'], # FibCap
        #     'modifications': {
        #         'MuA': '2.94',
        #         'Alpha2': '82.773',
        #         'Alpha1': '6.285',
        #         'Alpha4': '7.445'
        #     }
        # }
    ]

    # Call the function with the provided arguments
    modify_parameters(args.file, modification_sets)