import re
import csv

def parse_md_to_csv(input_md_path, en_csv_path):
    """
    Parses a structured markdown file and outputs a CSV file.
    """
    functions = []
    current_directory = ""
    current_filepath = ""
    current_function = None
    parsing_parameters = False

    with open(input_md_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()

            dir_match = re.match(r'^## Directory: `(.+?)`', line)
            if dir_match:
                current_directory = dir_match.group(1)
                parsing_parameters = False
                continue

            file_match = re.match(r'^### File: (.+)', line)
            if file_match:
                current_filepath = file_match.group(1)
                parsing_parameters = False
                continue

            func_match = re.match(r'^#### Function: (.+)', line)
            if func_match:
                if current_function:
                    functions.append(current_function)
                
                current_function = {
                    "id": len(functions) + 1,
                    "directory": current_directory,
                    "filepath": current_filepath,
                    "function_name": func_match.group(1),
                    "description": "",
                    "input_parameters": [],
                    "returns": ""
                }
                parsing_parameters = False
                continue

            if current_function:
                desc_match = re.match(r'^- \*\*Description:\*\* (.+)', line)
                if desc_match:
                    current_function["description"] = desc_match.group(1)
                    parsing_parameters = False
                    continue

                params_header_match = re.match(r'^- \*\*Parameters:\*\*', line)
                if params_header_match:
                    parsing_parameters = True
                    continue

                returns_match = re.match(r'^- \*\*Returns:\*\* (.+)', line)
                if returns_match:
                    current_function["returns"] = returns_match.group(1)
                    parsing_parameters = False
                    continue
                
                if parsing_parameters and (line.startswith('- ') or line.startswith('* ')):
                    # Clean up the parameter line
                    param_text = re.sub(r'^\s*[-*]\s*', '', line).strip()
                    current_function["input_parameters"].append(param_text)

    if current_function:
        functions.append(current_function)

    # Process collected parameters into a single string
    for func in functions:
        func['input_parameters'] = ", ".join(func['input_parameters'])

    # Write to CSV files
    header = ['id', 'directory', 'filepath', 'function_name', 'description', 'input_parameters', 'returns']
    
    with open(en_csv_path, 'w', newline='', encoding='utf-8') as f_en:
        writer_en = csv.DictWriter(f_en, fieldnames=header)
        writer_en.writeheader()
        writer_en.writerows(functions)

    print(f"Successfully created {en_csv_path}")


if __name__ == "__main__":
    MARKDOWN_FILE = 'myR/function_analysis.md'
    EN_CSV_FILE = 'function_list_en.csv'
    
    print("Starting markdown parsing and conversion to CSV...")
    parse_md_to_csv(MARKDOWN_FILE, EN_CSV_FILE)
    print("All tasks completed.")
