import csv
import deepl
import os
from dotenv import load_dotenv

def translate_csv_with_deepl(api_key, input_csv_path, output_csv_path):
    """
    Translates specified columns of a CSV file using the DeepL API.
    """
    try:
        translator = deepl.Translator(api_key)
        usage = translator.get_usage()
        if usage.character.limit is not None and usage.character.count >= usage.character.limit:
            print("DeepL character limit reached. Cannot proceed with translation.")
            return
        print(f"DeepL API usage: {usage.character.count} / {usage.character.limit or 'Unlimited'}")
    except Exception as e:
        print(f"Failed to initialize DeepL translator. Please check your API key. Error: {e}")
        return

    translated_rows = []
    
    with open(input_csv_path, 'r', encoding='utf-8') as f_in:
        reader = csv.DictReader(f_in)
        header = reader.fieldnames
        
        rows_to_translate = list(reader)
        total_rows = len(rows_to_translate)

        for i, row in enumerate(rows_to_translate):
            print(f"Translating row {i+1}/{total_rows}...")
            translated_row = row.copy()
            
            # Translate the required fields
            try:
                description = row.get('description', '')
                if description and description.strip():
                    translated_row['description'] = str(translator.translate_text(description, target_lang="KO"))
                
                params = row.get('input_parameters', '')
                if params and params.strip():
                    translated_row['input_parameters'] = str(translator.translate_text(params, target_lang="KO"))
                
                returns = row.get('returns', '')
                if returns and returns.strip():
                    translated_row['returns'] = str(translator.translate_text(returns, target_lang="KO"))

            except deepl.DeepLException as e:
                print(f"An error occurred during translation for row {i+1}: {e}")
                # Keep original text if translation fails
                pass

            translated_rows.append(translated_row)

    with open(output_csv_path, 'w', newline='', encoding='utf-8-sig') as f_out:
        writer = csv.DictWriter(f_out, fieldnames=header)
        writer.writeheader()
        writer.writerows(translated_rows)
        
    print(f"\nSuccessfully created translated file: {output_csv_path}")

if __name__ == "__main__":
    # Load environment variables from .env file
    load_dotenv()
    
    DEEPL_API_KEY = os.getenv("DEEPL_API_KEY")
    INPUT_CSV = 'function_list_en.csv'
    OUTPUT_CSV = 'function_list_ko.csv'
    
    if not DEEPL_API_KEY:
        print("Error: DEEPL_API_KEY is not set in the .env file.")
    else:
        translate_csv_with_deepl(DEEPL_API_KEY, INPUT_CSV, OUTPUT_CSV)
