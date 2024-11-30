import os
import json


def update_benefit_key(folder_path):
    """
    Reads all JSON files in a folder, replaces the 'Benefit' key with 'UrbanChallenge',
    and saves the modified JSON back to the same folder.

    Args:
        folder_path (str): Path to the folder containing the JSON files.
    """
    # List all files in the folder
    for file_name in os.listdir(folder_path):
        file_path = os.path.join(folder_path, file_name)

        # Process only JSON files
        if file_name.endswith(".json") and os.path.isfile(file_path):
            with open(file_path, 'r', encoding='utf-8') as file:
                try:
                    data = json.load(file)  # Load the JSON data
                except json.JSONDecodeError as e:
                    print(f"Error decoding JSON in file {file_name}: {e}")
                    continue

            # Check if 'Benefit' exists and replace it with 'UrbanChallenge'
            if "Benefit" in data:
                data["UrbanChallenge"] = data.pop("Benefit")

                # Save the updated JSON back to the file
                with open(file_path, 'w', encoding='utf-8') as file:
                    json.dump(data, file, indent=4)
                print(f"Updated 'Benefit' key to 'UrbanChallenge' in file: {file_name}")
            else:
                print(f"No 'Benefit' key found in file: {file_name}")


folder_path = "/mnt/DATA/Ricerca/5 - Paper/11 - Greening/Code/OptimalNBS/Instances/"
cities = ['Bologna', 'Catania', 'Milan', 'Naples', 'Rome', 'Turin']
splitLabel = ["XS", "S", "M", "L"]

for c in cities:
    for s in splitLabel:
        path = folder_path + c + "/" + s + "/"
        update_benefit_key(path)
