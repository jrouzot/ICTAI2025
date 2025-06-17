import os
import json

# Folder containing result files
result_folder = "results/generated"
output_json = "resultsInterruptions.json"

# Dictionary to store results
results = {}

# Process each file in the result folder
for filename in os.listdir(result_folder):
    file_path = os.path.join(result_folder, filename)

    # Skip if not a file
    if not os.path.isfile(file_path):
        continue

    # Initialize storage for ub and solutions
    lb = None
    timePL = None
    solutions = []

    # Read the file and extract ub and solutions
    with open(file_path, "r") as file:
        currentSolution = []
        for line in file:
            line = line.strip()
            if line.startswith("rmaxPL"):
                parts = line.strip().split("=")
                if len(parts) == 2:
                    lb = float(parts[1])  # Extract the ub value
                else:
                    print("Error in result file parsing: rmaxPL")
            elif line.startswith("timePL"):
                parts = line.strip().split("=")
                if len(parts) == 2:
                    # Append (objective, time) tuple
                    timePL = float(parts[1])
                else:
                    print("Error in result file parsing: timePL")
            elif line.startswith("time"):
                parts = line.strip().split("=")
                if len(parts) == 2:
                    # Append (objective, time) tuple
                    currentSolution.append(float(parts[1]))
                else:
                    print("Error in result file parsing: time")
            elif line.startswith("rmax"):
                parts = line.strip().split("=")
                if len(parts) == 2 and len(currentSolution) == 1:
                    # Append (objective, time) tuple
                    currentSolution.insert(0, float(parts[1]))
                    solutions.append(currentSolution)
                    currentSolution = []
                else:
                    print("Error in result file parsing: rmax")

        # Store ub and solutions in the dictionary
        results[filename] = {
            "lb": lb,
            "timePL": timePL,
            "solutions": solutions
        }

# Write the results to a JSON file
with open(output_json, "w") as json_file:
    json.dump(results, json_file, indent=4)

print(f"Solutions and upper bounds have been saved to {output_json}")
