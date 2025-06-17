import os
import subprocess
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process instance files with omdpi solver.")
parser.add_argument("instance_folder", type=str, help="Path to the folder containing instance files.")
parser.add_argument("results_folder", type=str, help="Path to the folder where results will be stored.")
args = parser.parse_args()

# Extract folder paths from arguments
instance_folder = args.instance_folder
results_folder = args.results_folder
executable = "build/omdpi"

# Ensure the results folder exists
os.makedirs(results_folder, exist_ok=True)

# Iterate over all files in the instance folder
for instance_file in os.listdir(instance_folder):
    instance_path = os.path.join(instance_folder, instance_file)

    # Skip if it's not a file
    if not os.path.isfile(instance_path):
        continue

    # Prepare the output file path
    result_file = os.path.join(results_folder, f"{os.path.splitext(instance_file)[0]}.txt")

    # Run the executable and save output
    with open(result_file, "w") as output_file:
        try:
            subprocess.run([executable, instance_path], stdout=output_file, stderr=subprocess.STDOUT, check=True)
            print(f"Processed: {instance_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error processing {instance_file}: {e}")
