import os

def convert_to_first_format(second_file, first_file):
    # Function to generate alphabetical names for instruments
    def get_instrument_name(index):
        name = ""
        while index >= 0:
            name = chr(index % 26 + ord('A')) + name
            index = index // 26 - 1
        return name

    # Read the second file
    with open(second_file, 'r') as f:
        lines = f.readlines()

    # Initialize storage for the converted content
    instruments = []
    downlinks = []
    events = {}
    buffer_count = 0
    downlink_count = 0

    # Parse the second file
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if "buffers" in line:
            buffer_count = int(line.split()[0])
            for b in range(buffer_count):
                buffer_data = lines[i + 1 + b].strip().split()
                instrument_name = get_instrument_name(b)
                events[instrument_name] = []  # Prepare event storage for this instrument
                instruments.append(f"{instrument_name} 10 100 {buffer_data[1]} {buffer_data[0]}")
            i += buffer_count + 1
        elif "downlinks" in line:
            downlink_count = int(line.split()[0])
            for d in range(downlink_count):
                downlinks.append(lines[i + 1 + d].strip())
            i += downlink_count + 1
        elif "events for" in line:
            buffer_id = int(line.split()[-1])  # Convert to 0-based index
            instrument_name = get_instrument_name(buffer_id)
            event_count = int(line.split()[0])
            for e in range(event_count):
                events[instrument_name].append(lines[i + 1 + e].strip())
            i += event_count + 1
        else:
            i += 1

    # Write the first file
    with open(first_file, 'w') as f:
        f.write(f"{buffer_count} instruments\n")
        f.writelines("\n".join(instruments) + "\n")
        f.write(f"{downlink_count} downlinks\n")
        f.writelines("\n".join(downlinks) + "\n")
        for instrument in events.keys():
            f.write(f"0 opportunities for {instrument}\n")
        for instrument, event_list in events.items():
            f.write(f"{len(event_list)} events for {instrument}\n")
            f.writelines("\n".join(event_list) + "\n")

def process_folder_recursively(in_folder, out_folder):
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    for root, dirs, files in os.walk(in_folder):
        # Recreate the folder structure in the output folder
        relative_path = os.path.relpath(root, in_folder)
        target_folder = os.path.join(out_folder, relative_path)
        if not os.path.exists(target_folder):
            os.makedirs(target_folder)

        for file in files:
            # Process only files with .txt extension
            if file.endswith(".txt"):
                in_file_path = os.path.join(root, file)
                out_file_path = os.path.join(target_folder, file)  # Use the same name in the out folder
                convert_to_first_format(in_file_path, out_file_path)
                print(f"Converted: {in_file_path} -> {out_file_path}")

# Example usage
process_folder_recursively("generatedTmp", "generated")
