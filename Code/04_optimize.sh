#!/bin/bash

# Directory containing the molecule files
molecule_dir="/home/student/phoeper/Projekt/data/xyz_sub"

# Output directory for optimized molecules
output_dir="/home/student/phoeper/Projekt/data/opt_sub"

# Log file to keep track of the optimization process
log_file="optimize.log"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Start logging
echo "Optimization started at $(date)" > "$log_file"

# Loop through each molecule file in the directory
for molecule_ in "$molecule_dir"/*.xyz; do
    # Get the base name of the molecule file (without directory and extension)
    base_name=$(basename "$molecule_" .xyz)
    
    echo "Optimizing $molecule_..." | tee -a "$log_file"

    # Run the xtb optimization
    if xtb "$molecule_" --opt --namespace "$base_name" &>> "$log_file"; then
        # Rename the xtbopt.xyz file to include the base name
        mv "$base_name.xtbopt.xyz" "$output_dir/${base_name}_opt.xyz"
        echo "Optimized $molecule_ to $output_dir/${base_name}_opt.xyz" | tee -a "$log_file"
    else
        echo "Failed to optimize $molecule_" | tee -a "$log_file"
    fi

    # Clean up any other generated files by xtb
    rm -f "$base_name."*
done

echo "All molecules optimized!" | tee -a "$log_file"

#TODO: 'find -name "*.xtboptok" -type f -delete' in Projekt
