import re
import subprocess
import pandas as pd
import os
import ast
from morfeus import XTB, BuriedVolume, read_xyz

def prepare_molecule_data(csv_file, xyz_folder, max_iterations=None):
    # Load the CSV file
    indices_df = pd.read_csv(csv_file)

    # Prepare a list to store the tuples of (xyz_file_path, hydrogenated_atoms)
    molecule_data = []

    # Process each row in the CSV file
    for i, row in indices_df.iterrows():
        if max_iterations is not None and i >= max_iterations:
            break

        # Extract hydrogenated atoms indices and parse the string as a list of integers
        hydrogenated_atoms = ast.literal_eval(row['hydrogenated_atoms'])
        
        # Corresponding XYZ file
        xyz_file = os.path.join(xyz_folder, f"molecule_{i+1}_opt.xyz")
        
        # Append the tuple (xyz_file_path, hydrogenated_atoms) to the list
        molecule_data.append((xyz_file, hydrogenated_atoms))

    return molecule_data

def run_xtb(filepath):
    command = ['xtb', filepath, '--gfn', '1', '--pop']
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error running xtb on {filepath}:\n{result.stderr}")
        return result.stdout
    return result.stdout

def extract_mulliken_charges(output, atom):
    charges = []
    in_mulliken_section = False

    for line in output.splitlines():
        if 'Mulliken/CM5 charges' in line:
            in_mulliken_section = True
            continue  # Skip the header line
        elif in_mulliken_section:
            if line.strip() == '':
                break  # End of the Mulliken section
            # Updated regex pattern
            match = re.match(r'\s*(\d+)([A-Za-z]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)\s+([-.\d]+)', line)
            if match:
                atom_index = int(match.group(1))
                if atom_index is atom:
                    charges.append((float(match.group(3)), float(match.group(4)), float(match.group(5)), float(match.group(6))))
            else:
                print(f"Debug: No match for line: '{line}'")

    return charges

def process_all_molecules(xyz_file, atom):
    if os.path.exists(xyz_file):
        output = run_xtb(xyz_file)
        charges = extract_mulliken_charges(output, atom)
        return charges
    else:
        print(f'File {xyz_file} does not exist.')

################################################################################

# Parameters
max_iterations = 1 # lower for testing
xyz_folder = "/home/student/phoeper/Projekt/data/optimized"
csv_file = "/home/student/phoeper/Projekt/data/mapped_hydro_checked.csv"


# Prepare the molecule data
molecule_data = prepare_molecule_data(csv_file, xyz_folder, max_iterations)
df = pd.read_csv(csv_file)
substrates= df['substrate'].tolist()[:max_iterations]

# Create a dictionary to store the calculated descriptors
descriptors_dict = {
    'atom': [],
    'substrate': [],
    'buried_volume': [],
    'distal_volume': [],
    'fraction_buried_volume': [],
    'free_volume': [],
    'electrophilicity': [],
    'nucleophilicity': [],
    'radical_attack': [],
    'dual_descriptor': [],
    'mulliken charge, n(s), n(p), n(d)': []
}

# Iterate over the molecule data and add the descriptors to the dictionary
for xyz_file, hydrogenated_atoms in molecule_data:
    elements, coordinates = read_xyz(xyz_file)

    xtb = XTB(elements, coordinates)
    electron_addition = xtb.get_fukui('electrophilicity')
    electron_removal = xtb.get_fukui('nucleophilicity')
    radical_attack = xtb.get_fukui('radical')
    dual_descriptor = xtb.get_fukui('dual')

    # Calculate the descriptors for each hydrogenated atom
    for atom in hydrogenated_atoms:

        # Calculate the buried volume
        bv = BuriedVolume(elements, coordinates, atom)
        buried_volume = bv.buried_volume
        distal_volume = bv.compute_distal_volume().distal_volume
        fraction_buried_volume = bv.fraction_buried_volume
        free_volume = bv.free_volume

        # Get the Fukui indices
        fukui_er = electron_removal[atom]
        fukui_ea = electron_addition[atom]
        fukui_ra = radical_attack[atom]
        fukui_dd = dual_descriptor[atom]

        # Get the Mulliken charges
        charges = process_all_molecules(xyz_file, atom)

        # Append the calculated to the dictionary
        descriptors_dict['atom'].append(atom)
        descriptors_dict['substrate'].append(substrates[molecule_data.index((xyz_file, hydrogenated_atoms))])
        descriptors_dict['buried_volume'].append(buried_volume)
        descriptors_dict['distal_volume'].append(distal_volume)
        descriptors_dict['fraction_buried_volume'].append(fraction_buried_volume)
        descriptors_dict['free_volume'].append(free_volume)

        descriptors_dict['electrophilicity'].append(fukui_er)
        descriptors_dict['nucleophilicity'].append(fukui_ea)
        descriptors_dict['radical_attack'].append(fukui_ra)
        descriptors_dict['dual_descriptor'].append(fukui_dd)

        descriptors_dict['mulliken charge, n(s), n(p), n(d)'].append(charges)

        print(f"Processed: {xyz_file}, Atom: {atom}")

# Create a DataFrame from the dictionary
output_df = pd.DataFrame(descriptors_dict)

# Save the DataFrame to a CSV file
output_df.to_csv("/home/student/phoeper/Projekt/data/descriptors/all_descriptors_single.csv", index=False)
