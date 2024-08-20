import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from xtb_func import create_optimized_3D_morfeus

# Parameters
input_csv = '/home/student/phoeper/Projekt/data/hydro_mapped.csv'   # Input CSV file path
output_dir = '/home/student/phoeper/Projekt/data/xyz_prod/'  # Output directory

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Read input CSV file
df = pd.read_csv(input_csv)

# Function to write XYZ file
def write_xyz(filename, elements, coordinates):
    with open(filename, 'w') as f:
        f.write(f"{len(elements)}\n\n")  # Number of atoms
        for element, coord in zip(elements, coordinates):
            x, y, z = coord
            f.write(f"{element} {x:.6f} {y:.6f} {z:.6f}\n")

# Process rows in the DataFrame
cap = 365
for index, row in enumerate(df.iterrows()):
    if index >= cap:
        break
    row = row[1]
    smiles = row['product'] #'substrate' or 'product'
    
    try:
        elements, coordinates = create_optimized_3D_morfeus(smiles, optimisation='MMFF94')
    except Exception as e:
        print(f"Error processing {smiles}: {e}")
        continue

    print(f"Processing molecule {index + 1}")

    # Define the filename for the XYZ file
    filename = os.path.join(output_dir, f"molecule_{index + 1}.xyz")
    
    # Write the XYZ data to the file
    write_xyz(filename, elements, coordinates)

print(f"Substrates and their XYZ coordinates saved to: {output_dir}")
