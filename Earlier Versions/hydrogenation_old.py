import csv
import rdkit.Chem as Chem
from localmapper import localmapper

# CSV-Reader
file = '/Users/peter/Documents/Code/Bachelorarbeit/Einführung/data/hydrierung.csv'
with open(file, newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

# CSV anpassen
concatenated_entries = [item[0] + '>>' + item[1] for item in data[1:] if len(item) > 1]
rxns = concatenated_entries

# Local-Mapper
mapper = localmapper()
results = mapper.get_atom_map(rxns, return_dict=True)

# Filter nach Confidence
confidence = False
filtered_results = [entry for entry in results if isinstance(entry, dict)] #and entry.get('confident') is confidence

# Extract templates
mapped_rxns = [entry['mapped_rxn'] for entry in filtered_results]

#------------------------------------------------------------------------------------------------------------------------

# Get the number of hydrogens of an atom in a molecule.
def get_atom_hydrogens(mol, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    return atom.GetTotalNumHs()

# Identify hydrogenated atoms in a mapped reaction SMILES string.
def identify_hydrogenated_atoms(mapped_reaction_smiles_list):
    results = []
    
    for mapped_reaction_smiles in mapped_reaction_smiles_list:
        # Split mapped reaction SMILES into reactant and product parts
        reactant_smiles, product_smiles = mapped_reaction_smiles.split(">>")
        
        # Convert SMILES to RDKit Mol objects
        reactant_mol = Chem.MolFromSmiles(reactant_smiles)
        product_mol = Chem.MolFromSmiles(product_smiles)
        
        if reactant_mol is None or product_mol is None:
            raise ValueError(f"Invalid SMILES provided: {mapped_reaction_smiles}")
        
        # Initialize list to store hydrogenated atoms
        hydrogenated_atoms = []
        
        # Get atom mapping from reactant to product
        reactant_atom_map = {atom.GetProp("molAtomMapNumber"): atom.GetIdx() for atom in reactant_mol.GetAtoms() if atom.HasProp("molAtomMapNumber")}
        product_atom_map = {atom.GetProp("molAtomMapNumber"): atom.GetIdx() for atom in product_mol.GetAtoms() if atom.HasProp("molAtomMapNumber")}
        
        # Iterate through the atom map
        for map_num, reactant_atom_idx in reactant_atom_map.items():
            if map_num not in product_atom_map:
                continue
            
            product_atom_idx = product_atom_map[map_num]
            
            # Get number of hydrogens
            reactant_h_count = get_atom_hydrogens(reactant_mol, reactant_atom_idx)
            product_h_count = get_atom_hydrogens(product_mol, product_atom_idx)
            
            # Compare hydrogen counts to detect hydrogenation
            if product_h_count > reactant_h_count:
                hydrogenated_atoms.append(int(map_num))
        
        results.append(hydrogenated_atoms)
    
    return results

# Usage:
mapped_reaction_smiles_list = mapped_rxns
hydrogenated_atoms_list = identify_hydrogenated_atoms(mapped_reaction_smiles_list)

#------------------------------------------------------------------------------------------------------------------------

# Append mapped reactions and hydrogenated atom indices to the original CSV data
updated_data = [data[0] + ['Mapped_Reaction', 'Hydrogenated_Atoms']]  # Adding headers for new columns

for i, row in enumerate(data[1:]):
    if i < len(filtered_results):
        mapped_reaction = filtered_results[i]['mapped_rxn']
        hydrogenated_atoms = hydrogenated_atoms_list[i]
        row.append(mapped_reaction)
        row.append(str(hydrogenated_atoms).strip('[]'))  # Convert list to string without square brackets
    updated_data.append(row)

# Remove quotation marks from entries
cleaned_data = [[str(cell).replace('"', '') for cell in row] for row in updated_data]

# Write the updated data to a new CSV file
output_file = '/Users/peter/Documents/Code/Bachelorarbeit/Einführung/data/mapped_hydro.csv'
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(cleaned_data)

print("Updated CSV file has been saved.")