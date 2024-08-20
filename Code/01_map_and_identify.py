import csv
import rdkit.Chem as Chem
from localmapper import localmapper

# Function to read CSV file
def read_csv(file_path):
    with open(file_path, newline='') as f:
        reader = csv.reader(f)
        return list(reader)

# Function to write CSV file
def write_csv(file_path, data):
    with open(file_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(data)
    print(f"Updated CSV file has been saved to {file_path}")

# Function to concatenate entries for reaction SMILES
def concatenate_entries(data):
    return [item[0] + '>>' + item[1] for item in data[1:] if len(item) > 1]

# Function to get mapped reactions using localmapper
def get_mapped_reactions(rxns):
    mapper = localmapper()
    results = []
    for idx, rxn in enumerate(rxns):
        print(f"Mapping reaction {idx + 1}/{len(rxns)}: {rxn}")
        try:
            result = mapper.get_atom_map([rxn], return_dict=True)
            if isinstance(result[0], dict):
                results.append(result[0])
            else:
                print(f"Unexpected result format for reaction {rxn}: {result}")
        except Exception as e:
            print(f"Error mapping reaction at index {idx}: {rxn}")
            print(e)
            break
    return results

# Function to get the number of hydrogens on an atom
def get_atom_hydrogens(mol, atom_idx):
    atom = mol.GetAtomWithIdx(atom_idx)
    return atom.GetTotalNumHs()

# Function to identify hydrogenated atoms in mapped reaction SMILES strings
def identify_hydrogenated_atoms(mapped_reaction_smiles_list):
    results = []
    for idx, mapped_reaction_smiles in enumerate(mapped_reaction_smiles_list):
        #print(f"Processing reaction {idx + 1}/{len(mapped_reaction_smiles_list)}: {mapped_reaction_smiles}")
        try:
            reactant_smiles, product_smiles = mapped_reaction_smiles.split(">>")
            reactant_mol = Chem.MolFromSmiles(reactant_smiles)
            product_mol = Chem.MolFromSmiles(product_smiles)

            if reactant_mol is None or product_mol is None:
                raise ValueError(f"Invalid SMILES provided at index {idx}: {mapped_reaction_smiles}")

            hydrogenated_atoms = []
            reactant_atom_map = {atom.GetProp("molAtomMapNumber"): atom.GetIdx() for atom in reactant_mol.GetAtoms() if atom.HasProp("molAtomMapNumber")}
            product_atom_map = {atom.GetProp("molAtomMapNumber"): atom.GetIdx() for atom in product_mol.GetAtoms() if atom.HasProp("molAtomMapNumber")}

            for map_num, reactant_atom_idx in reactant_atom_map.items():
                if map_num not in product_atom_map:
                    continue
                product_atom_idx = product_atom_map[map_num]
                reactant_h_count = get_atom_hydrogens(reactant_mol, reactant_atom_idx)
                product_h_count = get_atom_hydrogens(product_mol, product_atom_idx)
                if product_h_count > reactant_h_count:
                    hydrogenated_atoms.append(int(map_num))

            results.append(hydrogenated_atoms)
        except Exception as e:
            print(f"Error processing reaction at index {idx}: {mapped_reaction_smiles}")
            print(e)
            continue
    return results

# Function to update CSV data with mapped reactions and hydrogenated atom indices
def update_csv_data(data, mapped_reactions, hydrogenated_atoms_list):
    updated_data = [data[0] + ['mapped_reaction', 'hydrogenated_atoms']]
    for i, row in enumerate(data[1:]):
        if i < len(mapped_reactions):
            mapped_reaction = mapped_reactions[i]
            hydrogenated_atoms = hydrogenated_atoms_list[i]
            row.append(mapped_reaction)
            row.append(hydrogenated_atoms)
        updated_data.append(row)
    return updated_data

# Main function to run the process
def main():
    input_file = '/Users/peter/Documents/Code/Bachelorarbeit/Einführung/data/hydrierung_v2.csv'
    output_file = '/Users/peter/Documents/Code/Bachelorarbeit/Einführung/data/mapped_hydro_v2.csv'

    print("Reading CSV file...")
    data = read_csv(input_file)
    print("CSV file read successfully.")
    print("Concatenating entries for reaction SMILES...")
    concatenated_entries = concatenate_entries(data)
    print("Entries concatenated successfully.")
    print("Getting mapped reactions using localmapper...")
    filtered_results = get_mapped_reactions(concatenated_entries)
    print("Mapped reactions obtained successfully.")

    mapped_reactions = [entry['mapped_rxn'] for entry in filtered_results]
    hydrogenated_atoms_list = identify_hydrogenated_atoms(mapped_reactions)

    updated_data = update_csv_data(data, mapped_reactions, hydrogenated_atoms_list)
    write_csv(output_file, updated_data)

if __name__ == "__main__":
    main()
