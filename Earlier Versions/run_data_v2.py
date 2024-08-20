import csv
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

# Output templates
for mapped_rxn in mapped_rxns:
    print(mapped_rxns)

# # Example usage with the identify_hydrogenated_atoms function
# hydrogenated_atoms_list = identify_hydrogenated_atoms(templates)
# print("Hydrogenated atoms list:", hydrogenated_atoms_list)


'''
import csv
from localmapper import localmapper

#CSV-Reader
file = '/Users/peter/Documents/Code/Bachelorarbeit/Einführung/data/scopes.csv'
with open (file, newline='') as f:
    reader = csv.reader(f)
    data = list(reader)

#CSV anpassen
concatenated_entries = [item[0] + '>>' + item[8] for item in data[1:] if len(item) > 8]
rxns = concatenated_entries

#Local-Mapper
mapper = localmapper()
results = mapper.get_atom_map(rxns, return_dict=True)

#Filter nach Confidence
confidence = False
filtered_results = [entry for entry in results if isinstance(entry, dict) ] #and entry.get('confident') is confidence

#Output formatieren
formatted_entries = [
    f"{entry['mapped_rxn']} {index+1}\n{entry['template']} {index+1}\n" 
    for index, entry in enumerate(filtered_results)
]

for entry in formatted_entries:
    print(entry)
'''