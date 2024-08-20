import pandas as pd

# Load the CSV file
file_path = 'Einführung/data/hydro_mapped.csv'  # Replace with the actual file path
df = pd.read_csv(file_path)

# Function to process the mapped_reaction and extract substrates and products
def map_reaction(row):
    mapped_reaction = row['mapped_reaction']
    substrate, product = mapped_reaction.split('>>')
    return pd.Series([substrate, product])

# Apply the function to each row and update the first two columns
df[['substrate', 'product']] = df.apply(map_reaction, axis=1)

# Save the updated DataFrame back to a CSV file
output_file_path = 'hydro_mapped.csv'  # Replace with the desired output file path
df.to_csv(output_file_path, index=False)

print("CSV file has been updated and saved as:", output_file_path)
