import os
from pathlib import Path
from morfeus import read_xyz, BuriedVolume
import pandas as pd
from xtb_func import calc_SASA, calc_dispersion_descriptors, calc_xtb, create_sterimol

directory = "/home/student/phoeper/Projekt/data/optimized"

num_molecules_to_process = 3  # Number of molecules to process

# Initialize lists to match the number of molecules to process
num_rows = num_molecules_to_process+1
sterimol_L = [None] * num_rows
sterimol_B1 = [None] * num_rows
sterimol_B5 = [None] * num_rows
sasa_area = [None] * num_rows
sasa_volume = [None] * num_rows
disp_int_list = [None] * num_rows
disp_area_list = [None] * num_rows
disp_volume_list = [None] * num_rows
xtb_ip_list = [None] * num_rows
xtb_ea_list = [None] * num_rows
xtb_homo_list = [None] * num_rows
xtb_lumo_list = [None] * num_rows
xtb_dipole_list_1 = [None] * num_rows
xtb_dipole_list_2 = [None] * num_rows
xtb_dipole_list_3 = [None] * num_rows
xtb_ephilicity_list = [None] * num_rows
xtb_nphilicity_list = [None] * num_rows

for i in range(1, num_rows):
    filename = f"molecule_{i}_opt.xyz"
    filepath = os.path.join(directory, filename)
    if os.path.exists(filepath):
        with open(filepath, 'r') as file:
            elements, coordinates = read_xyz(filepath)
            print(f"Processing {filename}...")

        # Molecular calculations
        try:
            l, b_1, b_5 = create_sterimol(elements, coordinates)
            sterimol_L[i] = l
            sterimol_B1[i] = b_1
            sterimol_B5[i] = b_5
            area, volume = calc_SASA(elements, coordinates)
            sasa_area[i] = area
            sasa_volume[i] = volume
            disp_int, disp_area, disp_volume = calc_dispersion_descriptors(elements, coordinates)
            disp_int_list[i] = disp_int
            disp_area_list[i] = disp_area
            disp_volume_list[i] = disp_volume
            
            # xTB calculations
            try:
                xtb_ip, xtb_ea, xtb_homo, xtb_lumo, xtb_dipole_1, xtb_dipole_2, xtb_dipole_3, xtb_ephilicity, xtb_nphilicity = calc_xtb(elements, coordinates)
                xtb_ip_list[i] = xtb_ip
                xtb_ea_list[i] = xtb_ea
                xtb_homo_list[i] = xtb_homo
                xtb_lumo_list[i] = xtb_lumo
                xtb_dipole_list_1[i] = xtb_dipole_1
                xtb_dipole_list_2[i] = xtb_dipole_2
                xtb_dipole_list_3[i] = xtb_dipole_3
                xtb_ephilicity_list[i] = xtb_ephilicity
                xtb_nphilicity_list[i] = xtb_nphilicity

            except Exception as e:
                print(f"Failed xTB calculation for {filename} due to: {e}")

        except Exception as e:
            print(f"Failed descriptor calculation for {filename} due to: {e}")
    
        print(f"Successfully processed {filename} / {num_rows-1}")

    else:
        print(f"{filename} does not exist.")

# Create a new DataFrame with the new descriptors
output_df = pd.DataFrame({
    'Sterimol_L': sterimol_L[1:],
    'Sterimol_B1': sterimol_B1[1:],
    'Sterimol_B5': sterimol_B5[1:],
    'SASA_Area': sasa_area[1:],
    'SASA_Volume': sasa_volume[1:],
    'Disp_int': disp_int_list[1:],
    'Disp_Volume': disp_volume_list[1:],
    'Disp_Area': disp_area_list[1:],
    'XTB_ip': xtb_ip_list[1:],
    'XTB_ea': xtb_ea_list[1:],
    'XTB_HOMO': xtb_homo_list[1:],
    'XTB_LUMO': xtb_lumo_list[1:],
    'XTB_Dipol_1': xtb_dipole_list_1[1:],
    'XTB_Dipole_2': xtb_dipole_list_2[1:],
    'XTB_Dipole_3': xtb_dipole_list_3[1:],
    'XTB_Electrophilicity': xtb_ephilicity_list[1:],
    'XTB_Nucleophilicity': xtb_nphilicity_list[1:]
})

# Save the updated DataFrame to a CSV file
path_to_save = Path('/home/student/phoeper/Projekt/data')
filename_to_save = Path('descriptors_opt.csv')
combined_path = path_to_save / filename_to_save

output_df.to_csv(combined_path, index=False)

