from rdkit import Chem
from rdkit.Chem import AllChem
from morfeus import read_xyz, Sterimol, SASA, Dispersion, XTB
from morfeus.conformer import ConformerEnsemble

# SMILES -> XYZ Structure using Morfeus
def create_optimized_3D_morfeus(smiles_substrate: str, optimisation: str):
    ce = ConformerEnsemble.from_rdkit(smiles_substrate, optimize=optimisation)
    # Forcefield-Optimisations: 'MMFF94', 'MMFF94s', 'UFF'
    ce.prune_rmsd()
    ce.sort()
    elements = ce.conformers[0].elements
    coordinates = ce.conformers[0].coordinates
    return elements, coordinates

# SMILES -> XYZ Structure using RDkit
def create_optimized_3D_rdkit (smiles_substrate: str, nr_confs: int, nr_threads: int):

    mol = Chem.MolFromSmarts(smiles_substrate)
    mol_with_H = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(mol_with_H, numConfs=nr_confs)
    all_confs_optimized = AllChem.MMFFOptimizeMoleculeConfs(mol_with_H, numThreads=nr_threads)
    min_mol = min(all_confs_optimized, key=lambda t: t[1])
    y_value = min_mol[1]
    position = [x for x, y in enumerate(all_confs_optimized) if y[1] == y_value]
    lowest_index_number = position[0]
    Chem.MolToXYZFile(mol_with_H, 'xyz', confId=lowest_index_number)
    elements, coordinates= read_xyz('xyz')
    return elements, coordinates

# Steric Substituent Parameters
def create_sterimol (elements, coordinates):
    sterimol = Sterimol(elements, coordinates, 1, 2)
    L = sterimol.L_value
    b_1 = sterimol.B_1_value
    b_5 = sterimol.B_5_value
    return L, b_1, b_5

# Solvent accesible surface area
def calc_SASA (elements, coordinates):
    sasa = SASA(elements, coordinates)
    y = sasa.area
    z = sasa.volume
    return y, z

# Dispersion descriptors
def calc_dispersion_descriptors (elements, coordinates):
    disp = Dispersion(elements, coordinates)
    int = disp.p_int
    area = disp.area
    volume = disp.volume
    return int, area, volume

# Extended tight-binding descriptors
def calc_xtb (elements, coordinates):
    xtb = XTB(elements, coordinates)
    a = xtb.get_ip()        # Ionization-Potential
    b = xtb.get_ea()        # Electron-Affinity
    c = xtb.get_homo()      # Highest occupied molecular orbit energy
    d = xtb.get_lumo()      # Lowest unoccupied molecular orbit energy
    e = xtb.get_dipole()    # Dipole momentum
    f = xtb.get_global_descriptor("electrophilicity", corrected=True)   #Global Electrophilicity
    g = xtb.get_global_descriptor("nucleophilicity", corrected=True)    #Global Nucleophilicity
    return a, b, c, d, e[0], e[1], e[2], f, g



