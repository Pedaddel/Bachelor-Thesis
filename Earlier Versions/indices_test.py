from rdkit import Chem
from rdkit.Chem import AllChem

def get_reacting_atoms(reaction_smarts):
    # Parse the reaction SMARTS
    rxn = AllChem.ReactionFromSmarts(reaction_smarts)
    
    # Generate reactant and product molecules
    reactant_mol = rxn.GetReactants()[0]
    product_mol = rxn.GetProducts()[0]
    
    # Get the atom mapping from reactant to product
    atom_map = {}
    for atom in reactant_mol.GetAtoms():
        reactant_idx = atom.GetIdx()
        product_idx = atom.GetAtomMapNum()
        if product_idx != 0:  # Only consider mapped atoms
            atom_map[reactant_idx] = product_idx - 1  # Convert to zero-based index
    
    # Find the reacting atoms by comparing bonds
    reacting_atoms = set()
    for bond in reactant_mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        bond_type = bond.GetBondType()
        
        # Check if the bond exists in the product
        mapped_begin_idx = atom_map.get(begin_idx, None)
        mapped_end_idx = atom_map.get(end_idx, None)
        if mapped_begin_idx is not None and mapped_end_idx is not None:
            product_bond = product_mol.GetBondBetweenAtoms(mapped_begin_idx, mapped_end_idx)
            if product_bond is None or product_bond.GetBondType() != bond_type:
                reacting_atoms.add(begin_idx)
                reacting_atoms.add(end_idx)
        else:
            reacting_atoms.add(begin_idx)
            reacting_atoms.add(end_idx)
    
    # Convert the indices to SMARTS numeration (starting from 1)
    reacting_atoms = [idx + 1 for idx in reacting_atoms]
    
    return reacting_atoms

# Example usage
reaction_smarts = "[c:1]12[cH:2][cH:3][cH:4][cH:5][c:6]1[cH:7][c:8](-[c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[o:15]2>>[c:1]12[cH:2][cH:3][cH:4][cH:5][c:6]1[CH2:7][C@@H:8]([c:9]1[cH:10][cH:11][cH:12][cH:13][cH:14]1)[O:15]2"
reacting_atoms = get_reacting_atoms(reaction_smarts)
print("Reacting atoms indices in reactant:", reacting_atoms)

#Problem: Aromatic bonds turn into single or double bonds when the aromat gets hydrogenated making all of them "reacting"