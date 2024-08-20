from rdkit import Chem

def filter_reactants_with_template(reactant_smiles_list, template_smarts):
    template = Chem.MolFromSmarts(template_smarts)
    if template is None:
        raise ValueError(f"Invalid template SMARTS: {template_smarts}")

    matching_reactants = []
    for reactant_smiles in reactant_smiles_list:
        reactant = Chem.MolFromSmiles(reactant_smiles)
        if reactant is not None:
            if reactant.HasSubstructMatch(template):
                matching_reactants.append(reactant_smiles)
        else:
            print(f"Invalid SMILES: {reactant_smiles}")

    return matching_reactants

# Example usage
template_smarts = '[c:1]1:[c:10]:[n:9]:[c:8]2:[c:7]:[c:6]:[c:5]:[c:4]:[c:3]:2:[n:2]:1'

reactant_smiles_list = ['CC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'FC(C1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1)(F)F', 'CCC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'CCCC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'CCCCC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'CCCCCCCCCCC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'CC(C)CC1=CC=C(N=C(C2=CC=CC=C2)C(C3=CC=CC=C3)=N4)C4=C1', 'C1(C2=CC=CC=C2)=CC=C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N5)C5=C1', 'C1(CC2=CC=CC=C2)=CC=C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N5)C5=C1', 'C1(/C=C/C2=CC=CC=C2)=CC=C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N5)C5=C1', 'C1(CCC2=CC=CC=C2)=CC=C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N5)C5=C1', 'CC([Si](C1=CC=CC=C1)(OCC2=CC=C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N5)C5=C2)C6=CC=CC=C6)(C)C', 'CC1=C2C(N=C(C3=CC=CC=C3)C(C4=CC=CC=C4)=N2)=CC=C1', 'C12=CC=CC=C1C=C(C3=CC=CC=C3)O2', 'CC(C=CC=C1)=C1C2=CC3=CC=CC=C3O2', 'CC1=CC=CC(C2=CC3=CC=CC=C3O2)=C1', 'CC(C=C1)=CC=C1C2=CC3=CC=CC=C3O2', 'COC(C=C1)=CC=C1C2=CC3=CC=CC=C3O2', 'FC(C(C=C1)=CC=C1C2=CC3=CC=CC=C3O2)(F)F', 'FC(C=C1)=CC=C1C2=CC3=CC=CC=C3O2', 'CC1=CC2=CC=CC=C2O1', 'CCCCC1=CC2=CC=CC=C2O1', 'CCCCCCCCCCC1=CC2=CC=CC=C2O1', 'CC(C)C1=CC2=CC=CC=C2O1', 'CC(C)(C)C1=CC2=CC=CC=C2O1', 'C12=CC=CC=C1C=C(CC3=CC=CC=C3)O2', 'C12=CC=CC=C1C=C(C3=CC=CC=C3)O2', 'CC1=COC2=CC=CC=C21', 'C12=CC=CC=C1C=C(C3=NC=CC=C3)O2', 'C12=CC=CC=C1C=C(C3=CN=CC=C3)O2', 'CC1=CC=C(C2=CC=C(F)C=C2)O1', 'CC1=CC=C(C2=CC=C(OC)C=C2)O1', 'CCCCC1=CC=C(C2=CC=C(OC)C=C2)O1', 'CC1=CC=C(C2=CC(OC)=CC(OC)=C2)O1', 'CC1=CC=C(C2=CC(OC)=C(OC)C(OC)=C2)O1', 'CC1=CC=C(C2=CC=C(N(C)C)C=C2)O1', 'CC1=CC=C(C2=CC=CC(C)=C2)O1', 'CC1=CC=C(C2=CC=C(C)C=C2)O1', 'CC1=CC=C(C2=CC=CC(F)=C2)O1', 'CC1=CC=C(C2=CC=C3C(C=CC=C3)=C2)O1', 'CC1=CC=C(C2=CC=CC=C2)O1', 'CC1=CC=C(C2=CC=C(C(F)(F)F)C=C2)O1', 'CC1=CC=C(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)O1', 'CC1=CC(C2=CC=C(F)C=C2)=CO1', 'CC1=CC(C2=C(F)C=CC=C2)=CO1', 'CC1=CC(C2=CC=C(C)C=C2)=CO1', 'CC1=CC(C2=CC=CC=C2)=CO1', 'CC1=CC(C2=CC=C(OC)C=C2)=CO1', 'CC1=CC(C2=CC3=C(C=CC=C3)C=C2)=CO1', 'C12=CC=CC=C1C=CS2', 'CC1=CC2=CC=CC=C2S1', 'CCC1=CC2=CC=CC=C2S1', 'CCCC1=CC2=CC=CC=C2S1', 'CCCCC1=CC2=CC=CC=C2S1', 'CCCCCCCCCC(C1=CC2=CC=CC=C2S1)=O', 'CC(C)CC1=CC2=CC=CC=C2S1', 'C12=CC=CC=C1C=C(CC3=CC=CC=C3)S2', 'CC1=CSC2=CC=CC=C21', 'C1=CC=CS1', 'CCC1=CC=CS1', 'CC1=CSC=C1', 'FC(C=C1)=CC=C1C2=CC=CS2', 'CC1=CC=C(C2=CC=CC=C2)S1', 'FC(C=C1)=CC=C1C2=CC=C(C)S2', 'CC1=CC=C(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)S1', 'FC(C1=CC(C2=CC=C(CC)S2)=CC(C(F)(F)F)=C1)(F)F', 'FC(C=C1)=CC=C1C2=CSC=C2C', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=CC=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=CC(C)=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=C(C)C=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC(C)=CC=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=C(C)C=CC=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=C(OC)C=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=C([Si](C)(C)C)C=C21', 'CC1=CN(C(OC(C)(C)C)=O)C2=CC=C(C3CCCCC3)C=C21', 'CC(C)C1=CN(C(OC(C)(C)C)=O)C2=CC=CC=C21', 'CC(OC(N1C2=CC=CC=C2C(C3CCCC3)C1)=O)(C)C', 'CCCCCCC1=CN(C(OC(C)(C)C)=O)C2=CC=CC=C21', 'CCCCCCCCCCC1=CN(C(OC(C)(C)C)=O)C2=CC=CC=C21', 'CC1=CN(C(OC)=O)C2=CC=CC=C21', 'CC1=CC2=CC=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=C(C)C=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC(C)=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC=C(C)C=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC=CC(C)=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC(OC)=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC(C(C)(C)C)=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC(CCCC)=CC=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC(C)=CC(C)=C2N1C(OC(C)(C)C)=O', 'CC1=CC2=CC=CC=C2N1C(OC)=O', 'CC(OC(N1C2=CC=CC=C2C=C1CO[Si](C)(C)C(C)(C)C)=O)(C)C', 'CC1=CC2=CC=CC(C)=C2O1', 'CC1=CC2=CC=CC(OC)=C2O1', 'CC1=CC2=CC=CC(C(C)(C)C)=C2O1', 'CC1=CC2=CC(B3OC(C)(C)C(C)(C)O3)=CC=C2O1', 'CC1=CC2=CC(OC)=CC=C2O1', 'CC1=CC2=CC(C(C)(C)C)=CC=C2O1', 'CC1=CC2=CC(C(C)C)=CC=C2O1', 'CCCCC1=CC2=CC=CC=C2O1', '[H]C1=CC=C(OC(C)=C2)C2=C1', 'CC1=CC=C(OC(C)=C2)C2=C1', 'CC1=CC2=CC(CC)=CC=C2O1', 'CC1=CC2=CC(C(C)C)=CC=C2O1', 'CC1=CC2=CC(C(C)(C)C)=CC=C2O1', 'CC1=CC2=CC=CC(C)=C2O1', 'CC1=CC2=CC=CC(C(C)(C)C)=C2O1', 'CC1=CC2=CC=C(C)C(C)=C2O1', 'CC1=CC2=CC(C3=CC=CC=C3)=CC=C2O1', 'CC1=CC2=CC(OC)=CC=C2O1', 'CC1=CC2=CC=CC(OC)=C2O1', 'CC1=CC2=CC(C(F)(F)F)=CC=C2O1', 'CC1=CC2=CC(C(OC)=O)=CC=C2O1', 'CC1=CC2=CC(N3CCOCC3)=CC=C2O1', 'CC1=CC2=CC(N(CC)CC)=CC=C2O1', 'CC1=CC2=CC(NC(C)=O)=CC=C2O1', 'CC1=CC2=CC(N)=CC=C2O1', 'CC1=CC2=CC(B3OC(C)(C)C(C)(C)O3)=CC=C2O1', 'C12=CC=CC=C1C=C(C3=CC=CC=C3)O2', 'CCC1=CC2=CC=CC=C2O1', 'CCCCC1=CC2=CC=CC=C2O1', 'CC1=COC2=CC=CC=C21', 'CC1=CC2=CC=CC(F)=C2O1', 'CC1=CC2=CC(F)=CC=C2O1', 'CC1=CC2=CC=C(F)C(F)=C2O1', 'CC1=CC2=CC(F)=CC(F)=C2O1', 'CC1=CC2=C(F)C=CC(F)=C2O1', 'CC1=CC2=CC(F)=C(F)C(F)=C2O1', 'CC1=CC2=C(F)C(F)=C(F)C=C2O1']

matching_reactants = filter_reactants_with_template(reactant_smiles_list, template_smarts)
print("Reactants containing the template:")
for reactant in matching_reactants:
    print(reactant)



#Möglichkeit für breitere Suche Tanimoto Similarities 