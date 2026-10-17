import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors 
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

mol = Chem.MolFromSmiles("c1ccccc1")

img = Draw.MolToImage(mol)
img.save("benzene.png")

mol_h = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol_h, randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol_h) 

with open("first_benzene.xyz", "w") as f:
	f.write(Chem.MolToXYZBlock(mol_h))
