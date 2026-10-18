from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

mols = []
names = []
weights = []
legends = []


with open("molecules.txt") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        name, smiles = line.split(maxsplit=1)
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Could not parse {name}, skipping")
            continue
        weight = Descriptors.ExactMolWt(mol)
        legend = f"{name} ({weight:.2f})"
        mols.append(mol)
        names.append(name)
        weights.append(weight)
        legends.append(legend)

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(1000, 1000), legends=legends, legendFontSize=60)
img.save("molecule_images.pdf")
