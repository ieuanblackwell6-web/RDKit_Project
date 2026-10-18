from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, QED, Crippen
from rdkit.Chem.Draw import rdMolDraw2D

mols = []
names = []
weights = []
legends = []
complexities = []
drug_likenesses = []
logps = []

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
        complexity = Descriptors.BertzCT(mol)
        drug_likeness = QED.qed(mol) 
        logp = Crippen.MolLogP(mol)
        legend = f"{name} ({weight:.2f}) \n Bertz Complexity = {complexity:.2f} \n Drug Likeness = {drug_likeness:.2f} \n LogP = {logp:.2f}"
        mols.append(mol)
        names.append(name)
        weights.append(weight)
        complexities.append(complexity)
        drug_likenesses.append(drug_likeness)
        logps.append(logp)
        legends.append(legend)

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(1000, 1000), legends=legends, legendFontSize=60)
img.save("molecule_images.pdf")
