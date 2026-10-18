from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, QED, Crippen

mols = []
names = []
legends = []

def parse_text():
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
            mols.append(mol)
            names.append(name)

def build_sheet():
    for mol, name in zip(mols, names):
        weight = Descriptors.ExactMolWt(mol)
        complexity = Descriptors.BertzCT(mol)
        drug_likeness = QED.qed(mol)
        logp = Crippen.MolLogP(mol)
        legend = (f"{name} ({weight:.2f})\n"
                  f"Bertz Complexity = {complexity:.2f}\n"
                  f"Drug Likeness = {drug_likeness:.2f}\n"
                  f"LogP = {logp:.2f}")
        legends.append(legend)

parse_text()
build_sheet()

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(1000, 1000), legends=legends, legendFontSize=60)
img.save("molecule_images.pdf")
