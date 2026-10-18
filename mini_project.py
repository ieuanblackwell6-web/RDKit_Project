from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, QED, Crippen

def parse_text(path: str) -> tuple[list[str], list[Chem.Mol]]:
    names, mols = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            name, smiles = line.split(maxsplit=1)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Could not parse {name}, skipping")
                continue
            names.append(name)
            mols.append(mol)
    return names, mols

def make_legend(name: str, mol: Chem.Mol) -> str:
    weight = Descriptors.ExactMolWt(mol)
    complexity = Descriptors.BertzCT(mol)
    drug_likeness = QED.qed(mol)
    logp = Crippen.MolLogP(mol)
    return (f"{name} ({weight:.2f})\n"
            f"Bertz Complexity = {complexity:.2f}\n"
            f"Drug Likeness = {drug_likeness:.2f}\n"
            f"LogP = {logp:.2f}")

names, mols = parse_text("molecules.txt")
print(f"Parsed {len(mols)} molecules")  # quick check

legends = [make_legend(n, m) for n, m in zip(names, mols)]

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(1000, 1000), legends=legends, legendFontSize=60)
img.save("molecule_images.pdf")
