from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, QED, Crippen, AllChem

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

def mmff_energy(mol: Chem.Mol) -> float | None:
    """Return the MMFF94 energy in kcal/mol, or None if it can't be calculated."""
    mol_h = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol_h, randomSeed=42) < 0:
        return None
    if not AllChem.MMFFHasAllMoleculeParams(mol_h):
        return None
    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=2000)
    props = AllChem.MMFFGetMoleculeProperties(mol_h)
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
    return ff.CalcEnergy()

def make_legend(name: str, mol: Chem.Mol, energy: float | None) -> str:
    weight = Descriptors.ExactMolWt(mol)
    complexity = Descriptors.BertzCT(mol)
    drug_likeness = QED.qed(mol)
    logp = Crippen.MolLogP(mol)
    energy_text = f"{energy:.2f} kcal/mol" if energy is not None else "n/a"
    return (f"{name} ({weight:.2f})\n"
            f"Bertz Complexity = {complexity:.2f}\n"
            f"Drug Likeness = {drug_likeness:.2f}\n"
            f"LogP = {logp:.2f}\n"
            f"MMFF energy = {energy_text}")

names, mols = parse_text("molecules.txt")
print(f"Parsed {len(mols)} molecules")  # quick check

energies = [mmff_energy(m) for m in mols]
legends = [make_legend(n, m, e) for n, m, e in zip(names, mols, energies)]

img = Draw.MolsToGridImage(mols, molsPerRow=2, subImgSize=(600, 450), legends=legends, legendFontSize=50)
img.save("molecule_images.pdf")
