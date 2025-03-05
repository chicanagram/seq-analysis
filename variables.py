address_dict = {
    'examples': '../PIPS/PIPS-tools/examples/',
    'PIPS': '../PIPS/PIPS-GalOx-data/',
    'PIPS2': '../PIPS/PIPS2-UPOs-data/',
    'influenza-resistance': '../PIPS/influenza-resistance/',
    'ProtSolM': '../solubility-data/ProtSolM/',
    'SoluProtMut': '../solubility-data/SoluProtMut/',
    'SoluProtMut_PKS': '../solubility-data/SoluProtMut/',
    'SoluProtMut_LGK': '../solubility-data/SoluProtMut/',
    'SoluProtMut_bLac': '../solubility-data/SoluProtMut/',
    'SoluProtMut_PKS-bLac': '../solubility-data/SoluProtMut/',
    'SoluProtMut_PKS-bLac-LGK': '../solubility-data/SoluProtMut/',
    'PON-Sol2': '../solubility-data/PON-Sol2/',
    'ECOHARVEST': '../ECOHARVEST/',
    'pips-insilico': '../pips-insilico/data/',
    'databases': '../seq-db/'
}

subfolders = {
    'sequences': 'sequences/',
    'msa': 'msa/',
    'blast': 'blast/',
    'hmm': 'hmm/',
    'conservation_analysis': 'conservation_analysis/',
    'aggregation': 'aggregation/',
    'stability': 'stability/',
    'ml_prediction': 'ml_prediction/',
    'yasara': 'yasara/',
    'pdb': 'pdb/',
    'sce': 'sce/',
    'patents': 'patents/',
    'seqsearch': 'seqsearch/',
    'protein_embeddings': 'protein_embeddings/',
    'expdata': 'expdata/',
    'solubility': 'solubility/',
    'camsol': 'camsol/',
    'deepsolue': 'deepsolue/',
    'ponsol2': 'ponsol2/',
    'netsolp': 'netsolp/',
    'ohe': 'ohe/',
}


mapping = {
    'A': 'Ala',
    'H': 'His',
    'Y': 'Tyr',
    'R': 'Arg',
    'T': 'Thr',
    'K': 'Lys',
    'M': 'Met',
    'D': 'Asp',
    'N': 'Asn',
    'C': 'Cys',
    'Q': 'Gln',
    'E': 'Glu',
    'G': 'Gly',
    'I': 'Ile',
    'L': 'Leu',
    'F': 'Phe',
    'P': 'Pro',
    'S': 'Ser',
    'W': 'Trp',
    'V': 'Val'
    }

aaList = list("ACDEFGHIKLMNPQRSTVWY")
aaList_with_X = list("ACDEFGHIKLMNPQRSTVWYX")
# aaList = ['A','H','Y','R','T','K','M','D','N','C','Q','E','G','I','L','F','P','S','W','V']
aa2idx = {aa: i for i, aa in enumerate(aaList)}

element_mapping = {
    1: "H",  # Hydrogen
    6: "C",  # Carbon
    7: "N",  # Nitrogen
    8: "O",  # Oxygen
    15: "P", # Phosphorus
    16: "S", # Sulfur
    11: "Na",# Sodium
    12: "Mg",# Magnesium
    19: "K", # Potassium
    20: "Ca",# Calcium
    17: "Cl",# Chlorine
    9: "F",  # Fluorine
    35: "Br",# Bromine
    53: "I"  # Iodine
}