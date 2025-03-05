import numpy as np
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
# from Bio.Phylo.TreeConstruction import DistanceCalculator as calculator
import matplotlib.pyplot as plt
from variables import address_dict, subfolders
from utils import fetch_sequences_from_fasta, write_sequence_to_fasta

models= ['identity', 'blastn', 'trans', 'benner6', 'benner22', 'benner74', 'blosum100', 'blosum30', 'blosum35', 'blosum40', 'blosum45', 'blosum50',
    'blosum55', 'blosum60', 'blosum62', 'blosum65', 'blosum70', 'blosum75', 'blosum80', 'blosum85', 'blosum90', 'blosum95', 'feng', 'fitch', 'genetic',
    'gonnet', 'grant', 'ident', 'johnson', 'levin', 'mclach', 'miyata', 'nwsgappep', 'pam120', 'pam180', 'pam250', 'pam30', 'pam300', 'pam60', 'pam90',
    'rao', 'risler', 'structure']

# Input and output file paths
data_folder = address_dict['ECOHARVEST']
msa_fname = 'lipases_initialMSA_mafft.fasta' # 'CALB_phmmer_uniprot_trembl_incE=1e-03_E=1e-03_mafft.fasta'
tree_model = 'blosum62' # 'identity' #
# Load the MSA from a FASTA file
alignment = AlignIO.read(data_folder + subfolders['msa'] + msa_fname, "fasta")
print(len(alignment))

# Calculate pairwise distances
print('Calculating distances...')
calculator = DistanceCalculator(tree_model)
distance_matrix = calculator.get_distance(alignment)
np.save(data_folder + msa_fname.replace('.fasta', '_distmat.npy'), distance_matrix)
print(distance_matrix)

# Construct the phylogenetic tree using Neighbor-Joining (NJ)
print('Constructing tree...')
constructor = DistanceTreeConstructor()
tree = constructor.nj(distance_matrix)
# constructor = DistanceTreeConstructor(calculator)
# tree = constructor.build_tree(alignment)
print(tree)

# # Visualize the phylogenetic tree
fig = plt.figure(figsize=(6, 6), dpi=100) # create figure & set the size
plt.rc('font', size=12)             # fontsize of the leaf and node labels
plt.rc('xtick', labelsize=10)       # fontsize of the tick labels
plt.rc('ytick', labelsize=10)       # fontsize of the tick labels
#turtle_tree.ladderize()		   # optional way to reformat your tree
axes = fig.add_subplot(1, 1, 1)
Phylo.draw(tree, axes=axes)
fig.savefig(data_folder + subfolders['msa'] + msa_fname.replace('.fasta', '_phylo.png'), bbox_inches='tight')
plt.show()