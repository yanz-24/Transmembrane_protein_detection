'''
calculate the center of mass (COM) of a protein
'''
from Bio.PDB import PDBParser

# create a PDBParser object
parser = PDBParser()  
# create a structure object from a PDB file
structure = parser.get_structure('name', '../projet_python/1uaz.pdb')  
com = structure.center_of_mass()
print(com)

# For each chain
'''
for chain in structure.get_chains():
    print(chain.center_of_mass())
'''
# https://www.biostars.org/p/416228/