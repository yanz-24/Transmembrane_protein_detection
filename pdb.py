
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# Initialize a PandasPdb object + fetch PDB file
#NB :add a user input to chose the file
ppdb = PandasPdb().fetch_pdb('1uaz')
#print(ppdb.df['ATOM'].values[1][1])


#print(ppdb.df['ATOM'][ppdb.df['ATOM']['chain_id']])

p = PDBParser()
structure = p.get_structure("1UAZ", "1uaz.pdb")
model = structure[0]
dssp = DSSP(model, "1uaz.pdb" ,dssp='mkdssp')
# DSSP data is accessed by a tuple (chain_id, res_id)
a_key = list(dssp.keys())[2]
print(dssp[a_key])
