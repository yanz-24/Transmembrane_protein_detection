import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# Initialize a PandasPdb object + fetch PDB file
# NB :add a user input to chose the file
ppdb = PandasPdb().fetch_pdb('1uaz')
df_pdb = ppdb.df['ATOM']

# Save DSSP in dataframe
p = PDBParser()
structure = p.get_structure("1UAZ", "1uaz.pdb")
model = structure[0]
dssp = DSSP(model, "1uaz.pdb", dssp='mkdssp')

# keep ACC and residue number
df_dssp = pd.DataFrame(data=dssp)
df_dssp_acc = df_dssp.iloc[:,[0,3]]

# merge pdb and dssp
df_dssp_acc.iloc[:,0]= df_dssp_acc.iloc[:,0].astype(int)
df_pdb['residue_number']= df_pdb['residue_number'].astype(int)
df_pdb_acc = df_dssp_acc.merge(right=df_pdb, how='left', left_on=0, right_on='residue_number')
print(df_pdb_acc)

# remove all the atoms with ACC < 0
# df_pdb_acc[3]= df_pdb_acc[3].astype(int)
df_pdb_acc_clean = df_pdb_acc[df_pdb_acc[3] > 0]
print(df_pdb_acc_clean)