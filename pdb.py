import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

# writes the content of a list in a file according to an opening mode
def list_to_file(file_name, list_name, opening_mode):
    f = open(file_name, opening_mode)
    for element in list_name:
        f.write(str(element))
    f.close()


# Initialize a PandasPdb object + fetch PDB file
# NB :add a user input to chose the file
ppdb = PandasPdb().fetch_pdb('1uaz')

# DSSP data is accessed by a tuple (chain_id, res_id)
p = PDBParser()
structure = p.get_structure("1UAZ", "1uaz.pdb")
model = structure[0]
dssp = DSSP(model, "1uaz.pdb", dssp='mkdssp')

atoms_out_list = []
for i in range(len(ppdb.df['ATOM'].values)):
    #to avoid list index out of range
    if ppdb.df['ATOM'].values[i][8] < len(list(dssp.keys())):
        a_key = list(dssp.keys())[ppdb.df['ATOM'].values[i][8] - 1]
        if dssp[a_key][3] > 0.0:
            atoms_out_list.append(ppdb.df['ATOM'].values[i])


df = pd.DataFrame(data= atoms_out_list)
df = df.drop(df.columns[[19,20]], axis=1)  #remove the last two colums of the df

df_dssp = pd.DataFrame(data=dssp)
df_dssp1 = df_dssp.iloc[:,[0,3]]

#df_dssp1['0']= df_dssp1['0'].astype(int)
#df['8']= df['8'].astype(int)

#my_df = df_dssp1.merge(right= df, how= 'left', left_on= 0, right_on=5)

#print(my_df)
print(df_dssp1)
