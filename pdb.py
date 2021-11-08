import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import math

#generates an even sphere distribution of points around the center of mass
def fibonacci_sphere(samples=1):

    points = []
    phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians

    for i in range(samples):
        y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        radius = math.sqrt(1 - y * y)  # radius at y

        theta = phi * i  # golden angle increment

        x = math.cos(theta) * radius
        z = math.sin(theta) * radius

        points.append((com[0]+x, com[1]+y, com[2]+z))

    return points

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

# remove all the atoms with ACC < 0
# df_pdb_acc[3]= df_pdb_acc[3].astype(int)
df_pdb_acc_clean = df_pdb_acc[df_pdb_acc[3] > 0]

# select Ca from pdb
df_pdb_CA = df_pdb_acc_clean[df_pdb_acc_clean['atom_name'] == 'CA']

#compute the center of mass

# create a PDBParser object
parser = PDBParser()
# create a structure object from a PDB file
structure = parser.get_structure('name', '1uaz.pdb')
com = structure.center_of_mass()

#Determination of the lines passing through thecenter of mass

points_com = fibonacci_sphere(10000)
