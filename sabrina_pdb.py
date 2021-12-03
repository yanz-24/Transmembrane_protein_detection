import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import math


'''
keeping CA atoms and computing ACC
returns a df 
'''

def prepare_pdb( pdb_file, pdb_id ) :

    # Initialize a PandasPdb object + fetch PDB file
    ppdb = PandasPdb().fetch_pdb(pdb_id)
    df_pdb = ppdb.df['ATOM']

    # Save DSSP in dataframe
    p = PDBParser()
    structure = p.get_structure(pdb_id, pdb_file)

    model = structure[0]
    dssp = DSSP(model, pdb_file, dssp='mkdssp')

    # keep ACC and residue number
    df_dssp = pd.DataFrame(data=dssp)
    df_dssp_acc = df_dssp.iloc[:,[0,3]]

    # merge pdb and dssp
    df_dssp_acc.iloc[:,0]= df_dssp_acc.iloc[:,0].astype(int)
    df_pdb['residue_number']= df_pdb['residue_number'].astype(int)
    df_pdb_acc = df_dssp_acc.merge(right=df_pdb, how='left', left_on=0, right_on='residue_number')

    # remove all the atoms with ACC < 0
    df_pdb_acc_clean = df_pdb_acc[df_pdb_acc[3] > 0]

    # select Ca from pdb
    df_pdb_CA = df_pdb_acc_clean[df_pdb_acc_clean['atom_name'] == 'CA']
    df_pdb_CA = df_pdb_CA.rename(columns={3: 'acc'})

    return df_pdb_CA


'''
Compute center of mass
'''

def mass_center(pdb_file) :
    # create a PDBParser object
    parser = PDBParser()
    # create a structure object from a PDB file
    structure = parser.get_structure('name', pdb_file)
    com = structure.center_of_mass()

    return com


'''
generating an even sphere distribution of points 
around the center of mass
'''

def fibonacci_sphere(com, samples=1):

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


'''
The hydrophobic factor of the objective function is
the relative hydrophobic membrane-exposedsurface area
(hydrophobic area divided by all surface area).
'''

def hydrophobic_factor(df, res_hydropho, res_hydrophi ):

    df_hydrophobic = df[df['residue_name'].isin(res_hydropho)]
    df_hydrophilic = df[df['residue_name'].isin(res_hydrophi)]
    sum_acc_hydrophobic = sum(df_hydrophobic['acc'])
    sum_acc_hydrophilic = sum(df_hydrophilic['acc'])
    factor = sum_acc_hydrophobic / (sum_acc_hydrophilic + sum_acc_hydrophobic)

    return factor




#Determination of the lines passing through thecenter of mass
#points_com = fibonacci_sphere(10000)








