import argparse
import pandas as pd
import math
from pathlib import Path
import os

from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB import PDBList 
from Bio.PDB.DSSP import DSSP


import read_pdb as pdb
import fibonacci_sphere as fibo
import objective_function as obj


parser = argparse.ArgumentParser()
parser.add_argument("-i", default=False, help="input file (.pdb)")
parser.add_argument("-id", default=False, help="pdb code")
parser.add_argument("-pt", default=1000, help="number of points in fibonnaci sphere <int> (default=1000)")
args = parser.parse_args()
#use example :  python3 main.py -i 1uaz.pdb -id 1UAZ

# show every row
pd.set_option('display.max_columns', None)
# show every column
pd.set_option('display.max_rows', None)

if __name__ == "__main__":
    # get calpha of pdb in dataframe and calculate the com
    if args.id != False: # if user provide pdb id (preferred)
        ppdb = PandasPdb().fetch_pdb(args.id)

        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(args.id, file_format="pdb")
        # change pdb extension from ent to pdb
        pdb_file = Path(pdb_file)
        pdb_file = pdb_file.rename(pdb_file.with_suffix('.pdb'))

        df = pdb.prepare_pdb(ppdb, pdb_file)
        com = fibo.mass_center(pdb_file)
        os.remove(pdb_file)

    elif args.i != False: # if user provide input file
        ppdb = PandasPdb().read_pdb(args.i)
        df = pdb.prepare_pdb(ppdb, args.i)
        com = fibo.mass_center(args.i)

    else: # no input: Program termination
        print("No available pdb id or pdb file found.")
        os._exit()
    
    # fibo distribution
    fibo_sphere = fibo.fibonacci_sphere(com, args.pt)

    # iterate on normal vector and calculate Q-value
    for vector in fibo_sphere:
        # add the slice number of each residue
        df_current_vector = fibo.slicing(df, vector, com)
        # check if the residues are straight
        obj.add_is_straight(df_current_vector).to_csv('test.csv')
        
        # calculate the Q-value of each slice

        df_grouped = df_current_vector.groupby(['slice'])
        for name, group in df_grouped:
            # calculate the hydrophobic factor
            hydrophobic_factor = obj.hydrophobic_factor(group)
            #print(hydrophobic_factor)

            # calculate the straightness factor
            straightness_factor = obj.straightness_factor(group)
            print(straightness_factor)
        
        


