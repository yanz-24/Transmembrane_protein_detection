import argparse
import pandas as pd
from pathlib import Path
import os
import numpy as np

from biopandas.pdb import PandasPdb
from Bio.PDB import PDBList

import read_pdb as pdb
import fibonacci_sphere as fibo
import objective_function as obj
import coordinates_transformation as trans


parser = argparse.ArgumentParser()
parser.add_argument("-i", default=False, help="input file (.pdb)")
parser.add_argument("-id", default=False, help="pdb code")
parser.add_argument("-pt", default=10, type=int, help="number of points in fibonnaci sphere <int> (default=10)")
args = parser.parse_args()

# show every row
pd.set_option('display.max_columns', None)
# show every column
pd.set_option('display.max_rows', None)

if __name__ == "__main__":
    # get calpha of pdb in dataframe and calculate the com
    if args.id: # if user provide pdb id (preferred)
        ppdb = PandasPdb().fetch_pdb(args.id)

        pdbl = PDBList()
        pdb_file = pdbl.retrieve_pdb_file(args.id, file_format="pdb")
        # change pdb extension from ent to pdb
        pdb_file = Path(pdb_file)
        pdb_file = pdb_file.rename(pdb_file.with_suffix('.pdb'))

        df = pdb.prepare_pdb(ppdb, pdb_file)
        com = fibo.mass_center(pdb_file)
        os.remove(pdb_file)

    elif args.i: # if user provide input file
        ppdb = PandasPdb().read_pdb(args.i)
        df = pdb.prepare_pdb(ppdb, args.i)
        com = fibo.mass_center(args.i)

    else: # no input: Program termination
        print("No available pdb id or pdb file found.")
        os._exit()

    # fibo distribution
    fibo_sphere = fibo.fibonacci_sphere(com, args.pt)

    best_Qvalue, best_vector = obj.get_best_vector(df, com, fibo_sphere)
    
    # Classification of protein
    '''
    1. Q-value < lower limit: globular protein
    2. lower limit < Q-value < upper limit (Swissprot): the globular fragment 
        of a transmembrane protein (rare case, should be checked manually)
    3. upper limit < Q-value: transmembrane protein (alpha, beta, coil according to the DSSP algorithm)
    
    use the threshold in article? test 10+ pdb and decide arbitrarily? 
    Can we ignore case 2 to make a binary classification?
    

    # threshold in article (Fig. 1.)
    lower_limit = 40
    upper_limit = 50

    if best_Qvalue > upper_limit:
        print("transmembrane protein")
    else:
        print("globular protein")
        os._exit(1)
    '''

    # Membrane positioning
    '''
    method defined in https://doi.org/10.1093/protein/gzv063 

    two variables to be taken into consideration:
    1. iterate on membrane thickness (th) from 2.5 nm to 10 nm
    2. iterate on membrane center (Cmb) from -5 nm to 5 nm around COM
    condition: at least one atom between membrane

    step 1: Transform the coordinates so that the COM becomes the origin (translate)
            and the normal vector becomes the Z-axis (rotate).
    step 2: iterate on th and Cmb (complexity: ~7500)
    step 3: calculate C in each loop 
            (determine whether atom out of mb by z coordinate)
    step 4: draw two plots: a) C vs th; b) C vs Cmb
    '''
    
    z_axis = np.array([0, 0, 1])

    rotation_matrix = trans.get_rotation_matrix(vec1=best_vector, vec2=z_axis)

    print(best_vector, com)

    # view molecule