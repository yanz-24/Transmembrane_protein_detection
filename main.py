import argparse
from numpy.core.defchararray import lower
import pandas as pd
from pathlib import Path
import os
import numpy as np
# import matplotlib.pyplot as plt

from biopandas.pdb import PandasPdb
from Bio.PDB import PDBList

import read_pdb as pdb
import fibonacci_sphere as fibo
import objective_function as obj
import coordinates_transformation as trans
import draw_mb 

parser = argparse.ArgumentParser()
parser.add_argument("-i", default=False, help="input file (.pdb)")
parser.add_argument("-id", default=False, help="pdb code")
parser.add_argument("-o", help="output file (.pdb)")
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
        os._exit(0)

    # fibo distribution
    fibo_sphere = fibo.fibonacci_sphere(com, args.pt)

    # determine the best normal vector of membrane
    best_Qvalue, best_vector = obj.get_best_vector(df, com, fibo_sphere)

    # Classification of protein
    '''
    1. Q-value < lower limit: globular protein
    2. lower limit < Q-value < upper limit (Swissprot): the globular fragment 
        of a transmembrane protein (rare case, should be checked manually)
    3. upper limit < Q-value: transmembrane protein (alpha, beta, coil according to the DSSP algorithm)
    
    ignore case 2 to make a binary classification
    '''

    # threshold in article (Fig. 1.)
    # lower_limit = 40
    upper_limit = 35

    if best_Qvalue > upper_limit:
        print("This is a transmembrane protein.")
    else:
        print("This is a globular protein.")
        os._exit(1)
    

    # Membrane positioning
    '''
    method defined in https://doi.org/10.1093/protein/gzv063 

    two variables to be taken into consideration:
    1. iterate on membrane thickness (tk) from 2.5 nm to 10 nm
    2. iterate on membrane center (cmb) from -5 nm to 5 nm around COM
    condition: at least one atom between membrane
    '''
    # step 1: Transform the coordinates so that the COM becomes the origin (translate)
    df = trans.translate_df(df, com)
    # df.drop(df.columns[-3:], axis=1, inplace=True) # Remove the last three columns 
    # step 1: and the normal vector becomes the Z-axis (rotate).
    z_axis = np.array([0, 0, 1])
    rotation_matrix = trans.get_rotation_matrix(vec1=(best_vector-com), vec2=z_axis)
    df = trans.rotate_df(df, rotation_matrix)

    # step 2: iterate on tk and cmb (complexity: ~7500)
    M_residue = ['PHE', 'MET', 'GLY', 'ILE', 'LEU', 'TRP', 'VAL', 'CYS', 'SER', 'ALA', 'HIS']
    S_residue = ['ASP', 'GLU', 'LYS', 'ASN', 'PRO', 'GLN', 'ARG', 'THR', 'TYR']
    
    best_C = 0
    for tk in range(25, 100):
        for cmb in range(-50, 50):
            upper_mb = cmb + tk/2
            lower_mb = cmb - tk/2

            if df['z_coord'].min() > upper_mb or df['z_coord'].max() < upper_mb:
                continue # jump if mb does not include protein at all

            else:
                # inside (i) or outside (e) mb:
                ei = pd.cut(df['z_coord'], bins = [lower_mb-100, lower_mb, upper_mb, upper_mb+100], 
                    labels = ["e", "i", "e"], ordered=False)
                # if residue is M (bool):
                if_M = df['residue_name'].apply(lambda x: any([k in x for k in M_residue]))
                # cbind ei and if_M
                df_MSie  = pd.concat([ei.reset_index(drop=True), if_M.reset_index(drop=True)], axis=1)

                # count occurrences
                df_confusion = pd.crosstab(ei.reset_index(drop=True), if_M.reset_index(drop=True))
                Me = df_confusion.loc['e', True]
                Se = df_confusion.loc['e', False]
                Mi = df_confusion.loc['i', True]
                Si = df_confusion.loc['i', False]
               
                # step 3: calculate C in each loop (determine whether atom out of mb by z coordinate)
                C_value = (Mi*Se - Si*Me)/((Mi+Si)*(Mi+Me)*(Si+Se)*(Se+Me))**0.5

                if C_value > best_C: # if find a higher C
                    best_C = C_value
                    best_tk = tk
                    best_cmb = cmb
    
    
    # output: pdb file
    if args.o:
        df_A = ppdb.df['ATOM']
        df_A = df_A[df_A['chain_id'] == 'A']
        df_A = trans.translate_df(df_A, com)
        df_A = trans.rotate_df(df_A, rotation_matrix)

        ppdb.df['ATOM'] = df_A
        ppdb.to_pdb(path=args.o, 
                records=None, 
                gz=False, 
                append_newline=True)
        
        draw_mb.draw_mb(best_cmb, best_tk, args.o)
