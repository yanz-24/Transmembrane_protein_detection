import argparse
import pandas as pd
from pathlib import Path
import os

from biopandas.pdb import PandasPdb
from Bio.PDB import PDBList
import read_pdb as pdb
import fibonacci_sphere as fibo
import objective_function as obj


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

    # iterate on normal vector and calculate Q-value
    best_vector = fibo_sphere[0]
    best_Qvalue = 0
    for vector in fibo_sphere:
        # add the slice number of each residue
        df_current_vector = fibo.slicing(df, vector, com)
        # check if the residues are straight
        obj.add_is_straight(df_current_vector).to_csv('test.csv')
        
        # calculate the Q-value of each slice of a vector
        df_grouped = df_current_vector.groupby(['slice'])
        list_qvalue_1a = list() # list of Q-value of each slice 
        for name, group in df_grouped:
            qvalue_1a = obj.calculate_Qvalue(group)
            list_qvalue_1a.append(qvalue_1a)
        
        # calculate the Q-value of a 22A window of a vector
        for i in range(len(list_qvalue_1a)-22):
            qvalue_22a = sum(list_qvalue_1a[i:(i+22)])
            if qvalue_22a > best_Qvalue: # if find a better Qvalue
                best_Qvalue = qvalue_22a
                best_vector = vector

        

#print(pdb.prepare_pdb(ppdb))

# view molecule

#prendre le lineplot des carbon alpha en output
#renvoyer la prot aligner sur z en changant le beta factor de la partie dans et a l'ext et colorer en 2 couleurs
