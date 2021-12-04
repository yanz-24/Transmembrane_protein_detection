import pandas as pd
from biopandas.pdb import PandasPdb
from Bio.PDB import PDBParser
from Bio.PDB import PDBList 
from Bio.PDB.DSSP import DSSP
import math
from pathlib import Path


def prepare_pdb(pdb_file=False, pdb_id=False) :

    if pdb_file != False: # Custom files in .pdb format are preferred
        # Initialize a PandasPdb object + fetch PDB file
        ppdb = PandasPdb().read_pdb(pdb_file)
        df_pdb = ppdb.df['ATOM']

        # create a structure object from a PDB file?
        p = PDBParser() 
        structure = p.get_structure(pdb_id, pdb_file)
        model = structure[0]
        dssp = DSSP(model, pdb_file, dssp='mkdssp')


    elif pdb_id != False:
        ppdb = PandasPdb().fetch_pdb(pdb_id)
        df_pdb = ppdb.df['ATOM']

        # create a structure object from a PDB file?
        pdbl = PDBList()
        filename = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb")

        # change pdb extension from ent to pdb
        filename = Path(filename)
        filename = filename.rename(filename.with_suffix('.pdb'))

        p = PDBParser() 
        structure = p.get_structure(pdb_id, filename)
        model = structure[0]
        dssp = DSSP(model, filename, dssp='mkdssp')


    else:
        return("No available pdb id or pdb file found.")
        

    # Save DSSP in dataframe
    df_dssp = pd.DataFrame(data=dssp)

    # keep ACC and residue number
    df_dssp_acc = df_dssp.iloc[:,[0,3]]

    # merge pdb and dssp
    df_dssp_acc.iloc[:,0] = df_dssp_acc.iloc[:,0].astype(int)
    df_pdb['residue_number'] = df_pdb['residue_number'].astype(int)
    df_pdb_acc = df_dssp_acc.merge(right=df_pdb, how='left', left_on=0, right_on='residue_number')

    # remove all the atoms with ACC < 0
    df_pdb_acc_clean = df_pdb_acc[df_pdb_acc[3] > 0]

    # select Ca from pdb
    df_pdb_CA = df_pdb_acc_clean[df_pdb_acc_clean['atom_name'] == 'CA']
    df_pdb_CA = df_pdb_CA.rename(columns={3: 'acc'}).reset_index(drop=True)

    return df_pdb_CA
    
