import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

def prepare_pdb(ppdb, pdb_file=False) :

    # Initialize a PandasPdb object + fetch PDB file
    print(ppdb.df.head())
    df_pdb = ppdb.df['ATOM']
    df_pdb = df_pdb[df_pdb['chain_id'] == 'A']

    # Save DSSP in dataframe
    p = PDBParser()
    structure = p.get_structure('name', pdb_file)

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
    df_pdb_CA = df_pdb_CA.rename(columns={3: 'acc'}).reset_index(drop=True)

    return df_pdb_CA
    
#package free sassa pour la relative acc
