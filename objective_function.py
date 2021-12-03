hydrophobic = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'TRP', 'VAL', 'TYR']
hydrophilic = ['ALA', 'CYS', 'ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR']

def hydrophobic_factor(df, res_hydropho, res_hydrophi ):
    '''
    The hydrophobic factor of the objective func-tion is
    the relative hydrophobic membrane-exposedsurface area
    (hydrophobic area divided by all surface area).
    '''
    #create separate df for hydrophobic and hydrophilic residues
    df_hydrophobic = df[df['residue_name'].isin(res_hydropho)]
    df_hydrophilic = df[df['residue_name'].isin(res_hydrophi)]
    sum_acc_hydrophobic = sum(df_hydrophobic['acc'])
    sum_acc_hydrophilic = sum(df_hydrophilic['acc'])
    factor = sum_acc_hydrophobic / (sum_acc_hydrophilic + sum_acc_hydrophobic)

    return factor

#example
#print(hydrophobic_factor(df_pdb_CA, hydrophobic, hydrophilic))
