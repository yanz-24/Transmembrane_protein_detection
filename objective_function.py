import numpy as np
import fibonacci_sphere as fibo


def hydrophobic_factor(df):
    '''
    Calculate the hydrophobic factor of a slice of 1A.

    The hydrophobic factor of the objective function is
    the relative hydrophobic membrane-exposedsurface area
    (hydrophobic area divided by all surface area).
    '''
    hydrophobic = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'TRP', 'VAL', 'TYR']
    hydrophilic = ['ALA', 'CYS', 'ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR']

    # create separate df for hydrophobic and hydrophilic residues
    df_hydrophobic = df[df['residue_name'].isin(hydrophobic)]
    df_hydrophilic = df[df['residue_name'].isin(hydrophilic)]
    sum_acc_hydrophobic = sum(df_hydrophobic['acc'])
    sum_acc_hydrophilic = sum(df_hydrophilic['acc'])
    factor = sum_acc_hydrophobic / (sum_acc_hydrophilic + sum_acc_hydrophobic)

    return factor


def monotonic(x):
    """
    check list monotonicity and return a boolean value
    """
    dx = np.diff(x)
    return np.all(dx <= 0) or np.all(dx >= 0)


def add_is_straight(df):
    """
    Check if the residues are straight and add this information in dataframe

    "The i-th residue in a protein chain is part of a straight triplet if the projection of
    the Cα atoms of the previous third (i − 3) residue, itself (i)
    and the next third (i +3) residue onto a predefined vector (the
    normal vector of membrane planes, see below) are in a monotone decreasing
    or increasing order."

    input:
    - df: pdb dataframe of Calpha
    """
    df['straight'] = False

    for i, row in df.iterrows():
        if_straight = False

        if i > 3 and i < (len(df) - 4):
            # the projection of triplets around i in form of list
            projection_list = list(df.projection)[(i - 3):(i + 3)]
            if_straight = monotonic(projection_list)

        df.at[i, 'straight'] = if_straight
    return (df)


def straightness_factor(df_slice):
    '''
    calculate the straightness factor and return
    "The straightness factor is defined as 
    the relative frequency of ‘straight’ residues in a given protein slice."

    input:
    - df_slice: pdb dataframe of Calpha with slice column (result returned by add_is_straight)
    '''
    # number of straight residues (column 'straight' is True)
    number_straight = df_slice['straight'].values.sum()
    return number_straight / len(df_slice)


def turn_factor(df_slice):
    """
    The turn factor is defined as one minus the relative frequency of ‘turn’ residues in a given slice.
    Turn triplets have a similar definition as the straight triplets: the i-th residue in a protein
    chain is the center of a turn if the projection of the Cα atoms of the previous third (i − 3) residue,
    itself (i) and the next third (i + 3) residue onto the predefined vector
    are not in a  monotone decreasing or increasing order.
    
    - df_slice: pdb dataframe of Calpha with slice column (result returned by add_is_straight)

    """

    number_turn = (~df_slice['straight'].values.sum())
    return 1 - (number_turn / len(df_slice))

'''
def end_chain(df_slice):
    """
    The end-chain factor is one minus the relative frequency of
    chain end residues in a given slice
    """
    # 1 - freq(ter res in slice )
'''

def calculate_Qvalue(df_slice):
    """
    calculate the Q-value of a PDB dataframe 
    - df_slice: pdb dataframe of Calpha with slice column (result returned by add_is_straight)

    """
    qvalue = hydrophobic_factor(df_slice) * straightness_factor(df_slice) * turn_factor(df_slice)
    return qvalue


def get_best_vector(df, com, fibo_sphere):
    """
    return the best Q value and the best normal vector of a protein

    input:
    - df: pdb dataframe
    - com: center of mass
    - fibo_sphere: a list of 3D vectors returned by fibo.fibonacci_sphere function
    """
    best_vector = fibo_sphere[0]
    best_Qvalue = 0
    for vector in fibo_sphere:
        # add the slice number of each residue
        df_current_vector = fibo.slicing(df, vector, com)
        # check if the residues are straight
        add_is_straight(df_current_vector)
        
        # calculate the Q-value of each slice of a vector
        df_grouped = df_current_vector.groupby(['slice'])
        list_qvalue_1a = list() # list of Q-value of each slice 
        for name, group in df_grouped:
            qvalue_1a = calculate_Qvalue(group)
            list_qvalue_1a.append(qvalue_1a)
        
        # calculate the Q-value of a 22A window of a vector
        for i in range(len(list_qvalue_1a)-22):
            qvalue_22a = sum(list_qvalue_1a[i:(i+22)])
            if qvalue_22a > best_Qvalue: # if find a better Qvalue
                best_Qvalue = qvalue_22a
                best_vector = vector
    return [best_Qvalue, best_vector]
