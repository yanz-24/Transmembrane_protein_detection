import sabrina_pdb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="input file (.pdb)")
parser.add_argument("-id", required=True, help="pdb code")
parser.add_argument("-pt", default=1000, help="number of points in fibonnaci sphere <int> (default=1000)")
args = parser.parse_args()
#use example :  python3 main.py -i 1uaz.pdb -id 1UAZ

if __name__ == "__main__":
    #get the pdb
    pdb_prepared = sabrina_pdb.prepare_pdb( args.i, args.id)

    #compute the com
    com= sabrina_pdb.mass_center(args.i)

    #fibo distribution
    fibo_sphere = sabrina_pdb.fibonacci_sphere(com, args.pt)

    #lists of hydrophobic and hydrophilic residues
    hydrophobic = ['PHE', 'GLY', 'ILE', 'LEU', 'MET', 'TRP', 'VAL', 'TYR']
    hydrophilic = ['ALA', 'CYS', 'ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR']

    #get the hydrophobic factor
    hydorphobic_factor = sabrina_pdb.hydrophobic_factor(pdb_prepared, hydrophobic, hydrophilic)


