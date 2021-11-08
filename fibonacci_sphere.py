import math
from Bio.PDB import PDBParser

#compute the center of mass

# create a PDBParser object
parser = PDBParser()
# create a structure object from a PDB file
structure = parser.get_structure('name', '1uaz.pdb')
com = structure.center_of_mass()
#print(com)

def fibonacci_sphere(samples=1):

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

fibonacci_sphere(10000)
