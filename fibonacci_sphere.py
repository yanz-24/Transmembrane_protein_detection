import math
import numpy as np
from Bio.PDB import PDBParser

def mass_center(pdb_file) :
    # create a PDBParser object
    parser = PDBParser()
    # create a structure object from a PDB file
    structure = parser.get_structure('name', pdb_file)
    com = structure.center_of_mass()

    return com

def fibonacci_sphere(com, samples=1):
    '''
    evenly distribute points on a sphere and return all the points in a list

    input:
    - com: center of mass
    - samples: number of points to generate (default=1)
    '''
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

def point2vector(vector, point, center):
    '''
    project a point to a vector in a 3D space and return the reletive
    position of point as a ratio of distance relative to the center of mass

    input:
    - vector: normal vector calculated by fibonacci_sphere
    - point: Coordinates of Calpha residue
    - center: Coordinates of the coordinates origin (center fo mass) 
            because we translate the spatial coordinate system so that 
            the  centre of mass becomes the origin of the coordinates
    '''
    # convert into numpy array
    vector = np.array(vector)
    point = np.array(point)
    center = np.array(center)

    normal_vector = vector - center
    point_vector = point - center

    return (np.dot(normal_vector, point_vector)/np.dot(normal_vector, normal_vector))


def slicing(df, vector, center):
    '''
    cut the protein into 1 Å wide slices along a normal vector
    and add slice number to dataframe

    input:
    - df: pdb dataframe of Calpha
    - vector: a normal vector calculated by fibonacci_sphere
    - center: Coordinates of the coordinates origin (center fo mass) 
            because we translate the spatial coordinate system so that 
            the  centre of mass becomes the origin of the coordinates
    '''
    for index, row in df.iterrows():
        point = np.array([row['x_coord'], row['y_coord'], row['z_coord']])
        projection = point2vector(vector, point, center)
        df.loc[index, ["projection"]] = projection
        df.loc[index, ["slice"]] = int(math.floor(projection))

    return(df)



    








