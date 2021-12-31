import random
import math

def random_point_in_circle(r=5):
    # random angle
    alpha = 2 * math.pi * random.random()
    # calculating coordinates
    x = r * math.cos(alpha)
    y = r * math.sin(alpha)
    return x,y


def draw_mb(cmb, tk, filepath):   
    z_upper = cmb + tk / 2
    z_lower = cmb + tk / 2
    with open(filepath, 'a') as outfile:
        for i in range(5000, 6000):
            x,y = random_point_in_circle()
            line = f'HETATM {i}  O   DUM  {i}      {x:8.3f}{y:8.3f}{z_upper:8.3f}                          \n'
            outfile.write(line)
            line = f'HETATM {i}  N   DUM  {i}      {x:8.3f}{y:8.3f}{z_lower:8.3f}                          \n'
            outfile.write(line)

draw_mb(1, 38, '/mnt/d/M2/python2_project/bqf.pdb')