import random
import numpy as np
import math

def random_point_in_circle():
    length = np.sqrt(np.random.uniform(0, 400))
    angle = np.pi * np.random.uniform(0, 2)
    x = length * np.cos(angle)
    y = length * np.sin(angle)
    return x,y


def draw_mb(cmb, tk, filepath):   
    z_upper = cmb + tk / 2
    z_lower = cmb - tk / 2
    with open(filepath, 'a') as outfile:
        for i in range(5000, 5500):
            x,y = random_point_in_circle()
            line = f'\nHETATM{i:>5d}  O   DUM A5000    {x:8.3f}{y:8.3f}{z_lower:8.3f}'
            outfile.write(line)
            line = f'\nHETATM{i:>5d}  N   DUM A5000    {x:8.3f}{y:8.3f}{z_upper:8.3f}'
            outfile.write(line)

draw_mb(1, 38, '/mnt/d/M2/python2_project/bad.pdb')