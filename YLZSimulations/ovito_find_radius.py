#!/usr/bin/env python
import os
import sys
import copy
import math
import numpy as np
from scipy.optimize import minimize
import ovito

global_stride = 1

def circle_equation(params, x, y):
    cx, cy, r = params
    return (x - cx)**2 + (y - cy)**2 - r**2

def objective(params, x, y):
    return np.sum(circle_equation(params, x, y)**2)

def fit_circle_to_points(x, y):
    initial_guess = [np.mean(x), np.mean(y), np.std(x + y)]
    result = minimize(objective, initial_guess, args=(x, y))
    cx, cy, r = result.x
    return cx, cy, r

def generate_circle_points(cx, cy, r, num_points=1000):
    theta = np.linspace(0, 2*np.pi, num_points)
    x = cx + r * np.cos(theta)
    y = cy + r * np.sin(theta)
    return x,y


def get_output_file(filepath,suffix="",prefix=""):
    # Input file name
    path, infile = os.path.split(filepath)
    # Getting the outfile including the seed
    removable = [".xyz",".gz",".dat","lmpout_xyz_"]
    root_seeded = infile
    for item in removable:
        root_seeded = root_seeded.replace(item,"")

    bits = root_seeded.split("_")
    seedword = "seed"
    for bit in bits:
        if seedword==bit[:len(seedword)]:
            seed = int(bit.replace(seedword,""))
            bits.remove(bit)
            break
    root = "_".join(bits)

    return os.path.join(path,prefix+root_seeded+suffix), root, seed
     

def find_smallest_radius(filename,stride=1):
    # Given a filename, this function extracts the progress of the radius with time
    # Based on the name, it decides on whether it takes the endpoints of the middle of the simulation box
    pipeline = ovito.io.import_file(filename, columns = ["ID", "Particle Type", "Position.X", "Position.Y", "Position.Z", "shapex", "shapey", "shapez","quatw","quati","quatj","quatk","Transparency"])
    numFrames = pipeline.source.num_frames
    radii = []

    def rangefun(x,halfbox):
        # Taking data from the middle of the box (arbitrary, can be changed)
        return abs(x)<2

    for frame in range(0,numFrames,stride):
        data = pipeline.compute(frame)
        halfbox = np.array(data.cell)[0,0]/2
        print("frame:",frame, " out of ", numFrames)
        timestep = data.attributes['Timestep']
        print("Timestep: ", timestep)
        selection_y = []
        selection_z = []
        for part in range(data.particles.count):
            if rangefun(data.particles.positions[part,0],halfbox):
                selection_y.append(data.particles.positions[part,1])
                selection_z.append(data.particles.positions[part,2])
        if len(selection_y)>20:
            selection_y = np.array(selection_y)
            selection_z = np.array(selection_z)
            # Fitting a circle
            cy, cz, r = fit_circle_to_points(selection_y, selection_z)
            radii.append([timestep,abs(r)])
        else:
            radii.append([timestep,0])
    radii = np.array(radii)
    return radii

if __name__=="__main__":
    arguments = sys.argv
    for filename in arguments[1:]:
        print("Processing file:", filename)
        outfile,root,seed = get_output_file(filename,".dat","radii_")
        radii = find_smallest_radius(filename,stride=global_stride)
        np.savetxt(outfile,radii)
