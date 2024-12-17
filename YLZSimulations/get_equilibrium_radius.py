#!/usr/bin/env python
import os
import sys
import copy
import math
import numpy as np
from scipy.optimize import minimize

def get_output_file(filepath,suffix="",prefix=""):
    # Input file name
    path, infile = os.path.split(filepath)
    # Getting the outfile including the seed
    removable = [".xyz",".gz","radii_",".dat","lmpout_xyz_"]
    root_seeded = infile
    for item in removable:
        root_seeded = root_seeded.replace(item,"")

    bits = root_seeded.split("_")
    seedword = "P"
    seedword2 = "mu"
    for bit in bits[:]:
        if seedword2==bit[:len(seedword2)]:
            seed2 = bit.replace(seedword2,"")
            bits.remove(bit)
        if seedword==bit[:len(seedword)]:
            seed = bit.replace(seedword,"")
            bits.remove(bit)
            # This always comes last
            break
    root = "_".join(bits)

    return os.path.join(path,prefix+root_seeded+suffix), root, seed, seed2
     

if __name__=="__main__":
    arguments = sys.argv
    mydict = {}
    # have a dict:
    # keys are "seeds"=list, "values"=numpy array, "stdev"=numpy array + "samples" = int
    for filename in arguments[1:]:
        outfile,root,P,mu = get_output_file(filename,".dat","radii_")
        if not (mu,P) in mydict:
            mydict[(mu,P)] = [filename]
        else:
            mydict[(mu,P)].append(filename)
    for key,files in mydict.items():
        mu,P = key
        name = "radius_eq_mu"+mu+"_P"+P
        if os.path.isfile(name+".dat"):
            print("Skipping preexisting "+name) 
        else:
            print("Analysing mu = "+mu+" P = "+P)
            # Ask for limits
            data = np.genfromtxt(files[0])
            num1 = int(float(input("Enter the first timestep to be considered\n")))
            num2 = int(float(input("Enter the last timestep to be considered\n")))
            # normalise limits
            ind1 = np.searchsorted(data[:,0], num1)-1
            ind2 = np.searchsorted(data[:,0], num2)
            # Get average value for each file
            myavg = 0
            mystd = 0
            print("Computing averages:")
            for myfile in files:
                data = np.genfromtxt(myfile)
                x = np.mean(data[ind1:ind2,1])
                print(x)
                myavg += x
                mystd += x**2
            N = len(files)
            print("N files = ", N)
            myavg = myavg/N
            mystd = mystd/N
            P_num = float(P)
            a = np.array([P_num,myavg])
            np.savetxt(name+".dat",a.reshape(1, a.shape[0]))
            if N>1:
                # var = <x^2> - <x>^2
                # unbiased var = (N/(N-1))*var
                # var of mean = var/N
                mystd = (N/(N-1))*(mystd-myavg**2)
                a = np.array([P_num,mystd])
                np.savetxt(name+".err",a.reshape(1, a.shape[0]))
        



