# import packages
import numpy as np, math, pylab, h5py
import verlet, initialization as init

#TODO needs to import data from .hdf5, ideally finding last written thing

RCache: # file to cache results of iteration in, write every 50 iterations to hdf5

def interation(n):
    
