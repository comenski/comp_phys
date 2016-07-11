# first MD simulation: create initial data, import verlet.py to compute, write data to array and visualize movement of particles with time

# import packages
import numpy as np

# declare constants, variables etc.
numParticles = 100

# initialList = create [[xi,yi]*100] # random distribution with abs(xi) < boxLength > abs(yi)

initVel = create [[dxi,dyi]*numParticles] # create array according to maxwell-boltzmann

AvgVel = [0,0]
for i < range(len(initVel)):
    AvgVel += initVel[i]

initVelNorm = initVel
for i < range(len(initVel)):
    initVelNorm[i] -= AvgVel/len(initVel)




# initVel= [initVel0
