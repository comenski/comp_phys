# first MD simulation: create initial data, import verlet.py to compute, write data to array and visualize movement of particles with time

# import packages
import numpy as np
import verlet
from scipy.constants import g

# declare constants, variables etc.
SIGMA=3.405e-10
EPSILON=1.651e-21
MASS=6.644e-26
numParticles = 100

initCoords = []
for i in range(numParticles):
    xList = np.random.sample(100)
    yList = np.random.sample(100)
    initCoords.append([xList[i],yList[i]]*length)

initVel0 = create [[dxi,dyi]*numParticles] # create array according to maxwell-boltzmann

AvgVel = [0,0]
for i in range(len(initVel0)):
    AvgVel += initVel[i]

initVel = initVel0
for i in range(len(initVel0)):
    initVel[i] -= AvgVel/len(initVel0)

#define interaction according to leonard-jones-potential
def forceLJ(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]

FG = [0,g*MASS]
def forceLJG(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]+FG
