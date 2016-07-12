# first MD simulation: create initial data, import verlet.py to compute, write data to array and visualize movement of particles with time

# import packages
import numpy as np
import verlet
from scipy.constants import g as GRAV, k as kB
import pylab

# declare constants, variables etc.
SIGMA=3.405e-10
EPSILON=1.651e-21
MASS=6.644e-26
numParticles = 100
boxLength = 5 # TODO get box length from volume taken by 100*10 particles as given by ideal gas theory OR 10*particle free range, whichever is bigger
TEMP = 288.15
VAR = (kB*TEMP/MASS)**(1/2)

initCoords = []
xList = np.random.sample(numParticles)
yList = np.random.sample(numParticles)
for i in range(numParticles):
    initCoords.append([xList[i],yList[i]]*boxLength)

# create velocity, x and y normal distributed -> [x,y] maxwell-boltzmann distributed
# sum over all velocities to get net momentum of all particles, substract net momentum/number of particles from each velocity to retain distribution and get zero net momentum

initVel0 = []
dxList = np.random.normal(0,VAR,100)
dyList = np.random.normal(0,VAR,100)
for i in range(numParticles):
    initVel0.append([dxList[i],dyList[i]])

netTotalVel = [0,0]
for i in range(len(initVel0)):
    for j in range(len(initVel0[i])):
        netTotalVel[j] += initVel0[i][j]

initVel = initVel0
for i in range(len(initVel0)):
    for j in range(len(netTotalVel)):
        initVel[i][j] -= netTotalVel[j]/len(initVel0)

#-------------------------------

#define interaction according to leonard-jones-potential for problems 5,7
def forceLJ(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]

FG = [0,GRAV*MASS]
def forceLJG(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]+FG
