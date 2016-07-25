import numpy as np, math
from scipy.constants import g as GRAV, k as kB
import time

# declare constants, variables etc. TODO change from SI to stupid argon units
SIGMA=3.405e-10
EPSILON=1.651e-21
MASS=6.644e-26
numParticles = 100
boxLength = 5 # TODO get box length from volume taken by 100*10 particles as given by ideal gas theory OR 10*particle free range, whichever is bigger
TEMP = 288.15
VAR = (kB*TEMP/MASS)**(1/2) # According to ideal gas theory
rMax = 0.5*boxLength #TODO completely arbitrary atm

time1 = time.clock()

# get initial coordinates, randomly distributed in space # TODO write to data file[0]
initCoords = np.random.sample((numParticles,2))*boxLength

time2 = time.clock()

#create velocity distribution, dx and dy normal distributed -> [dx,dy] maxwell-boltzmann distributed
#initVel0 = []
#dxList = np.random.normal(0,VAR,100)
#dyList = np.random.normal(0,VAR,100)
#for i in range(numParticles):
#    initVel0.append([dxList[i],dyList[i]])
#
## sum over all velocities to get net total velocity of all particles
#netTotalVel = [0,0]
#for i in range(len(initVel0)):
#    for j in range(len(initVel0[i])):
#        netTotalVel[j] += initVel0[i][j]
#
## substract from each velocity to retain distribution and get zero net momentum.
#initVel = initVel0
#for i in range(len(initVel0)):
#    for j in range(len(netTotalVel)):
#        initVel[i][j] -= netTotalVel[j]/len(initVel0)

#-------------------------------

# define interaction according to leonard-jones-potential for problems 5,7
# goes into respective module to be used by verlet.py or
# goes directly into verlet.py or gets provided from definitions.py to verlet.py and module takes an additional argument LJ/LJG/throw for force

# creating matrices for distances, masked for maximum range of interaction and added together to satisfy periodic boundary conditions
rMatrix = np.tile(initCoords,(numParticles,1)).reshape(numParticles,numParticles,2) # create 100x100x2 from 100x2
rMatrix -= np.transpose(rMatrix,(1,0,2)) # transpose and substract to create all possible pairs
rMatrix = np.ma.masked_outside(rMatrix, -rMax, rMax) # mask when outside of maximum range of computation
rMatrix = np.ma.masked_values(rMatrix, 0) # finally mask where r=0
#rMatrix = np.ma.filled(rMatrix, 0) # set masking value to 0 for summation of arrays

time3 = time.clock()
# same for particles created by boundary conditions
#rMatrixBound1 = np.tile(initCoords,(numParticles,1)).reshape(numParticles,numParticles,2)
#rMatrixBound1 -= np.transpose(rMatrixBound1+[0,boxLength],(1,0,2))
#rMatrixBound1 = np.ma.masked_outside(rMatrixBound1, -rMax, rMax)
#rMatrixBound1 = np.ma.filled(rMatrixBound1, 0)
#
#rMatrixBound2 = np.tile(initCoords,(numParticles,1)).reshape(numParticles,numParticles,2)
#rMatrixBound2 -= np.transpose(rMatrixBound2+[boxLength,0],(1,0,2))
#rMatrixBound2 = np.ma.masked_outside(rMatrixBound2, -rMax, rMax)
#rMatrixBound2 = np.ma.filled(rMatrixBound2, 0)
#
#rMatrixBound3 = np.tile(initCoords,(numParticles,1)).reshape(numParticles,numParticles,2)
#rMatrixBound3 -= np.transpose(rMatrixBound3-[0,boxLength],(1,0,2))
#rMatrixBound3 = np.ma.masked_outside(rMatrixBound3, -rMax, rMax)
#rMatrixBound3 = np.ma.filled(rMatrixBound1, 0)
#
#rMatrixBound4 = np.tile(initCoords,(numParticles,1)).reshape(numParticles,numParticles,2)
#rMatrixBound4 -= np.transpose(rMatrixBound4-[boxLength,0],(1,0,2))
#rMatrixBound4 = np.ma.masked_outside(rMatrixBound4, -rMax, rMax)
#rMatrixBound4 = np.ma.filled(rMatrixBound2, 0)
#
#rMatrix += rMatrixBound1 + rMatrixBound2 + rMatrixBound3 + rMatrixBound4 # add up all the effects
#rMatrix = np.ma.masked_values(rMatrix, 0) # finally mask where r=0
time4 = time.clock()
#print(rMatrix)

# define function that computes force from leonard-jones-potential, gets called by verlet.py
def forceLJ(r):
    rAbs = math.sqrt(r[0]**2+r[1]**2)
    return 4*EPSILON*(12*SIGMA**12/rAbs**14 - 6*SIGMA**6/rAbs**8)*r

testForceMatrix = np.ma.apply_along_axis(forceLJ, 2, rMatrix) # TODO just for testing purposes
time5 = time.clock()

print(np.sum(testForceMatrix,axis=0))
time6 = time.clock()

print(np.array([time2,time3,time4,time5,time6])-time1)

# TODO call with np.apply_along_axis(forceLJ, 2, rMatrix) in respective modules
FG = [0,GRAV*MASS]


#print(forceLJ(rMatrix))




#def forceLJ(r):
#    [xi,yi,xj,yj]=[r[i][0],r[i][1],r[j][0],r[j][1]]
#    rAbs = math.sqrt(xi**2-xj**2+yi**2-yj**2)
#    if rAbs < rMax:
#        return 4*EPSILON*(12*SIGMA**12/rAbs**14 - 6*SIGMA**6/rAbs**8)*[xj-xi,yj-yi]
#    xj += boxLength
#    if rAbs < rMax:
#        return 4*EPSILON*(12*SIGMA**12/r**14 - 6*SIGMA**6/r**8)*[xj-xi,yj-yi]
#    xj -= boxLength
#    yj += boxLength
#    if rAbs < rMax:
#        return 4*EPSILON*(12*SIGMA**12/r**14 - 6*SIGMA**6/r**8)*[xj-xi,yj-yi]
#    else:
#        return 0




#----------------------------------------------OLD STUF, PROBS NO LONGER APPLICABLE

# TODO two better ideas performance-wise for forceLJG:
# use forceLJ(r) for heated case as well, add after summation over all j has happened (doesn't need verlet b/c FG is constant):
# -------- goes to the end of iteration(n) ------------------
# drFG = integral(1/2*GRAV*dt**2)
# dvFG = GRAV*dt
# if j==numParticles+1:
#     rho+=drFG , nu+=dvFG , n
# this is equal to treating gravity as n+1th particle with constant interaction and should be warranted, given the
# small timesteps used for simulation and the relative weakness of gravitational pull
