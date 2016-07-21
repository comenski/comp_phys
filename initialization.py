import numpy as np, math
from scipy.constants import g as GRAV, k as kB

# declare constants, variables etc. TODO change from SI to stupid argon units
SIGMA=3.405e-10
EPSILON=1.651e-21
MASS=6.644e-26
numParticles = 100 
boxLength = 5 # TODO get box length from volume taken by 100*10 particles as given by ideal gas theory OR 10*particle free range, whichever is bigger
TEMP = 288.15
VAR = (kB*TEMP/MASS)**(1/2) # According to ideal gas theory
rMax = 0.5*boxLength


# get initial coordinates, randomly distributed in space
initCoords = [] #TODO if EFDL = []: write to EFDL in correct format: EFDL0=[ [[xi,yi],[dxi,dyi],0] ]
xList = np.random.sample(numParticles)
yList = np.random.sample(numParticles)
for i in range(numParticles):
    initCoords.append([xList[i],yList[i]]*boxLength)

# create velocity distribution, dx and dy normal distributed -> [dx,dy] maxwell-boltzmann distributed
initVel0 = []
dxList = np.random.normal(0,VAR,100)
dyList = np.random.normal(0,VAR,100)
for i in range(numParticles):
    initVel0.append([dxList[i],dyList[i]])

# sum over all velocities to get net total velocity of all particles
netTotalVel = [0,0]
for i in range(len(initVel0)):
    for j in range(len(initVel0[i])):
        netTotalVel[j] += initVel0[i][j]

# substract from each velocity to retain distribution and get zero net momentum.
initVel = initVel0
for i in range(len(initVel0)):
    for j in range(len(netTotalVel)):
        initVel[i][j] -= netTotalVel[j]/len(initVel0)

#------------------------------- 

# define interaction according to leonard-jones-potential for problems 5,7 
# goes into respective module to be used by verlet.py or 
# goes directly into verlet.py or gets provided from definitions.py to verlet.py and module takes an additional argument LJ/LJG/throw for force 
 
FG = [0,GRAV*MASS] 
 
 
#def forceLJ(r):
#   return 4*EPSILON*(12*SIGMA**12/r**14 - 6*SIGMA**6/r**8)*[xj-xi,yj-yi] #TODO correct nomenclature
 
# let's try with vector this time:
F= np.zeros(numParticles,numParticles) 

def position():
    for i in range(1, numParticles):
        for j in range(i-1):
            F[i][j]=forceLJ(r,i,j)
    for j in range(1,numParticles):
        for i in range(i-1):
            F[j][i]=-F[i][j]
position()
 
def forceLJ(r,i,j):
    [xi,yi,xj,yj]=[r[i][0],r[i][1],r[j][0],r[j][1]]
    rAbs = math.sqrt(xi**2-xj**2+yi**2-yj**2)
    if rAbs < rMax:
        return 4*EPSILON*(12*SIGMA**12/rAbs**14 - 6*SIGMA**6/rAbs**8)*[xj-xi,yj-yi]
    xj += boxLength
    if rAbs < rMax:
        return 4*EPSILON*(12*SIGMA**12/r**14 - 6*SIGMA**6/r**8)*[xj-xi,yj-yi]
    xj -= boxLength
    yj += boxLength
    if rAbs < rMax:
        return 4*EPSILON*(12*SIGMA**12/r**14 - 6*SIGMA**6/r**8)*[xj-xi,yj-yi]
    else:
        return 0




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
