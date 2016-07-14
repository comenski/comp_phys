# first MD simulation: create initial data, import verlet.py to compute, write data to array and visualize movement of particles with time

# import packages
import numpy, math, pylab
from scipy.constants import g as GRAV, k as kB
import verlet

# declare constants, variables etc. TODO change from SI to stupid argon units
SIGMA=3.405e-10
EPSILON=1.651e-21
MASS=6.644e-26
numParticles = 100
boxLength = 5 # TODO get box length from volume taken by 100*10 particles as given by ideal gas theory OR 10*particle free range, whichever is bigger
TEMP = 288.15
VAR = (kB*TEMP/MASS)**(1/2)

initCoords = [] #TODO if EFDL = []: write to EFDL in correct format: EFDL0=[ [[xi,yi],[dxi,dyi],0] ]
xList = np.random.sample(numParticles)
yList = np.random.sample(numParticles)
for i in range(numParticles):
    initCoords.append([xList[i],yList[i]]*boxLength)

# create velocity, x and y normal distributed -> [x,y] maxwell-boltzmann distributed
# sum over all velocities to get net total velocity of all particles, substract from each velocity to retain distribution and get zero net momentum

initVel0 = []
dxList = np.random.normal(0,VAR,100)
dyList = np.random.normal(0,VAR,100)
for i in range(numParticles):
    initVel0.append([dxList[i],dyList[i]])

netTotalVel = [0,0]
for i in range(len(initVel0)):
    for j in range(len(initVel0[i])):
        netTotalVel[j] += initVel0[i][j]

initVel = initVel0 #TODO append to external_file_data_list[0]
for i in range(len(initVel0)):
    for j in range(len(netTotalVel)):
        initVel[i][j] -= netTotalVel[j]/len(initVel0)

#-------------------------------

# define interaction according to leonard-jones-potential for problems 5,7
# goes into respective module to be used by verlet.py or
# goes directly into verlet.py or gets provided from definitions.py to verlet.py and module takes an additional argument LJ/LJG/throw for force
def forceLJ(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]

FG = [0,GRAV*MASS]
def forceLJG(r):
    R=SIGMA/r
    return 4*EPSILON*(12*R**13-6*R**7)*[xj-xi,yj-yi]+FG



#---------------------------------- ideas for iteration funtion TODO still need to implement the force at some point, possibly directly in iteration(n)

def interation(n):
    if n > EFDL[-1][-1]:
        print("done")
    else:
        for i in range(numParticles):
            rho = [[0,0]]*numParticles
            nu = [[0,0]]*numParticles
            Ri = EFDL[-1][i][0]
            dRi = EFDL[-1][i][1]
            for j in range(numParticles):
                Rj = EFDL[-1][j][0]
                dRj = EFDL[-1][j][1]
                if i==j:
                    skip #TODO find correct way to jump to next j (continue, break etc.)
                else:
                    rAbs = math.sqrt(r[0]**2+r[1]**2)
                    mu = verlet.LJ(Ri,dRi,Rj,dRj)
                    rho[i][0] += mu[0][0]
                    rho[i][1] += mu[0][1]
                    nu[i][0] += mu[1][0]
                    nu[i][1] += mu[1][1]
#        for i in range(numParticles): # TODO write rho, nu, n to last applicable line



        n += 1
#        if n % 50 == 0:
#            # TODO append EFDL[-1] to EFDL
