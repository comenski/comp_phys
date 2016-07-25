# this module defines all functions to be imported for later calculations

import initialization as init
import numpy as np

#verlet takes position und velocities from calling module to compute next positions and velocities, returns x1, v2
def verlet(v, rMatrix):    # possible depends on forceLJ/G
    
       
    forceMatrix = np.ma.apply_along_axis(init.forceLJ, 2, rMatrix) 
    force= np.sum(forceMatrix,axis=0)
    v1= v + force*(init.dt/2)
    x1 = rMatrix + v1 * init.dt
    
    forceMatrixT = np.ma.apply_along_axis(init.forceLJ, 2, x1)
    forceT = np.sum(forceMatrixT, axis =0)
    v2= v1 + (forceT) * (init.dt/2)
    
    return x1, v2 
   
    
    
    
     
    #v1= v0 + force *(init.dt/2)
    #x1 = x0 + v1 * init.dt
    #v2= v1 + (forceT) * (init.dt/2)
    
    #v(t+dt/2) = v(t) + F(t)/m * dt/2
    #x(t+dt) = x(t) + v(t+dt/2)* dt
    #v(t+dt) = v(t+dt/2)+F(t+dt)/m * dt/2
    
    
