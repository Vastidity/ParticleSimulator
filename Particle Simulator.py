"""
A simulation of a mutiple body system
Author: Nick Bilgri

21 Dec 2015
Updated 10 Jan 2016
Updated 8 Mar 2016

TODO:
        add graphics
        implement fractional timesteps in update function to allow for more accurate
            timing of collisions.
"""

import math
import random
import numpy as np

#initial set up:
#cube sides of 10000m
#cube density of 3500 kg/m^3 --> cube mass of 3.5*10**15
#space to populate: a sphere of radius 1000000m


elasticity = .800
#gooiness of particles. Found in resolveCollision
G = 6.67e-11
#gravitational constant

class Particle:
    """
    pos = an array of three floats. Initializes particle position. Defaults to [0,0,0]
    vel = an array of three floats. Initializes particle velocity. Defaults to [0,0,0]
    rad = a float. Radius of particle in meters.
    mass = a float. Mass of particle in kg.
    """
    def __init__(self, pos = np.zeros(3), vel = np.zeros(3), rad = 10000.0, mass = 3.5*10**15):
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.rad= rad
        self.mass = mass
    
    def getPos(self):
        return self.pos
        
    def getVel(self):
        return self.vel
        
    def getRad(self):
        return self.rad
    
    def getMass(self):
        return self.mass
        
    def setPos(self,newPos):
        self.pos = np.array(newPos)
    
    def setVel(self,newVel):
        self.vel = np.array(newVel)
    
    def getDist(self,otherParticle):
        dist = 0
        for i in range(3):
            dist += math.pow(self.getPos()[i]-otherParticle.getPos()[i],2)
        return math.sqrt(dist)
                      
    def isColliding(self,otherParticle):
        if self.getDist(otherParticle) <= self.getRad()+otherParticle.getRad():
            return True
        else:
            return False
            
    def __str__(self):
        return ("Particle position: " + str(self.pos)+ " Particle velocity: " + str(self.vel))



def getVector(vecOne,vecTwo):
    """
    Takes two length-three arrays of floats, which it treats as vectors [x,y,z]
    Returns one length-three array of floats equivalent to the vector
            vecTwo - vecOne with a length equal to one.
    """
    if np.array_equal(vecOne,vecTwo):
        return np.array(np.zeros(3))
    difVec = []
    for i in range(3):
         x = vecTwo[i] - vecOne[i]
         difVec += [x]
    difVec = np.array(difVec)
    mag = math.sqrt(difVec.dot(difVec))
    return difVec/mag
    
    
    
def resolveCollision(troublemakers):
    """
    A function to  resolve collisions between two Particles
    troublemakers: a list of two Particle objects
    Returns: a list of two Particle objects. These are the same Particle objects
        from the input list, with their velocities updated to simulate an inelastic collision.
    """
    tempPos1 = troublemakers[0].getPos()
    tempVel1 = troublemakers[0].getVel()
    tempPos2 = troublemakers[1].getPos()
    tempVel2 = troublemakers[1].getVel()
    
    cmPos = np.zeros(3)
    cmVel = np.zeros(3)
    cmPos1 = np.zeros(3)
    cmVel1 = np.zeros(3)
    cmPos2 = np.zeros(3)
    cmVel2 = np.zeros(3)
    #Convert to center of mass frame
    for i in range(3):
        cmPos[i] = .5*(tempPos1[i]+tempPos2[i])
        cmVel[i] = .5*(tempVel1[i]+tempVel2[i])
        cmPos1[i] = .5*(tempPos1[i]-tempPos2[i])
        cmVel1[i] = .5*(tempVel1[i]-tempVel2[i])
        cmPos2[i] = .5*(tempPos2[i]-tempPos1[i])
        cmVel2[i] = .5*(tempVel2[i]-tempVel1[i])
    
    #construct frame for collision
    w = getVector(cmPos1,cmPos2)
    w = w * (1/math.sqrt(w.dot(w)))
    
    v = np.cross(w,np.array([1,1,1]))
    v = v * (1/math.sqrt(v.dot(v)))
    
    u = np.cross(w,v)
    u = u * (1/math.sqrt(u.dot(u)))
    
    #create unit conversion matrix (x,y,z) --> (u,v,w)
    converter = np.matrix([u,v,w])
    try:
        inverse = converter.getI()
    except LinAlgError:
        print "Problem with singular matrix in collision handling"
    
    #convert to (u,v,w)
    uvwVel1 = converter * ((np.asmatrix(cmVel1)).transpose())
    uvwVel2 = converter * ((np.asmatrix(cmVel2)).transpose())
    
    #collide
    uvwVel1[2] = -elasticity * uvwVel1[2]
    uvwVel2[2] = -elasticity * uvwVel2[2]
    
    #convert back to (x,y,z)
    cmVel1 = inverse * (uvwVel1)
    cmVel2 = inverse * (uvwVel2)
    cmVel1 = (cmVel1.transpose().getA())[0]
    cmVel2 = (cmVel2.transpose().getA())[0]

    
    #leave cm frame
    tempPos1 = cmPos1 + cmPos
    tempVel1 = cmVel1 + cmVel
    tempPos2 = cmPos2 + cmPos
    tempVel2 = cmVel2 + cmVel
    
    #set new values
    troublemakers[0].setVel(tempVel1)
    troublemakers[1].setVel(tempVel2)



         
def runSimulation(num = 3, time = 500, scale = 100000):
    """
    num = an integer. Number of particles to simulate
    time = an integer. Number of steps to simulate
    scale = an integer. Particles are randomly generated with position and velocity
        vectors between 0 and 1. The scale factor is multiplied by these values
    returns None, but is currently set to print positions and vectors for each step in the
        simulation. This can be disabled by commenting out the last three lines of the 
        function.
    """
    particles = []
    for i in range(num):
        pos=np.array(np.zeros(3))
        vel= np.array(np.zeros(3))
        for j in range(3):
            pos[j] = scale * random.random()
            vel[j] = scale * random.random()
        particles.append(Particle(pos,vel))
        #testParticle1 = Particle(pos = np.array([0.,0.,0.]),vel = np.array([0.,0.,0.]))
        #testParticle2= Particle(pos = np.array([3001.,500.,0,]),vel = np.array([-.5,0.,0.]))
        #particles.append(testParticle1)
        #particles.append(testParticle2)
    for i in range(time):
        update(particles,1)
        print(i)
        for j in range(num):
            print(str(particles[j]))

def update(particles, timestep):
    """
    A helper function to runSimulaton. It updates the particles over the number
        of timesteps. Multiple timesteps not implemented yet, so this input is currently
        nonfunctional.
    particles: a list of Particle objects.
    timestep: an integer
    Returns: an updated particle list.
    """
    #create temporary lists to store new data. 
    #Use to update all particles at once at the end of update
    newPosList = []
    newVelList = []
    
    for rock in particles:
        #get force on particle in 3-element array
        Fg = np.array(np.zeros(3))
        for otherRock in particles:
            if otherRock != rock:
                forceDirection = getVector(rock.getPos(),otherRock.getPos())
                forceMag = G*rock.getMass()*otherRock.getMass()/(math.pow(rock.getDist(otherRock),2))
                Fg += forceMag * forceDirection
        
        #convert to acceleration
        accel = (Fg)/rock.getMass()
        #get new position
        newPos = np.array(np.zeros(3))
        for i in range(3):
            newPos[i] = rock.getPos()[i] + rock.getVel()[i]*timestep + .5*accel[i]*timestep*timestep
        newPosList.append(newPos)
        
        #get new velocity
        newVel = np.array(np.zeros(3))
        for i in range(3):
            newVel[i] = rock.getVel()[i] + accel[i]*timestep
        newVelList.append(newVel)
            
    # set new pos and vel data
    i=0
    for rock in particles:
        rock.setPos(newPosList[i])
        rock.setVel(newVelList[i])
        i+=1
        
    #check for collisions. Create a list to contain lists of pairs of colliding Particles.
    collisionsToResolve = []
    for rock in particles:
        for otherRock in particles:
            if otherRock != rock:
                if rock.isColliding(otherRock):
                    collisionsToResolve.append([rock,otherRock])
    
    while (len(collisionsToResolve) > 0):
        resolveCollision(collisionsToResolve[0])
        inverse = [collisionsToResolve[0][1],collisionsToResolve[0][0]]
        collisionsToResolve.remove(collisionsToResolve[0])
        collisionsToResolve.remove(inverse)

