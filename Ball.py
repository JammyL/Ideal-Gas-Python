# -*- coding: utf-8 -*-
"""
This file contains details of the ball class.
"""
import reactions as r
import numpy as np
from matplotlib.patches import Circle




class Ball(r.reactant):
    
    #Defines a ball id neccessary for keeping track of which ball is involved in a collision
    _id = 0
    
    #resets the id counter to allow simulations to be run on for loops for data collection
    def reset_id():
        Ball._id = 0
    
    def __init__(self, mass = 1.0, pos = [0,0], vel = [0,0], container = False, radius = 1.0, animate = True, reactant_type = 'reactant_A', distance_check = True):
        
        #Intialises a reactant with the appropriate type
        r.reactant.__init__(self, reactant_type = reactant_type, animate = animate)
        
        #Exceptions to prevent incorrect arguments being given
        if len(pos) != 2:
            raise Exception('r takes 2 arguments but {} were given'.format(len(pos)))
        if len(vel) != 2:    
            raise Exception('v takes 2 arguments but {} were given'.format(len(vel)))
        if type(container) != bool:
            raise Exception('container must be of type bool but type {} was given'.format(type(container)))
        
        #Sets variables. 
        self.__container = container
        self.__distance_check = distance_check
        if container == False:
            self.__m = float(mass)
            self.__r = np.array(pos, dtype = 'float')
            self.__v = np.array(vel, dtype = 'float')
            self.__rad = float(radius)
            self.__time_to_collision = float(np.Inf)
            self.__colliding_ball = self
            self.__animate = animate
            if self.__animate:
                self.__patch = Circle(self.__r, self.__rad, color = self.get_patch_color(), ec = 'k')
        else:
            self.__m = np.Inf
            self.__r = np.array([0,0], dtype = 'float')
            self.__v = np.array([0,0], dtype = 'float')
            self.__rad = float(radius)
            self.__time_to_collision = float(np.Inf)
            self.__colliding_ball = self
            self.__animate = animate
            if self.__animate:
                self.__patch = Circle(self.__r, self.__rad, color = 'b', fill = False)
        
        #Assigns the balls id and increments id counter
        self.__hash = Ball._id
        Ball._id += 1
   
    #Defines hash function to return the ball id
    def __hash__(self):
        return hash(self.__hash)
    
    #Defines eq function to check if ids are the same
    def __eq__(self, other):
        return self is other
    
    #Defines neq function to check ids are different
    def __neq__(self, other):
        return not self.__eq__(other)
      
    #Funcitons return properties of the ball              
    def mass(self):
        return self.__m
    def pos(self):
        return self.__r
    def vel(self):
        return self.__v
    def kinetic_energy(self):
        return 0.5*self.__m * np.dot(self.__v, self.__v)
    def radius(self):
        return self.__rad
    def container(self):
        return self.__container
    def get_time(self):
        return float(self.__time_to_collision)
    def get_colliding_ball(self):
        return hash(self.__colliding_ball)
    def get_patch(self):
        if self.__animate:
            return self.__patch
        else:
            raise Exception('Animation is disabled')
            
    '''
    Resets the time to collision of a ball to infinity
    This means any new calculation between it and another ball will
    result in that other ball becoming its new collision partner
    '''    
    def reset_time(self):
        self.__time_to_collision = np.Inf
    
    #Moves the ball at its current velocity over time dt
    def move(self,dt):
        dt = float(dt)
        self.__r += self.__v*dt
        self.__time_to_collision -= dt
        if self.__animate:
            self.__patch.centre = self.pos()
        return self.__r
        
    '''
    Calculates the time to collision between two balls.
    If this time is lower than the current time being stored by self,
    the ball replaces its time to collision and colliding ball with the values
    calculated here
    '''
    def time_to_collision(self,other):
        r = self.__r - other.__r
        v = self.__v - other.__v
        if self.__container:
            raise Exception('Cannot collide container with other balls. Container must be second ball')
        if other.__container:
            rad = self.__rad - other.__rad
            det = (np.dot(r,v)**2) - ((np.dot(r,r) - (rad**2))*np.dot(v,v))
            if det < 0:
                t_1 = np.Inf
            else:    
                t_1 = (- (np.dot(r,v)) + np.sqrt(det))/(np.dot(v,v))
                if t_1 <= 0:
                    t_1 = np.Inf

        else:
            rad = self.__rad + other.__rad
            det = (np.dot(r,v)**2) - ((np.dot(r,r) - (rad**2))*np.dot(v,v))
            if det < 0:
                t_1 = np.Inf
            else:
                t_1 = (- (np.dot(r,v)) - np.sqrt(det))/(np.dot(v,v))
                if t_1 <= 0:
                    t_1 = np.Inf
                    
        if t_1 < self.__time_to_collision:
            self.__time_to_collision = t_1
            self.__colliding_ball = other
            
    '''
    Calculates the new velocities two balls will have should they collide at
    this current moment in time
    '''
    def new_vel(self,other):
        if other.__container:
            v_1 = self.__v - ((2*np.dot(self.__v - other.__v, self.__r - other.__r) * (self.__r - other.__r))/(np.dot(self.__r - other.__r, self.__r - other.__r))) 
            self.__v = v_1
        else:
            v_1 = self.__v - ((2*other.__m*np.dot(self.__v - other.__v, self.__r - other.__r) * (self.__r - other.__r))/((self.__m + other.__m)*np.dot(self.__r - other.__r, self.__r - other.__r)))
            v_2 = other.__v - ((2*self.__m*np.dot(other.__v - self.__v, other.__r - self.__r) * (other.__r - self.__r))/((other.__m + self.__m)*np.dot(other.__r - self.__r, other.__r - self.__r)))
            self.__v = v_1
            other.__v = v_2
      
    '''
    Moves two balls by the time to collision of self. This should also be
    the time to collision of other. The new velocities are then calculated
    for these two balls and their times to collision are reset. 
    '''
    def collide(self, other):
        dt = self.__time_to_collision
        #Raises exception if particles will not collide.
        if dt == np.Inf:
            raise Exception('Particles will not collide')
        else:
            self.move(dt)
            other.move(dt)
            #Checks the particles are the right distance apart when colliding
            #Can be switched off.
            if self.__distance_check:
                displacement_apart = self.pos() - other.pos()
                distance_apart = np.sqrt(np.dot(displacement_apart, displacement_apart))
                if other.__container:
                    if abs(distance_apart - (other.__rad - self.__rad)) > 1e-4:
                        raise Exception('Particle and container are colliding but not touching')
                else:
                    if abs(distance_apart - (other.__rad + self.__rad)) > 1e-4:
                        raise Exception('Particle and particle are colliding but not touching')
                
            if self.__animate:
                #Sets the patch colour in the event of a reaction and animation is on
                self.__patch.set_facecolor(self.get_patch_color())
                other.__patch.set_facecolor(other.get_patch_color())
            self.reset_time()
            other.reset_time()
            self.new_vel(other)
