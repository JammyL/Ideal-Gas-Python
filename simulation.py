#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Simulation class for the reactions
"""
#import csv
#import time
import numpy as np
import Ball as b
import reactions as r
import matplotlib.pyplot as plt
from scipy.constants import golden
from random import sample
from random import randint
from tqdm import tqdm
from copy import deepcopy


        

class reaction_Simulation:
    
    
    def __init__(self ,ball_radius = 1.0, container_radius = 10.0, v_0 = 20, numbers = [1,0,0,0,0,0], animate = True, save_location = 'Data/'):
        
        #Exceptions to check the right length of array is given
        while len(numbers) < 6:
            numbers.append(0.)
        r.reactant.reset_values()  
        if len(numbers) > 6:
            raise Exception('List of number of balls must be length 6 or less but length {} was given'.format(len(numbers)))
        #Defines various variables
        self.__save_location = save_location
        #Stores current time in simulation
        self.__total_time = 0.
        #Stores the time at each iteration
        self.__time_list = []
        #Stores the total momentum change on the container
        self.__total_delta_p = 0.
        #Stores the total momentum change at each iteration
        self.__delta_p_list = []
        self.__ball_rad = float(ball_radius)
        self.__container_rad = float(container_radius)
        #Main list of balls ordered in time to collision
        self.__balls_list = []
        #List of balls in order of creation
        self.__reference_list = []
        #Stores the number of each type of each particle
        self.__reaction_values = []
        #Stores the starting values of each type of particle
        self.__starting_numbers = deepcopy(numbers)
        self.__current_frame = 0
        #Stores the current values of each type of yet to be created!!!!!!!!!!!!!
        self.__numbers = numbers
        #Stores the total number of particles to be created
        self.__n = sum(numbers)
        #Stores the average value of the gaussian from which speeds are calculated
        self.__v_0 = float(v_0)
        self.__animate = animate    
        #Creates a dictionary
        self.__d = {}

        '''
        Arranges the particles evenly using the golden ratio to obtain an 
        even distribution rotationally
        '''
        
        #Sets a base magnitude of 1
        x =  1/np.sqrt(2)
        y =  1/np.sqrt(2)
        phi = 0
        
        #Sets particle 1 out by its radius from the origin
        rho = self.__ball_rad
        #Creates a list of magnitudes and angles 
        rho_list = []
        phi_list = []
        
        #Creates two lists of velocities
        v_x = []
        v_y = []
        
        #If number of particles is odd, adds one extra value.
        if self.__n % 2 != 0:
            v_x.append(v_0)
            v_y.append(v_0)
            
        #Generates n/2 random positive velocities sampled from a normal distriubtion
        for i in range(self.__n//2):
            v_x_i = abs(np.random.normal(v_0/np.sqrt(2), v_0/4))
            v_y_i = abs(np.random.normal(v_0/np.sqrt(2), v_0/4))
            
            #Appends both a +ve and -ve value of velocity for every value generated
            #This means that for even numbers, mean velocity is 0
            v_x.append(v_x_i)
            v_y.append(v_y_i)
            v_x.append(-v_x_i)
            v_y.append(-v_y_i)
            
        #Shuffles the velocities in the two arrays. 
        v_x = sample(v_x, len(v_x))
        v_y = sample(v_y, len(v_y))
        
        #For loop to generate particles
        for i in range(self.__n+1):   
            #Prevents particles being outside container
            if rho > self.__container_rad - self.__ball_rad:
                raise Exception('Ball is outside of the container')
            rho_list.append(rho)
            phi_list.append(phi)
            
            #Rotates round by a full rotation * golden ratio
            phi += 2*np.pi*golden
            
            #Increments distance with 1/rho to account for area relationship
            rho += (self.__ball_rad**2/rho)
            
        #Moves all particles out to fill the space of the container
        alpha = 0.999*(self.__container_rad - self.__ball_rad)/rho_list[-1]
        rho_list = np.array(rho_list)*alpha
        
        print('Creating Reference list')
        reference_progress = tqdm(total = self.__n, postfix = None) 
        
        #Creates balls using generated positions and velocities
        ball_number = 0
        while sum(numbers) > 0:
            #Selects random type of ball to be created.
            type_pos = randint(0,5)
            if numbers[type_pos] == 0:
                continue
            numbers[type_pos] -= 1
            
            #Determines position of ball
            x_i = (rho_list[ball_number] * x * np.cos(phi_list[ball_number])) - (rho_list[ball_number] * y * np.sin(phi_list[ball_number]))
            y_i = (rho_list[ball_number] * x * np.sin(phi_list[ball_number])) + (rho_list[ball_number] * y * np.cos(phi_list[ball_number]))
           
            #Determines type of reactant
            if type_pos == 0:
                ball_type = 'reactant_A'
            elif type_pos == 1:
                ball_type = 'reactant_B'
            elif type_pos == 2:
                ball_type = 'product_C'
            elif type_pos == 3:
                ball_type = 'product_D'
            elif type_pos == 4:
                ball_type = 'catalyst_E'
            elif type_pos == 5:
                ball_type = 'catalyst_F'
                
            ball = b.Ball(radius = self.__ball_rad, vel = [v_x[ball_number],v_y[ball_number]], pos = [x_i,y_i], animate = self.__animate, reactant_type = ball_type)
            self.__reference_list.append(ball)
            reference_progress.update(n=1)
            ball_number += 1
            
        reference_progress.close()
        print('Done')
        
        #Creates the container
        self.__container = b.Ball(container = True, radius = self.__container_rad, reactant_type = 'container')
        self.__reference_list.append(self.__container)
        
        
        print('Calculating time to collisions')
        #Calculates the number of calculations for progress bar
        number_of_calcs = 0
        for i in range(0, self.__n):
            number_of_calcs += i 
        ball_time_progress = tqdm(total = number_of_calcs, postfix = None) 
        
        #Calculates the time to collision for all balls in the reference list
        for i in range(0, self.__n):
            for j in range(i+1, self.__n+1):
               self.__reference_list[i].time_to_collision(self.__reference_list[j])
               ball_time_progress.update(n = 1)
        
        #Puts the container in the main balls list
        self.__balls_list.append(self.__reference_list[-1])
        
        ball_time_progress.close()
        print('Done')
        print('Placing Balls')
        ball_place_progress = tqdm(total = self.__n, postfix = None) 
        
        #Places ball in main ordered list by time to collision
        for i in range(0,self.__n):
            place = 0
            for el in self.__balls_list:
                if self.__reference_list[i].get_time() < el.get_time():
                    break
                else:
                    place += 1     
            self.__balls_list.insert(place,self.__reference_list[i])
            ball_place_progress.update(n=1)
        
        ball_place_progress.close()
        
        print('Done')
        print('Creating dictionary')
        
        hash_progress = tqdm(total = self.__n, postfix = None) 
        
        #Creates dictionary of each ball
        for el in self.__reference_list:
            self.__d[hash(el)] = set()
            hash_progress.update(n=1)
        hash_progress.close()
        print('Done')    
        print('Filling dictionary')
        
        fill_progress = tqdm(total = self.__n, postfix = None) 
        
        #Fills dictionaries with corresponding balls to collide
        for el in self.__reference_list[:-1]:
            self.__d[el.get_colliding_ball()].add(el)
            fill_progress.update(n=1)
        fill_progress.close()
        print('Done')
        
        #Calculates area of contianer
        self.__A = np.pi * self.__reference_list[-1].radius()**2
        
        #Resets id counter for multiple simulations
        b.Ball.reset_id()
                
    
    #Function which puts an item into the main ball list by time to collision
    def put(self, item, list_A):
        place = len(list_A)
        for el in reversed(list_A):
            if item.get_time() < el.get_time():
                place -= 1
            else:
                break
        list_A.insert(place, item)
        
    #Runs the next collision in the simulation
    def next_collision(self):
        
        self.__current_frame += 1
        
        #Removes first item from main balls list and finds corresponding ball to collide with
        colliding_ball_1 = self.__balls_list.pop(0)
        colliding_ball_2 = self.__reference_list[colliding_ball_1.get_colliding_ball()]
        
        #Removes corresponding ball from main list
        self.__balls_list.remove(colliding_ball_2)
        
        #Removes balls from each others dictionaries
        if colliding_ball_1 in self.__d[hash(colliding_ball_2)]:
            self.__d[hash(colliding_ball_2)].remove(colliding_ball_1)
        if colliding_ball_2 in self.__d[hash(colliding_ball_1)]:
            self.__d[hash(colliding_ball_1)].remove(colliding_ball_2)
    
        #Calculates the time moves from this collision and updates relevant values
        t_min = colliding_ball_1.get_time()
        
        self.__total_time += t_min
        self.__time_list.append(self.__total_time)
        
        #Collects data for Maxwell distribution
        if self.__collect_data:
            if self.__current_frame % self.__frame_distance == 0:  
                velocities = []
                for el in self.__reference_list[:-1]:
                    mag_v = np.sqrt(np.dot(el.vel(), el.vel()))
                    velocities.append(mag_v)
                self.__velocities_list.append(velocities)  
        
        #Collects data for pressures
        if self.__collect_data:
            if colliding_ball_2.container():
                v_i = colliding_ball_1.vel()
            else:
                v_i = 0
        
        #Reacts and collides the two balls in question
        colliding_ball_1.react(colliding_ball_2)
        colliding_ball_1.collide(colliding_ball_2)
        
        #Moves all other balls by the timestep dt
        for el in self.__balls_list:
            el.move(t_min)
                
    
        #Recalculates the time to collision of colliding ball 1 if ball 2 is the container
        if colliding_ball_2.container():
            colliding_ball_1.time_to_collision(colliding_ball_2)

        
        #Recalculates the time to collision of ball 1 and ball 2 with all other balls
        for el in self.__balls_list:
            if colliding_ball_2.container():
                colliding_ball_1.time_to_collision(el)

            else:
                colliding_ball_1.time_to_collision(el)
                colliding_ball_2.time_to_collision(el)

        #Creates temporary sets for new dictionaries to be created
        dummy_set_1 = set()
        dummy_set_2 = set()
        
        '''
        Systematically removes every element from ball 1's dictionary and
        main ball list.
        Recalculates its time to collision with all other balls.
        Reassigns the ball to a new dictionary and replaces it in the main list
        '''
        while self.__d[hash(colliding_ball_1)]:
            element = self.__d[hash(colliding_ball_1)].pop()
            self.__balls_list.remove(element)
            element.reset_time()
            for balls in self.__reference_list:
                if hash(balls) == hash(element):
                    continue
                element.time_to_collision(balls)
            if element.get_colliding_ball() == hash(colliding_ball_1):
                dummy_set_1.add(element)
            elif element.get_colliding_ball() == hash(colliding_ball_2):
                dummy_set_2.add(element)
            else:
                self.__d[element.get_colliding_ball()].add(element)
            self.put(element, self.__balls_list)
        
        '''
        Does the same as above if colliding ball 2 is not the container
        '''
        if colliding_ball_2.container() != True:
            while self.__d[hash(colliding_ball_2)]:
                element = self.__d[hash(colliding_ball_2)].pop()
                self.__balls_list.remove(element)
                element.reset_time()
                for balls in self.__reference_list:
                    if hash(balls) == hash(element):
                        continue
                    element.time_to_collision(balls)
                if element.get_colliding_ball() == hash(colliding_ball_1):
                    dummy_set_1.add(element)
                elif element.get_colliding_ball() == hash(colliding_ball_2):
                    dummy_set_2.add(element)
                else:
                    self.__d[element.get_colliding_ball()].add(element)
                self.put(element, self.__balls_list)

        #restores the dictionarys for ball 1 and 2
        self.__d[hash(colliding_ball_1)] = dummy_set_1       
        self.__d[hash(colliding_ball_2)] = dummy_set_2   
        
        #Puts ball 1 and 2 into the dictionaries of the balls they will now collide with
        if colliding_ball_2.container():
            self.__d[hash(colliding_ball_1.get_colliding_ball())].add(colliding_ball_1)
        else:
            self.__d[hash(colliding_ball_1.get_colliding_ball())].add(colliding_ball_1)
            self.__d[hash(colliding_ball_2.get_colliding_ball())].add(colliding_ball_2)

        #Puts ball 1 and 2 back into the main balls list
        self.put(colliding_ball_1, self.__balls_list)
        self.put(colliding_ball_2, self.__balls_list)
          
        #Collects data about that collision
        if self.__collect_data:         
            if colliding_ball_2.container():
                v_f = colliding_ball_1.vel()
            else:
                v_f = 0
            
            v = v_f - v_i
            delta_p = colliding_ball_1.mass() * np.sqrt(np.dot(v,v))                        
            self.__total_delta_p += delta_p
            self.__delta_p_list.append(self.__total_delta_p)
        
            values = r.reactant.get_values()
            values.append(self.__total_time)
            self.__reaction_values.append(values)

        
    def run(self, num_frames, collect_data = False):
        
        #Defines variables for data collection
        self.__collect_data = collect_data
        if collect_data:
            self.__num_frames  = num_frames
            self.__frame_distance = num_frames/5
            self.__velocities_list = []
            
            
        print('Running')
        
        #Creates the scene to animate the simulation
        if self.__animate:
            fig, ax = plt.subplots(subplot_kw={'aspect': 'equal'})
            plt.xlim((-self.__container_rad, self.__container_rad))
            plt.ylim((-self.__container_rad, self.__container_rad))
                        
            for el in self.__balls_list:
                ax.add_patch(el.get_patch())
                            
        progress = tqdm(total = num_frames, postfix = None)
        
        # Simulates n frames by repeatedly calling next collision
        for frame in range(num_frames):
            progress.update(n = 1)            
            self.next_collision()
            if self.__animate:
                plt.pause(0.001)
        
        progress.close()  
           
        #Organises and writes data into txt files
        
        if collect_data:
            
            data_r = self.__reaction_values
            
            #Defines names and locations for txt files
            save_string = self.__save_location + "Data_A_{}_N_{}_v_{}_rad_{}.txt"
            save_string_maxwell =  self.__save_location + "Data_A_{}_N_{}_v_{}_rad_{}_Maxwell.txt"
            save_string_reactions = self.__save_location + "Reactions_data_A_{}_B_{}_C_{}_D_{}_E_{}_F_{}.txt"
            
            #Saves data for reactions
            with open(save_string_reactions.format(self.__starting_numbers[0], self.__starting_numbers[1],self.__starting_numbers[2], self.__starting_numbers[3],self.__starting_numbers[4], self.__starting_numbers[5]), mode = 'w') as f:
                for el in data_r:
                    f.write('{} {} {} {} {} {} {}\n'.format(el[0], el[1], el[2], el[3], el[4], el[5], el[6]))
            
            #Creates a list of kinetic energies. 
            #Since kinetic energy is consereved this can just be repeated.
            kinetic_energy_list = []
            for frame in range(num_frames):
                total_K_E = 0
                for el in self.__balls_list[:-1]:
                    total_K_E += el.kinetic_energy() 
                kinetic_energy_list.append(total_K_E)
                            
            pressure_list = []
            container_radius = self.__balls_list[-1].radius()
            
            #Calculates the pressure over times dt
            for i in range(num_frames):
                pressure = self.__delta_p_list[i]/(self.__time_list[i] * 2 * np.pi * container_radius)
                pressure_list.append(pressure)

            l = [self.__time_list,pressure_list, kinetic_energy_list]
            data = zip(*l)
            data_maxwell = self.__velocities_list
            
            #Writes data for pressures and kinetic energy to txt file
            with open(save_string.format(self.__R, self.__n, self.__v_0, self.__ball_rad), mode = 'w') as f:
                f.write('Rad:{}, Number:{}, v_0:{}\n'.format(self.__A, self.__n, self.__v_0))

                for el in data:
                    f.write('{} {} {} \n'.format(el[0], el[1], el[2]))
            
            #Writes data for maxwell-boltzmann to txt file
            with open(save_string_maxwell.format(self.__container_rad, self.__n, self.__v_0, self.__ball_rad), mode = 'w') as f_maxwell:
                f_maxwell.write('Rad:{}, Number:{}, v_0:{}\n'.format(self.__container_rad, self.__n, self.__v_0))
                for i in range(0,5):
                    for el in data_maxwell[i]:
                        f_maxwell.write('{} '.format(el))
                    f_maxwell.write('\n')
                    
        if self.__animate:
            plt.show()
            

