#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 12:33:08 2019

@author: jamie
"""
import numpy as np

class reactant:
    
    _number_A = 0
    _number_B = 0
    _number_C = 0
    _number_D = 0
    _number_E = 0
    _number_F = 0
    _regular_activation_energy = 80000
    _catalyst_activation_energy = 20000
    
    def get_values():
        return [reactant._number_A, reactant._number_B, reactant._number_C, reactant._number_D, reactant._number_E, reactant._number_F]
    
    def reset_values():
        reactant._number_A = 0
        reactant._number_B = 0
        reactant._number_C = 0
        reactant._number_D = 0
        reactant._number_E = 0
        reactant._number_F = 0
        
    def __init__(self, reactant_type = 'reactant_A', animate = True):
        self.__animate = animate
        self.__type = reactant_type
        if self.__type == 'reactant_A':
            reactant._number_A += 1
            if animate:
                self.__patch_color = 'r'
        elif self.__type == 'reactant_B':
            reactant._number_B += 1
            if animate:
                self.__patch_color = 'b'
        elif self.__type == 'product_C':
            reactant._number_C += 1
            if animate:
                self.__patch_color = 'g'
        elif self.__type == 'product_D':
            reactant._number_D += 1
            if animate:
                self.__patch_color = 'y'
        elif self.__type == 'catalyst_E':
            reactant._number_E += 1
            if animate:
                self.__patch_color = 'k'
        elif self.__type == 'catalyst_F':
            reactant._number_F += 1
            if animate:
                self.__patch_color = 'w'
        elif self.__type == 'container':
            if animate:
                self.__patch_color = 'm'
    
    def get_patch_color(self):
        return self.__patch_color
        
    def get_type(self):
        return self.__type
        
    def regular_energy_check(self,other):
        if np.dot(self.vel(), other.vel()) < -reactant._regular_activation_energy:
            return True
        else:
            return False
        
    def catalyst_energy_check(self,other):
        if np.dot(self.vel(), other.vel()) < -reactant._catalyst_activation_energy:
            return True
        else:
            return False
            
    def reaction_A_B(self, other):
        self.__type = 'product_C'
        other.__type = 'product_D'
        if self.__animate:
            self.__patch_color = 'g'
            other.__patch_color = 'y'
        reactant._number_A -= 1
        reactant._number_B -= 1
        reactant._number_C += 1
        reactant._number_D += 1
        
    def reaction_C_D(self, other):
        self.__type = 'reactant_A'
        other.__type = 'reactant_B'
        if self.__animate:
            self.__patch_color = 'r'
            other.__patch_color = 'b'
        reactant._number_C -= 1
        reactant._number_D -= 1
        reactant._number_A += 1
        reactant._number_B += 1
            
    def reaction_A_E(self,other):
        self.__type = 'product_C'
        other.__type = 'catalyst_F'
        if self.__animate:
            self.__patch_color = 'g'
            other.__patch_color = 'w'
        reactant._number_A -= 1
        reactant._number_E -= 1
        reactant._number_C += 1
        reactant._number_F += 1
        
    def reaction_B_F(self,other):
        self.__type = 'product_D'
        other.__type = 'catalyst_E'
        if self.__animate:
            self.__patch_color = 'y'
            other.__patch_color = 'k'
        reactant._number_B -= 1
        reactant._number_F -= 1
        reactant._number_D += 1
        reactant._number_E += 1
                
    
    def react(self, other):
        if other.__type == 'container':
            return
        if self.__type == 'reactant_A':
            if other.__type == 'reactant_B':
                if self.regular_energy_check(other):
                    self.reaction_A_B(other)
            elif other.__type == 'catalyst_E':
                if self.catalyst_energy_check(other):
                    self.reaction_A_E(other)
                    
        elif self.__type == 'reactant_B':
            if other.__type == 'reactant_A':
                if self.regular_energy_check(other):
                    other.reaction_A_B(self)
            elif other.__type == 'catalyst_F':
                if self.catalyst_energy_check(other):
                    self.reaction_B_F(other)
            
        elif self.__type == 'product_C':
            if other.__type == 'product_D':
                if self.regular_energy_check(other):
                    self.reaction_C_D(other)
                    
        elif self.__type == 'product_D':
            if other.__type == 'product_C':
                if self.regular_energy_check(other):
                    other.reaction_C_D(self)
            
        elif self.__type == 'catalyst_E':
            if other.__type == 'reactant_A':
                if self.catalyst_energy_check(other):
                    other.reaction_A_E(self)
        
        elif self.__type == 'catalyst_F':
            if other.__type == 'reactant_B':
                if self.catalyst_energy_check(other):
                    other.reaction_B_F(self)
            
                    
            
                    


