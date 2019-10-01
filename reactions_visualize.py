#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 07:03:37 2019

@author: jamie
"""
import simulation as sim

numbers = [20,20,0,0,10,10]
a = sim.reaction_Simulation(ball_radius= 0.3, container_radius = 10, v_0 = 1000, animate = True, numbers = numbers)
a.run(num_frames = 500, collect_data = False)   