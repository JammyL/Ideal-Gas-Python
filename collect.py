#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 07:07:28 2019

@author: jamie
"""

import simulation as sim

numbers = [500,500,0,0,50,50]
a = sim.reaction_Simulation(ball_radius= 0.01, container_radius = 10, v_0 = 200, animate = False, numbers = numbers, save_location = 'example_data/')
a.run(num_frames = 300, collect_data = True)   