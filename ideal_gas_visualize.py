#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 07:00:57 2019

@author: jamie
"""

import simulation as sim

numbers = [100,0,0,0,0,0]
a = sim.reaction_Simulation(ball_radius= 0.1, container_radius = 10, v_0 = 500, animate = True, numbers = numbers)
a.run(num_frames = 1000, collect_data = False)   