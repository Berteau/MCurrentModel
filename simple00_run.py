# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 15:58:58 2019

@author: yohan
"""

from simple00_Simulation import *

params = {}

params["maxTime"] = 400 #milliseconds
params["tau"] = 0.1 

params["popCount"] = 10 #Number of neurons

sim = simple00_Simulation(params)

sim.run()

sim.plotColumns()
