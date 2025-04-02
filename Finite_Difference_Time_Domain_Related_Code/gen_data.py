#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 21:55:34 2022

@author: nanophotonics
"""
import numpy as np
import random
import functions as f2d
import importlib.util


spec_lin = importlib.util.spec_from_file_location('lumapi', "/opt/lumerical/v231/api/python/lumapi.py")
 
lumapi = importlib.util.module_from_spec(spec_lin) 
spec_lin.loader.exec_module(lumapi)

fdtd = lumapi.FDTD(hide = True)

fdtd.addfdtd()
FDTD = f2d.LUM(fdtd)

lum_para_file = "/simulations/2D/parameters/input/"
file = 'the_parameter_file.dat'

FDTD.parameters(lum_para_file, file, 4, 4000)
