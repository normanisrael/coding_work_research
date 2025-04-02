#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 08:57:05 2022

@author: nanophotonics
"""

import importlib.util
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random 
import os

class LUM:
    def __init__(self, fdtd):
        self.fdtd = fdtd
    #fdtd = lumapi.FDTD(hide = True)
    #fdtd = lumapi.FDTD()
    #fdtd.addfdtd()
    #input_var = [label, material, isidx, ux, uxspan, uxmin, uxmax, uy, uyspan,uymin, uymax, uz, uzspan, uzmin, uzmax]
    def block(self, fdtd, input_var):
        #self.fdtd = fdtd
        fdtd.addrect(name=input_var[0])
        #Geometry
        fdtd.set("x", input_var[2])
        fdtd.set("x span", input_var[3])
        fdtd.set("x min", input_var[4])
        fdtd.set("x max", input_var[5])
        
        fdtd.set("y", input_var[6])
        fdtd.set("y span", input_var[7])
        fdtd.set("y min", input_var[8])
        fdtd.set("y max", input_var[9])
        
        fdtd.set("z", input_var[10])
        fdtd.set("z span", input_var[11])
        fdtd.set("z min", input_var[12])
        fdtd.set("z max", input_var[13])
        
        if type(input_var[1]) == str:
            lumerical_material = 1
            
        elif type(input_var[1]) != str:
            lumerical_material = 0
            
        if lumerical_material == 1:
            #Material
            fdtd.set("material", input_var[1])
    
        if lumerical_material == 0:
            #Material
            fdtd.set("index", input_var[1])
            
    
    """
    FDTD
    """
    #inputs = [dimension, xin, xspan, xmin, xmax, yin, yspan, ymin, ymax, zin, mesh type, mesh refinement, x-step size, y-step size
    # xmin bc, xmax bc, ymin bc, ymax bc, auto shutoff min]
    def simulation_space(self, fdtd, inputs):
        #General
        fdtd.set("dimension", inputs[0])
        
        #Geometry
        fdtd.set("x", inputs[1])
        fdtd.set("x span", inputs[2])
        fdtd.set("x min", inputs[3])
        fdtd.set("x max", inputs[4])
        
        fdtd.set("y", inputs[5])
        fdtd.set("y span", inputs[6])
        fdtd.set("y min", inputs[7])
        fdtd.set("y max", inputs[8])
        
        fdtd.set("z", inputs[9])
        
        #Mesh settings
        fdtd.set("mesh type", inputs[10])
        fdtd.set("mesh refinement", inputs[11])
        
        fdtd.set("dx", inputs[12])
        fdtd.set("dy", inputs[13])
        
        #Boundary conditions
        fdtd.set("x min bc", inputs[14])
        fdtd.set("x max bc", inputs[15])
        fdtd.set("y min bc", inputs[16])
        fdtd.set("y max bc", inputs[17])
        
        #Advance options
        fdtd.set("auto shutoff min", inputs[18])
    
    def source(self, fdtd, source_input):
        fdtd.addmode()
        
        #General
        fdtd.set("injection axis", source_input[0])
        fdtd.set("direction", source_input[1])
        fdtd.set("mode selection", source_input[2])
        fdtd.set("number of trial modes", source_input[3])
        fdtd.set("frequency dependent profile", source_input[4]) #1 0
        if source_input[4] == 1:
            fdtd.set("number of field profile samples", source_input[5]) #10 1000
        
        #Geometry
        fdtd.set("x", source_input[6])
        
        fdtd.set("y", source_input[7])
        fdtd.set("y span", source_input[8])
        fdtd.set("y min", source_input[9])
        fdtd.set("y max", source_input[10])
        
        fdtd.set("z", source_input[11])
        fdtd.set("z span", source_input[12])
        fdtd.set("z min", source_input[13])
        fdtd.set("z max", source_input[14])
        
        #Frequency/Wavelength
        fdtd.set("wavelength start", source_input[15])
        fdtd.set("wavelength stop", source_input[16])
        fdtd.setglobalmonitor('frequency points', source_input[17])
        fdtd.setglobalsource('wavelength stop', source_input[18])
    
    def transmissionX(self, fdtd, transX_input):
        fdtd.addpower() 
        #fdtd.nonorm
        
        #Geometry
        fdtd.set("name", transX_input[0])
        fdtd.set("monitor type", transX_input[1])
        
        fdtd.set("x", transX_input[2])
        
        fdtd.set("y", transX_input[3])
        fdtd.set("y span", transX_input[4])
        fdtd.set("y min", transX_input[5])
        fdtd.set("y max", transX_input[6])
        
        fdtd.set("z", transX_input[7])
        
    def transmissionY(self, fdtd, transY_input):
        fdtd.addpower() 
        
        #Geometry
        fdtd.set("name", transY_input[0])
        fdtd.set("monitor type", transY_input[1])
        
        fdtd.set("x", transY_input[2])
        fdtd.set("x span", transY_input[3])
        fdtd.set("x min", transY_input[4])
        fdtd.set("x max", transY_input[5])
        
        fdtd.set("y", transY_input[6])
        
        fdtd.set("z", transY_input[7])
    
    def movie(self, fdtd, movie_input):
        fdtd.addmovie() #movie monitor 
        
        #z-normal
        #Geometry
        fdtd.set("name", movie_input[0])
        fdtd.set("monitor type", movie_input[1])
        
        fdtd.set("x", movie_input[2])
        fdtd.set("x span", movie_input[3])
        fdtd.set("x min", movie_input[4])
        fdtd.set("x max", movie_input[5])
        
        fdtd.set("y", movie_input[6])
        fdtd.set("y span", movie_input[7])
        fdtd.set("y min", movie_input[8])
        fdtd.set("y max", movie_input[9])
        
        fdtd.set("z", movie_input[10])
        
    def lumerical_parameters(self, lum_input, fpath, fname):
        input_var = ['label', 'material', 'x', 'xspan', 'xmin', 'xmax', 'y', 'yspan', 'ymin', 'ymax', 'z', 'zspan', 'zmin', 'zmax']
        block_info = pd.DataFrame.from_dict(lum_input, orient='index', columns=input_var)
        print(block_info.to_markdown())
        block_info.columns = [column_header.ljust(10, ' ') for column_header in block_info.columns]
        block_info.to_csv(fpath+fname, float_format='%.3e', header = block_info.columns, index=None, sep='\t')
        
    def output_parameters(self, lum_input, fpath, fname):
        input_var = ['label', 'material', 'xspan', 'yspan', 'L1', 'L2', 'L3', 'L4', 'Lam', 'P', 'R', 'Td', 'Tf', 'Tu']
        block_info = pd.DataFrame.from_dict(lum_input, orient='index', columns=input_var)
        print(block_info.to_markdown())
        block_info.columns = [column_header.ljust(8, ' ') for column_header in block_info.columns]
        block_info.to_csv(fpath+fname, float_format='%.3e ', header = block_info.columns, index=None, sep='\t')
    
    def results(self, fdtd, labels):
        Monitor_R = fdtd.getresult(labels[0], "T")
        Monitor_T = fdtd.getresult(labels[1], "T")
        Monitor_Tup = fdtd.getresult(labels[2], "T")
        Monitor_Tdown = fdtd.getresult(labels[3], "T")
        the_source = fdtd.getresult("source","mode profile")
        #fdtd.close()
        
        wavelength = Monitor_T['lambda']
        transmit = Monitor_T['T']
        transmissionUp = Monitor_Tup['T']
        transmissionDown = Monitor_Tdown['T']
        reflect = Monitor_R['T']
        total = transmit - reflect + transmissionUp - transmissionDown
        error = abs(1.0-total)
        return wavelength, error, transmit, transmissionUp, transmissionDown, reflect
    '''
    def results(self, fdtd, labels, positions): #[source_location, reflection_monitor_pos, transmission_pos, trans_up_pos, trans_down_pos, min_pos_r, max_pos_t, min_pos_tdown, max_pos_tup]
        Monitor_R = fdtd.getresult(labels[0], "E")
        Monitor_T = fdtd.getresult(labels[1], "E")
        Monitor_Tup = fdtd.getresult(labels[2], "E")
        Monitor_Tdown = fdtd.getresult(labels[3], "E")
        
        frequencies = Monitor_R['f'] 
        
        wavelength = np.ravel(Monitor_T['lambda'])
        reflect = Monitor_R['E']
        transmit = Monitor_T['E']
        transmissionUp = Monitor_Tup['E']
        transmissionDown = Monitor_Tdown['E']
        
        source_normalization = np.ravel(fdtd.sourcenorm(frequencies)) 
        source_normalization_rescale = 1.0 #source_normalization*np.exp(1j*2*np.pi*np.abs(positions[0] - positions[5])/wavelength)
        #source_normalization_rescale = source_normalization_rescale/np.max(source_normalization_rescale)#*np.exp(1j*2*np.pi*np.abs(positions[0] - positions[5])/wavelength)
        #plt.plot(frequencies, abs(source_normalization_rescale))
        #print(1.0/frequencies)
       # reflection = reflect[0,:,0,:,1],axi/source_normalization_rescale#*np.exp(-1j*2*np.pi*np.abs(positions[1] - positions[5])/wavelength)
        transmission_loc = int((transmit[0,:,0,:,1].shape[0]-1)/2)
        transmission_loc_up = int((transmissionUp[0,:,0,:,1].shape[0]-1)/2)
        transmission_loc_down = int((transmissionDown[0,:,0,:,1].shape[0]-1)/2)
        print(transmit[0,:,0,:,1].shape)
        print(transmission_loc)
        print(transmit[0,:,0,:,1])
        transmission = (transmit[0,transmission_loc,0,:,1]/source_normalization_rescale)#*np.exp(-1j*2*np.pi*np.abs(positions[2] - positions[6])/wavelength)
        transmission_up = (transmissionUp[0,transmission_loc_up,0,:,1])/source_normalization_rescale#*np.exp(-1j*2*np.pi*np.abs(positions[3] - positions[7])/wavelength)
        transmission_down = (np.median(transmissionDown[0,:,0,:,1],axis = 0)/source_normalization_rescale)#*np.exp(-1j*2*np.pi*np.abs(positions[4] - positions[8])/wavelength)
       
        
        #print(reflect[0,:,0,:,1])
       # total = transmission - reflection + transmission_up - transmission_down
        print(np.abs(transmission)+np.abs(transmission_up))
        plt.plot(wavelength,np.abs(transmission)+np.abs(transmission_up))
        #error = abs(1.0-total)
        return wavelength, transmission, transmission_up, transmission_down#, reflection
        #return wavelength, error, transmit, transmissionUp, transmissionDown, reflect
    '''
    def error_results(self, fdtd, labels):
        Monitor_R = fdtd.getresult(labels[0], "T")
        Monitor_T = fdtd.getresult(labels[1], "T")
        Monitor_Tup = fdtd.getresult(labels[2], "T")
        Monitor_Tdown = fdtd.getresult(labels[3], "T")
        the_source = fdtd.getresult("source","mode profile")
        
        wavelength = Monitor_T['lambda']
        transmit = Monitor_T['T']
        transmissionUp = Monitor_Tup['T']
        transmissionDown = Monitor_Tdown['T']
        reflect = Monitor_R['T']
        total = transmit - reflect + transmissionUp - transmissionDown
        error = abs(1.0-total)
        return wavelength, error
    
    def layout(self, fdtd, par):
        fdtd.switchtolayout()
        fdtd.setnamed("FDTD","dx",par)
        fdtd.setnamed("FDTD","dy",par)
        
    def save_sim_file(self, fdtd, fpath, fname):
        fdtd.save(fpath+fname)
        fdtd.run()
    
    def save_data(self, data_array, fpath, fname):
        np.savetxt(fpath+fname, np.transpose(data_array),fmt='%1.6e',delimiter='  ')
        
    def save_data_one_by_one(self, data_array, fpath, data_index, fname):
        print("Loop: ", data_index)
        with open(fpath+fname, 'ab') as f:
            np.savetxt(f, data_array,fmt='%1.6e',delimiter='  ')
    
    def read_data(self, fpath, fname, row_skip):
        file = np.loadtxt(fpath+fname, delimiter='\t', unpack=True, skiprows=row_skip) #pd.read_csv(fpath+fname, sep=" ") #np.loadtxt(fpath+fname, delimiter='  ', unpack=True, skiprows=1)
        return file
    
    def number_of_lines(self, path, name):
        with open(path+name, 'r') as fp:
            lines = len(fp.readlines())
        return lines
    
    def parameters(self, fpath, fname, nunits, npar):
        titles = ['L1', 'L2', 'L3', 'L4', 'nswg', 'Si length', 'P', 'Lam']
        def gen_par(no_units):
            Lam = np.random.uniform(400,1000)
            pick = [0]
            nswg = np.random.uniform(1.6,3.0)
            Si_Length = 100.0 
            while True:
                parameters = np.random.uniform(40,Lam,no_units-1)
                parameters = list(parameters)
                parameters[0] = 0.0
                Lf = Lam - np.sum(parameters)
                if Lf > 40.0:
                    parameters.append(Lf)
                    break
            
            parameters.append(nswg)
            parameters.append(Si_Length)
            parameters.append(random.choice(pick))
            parameters.append(Lam)
            return parameters
        
        inputs = []
       
        for idx in range(0, npar):
            vals = gen_par(nunits)
            inputs.append(vals)
        
        df = pd.DataFrame(inputs, columns = titles)
        df.columns = [column_header.ljust(10, ' ') for column_header in df.columns]
        df.to_csv(fpath+fname, float_format='%.5e', header = df.columns, index=None, sep='\t')
    
    def find_files(self, filename, search_path):
        result = []
        
        for root, dir, files in os.walk(search_path):
            if filename in files:
                result.append(os.path.join(root, filename))
        return result 
    
    def plot(self, min_val, total, blocks, vis):
        max_blocks = 2*total 
        points = {
                }
        for j in range(min_val, max_blocks):
            points.update({str(j):[[blocks[str(j)][4], blocks[str(j)][8]], [blocks[str(j)][5], blocks[str(j)][8]], 
                       [blocks[str(j)][5], blocks[str(j)][9]],[blocks[str(j)][4], blocks[str(j)][9]]]})
        
        #plots 
        for i in range(min_val, max_blocks):
            fc = 'r'
            if i == vis:
                fc = 'g'
                
            if i == 0 or i == 29:
                fc = 'b'
            polygon = plt.Polygon(points[str(i)],fc=fc)
            plt.gca().add_patch(polygon)
        
        plt.axis('scaled')
        plt.show()    
