#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 09:04:08 2022

@author: nanophotonics
"""

import numpy as np
import functions as sf
import importlib.util
#for step in range(40,44,2):
    
spec_lin = importlib.util.spec_from_file_location('lumapi', "/opt/lumerical/v231/api/python/lumapi.py")
 
#Functions that perform the actual loading
lumapi = importlib.util.module_from_spec(spec_lin) # 
spec_lin.loader.exec_module(lumapi)

fdtd = lumapi.FDTD(hide = False)
#fdtd = lumapi.FDTD()

#fdtd.addfdtd()
FDTD = sf.LUM(fdtd)

"""
FILE PATHS
"""
sim_file = "/home/nanophotonics/Research/Project_3swg/simulations/2D/movies/"
data_file = "/home/nanophotonics/Research/Project_3swg/simulations/2D/parameters/output/"
lum_para_file = "/home/nanophotonics/Research/Project_3swg/simulations/2D/parameters/"
input_para = "/home/nanophotonics/Research/Project_3swg/simulations/2D/parameters/input/"

#file = "ANT1SWG_pList_test_2.dat"
file = "ANT1SWG_pList_jun_6_comparison.dat"
output_file = "ANT1SWG_pList_jun_6_comparison_output_1.dat"
input_data = FDTD.read_data(input_para, file,1)

find_the_output_file = FDTD.find_files(output_file, data_file)
#read_output = 
#print(find_the_output_file)
if len(find_the_output_file) == 0:
    start = 0
    
else:
    no_of_lines = FDTD.number_of_lines(data_file, output_file)
    start = no_of_lines
    
full_output = []
#print(start)
stop = len(input_data[0,:])

total_no_parameters = len(input_data[0])-1

full_output = []

for idx in range(0,1):
    fdtd = lumapi.FDTD(hide = False)
    fdtd.addfdtd()
    FDTD = sf.LUM(fdtd)
    data = input_data[:,idx]
    """
    UNITS
    """
    um = 1e-6
    nm = 1e-9
     
    FULL = 1 #etch type, 0 is partial, 1 is full
     
    """
    CONSTANTS
    """
    etch = {220:[110, 70, 220],
            300:[110, 300]
           }
     
    etch_depth = etch[220][0]*nm
    Si_length = data[5]*nm
    center_height = [220.0*nm, 300.0*nm]
    bottom_height = center_height[0] - etch_depth
     
    """
    Polarization
    """
    Pmode = data[6]
    if Pmode == 0:
        P = "fundamental TM mode"
        
    elif Pmode == 1:
        P = "fundamental TE mode"
     
    """
    OTHER VARIABLES
    """
    x0 = 0.0*nm
    y0 = 0.0*nm
     
    nswg = data[4]
     
    L = data[0:4]*nm
     
    lam = np.sum(L)
    #print(lam)
    LTotal = 3*lam + 2*Si_length
    HSiO2 = 1055*nm
     
    blockHeight = etch[220][0]*nm
    ycmin = -0.25*center_height[0]
    ycmax = 0.75*center_height[0]
     
    #start here
    Lymin = ycmin
    Lymax = ycmin + bottom_height
    Uymin = Lymax 
    Uymax = ycmax
     
    USiO2ymin = ycmax
    USiO2ymax = ycmax + HSiO2
     
    LSiO2ymin = ycmin - HSiO2
    LSiO2ymax = Lymin
     
    shiftx = 0.0*nm #LTotal/2.0
    shifty = -50.0*nm
     
    zc = 0.0*nm
    zcspan = 10.0*nm
    zcmin = -5.0*nm
    zcmax = 5.0*nm
     
     
    """
    BLOCK XMIN & YMIN
    """
    bxmin0 = -0.5*Si_length - shiftx
    bxmin1 = -0.5*Si_length - shiftx
    bxmin2 = 0.5*Si_length - shiftx
    bxmin3 = 0.5*Si_length - shiftx + L[0]
    bxmin4 = 0.5*Si_length - shiftx + L[0] + L[1]
    bxmin5 = 0.5*Si_length - shiftx + L[0] + L[1] + L[2]
    bxmin6 = 0.5*Si_length - shiftx + L[0] + L[1] + L[2] + L[3]
    bxmin7 = 0.5*Si_length - shiftx + 2*L[0] + L[1] + L[2] + L[3]
    bxmin8 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + L[2] + L[3]
    bxmin9 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + 2*L[2] + L[3]
    bxmin10 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + 2*L[2] + 2*L[3]
    bxmin11 = 0.5*Si_length - shiftx + 3*L[0] + 2*L[1] + 2*L[2] + 2*L[3]
    bxmin12 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 2*L[2] + 2*L[3]
    bxmin13 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 3*L[2] + 2*L[3]
    bxmin14 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 3*L[2] + 3*L[3]
     
    Ubymin = Uymin - shifty
    Lbymin = Lymin - shifty
     
    USiO2bymin = USiO2ymin - shifty
    LSiO2bymin = LSiO2ymin - shifty 
     
     
    """
    BLOCK XMAX
    """
    bxmax0 = -0.5*Si_length - shiftx + LTotal
    bxmax1 = 0.5*Si_length - shiftx
    bxmax2 = 0.5*Si_length - shiftx + L[0]
    bxmax3 = 0.5*Si_length - shiftx + L[0] + L[1]
    bxmax4 = 0.5*Si_length - shiftx + L[0] + L[1] + L[2]
    bxmax5 = 0.5*Si_length - shiftx + L[0] + L[1] + L[2] + L[3]
    bxmax6 = 0.5*Si_length - shiftx + 2*L[0] + L[1] + L[2] + L[3]
    bxmax7 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + L[2] + L[3]
    bxmax8 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + 2*L[2] + L[3]
    bxmax9 = 0.5*Si_length - shiftx + 2*L[0] + 2*L[1] + 2*L[2] + 2*L[3]
    bxmax10 = 0.5*Si_length - shiftx + 3*L[0] + 2*L[1] + 2*L[2] + 2*L[3]
    bxmax11 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 2*L[2] + 2*L[3]
    bxmax12 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 3*L[2] + 2*L[3]
    bxmax13 = 0.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 3*L[2] + 3*L[3]
    bxmax14 = 1.5*Si_length - shiftx + 3*L[0] + 3*L[1] + 3*L[2] + 3*L[3]
     
    Ubymax = Uymax - shifty
    Lbymax = Lymax - shifty
     
    USiO2bymax = USiO2ymax - shifty
    LSiO2bymax = LSiO2ymax - shifty
     
    """
    BLOCK CENTERS
    """
    beginx = x0 - shiftx
     
    bcx0 = 0.5*(bxmax0 + bxmin0) 
    bcx1 = 0.5*(bxmax1 + bxmin1)  
    bcx2 = 0.5*(bxmax2 + bxmin2) 
    bcx3 = 0.5*(bxmax3 + bxmin3) 
    bcx4 = 0.5*(bxmax4 + bxmin4) 
    bcx5 = 0.5*(bxmax5 + bxmin5) 
    bcx6 = 0.5*(bxmax6 + bxmin6) 
    bcx7 = 0.5*(bxmax7 + bxmin7) 
    bcx8 = 0.5*(bxmax8 + bxmin8) 
    bcx9 = 0.5*(bxmax9 + bxmin9) 
    bcx10 = 0.5*(bxmax10 + bxmin10) 
    bcx11 = 0.5*(bxmax11 + bxmin11) 
    bcx12 = 0.5*(bxmax12 + bxmin12) 
    bcx13 = 0.5*(bxmax13 + bxmin13) 
    bcx14 = 0.5*(bxmax0 + bxmin0) 
     
    Ubyc = 0.5*(Ubymax + Ubymin)
    Lbyc = 0.5*(Lbymax + Lbymin)
     
    USiO2byc = 0.5*(USiO2bymax + USiO2bymin)
    LSiO2byc = 0.5*(LSiO2bymax + LSiO2bymin)
     
    nSWG = {
            '0':['SiO2 (Glass) - Palik'],
            '1':['Si (Silicon) - Palik'],
            '2':['Si (Silicon) - Palik'],
            '3':['Si (Silicon) - Palik'],
            '4':['Si (Silicon) - Palik'],
            '5':['Si (Silicon) - Palik'],
            '6':['Si (Silicon) - Palik'],
            '7':['Si (Silicon) - Palik'],
            '8':['Si (Silicon) - Palik'],
            '9':['Si (Silicon) - Palik'],
            '10':['Si (Silicon) - Palik'],
            '11':['Si (Silicon) - Palik'],
            '12':['Si (Silicon) - Palik'],
            '13':['Si (Silicon) - Palik'],
            '14':['Si (Silicon) - Palik'],
            '15':['Si (Silicon) - Palik'],
            '16':['Si (Silicon) - Palik'],
            '17':['Si (Silicon) - Palik'],
            '18':['Si (Silicon) - Palik'],
            '19':['Si (Silicon) - Palik'],
            '20':['Si (Silicon) - Palik'],
            '21':['Si (Silicon) - Palik'],
            '22':['Si (Silicon) - Palik'],
            '23':['Si (Silicon) - Palik'],
            '24':['Si (Silicon) - Palik'],
            '25':['Si (Silicon) - Palik'],
            '26':['Si (Silicon) - Palik'],
            '27':['Si (Silicon) - Palik'],
            '28':['Si (Silicon) - Palik'],
            '29':['SiO2 (Glass) - Palik']
            }
     
    if FULL == 1:
        Lbox = {
                'L1':['2','16','6','20','10','24'],
                'L2':['3','17','7','21','11','25'],
                'L3':['4','18','8','22','12','26'],
                'L4':['5','19','9','23','13','27']
                }
        
        swg_l4 = Lbox['L4'] 
        for swl4 in swg_l4:
            nSWG[swl4][0] = nSWG['1'][0] #Si
            
        swg_l3 = Lbox['L3'] 
        for swl3 in swg_l3:
            nSWG[swl3][0] = nSWG['0'][0] #SiO2
            
        swg_l2 = Lbox['L2'] 
        for swl2 in swg_l2:
            nSWG[swl2][0] = nswg #nswg
        
        swg_l1 = Lbox['L1'] 
        for swl1 in swg_l1:
            nSWG[swl1][0] = nSWG['0'][0] #Si
            
    if FULL == 0:
        Lbox = {
                'L1':['2','6','10'],
                'L2':['3','7','11'],
                'L3':['4','8','12'],
                'L4':['5','9','13']
                }
        
        change = Lbox['L2'] 
        for ch in change:
            nSWG[ch][0] = nSWG['0'][0]
        
        swg = Lbox['L1'] 
        for sw in swg:
            nSWG[sw][0] = nswg
    #print(nSWG['3'][0])
     
    blockdata = {
                 '0':['TOP', nSWG['0'][0], bcx0, LTotal, bxmin0, bxmax0, USiO2byc, HSiO2, USiO2bymin, USiO2bymax, zc, zcspan, zcmin, zcmax],
                 '1':['1', nSWG['1'][0], bcx1, Si_length, bxmin1, bxmax1, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '2':['2', nSWG['2'][0], bcx2, L[0], bxmin2, bxmax2, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '3':['3', nSWG['3'][0], bcx3, L[1], bxmin3, bxmax3, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '4':['4', nSWG['4'][0], bcx4, L[2], bxmin4, bxmax4, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '5':['5', nSWG['5'][0], bcx5, L[3], bxmin5, bxmax5, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '6':['6', nSWG['6'][0], bcx6, L[0], bxmin6, bxmax6, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '7':['7', nSWG['7'][0], bcx7, L[1], bxmin7, bxmax7, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '8':['8', nSWG['8'][0], bcx8, L[2], bxmin8, bxmax8, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '9':['9', nSWG['9'][0], bcx9, L[3], bxmin9, bxmax9, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '10':['10', nSWG['10'][0], bcx10, L[0], bxmin10, bxmax10, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '11':['11', nSWG['11'][0], bcx11, L[1], bxmin11, bxmax11, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '12':['12', nSWG['12'][0], bcx12, L[2], bxmin12, bxmax12, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '13':['13', nSWG['13'][0], bcx13, L[3], bxmin13, bxmax13, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '14':['14', nSWG['14'][0], bcx14, Si_length, bxmin14, bxmax14, Ubyc, blockHeight, Ubymin, Ubymax, zc, zcspan, zcmin, zcmax],
                 '15':['15', nSWG['15'][0], bcx1, Si_length, bxmin1, bxmax1, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '16':['16', nSWG['16'][0], bcx2, L[0], bxmin2, bxmax2, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '17':['17', nSWG['17'][0], bcx3, L[1], bxmin3, bxmax3, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '18':['18', nSWG['18'][0], bcx4, L[2], bxmin4, bxmax4, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '19':['19', nSWG['19'][0], bcx5, L[3], bxmin5, bxmax5, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '20':['20', nSWG['20'][0], bcx6, L[0], bxmin6, bxmax6, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '21':['21', nSWG['21'][0], bcx7, L[1], bxmin7, bxmax7, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '22':['22', nSWG['22'][0], bcx8, L[2], bxmin8, bxmax8, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '23':['23', nSWG['23'][0], bcx9, L[3], bxmin9, bxmax9, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '24':['24', nSWG['24'][0], bcx10, L[0], bxmin10, bxmax10, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '25':['25', nSWG['25'][0], bcx11, L[1], bxmin11, bxmax11, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '26':['26', nSWG['26'][0], bcx12, L[2], bxmin12, bxmax12, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '27':['27', nSWG['27'][0], bcx13, L[3], bxmin13, bxmax13, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '28':['28', nSWG['28'][0], bcx14, Si_length, bxmin14, bxmax14, Lbyc, blockHeight, Lbymin, Lbymax, zc, zcspan, zcmin, zcmax],
                '29':['BOTTOM', nSWG['29'][0], bcx14, LTotal, bxmin0, bxmax0, LSiO2byc, HSiO2, LSiO2bymin, LSiO2bymax, zc, zcspan, zcmin, zcmax]
                }
    #blockdata['3'][1] = nswg
    #print(blockdata['3'])
     
    """
    SIMULATION PARAMETER LIST
    """
    dimension = "2D"
    #sim_xspan = LTotal
    sim_xmin = blockdata['0'][4]
    sim_xmax = blockdata['0'][5]
    sim_xspan = sim_xmax - sim_xmin
    xcsim = 0.5*(sim_xmax + sim_xmin)
     
    #+ shifty
    #sim_yspan = 2*HSiO2 + blockdata['1'][7] + blockdata['15'][7]
    sim_ymin = blockdata['29'][8] #+ shifty
    sim_ymax = blockdata['0'][9] #+ shifty
    sim_yspan = sim_ymax - sim_ymin
    ycsim = 0.5*(sim_ymax + sim_ymin)
     
    zcsim = 0.0
     
    mesh_type = "uniform"
    mesh_refinement = "conformal variant 0"
     
    xstep = 12*nm
    ystep = 12*nm
     
    xminbc = "PML"
    xmaxbc = "PML"
    yminbc = "PML"
    ymaxbc = "PML"
     
    shutoffmin = 1e-11
     
     
    """
    SOURCE PARAMETER LIST
    """
    inject = "x-axis"
    direction = "Forward"
    mode = P
    no_trial_mode = 1000
    freq_dep_prof = 1 #100
    no_field_prof_samp = 10 #1000
     
    xcsrce = 0.75*blockdata['1'][3] + blockdata['1'][4]
     
    #ycsrce = y0 - shifty
    #srce_yspan = 2*HSiO2 + blockdata['1'][7] + blockdata['15'][7]
    srce_ymin = blockdata['29'][8] - shifty - bottom_height
    srce_ymax = blockdata['0'][9] - shifty
    srce_yspan = srce_ymax - srce_ymin
    ycsrce = 0.5*(srce_ymax + srce_ymin)
     
    srce_zc = zc
    srce_zspan = zcspan
    srce_zmin = zcmin 
    srce_zmax = zcmax 
     
    wvelen_start = 1450.0*nm
    wvelen_stop = 1650.0*nm
     
    globfreq_points = 100
    globwvelen_stop = 2000.0*nm
     
     
    """
    TRANSMISSION MONITOR PARAMETERS
    """
    #REFLECTION
    rfl_label = "Reflection Monitor"
    rfl_mon_type = "Linear Y"
     
    xcrfl = 0.5*blockdata['1'][3] + blockdata['1'][4]
     
    #ycrfl = blockdata['1'][6] 
    rfl_yspan = 2*HSiO2 + blockdata['1'][7] + blockdata['15'][7] - 60*nm
    rfl_ymin = blockdata['29'][8] + 30*nm 
    rfl_ymax = blockdata['0'][9] - 30*nm 
    ycrfl = 0.5*(rfl_ymax + rfl_ymin)
     
    zcrfl = 0.0
     
    #TRANSMISSION
    tran_label = "Transmission Monitor"
    tran_mon_type = "Linear Y"
     
    xctran = 0.5*blockdata['14'][3] + blockdata['14'][4]
     
    #yctran = blockdata['14'][6] - shifty
    tran_yspan = 2*HSiO2 + blockdata['14'][7] + blockdata['28'][7] - 60*nm
    tran_ymin = blockdata['29'][8] + 30*nm 
    tran_ymax = blockdata['0'][9] - 30*nm 
    yctran = 0.5*(tran_ymin + tran_ymax)
     
    zctran = 0.0
     
    #UPWARD TRANSMISSION
    tran_up_label = "Transmission Up Monitor"
    tran_up_mon_type = "Linear X"
     
    tran_up_xspan = 0.5*blockdata['14'][3] + blockdata['14'][4] - (0.5*blockdata['1'][3] + blockdata['1'][4])
    tran_up_xmin = 0.5*blockdata['1'][3] + blockdata['1'][4]
    tran_up_xmax = 0.5*blockdata['14'][3] + blockdata['14'][4] 
    xctranup = 0.5*(tran_up_xmax + tran_up_xmin)
     
    yctranup = blockdata['0'][9] - 30*nm 
    '''
    tran_up_yspan = 0.5*blockdata['14'][3] + blockdata['14'][4] - (0.5*blockdata['1'][3] + blockdata['1'][4])
    tran_up_ymin = 0.5*blockdata['1'][3] + blockdata['1'][4]
    tran_up_ymax = 0.5*blockdata['14'][3] + blockdata['14'][4] 
    '''
    #yctranup = 0.5*(tran_up_ymax + tran_up_ymin)
     
    zctranup = 0.0
     
    #DOWNWARD TRANSMISSION
    tran_down_label = "Transmission Down Monitor"
    tran_down_mon_type = "Linear X"
     
    tran_down_xspan = 0.5*blockdata['14'][3] + blockdata['14'][4] - (0.5*blockdata['1'][3] + blockdata['1'][4])
    tran_down_xmin = 0.5*blockdata['1'][3] + blockdata['1'][4] 
    tran_down_xmax = 0.5*blockdata['14'][3] + blockdata['14'][4] 
    xctrandown = 0.5*(tran_down_xmax + tran_down_xmin)
     
    yctrandown = blockdata['29'][8] + 30*nm 
    '''
    tran_down_yspan = 0.5*blockdata['14'][3] + blockdata['14'][4] - (0.5*blockdata['1'][3] + blockdata['1'][4])
    tran_down_ymin = 0.5*blockdata['1'][3] + blockdata['1'][4] 
    tran_down_ymax = 0.5*blockdata['14'][3] + blockdata['14'][4] 
    '''
    #yctrandown = 0.5*(tran_down_ymax + tran_down_ymin) 
     
    zctrandown = 0.0
     
    #MOVIE 
    mv_label = "Movie Monitor"
    mv_mon_type = "2D Z-normal"
     
    #xcmv = 0.5*(0.5*blockdata['14'][3] + blockdata['14'][4] - (0.5*blockdata['1'][3] + blockdata['1'][4]))
    mv_xspan = LTotal
    mv_xmin = blockdata['0'][4]
    mv_xmax = blockdata['0'][5]
    xcmv = 0.5*(mv_xmax + mv_xmin)
     
    #ycmv = 0.0 
    mv_yspan = 2*HSiO2 + blockdata['1'][7] + blockdata['15'][7]
    mv_ymin = blockdata['29'][8] 
    mv_ymax = blockdata['0'][9] 
    ycmv = 0.5*(mv_ymax + mv_ymin)
     
    zcmv = 0.0
     
     
    """
    PARAMETER LISTS AND MONITOR LABELS
    """
    sim_parameters = [dimension, xcsim, sim_xspan, sim_xmin, sim_xmax, ycsim, sim_yspan, sim_ymin, sim_ymax, zcsim, mesh_type, 
                      mesh_refinement, xstep, ystep, xminbc, xmaxbc, yminbc, ymaxbc, shutoffmin]
    source_parameters = [inject, direction, mode, no_trial_mode, freq_dep_prof, no_field_prof_samp, xcsrce, ycsrce, srce_yspan, srce_ymin, 
                         srce_ymax, srce_zc, srce_zspan, srce_zmin, srce_zmax, wvelen_start, wvelen_stop, globfreq_points, globwvelen_stop]
    reflection_parameters = [rfl_label, rfl_mon_type, xcrfl, ycrfl, rfl_yspan, rfl_ymin, rfl_ymax, zcrfl]
    transmission_parameters = [tran_label, tran_mon_type, xctran, yctran, tran_yspan, tran_ymin, tran_ymax, zctran]
    transmission_up_parameters = [tran_up_label, tran_up_mon_type, xctranup, tran_up_xspan, tran_up_xmin, tran_up_xmax, yctranup, zctranup]
    transmission_down_parameters = [tran_down_label, tran_down_mon_type, xctrandown, tran_down_xspan, tran_down_xmin, 
                                    tran_down_xmax, yctrandown, zctrandown]
    movie_parameters = [mv_label, mv_mon_type, xcmv, mv_xspan, mv_xmin, mv_xmax, ycmv, mv_yspan, mv_ymin, mv_ymax, zcmv]
     
    #scale_results = [xcsrce, xcrfl, xctran, xctranup, xctrandown, bxmin2, bxmax14, USiO2ymin - shifty, USiO2ymax - shifty]
     
    monitor_labels = ["Reflection Monitor", "Transmission Monitor", "Transmission Up Monitor", "Transmission Down Monitor"]
     
    #end here
    """
    SIMULATION HAPPENS HERE
    """
    def sim(step):
        fdtd.nonorm
        sim_parameters[12] = step*nm
        sim_parameters[13] = sim_parameters[12]
        FDTD.simulation_space(fdtd,sim_parameters)
        for bidx in range(0,len(blockdata)):
            bidx = str(bidx)
            #blockdata['3'][1] = float(nswg)
            FDTD.block(fdtd,blockdata[bidx])
            
        FDTD.source(fdtd,source_parameters)
        FDTD.transmissionX(fdtd,reflection_parameters)
        FDTD.transmissionX(fdtd,transmission_parameters)
        FDTD.transmissionY(fdtd,transmission_up_parameters)
        FDTD.transmissionY(fdtd,transmission_down_parameters)
        FDTD.movie(fdtd,movie_parameters)
        
        FDTD.save_sim_file(fdtd, sim_file, "ANT1SWG_test_TM_nov_22_22_test_purposes_3.fsp")
        
        wavelength, error, transmit, transmissionUp, transmissionDown, reflect = FDTD.results(fdtd,monitor_labels)
        #wavelength, transmit, transmissionUp, transmissionDown, reflect = FDTD.results(fdtd,monitor_labels, scale_results)
        return wavelength, transmit, transmissionUp, transmissionDown, reflect
     
    wavelength, transmit, transmissionUp, transmissionDown, reflect = sim(10)
    #wavelength, _, transmit, transmissionUp, transmissionDown, reflect = sim(10)
     
    entries = [np.squeeze(wavelength), transmit, transmissionUp, transmissionDown, reflect]
    
    outputs = []
    transmit = list(transmit)
    transmissionUp = list(transmissionUp)
    transmissionDown = list(transmissionDown)
    reflect = list(reflect)
    
    outputs.extend(transmit)
    outputs.extend(transmissionUp)
    outputs.extend(transmissionDown)
    outputs.extend(reflect)
    outputs = np.append(outputs, L[0])
    outputs = np.append(outputs, L[1])
    outputs = np.append(outputs, L[2])
    outputs = np.append(outputs, L[3])
    outputs = np.append(outputs, lam)
    outputs = np.append(outputs, nswg)
    outputs = np.append(outputs, center_height[0])
    outputs = np.append(outputs, Si_length)
    outputs = np.append(outputs, Pmode)
    outputs = np.array([outputs])
    
    
    FDTD.save_data_one_by_one(outputs, data_file, idx, output_file)
