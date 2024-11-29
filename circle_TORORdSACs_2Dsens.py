# -*- coding: utf-8 -*-
"""
Code for  "Spatiotemporal determinants of stretch-activated channel-induced 
re-entry in ventricular tissue: an in-silico study" 
by Melania Buonocunto, Aurore Lyon, Tammo Delhaas, Joost Lumens, and Jordi Heijman
contact: work: m.buonocunto@maastrichtuniversity.nl ; melania.buonocunto@medunigraz.at
personal: melania.buonoc@gmail.com  
"""

import numpy as np
import myokit
import matplotlib.pyplot as plt
import sys
import pandas as pd
import datetime
import time
import math

def model(j,data,q_inbetween, q_final):
    n, ns_ratio, ts , count = data
    ns = math.sqrt(n**2*ns_ratio) 
    log_interval = 1
    bcl = 1000
    tinterval=10 #change for continuous stretch
    
    #import model from myokit
    m, p, x = myokit.load('ToRORdSACs_2D_model.mmt')

   #simulate 2D
    p = myokit.pacing.blocktrain(1000, 1, 0, 5, 0)  # (period, duration, offset=0, level=1.0, limit=0) period = 1000 if i want it periodic, otherwise 0
    s = myokit.SimulationOpenCL(m, p, ncells=[n, n])
    s.set_paced_cells(1, n, 0, 0)  # set which cells I want to pace
    
    #set smulation parameters
    s.set_conductance(1.5,1.5)
    s.set_step_size(0.001)
    s.set_constant('SACs.tstretch', ts)
    s.set_constant('SACs.tinterval',tinterval)
    
    #define stretch amplitude
    ampl = 1.3 
    

    # Define the size of the matrix and the radius of the stretched circle (rs) and the inner filled circle (ri)
    
    ## total area of if it was a square (ns**2) = filled circle (ri) + gradient area (rs)
    rs = ns/(math.sqrt(math.pi))   #radius of the stretched circle
    ri = rs #0, rs/2 or rs #radius of inner circle (fully stretched!)

    # Create a matrix filled with ones
    lam = np.ones((n, n))
    
    # Determine the center of the circle
    #center = (ns // 2 , ns // 2) #upper left
    center = (n // 2 , n // 2) #center
    
    #FILLED
    for ii in range(n):
        for jj in range(n):
       
            if (ii - center[0]) ** 2 + (jj - center[1]) ** 2 <= (rs) ** 2:  
                lam[ii, jj] = ampl
               
               
    #GRADIENT
   # Iterate over the matrix indices and check if each point is within the circle
    # for ii in range(n):
    #     for jj in range(n):
       
    #         if (ii - center[0]) ** 2 + (jj - center[1]) ** 2 <= (ri) ** 2:  
    #             lam[ii, jj] = ampl
               
    #         elif (((ii - center[0]) ** 2 + (jj - center[1]) ** 2 > (ri) ** 2)\
    #               and ((ii - center[0]) ** 2 + (jj - center[1]) ** 2 <= (rs) ** 2)): 
    #                 rho = np.sqrt((ii - center[0]) ** 2 + (jj - center[1]) ** 2) #radial direction 
    #                 lam[ii, jj] = (1-ampl)/(rs-ri)*(rho-ri)+ampl

    
    # #Display the matrix
    # plt.imshow(lam, cmap='binary')
    # plt.show()
    

    ##set field in the raw model
    s.set_field('SACs.Lambda' , lam)
    ## Set up logging
    log = ['environment.time', 'membrane.v']  
    i = 0
    Vmm = 0
    t_run = [500, 200, 299, 31, 470]
    #Upload the prepaced state  
    pre_pace = np.load('Prepace_0D500000_2D10000_conductance_1.5_nostretch_120723.npy')
    s.set_state(pre_pace)
    
    for i in range(0,5): #5 #runs for around 3 beats maximum   #range (0,4) if fixed step size, otherwise use range (0,5)
        i=i+1 
        now = datetime.datetime.now()
        q_inbetween.put(str(now) + f', simulation {j} , parameters: n={n} ns_ratio={ns_ratio} ts={ts}, iteration: {i}')    
        time.sleep(np.random.random())
        
        #decrease step size at the onset of the new beat to avoid numerical errors (between t=999 and 1030)
        if i == 4:
            s.set_step_size(0.0001)
        else:
            s.set_step_size(0.001) 
                    
        log = s.run(t_run[i-1], log=log, log_interval=log_interval) 
        block = log.block2d()
        Vm = block.get2d('membrane.v')
        maxval = np.max(Vm[len(Vm)-1].flatten())
            
        if maxval < 0: #-35: for underthreshold activities #0: for afterdepolarizations
          print (f'simulation {j}: No front waves detected after {i} iterations')
          break
           
    # Store
    if i>=1: #3
        block = log.block2d()
        block.save('simulation_stretch_centre= %d.zip' %count)

    result = {
        'n': n,
        'ns_ratio': ns_ratio,
        'ts': ts,
        'count': count,
        'mv': maxval,
        'n_it': i
        }
    
    q_final.put((j, result))









