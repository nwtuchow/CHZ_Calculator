#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute Tuchow and Wright 2020 metrics for a grid of models
"""

from isochrones.mist import MISTEvolutionTrackGrid
import utils.hz_utils as hz
import numpy as np
#import matplotlib.pyplot as plt
import pandas as pd

#import MIST evolutionary track grid and get multi-index
track_grid = MISTEvolutionTrackGrid()
df=track_grid.df
index_names= df.index.names

fehs= df.index.levels[0]
masses= df.index.levels[1]

#create array and dataframe to store metric values
Bexp_arr= np.zeros(len(df))

B_df = pd.DataFrame({'Bexp':Bexp_arr})

B_df.index=df.index
B_df['star_age']=df['star_age']

#constant to use in exponential formulation of H 
b = 1 /30 #e-folding time for biosignature emergence in 1/Gyr
hab_start_age=1e7
import time

t1=time.perf_counter()
for feh in fehs:
    for mass in masses:
        track=df.xs((feh,mass),level=(0,1))
        for eep in track.index:
            try:
                evol= hz.HZ_evolution(track, eep)
                B_df.loc[(feh,mass,eep),'Bexp']= evol.obj_calc_B(H_form='exp',
                        Gamma_form='lna', b= b,hab_start_age=hab_start_age,
                        cold_starts=False, nd=1000,HZ_form='R18')
                
            except:
                print("Problem at FeH = %.1f , Mass = %.2f, EEP = %d" %(feh,mass,eep))
                B_df.loc[(feh,mass,eep),'Bexp']=np.nan
            
            

t2=time.perf_counter()
print("Took ", (t2-t1), " seconds")

B_df.to_csv("outputs/Bgrid_df.csv")