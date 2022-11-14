#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute CHZ fractions for MIST grid
"""

from isochrones.mist import MISTEvolutionTrackGrid
import utils.hz_utils as hz
import numpy as np
import pandas as pd
import time

#import MIST evolutionary track grid and get multi-index
track_grid = MISTEvolutionTrackGrid()
df=track_grid.df

index_names= df.index.names
fehs= df.index.levels[0]
masses= df.index.levels[1]

#create pandas dataframe to store CHZ fractions
fp_arr= np.zeros(len(df))
fd_arr= np.zeros(len(df))

CHZ_df = pd.DataFrame({'f_p':fp_arr, 'f_d':fd_arr})

CHZ_df.index=df.index
CHZ_df['star_age']=df['star_age']

CHZ_df2 = pd.DataFrame({'f_p':fp_arr, 'f_d':fd_arr})

CHZ_df2.index=df.index
CHZ_df2['star_age']=df['star_age']

#two starting ages for different scenerios for volatile delivery
hab_start_age=1e7 # rough age of protoplanetary disk
hab_start_age2=1e8 #~ late veneer start for habitability

t1=time.perf_counter()

#loop over all tracks in mist grid
for feh in fehs:
    for mass in masses:
        track=df.xs((feh,mass),level=(0,1))
        for eep in track.index:
            try:
                #generate hz evolution objects for both scenarios
                #calculate CHZ distance and planet fractions
                #i.e fraction of hz width, or fraction of hz planets in CHZ
                evol= hz.HZ_evolution(track, eep)
                evol.get_sustained_CHZ(CHZ_start_age=hab_start_age)
                CHZ_df.loc[(feh,mass,eep),'f_d']=evol.CHZ_dist_fraction(form="sustained")
                CHZ_df.loc[(feh,mass,eep),'f_p']=evol.CHZ_planet_fraction(form="sustained")
                
                evol.get_sustained_CHZ(CHZ_start_age=hab_start_age2)
                CHZ_df2.loc[(feh,mass,eep),'f_d']=evol.CHZ_dist_fraction(form="sustained")
                CHZ_df2.loc[(feh,mass,eep),'f_p']=evol.CHZ_planet_fraction(form="sustained")
                
            except:
                print("Problem at FeH = %.1f , Mass = %.2f, EEP = %d" %(feh,mass,eep))
                CHZ_df.loc[(feh,mass,eep),'f_p']=np.nan
                CHZ_df.loc[(feh,mass,eep),'f_d']=np.nan
                
                CHZ_df2.loc[(feh,mass,eep),'f_p']=np.nan
                CHZ_df2.loc[(feh,mass,eep),'f_d']=np.nan
            

t2=time.perf_counter()
print("Took ", (t2-t1), " seconds")

CHZ_df.to_csv("outputs/sCHZ_df_1e7.csv")
CHZ_df2.to_csv("outputs/sCHZ_df_1e8.csv")