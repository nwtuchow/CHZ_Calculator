#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute CHZ fractions for MIST grid
takes on order of 7 - 8 hrs to run
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
CHZ_df2 = pd.DataFrame({'f_p':fp_arr, 'f_d':fd_arr})

CHZ_df.index=df.index
CHZ_df['star_age']=df['star_age']

CHZ_df2.index=df.index
CHZ_df2['star_age']=df['star_age']

#two different durations for fixed duration CHZ
fixed_age1= 2e9
fixed_age2= 4e9



t1=time.perf_counter()
for feh in fehs:
    for mass in masses:
        track=df.xs((feh,mass),level=(0,1))
        for eep in track.index:
            try:
                evol= hz.HZ_evolution(track, eep)
                evol.obj_calc_tau(t_0=0.0,nd=500,mode='coarse')
                #else:
                #    CHZ_df.loc[(feh,mass,eep),'f_d']=0.0
                #    CHZ_df.loc[(feh,mass,eep),'f_p']=0.0
                #    CHZ_df2.loc[(feh,mass,eep),'f_d']=0.0
                #    CHZ_df2.loc[(feh,mass,eep),'f_p']=0.0
                #    continue
                
                #evol.get_CHZ(CHZ_start_age=hab_start_age)
                evol.get_fixed_age_CHZ(fixed_age=fixed_age1)
                CHZ_df.loc[(feh,mass,eep),'f_d']=evol.CHZ_dist_fraction(form="fixed")
                CHZ_df.loc[(feh,mass,eep),'f_p']=evol.CHZ_planet_fraction(form="fixed")
                
                evol.get_fixed_age_CHZ(fixed_age=fixed_age2)
                CHZ_df2.loc[(feh,mass,eep),'f_d']=evol.CHZ_dist_fraction(form="fixed")
                CHZ_df2.loc[(feh,mass,eep),'f_p']=evol.CHZ_planet_fraction(form="fixed")
            except:
                print("Problem at FeH = %.1f , Mass = %.2f, EEP = %d" %(feh,mass,eep))
                CHZ_df.loc[(feh,mass,eep),'f_p']=np.nan
                CHZ_df.loc[(feh,mass,eep),'f_d']=np.nan
                
                CHZ_df2.loc[(feh,mass,eep),'f_p']=np.nan
                CHZ_df2.loc[(feh,mass,eep),'f_d']=np.nan
            

t2=time.perf_counter()
print("Took ", (t2-t1), " seconds")

CHZ_df.to_csv("outputs/CHZ_df_2Gyr.csv")
CHZ_df2.to_csv("outputs/CHZ_df_4Gyr.csv")