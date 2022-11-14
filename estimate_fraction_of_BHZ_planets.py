# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 16:56:35 2022
Estimate BHZ fraction
@author: nwtuc
"""

from isochrones.mist import MISTEvolutionTrackGrid
from isochrones.mist import MIST_EvolutionTrack
from isochrones.interp import DFInterpolator
import utils.hz_utils as hz
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from scipy.interpolate import griddata,CloughTocher2DInterpolator,interp2d
import matplotlib.tri as tri
from scipy.ndimage import gaussian_filter

track_grid = MISTEvolutionTrackGrid()
mist_track = MIST_EvolutionTrack()

df=track_grid.df
index_names= df.index.names
fehs= df.index.levels[0]
masses= df.index.levels[1]

#CHZ_df= pd.read_csv("outputs/CHZ_df_late_start.csv")
#CHZ_df= pd.read_csv("outputs/CHZ_df.csv") #for linux
CHZ_df= pd.read_csv(hz.DATA_DIR+"/outputs/CHZ_df.csv")
#CHZ_df= pd.read_csv(hz.DATA_DIR+"/outputs/CHZ_df_late_start.csv")
CHZ_df = CHZ_df.set_index(index_names)

'''mass_arr=CHZ_df['initial_mass'].to_numpy()
age_arr= CHZ_df['star_age'].to_numpy()/1e9 #in Gyr
fd_arr= CHZ_df['f_d'].to_numpy()
fp_arr = CHZ_df['f_p'].to_numpy()

for i in range(len(fd_arr)):
    if np.isnan(fd_arr[i]):
        fd_arr[i]=0
    
    if np.isnan(fp_arr[i]):
        fp_arr[i]=0
'''

#%%

cond1=np.isnan(CHZ_df[["f_p","f_d"]])

CHZ_df = CHZ_df.where(~cond1,0)
#CHZ_df.loc[cond1,["f_p","f_d"]]

interp= DFInterpolator(CHZ_df)

def interpfunc(mass,age,feh):
    eep=mist_track.get_eep(mass, np.log10(age), feh,accurate=False)
    farr= interp([feh,mass,eep],['f_d','f_p'])
    return farr
        
#%%
starlist_df= pd.read_csv("target_lists/mission_stars_exocat_formatted.csv")

stmasses= starlist_df['st_mass'].to_numpy()
stmasses= stmasses[~np.isnan(stmasses)]

stfeh = starlist_df['st_metfe'].to_numpy()
stfeh= stfeh[~np.isnan(stfeh)]

nstar=10000
rand_masses=np.random.choice(stmasses,size=nstar,replace=True)
rand_feh = np.random.choice(stfeh, size= nstar, replace =True)
rand_ages= np.random.uniform(low=0.2, high=10.0,size=nstar) * (10**9)

fd_sample= np.empty(nstar)
fp_sample= np.empty(nstar)

for i in range(nstar):
    if i%50 ==0:
        print("Index =", i)
    temp_f= interpfunc(rand_masses[i],rand_ages[i],rand_feh[i])
    fd_sample[i]=temp_f[0]
    fp_sample[i]=temp_f[1]
    

fd_sample= np.where(~np.isnan(fd_sample),fd_sample,0.0)
fp_sample= np.where(~np.isnan(fp_sample),fp_sample,0.0)

#fp_sample= np.ma.filled(fp_sample,fill_value=0)

mean_fp= np.mean(fp_sample)
fBHZ=1-mean_fp
print("f_BHZ = ", fBHZ)

#for early start
#f_BHZ = 0.743

#for late start
#f_BHZ = 0.300