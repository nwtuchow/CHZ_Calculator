#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot fraction of HZ planets in the sustained CHZ
"""
from isochrones.mist import MISTEvolutionTrackGrid
import utils.hz_utils as hz
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from scipy.interpolate import griddata,CloughTocher2DInterpolator,interp2d
import matplotlib.tri as tri
from scipy.ndimage import gaussian_filter

#import mist tracks and store index names
track_grid = MISTEvolutionTrackGrid()

df=track_grid.df
index_names= df.index.names
fehs= df.index.levels[0]
masses= df.index.levels[1]

#import results of sustained CHZ calculations from "compute_sCHZ_grid.py"
CHZ_df= pd.read_csv("outputs/sCHZ_df_1e7.csv")
CHZ_df2= pd.read_csv("outputs/sCHZ_df_1e8.csv")

#select the slice of the grid that we are going to plot
CHZ_df=CHZ_df.query('(initial_mass <= 2) and (star_age < 1.4e10) and (initial_feh ==0.0)')
CHZ_df2=CHZ_df2.query('(initial_mass <= 2) and (star_age < 1.4e10) and (initial_feh ==0.0)')

mass_arr=CHZ_df['initial_mass'].to_numpy()
age_arr= CHZ_df['star_age'].to_numpy()/1e9 #in Gyr
fd_arr= CHZ_df['f_d'].to_numpy()
fp_arr = CHZ_df['f_p'].to_numpy()

mass_arr2=CHZ_df2['initial_mass'].to_numpy()
age_arr2= CHZ_df2['star_age'].to_numpy()/1e9 #in Gyr
fd_arr2= CHZ_df2['f_d'].to_numpy()
fp_arr2 = CHZ_df2['f_p'].to_numpy()

#set nans to zero to not confuse the interpolation scheme
for i in range(len(fd_arr)):
    if np.isnan(fd_arr[i]):
        fd_arr[i]=0
    
    if np.isnan(fp_arr[i]):
        fp_arr[i]=0
        
for j in range(len(fd_arr2)):
    if np.isnan(fd_arr2[j]):
        fd_arr2[j]=0
    
    if np.isnan(fp_arr2[j]):
        fp_arr2[j]=0

#%% Use LinearTriInterpolator to interpolate grid of irregularly spaced datapoints
npt=200
M= np.linspace(masses[0],1.45,npt)
A=np.linspace(0.0, 10.0,npt)

Mgrid, Agrid = np.meshgrid(M,A)

triang = tri.Triangulation(mass_arr, age_arr)
fd_interp= tri.LinearTriInterpolator(triang,fd_arr)
fp_interp= tri.LinearTriInterpolator(triang,fp_arr)
fd_interp2= tri.LinearTriInterpolator(triang,fd_arr2)
fp_interp2= tri.LinearTriInterpolator(triang,fp_arr2)


FD= fd_interp(Mgrid,Agrid)
FP= fp_interp(Mgrid,Agrid)
FD2= fd_interp2(Mgrid,Agrid)
FP2= fp_interp2(Mgrid,Agrid)

FD= np.ma.filled(FD, fill_value=0)
FP= np.ma.filled(FP, fill_value=0)
FD2= np.ma.filled(FD2, fill_value=0)
FP2= np.ma.filled(FP2, fill_value=0)

FD= gaussian_filter(FD, sigma=2)
FP= gaussian_filter(FP, sigma=2)
FD2= gaussian_filter(FD2, sigma=2)
FP2= gaussian_filter(FP2, sigma=2)


#%% get zams and tams contours

ZAMS_df= CHZ_df[CHZ_df["EEP"]==202]
TAMS_df= CHZ_df[CHZ_df["EEP"]==454]
ZAMS_df2= CHZ_df2[CHZ_df2["EEP"]==202]
TAMS_df2= CHZ_df2[CHZ_df2["EEP"]==454]
#%% plot fraction of planets in BHZ = 1 - fraction in sCHZ

nlvl=20
cmapname='Blues' #"Blues_r

fsize=16

fig, ax = plt.subplots()
cs=ax.contourf(M,A,1-FP, levels=nlvl,vmin=0.0, vmax=1.0,cmap=cmapname)

cbar = fig.colorbar(cs, ax=ax, shrink=0.95)
ax.plot(ZAMS_df["initial_mass"],ZAMS_df["star_age"]/1e9, color='white',ls='--')
ax.plot(TAMS_df["initial_mass"],TAMS_df["star_age"]/1e9, color='white',ls='--')
ax.scatter([1.0],[4.6],color="yellow",marker='*')
ax.set_xlabel(r"Mass ($M_{sun}$)", fontsize=fsize)
ax.set_ylabel("Age (Gyr)", fontsize=fsize)
ax.set_title(r"$t_0 = 10$ Myr", fontsize=fsize)
ax.set_ylim([0.0,10.0])
ax.set_xlim([masses[0], 1.45])
fig.tight_layout()



fig2, ax2 = plt.subplots()
cs2=ax2.contourf(M,A,1-FP2, levels=nlvl,vmin=0.0, vmax=1.0,cmap=cmapname)

cbar2 = fig2.colorbar(cs2, ax=ax2, shrink=0.95)
ax2.plot(ZAMS_df2["initial_mass"],ZAMS_df2["star_age"]/1e9, color='white',ls='--')
ax2.plot(TAMS_df2["initial_mass"],TAMS_df2["star_age"]/1e9, color='white',ls='--')
ax2.scatter([1.0],[4.6],color="yellow",marker='*')
ax2.set_xlabel(r"Mass ($M_{sun}$)", fontsize=fsize)
ax2.set_ylabel("Age (Gyr)", fontsize=fsize)
ax2.set_title(r"$t_0 = 100$ Myr", fontsize=fsize)
ax2.set_ylim([0.0,10.0])
ax2.set_xlim([masses[0], 1.45])
fig2.tight_layout()


