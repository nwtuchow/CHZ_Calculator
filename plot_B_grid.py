#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot the values for an example biosignature yield metric
"""
from isochrones.mist import MISTEvolutionTrackGrid
import utils.hz_utils as hz
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.tri as tri
from scipy.ndimage import gaussian_filter
import random as rand

#import mist tracks and store index names
track_grid = MISTEvolutionTrackGrid()

df=track_grid.df
index_names= df.index.names
fehs= df.index.levels[0]
masses= df.index.levels[1]

#import results of metric calculations from "compute_BHZ_grid.py"
#B_df= pd.read_csv("outputs/Bexp_df.csv")
B_df= pd.read_csv("outputs/Bgrid_df.csv")

#select the slice of the grid that we are going to plot
B_df=B_df.query('(initial_mass <= 1.5) and (star_age < 1.1e10) and (initial_feh ==0.0)')

mass_arr=B_df['initial_mass'].to_numpy()
age_arr= B_df['star_age'].to_numpy()/1e9 #in Gyr
B_arr = B_df['Bexp'].to_numpy()

 
#%%  Use LinearTriInterpolator to interpolate grid of irregularly spaced datapoints
npt=200
M= np.linspace(masses[0],1.45,npt)
A=np.linspace(0.0, 10.0,npt)

Mgrid, Agrid = np.meshgrid(M,A)

triang = tri.Triangulation(mass_arr, age_arr)
B_interp= tri.LinearTriInterpolator(triang,B_arr)

BEXP = B_interp(Mgrid,Agrid)

BEXP= np.ma.filled(BEXP, fill_value=0)

#smooth over rough edges caused by interpolation scheme
BEXP= gaussian_filter(BEXP, sigma=2)

B_sun= B_interp(1.0,4.57)

#normalize to solar values
BEXP=BEXP/B_sun
#%% get zams and tams contours

ZAMS_df= B_df[B_df["EEP"]==202]
TAMS_df= B_df[B_df["EEP"]==454]
#%% plot grid of biosignature yield metrics

nlvl=20
cmapname='Blues_r' #"Blues_r
fsize=16

fig, ax = plt.subplots()
cs=ax.contourf(M,A,BEXP, levels=nlvl,cmap=cmapname)
#cs=ax.tricontourf(mass_arr,age_arr,fd_arr, levels=nlvl,vmin=0.0, vmax=1.0,cmap=cmapname)
cbar = fig.colorbar(cs, ax=ax, shrink=0.95)
cbar.set_label("B (solar units)",fontsize=fsize)

ax.plot(ZAMS_df["initial_mass"],ZAMS_df["star_age"]/1e9, color='white',ls='--')
ax.plot(TAMS_df["initial_mass"],TAMS_df["star_age"]/1e9, color='white',ls='--')
ax.scatter([1.0],[4.57],color="yellow",marker='*')
ax.set_xlabel("Mass (Msun)", fontsize=fsize)
ax.set_ylabel("Age (Gyr)", fontsize=fsize)
#ax.set_title("B_exp")
ax.set_ylim([0.0,10.0])
ax.set_xlim([masses[0], 1.45])
