#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
example usage of CHZ functions
"""
from isochrones.mist import MISTEvolutionTrackGrid
import utils.hz_utils as hz
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#import stellar tracks from isochrones package
track_grid = MISTEvolutionTrackGrid()
df=track_grid.df

index_names= df.index.names

#arrays for isochrone grid feh and mass points
fehs= df.index.levels[0]
masses= df.index.levels[1]

#%% generate a habitable zone evolution object for an example MIST track
#this example is close to a solar analog
test_mass= masses[29] #1.0 M_sun
test_feh= fehs[-3] #0.0
test_eep=355 #corresponds to age = 4.58e9 yr
track=df.xs((test_feh,test_mass),level=(0,1))

#habitable zone evolution object using Kopparapu et al 2013 habitable zone
evol= hz.HZ_evolution(track, test_eep,HZ_form="K13")


#%%
#free parameters for different formulations of the CHZ
hab_start_age= 1e8 #yr
fixed_age= 2e9 #yr

#calculate duration spent in the habitable zone and plot it as a function of distance from star
tau= evol.obj_calc_tau(t_0=hab_start_age,nd=1000,mode='default')
d_arr= evol.d_range #corresponding distances for each tau value
tau_fig, tau_ax= evol.plot_tau()

#%%
#calculate different CHZ boundaries, store in evol object
evol.get_sustained_CHZ(CHZ_start_age=hab_start_age)
evol.get_fixed_age_CHZ(fixed_age= fixed_age)

#calculate fraction of habitable zone occupied by CHZ, fraction of hz planets in CHZ
fd=evol.CHZ_dist_fraction(form='sustained')
fp=evol.CHZ_planet_fraction(form='sustained')

fd2= evol.CHZ_dist_fraction(form='fixed age')
fp2=evol.CHZ_planet_fraction(form='fixed age')

#print hz boundaries and planet fractions
print("HZ: ", evol.current_i, " to ", evol.current_o)
print("Sustained CHZ: ", evol.sCHZ_i, " to ", evol.sCHZ_o)
print("2Gyr fixed duration CHZ: ", evol.fCHZ_i, " to ",evol.fCHZ_o)
print("Fraction of HZ distance in sustained CHZ: ", fd)
print("Fraction of HZ planets in sustained CHZ: ", fp)
print("Fraction of HZ distance in 2Gyr fixed duration CHZ: ", fd2)
print("Fraction of HZ planets in 2Gyr fixed duration CHZ: ", fp2)


#%% Calculate Tuchow and Wright 2020 biosignature yield metrics for different 
#assumptions about biosignature emergence
BBHZ = evol.obj_calc_B(H_form='linear',Gamma_form='lna',
                        cold_starts=True, nd=1000,hab_start_age=hab_start_age,HZ_form='R18')


b = 1 /30
A = 0.709
C0= 0.425


BBHZ2 = evol.obj_calc_B(H_form='exp',Gamma_form='lna', b= b, A=A, Gamma_norm= C0,
                        cold_starts=True, nd=1000,hab_start_age=hab_start_age,HZ_form='R18')


#%%
#plot habitable zone boundaries, dashed horizontal lines represent sustained CHZ
hz_fig, hz_ax= evol.plot_HZ(CHZ_start_age=hab_start_age,include_sCHZ=True, include_start_age=True)
hz_ax.set_xlim(5.0e7,4.8e9)
hz_ax.set_ylim(0.8,1.8)
