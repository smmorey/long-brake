# -*- coding: utf-8 -*-
"""
Created on Mon May 30th 22:21:07 2016
@author: Charlie

Version updated 2022-2023
@author: Charles M. Shobe (West Virginia University) and Susannah M. Morey (University of Washington)

Simple driver for longbrake model

This particular one allows variations in gamma and initial block size.
It also expands from the reach scale to the longitudinal profile scale and includes block deposition events.
"""
import os
import sys
import re
from subprocess import call
import brake
import numpy as np
import scipy.special as sp

incision = brake.Brake() #create an instance of the class 

#INPUT VARIABLES
n_cells = 2000
dx = 250 #m

chan_width = 200 #channel width at left-most node [m]
area_0 = 2e11 #[sq m]
k_w = 0.00001 #coefficient of channel width increase with downstream distance
hack_c = 0.01 #Hack's law coefficient
hack_exp = 0.6 #Hack's law exponent

delay_timescale = 10000.
run_time =  20000 #years
record_time_int =  500 #years
timestep = 1 #years

#uplift: currently imposes the same value everywhere. Just change the next line to have spatial variability.
#the key thing is that the array needs to have (n_cells - 1) values!
imposed_rock_uplift_rate = np.repeat(0.001, n_cells - 1) # uplift units: [m/ yr]

# scaled to Yarlung Siang River
bed_k = 1.6e-13 # modify this and 'imposed_rock_uplift_rate' to change total relief
hillslope_block_k = 1.6e-13 
fluvial_block_k = 1.6e-13
gamma = 0.0 # change to add hillslope-derived boulders
side_length = 5. 
z_0 = 0.1 #roughness height for the channel bed excluding boulders
tau_c = 0.0 #critical detachment shear stress

############IMPOSING BLOCK DEPOSITION EVENTS##########################
#block deposition events
depo_events_timing = np.array([1]) #array n_depo_events long that sets when depo events occur [yrs]
#to have a depositional event at the very beginning of the model run, the above array...
#should have only one element in it, and that value should be equal to the timestep.

#important: this section has to be modified such that the number of depositional events...
#provided below matches the number of entries in the timing array above.

#FIRST DEPOSITIONSLAL EVENT
#depo_distances_1 = np.array([]) #array n_blocks long that holds the position along the model domain of every block [m]
#depo_sizes_1 = np.array([]) #array n_blocks long that holds the size of every block [m]

# building boulder distribution

# load in boulder distribution files
depo_distances_1 = np.load('brake_boulder_locations_1.npy')
depo_sizes_1 = np.load('brake_boulder_sizes_1.npy')

#note: size and distance arrays go together, so the first element in size gets put at the position indicated by the first element in distance, and so on
#depo_distances_1 = np.concatenate((depo_distances_1a, depo_distances_1b, depo_distances_1c, depo_distances_1d), axis=None)
#depo_sizes_1 = np.concatenate((depo_sizes_1a, depo_sizes_1b, depo_sizes_1c, depo_sizes_1d), axis=None)

#SECOND DEPOSITIONAL EVENT
# depo_distances_2 = np.array([850, 850, 850]) #array n_blocks long that holds the position along the model domain of every block [m]
# depo_sizes_2 = np.array([2, 2, 2]) #array n_blocks long that holds the size of every block [m]

#depo_blocks_info has n_events*2 columns, one for distance and size of blocks from each event
depo_blocks_info = np.empty([len(depo_distances_1), 2]) # one event
#depo_blocks_info = np.empty([len(depo_distances_1), 4] # two events

#depo_blocks_info = np.empty((max(len(depo_distances_1), len(depo_distances_2)), 4,))
depo_blocks_info[:] = np.nan

depo_blocks_info[:len(depo_distances_1), 0] = depo_distances_1
depo_blocks_info[:len(depo_sizes_1), 1] = depo_sizes_1
#depo_blocks_info[:len(depo_distances_2), 2] = depo_distances_2
#depo_blocks_info[:len(depo_sizes_2), 3] = depo_sizes_2


####end block deposition events#######################################

#water discharge distribution
weibull_k = 0.5 # unitless, c in Rossi et al (2016) ; metric of variability of the flow; high = low variability discharge, you don't get a lot of big flows
meanq = 1000 # change this to change the mean discharge, which will calculate the weibull scale in the next line
weibull_scale = meanq/(sp.gamma(1+(1/weibull_k)))
#print(weibull_scale)
#weibull_scale = 5.0#calc'ed outside with sp.special.gamma() ; shifts distribution higher or lower discharge

suffix = 'test_1' #text string to append onto run data

print('Running brake: gamma ', str(gamma), ' and hillslope block size ', str(side_length))


#run the model
incision.run_model(run_time, timestep, record_time_int, imposed_rock_uplift_rate, delay_timescale, gamma, 
					side_length, suffix, n_cells, dx, bed_k, hillslope_block_k, fluvial_block_k, z_0, tau_c, 
					weibull_k, weibull_scale, depo_blocks_info, depo_events_timing, chan_width, k_w,
					hack_c, hack_exp, area_0)
print('DONE')
