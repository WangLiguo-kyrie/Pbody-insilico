#! /usr/bin/python

import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis import *
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate_small/'

u = mda.Universe(trajectory_dir+'production8-protein.tpr',trajectory_dir+'production8-Verlet-pbc.xtc')
condensate=u.select_atoms('all')
scaffold=u.select_atoms('index 0-68742')
client=u.select_atoms('index 68743-73663')

#discard initial 2us 
start = 4000
end = -1
#every 1ns
step = 2

#define log-space with sparse interval at small values
edges=np.log(np.linspace(1, np.power(1.05,70), 70))/np.log(1.05)*3
print(edges)

hist_collect=[]
for index, ts in enumerate(u.trajectory[start:end:step]): 
    dist = contacts.distance_array(condensate.center_of_geometry(), scaffold.positions,ts.dimensions)
    hist,edge=np.histogram(dist,bins=edges)

    bins=[(edge[i]+edge[i+1])/2 for i in range(69)]
    #print(bins)
    vol_shell=[4*np.pi*(edge[i+1]**3-edge[i]**3)/3 for i in range(69)]
    #print(vol_shell)
    bins=np.array(bins)
    density=[hist[i]/vol_shell[i] for i in range(69)]
    hist_collect.append(density)
hist_collect=np.array(hist_collect)
print(hist_collect.shape)
results=np.vstack((bins, hist_collect.mean(axis=0), hist_collect.std(axis=0)))
print(results)   
np.savetxt('droplet_scaffold_density.xvg',results,delimiter='	')


hist_collect=[]
for index, ts in enumerate(u.trajectory[start:end:step]): 
    dist = contacts.distance_array(condensate.center_of_geometry(), client.positions,ts.dimensions)
    hist,edge=np.histogram(dist,bins=edges)

    bins=[(edge[i]+edge[i+1])/2 for i in range(69)]
    #print(bins)
    vol_shell=[4*np.pi*(edge[i+1]**3-edge[i]**3)/3 for i in range(69)]
    #print(vol_shell)
    bins=np.array(bins)
    density=[hist[i]/vol_shell[i] for i in range(69)]
    hist_collect.append(density)
hist_collect=np.array(hist_collect)
print(hist_collect.shape)
results=np.vstack((bins, hist_collect.mean(axis=0), hist_collect.std(axis=0)))
print(results)   
np.savetxt('droplet_client_density.xvg',results,delimiter='	')



