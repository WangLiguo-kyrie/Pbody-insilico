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

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate-replica3/'

u = mda.Universe(trajectory_dir+'production9-protein.tpr',trajectory_dir+'production9-Verlet-pbc.xtc')
condensate=u.select_atoms('all')
ag = u.select_atoms('name BB and index 73063-73662')  #polyD10
print(len(ag.atoms))
fragments = ag.atoms.fragments

#discard initial 1us 
start = 2000
end = -1
#every 25ns
step = 50

dist_result=[]
for index, ts in enumerate(u.trajectory[start:end:step]):  
    dist_collect=[]
    dist_collect.append(ts.time)
    for frag in fragments:
        print(frag.atoms)
        print(len(frag.atoms))
        
        #distance based on folded domain COG to condensate COG to determine if belongs to interface region
        dist = contacts.distance_array(condensate.center_of_geometry(), frag.center_of_geometry(),box=ts.dimensions)
        dist_collect.append(dist[0][0])
    dist_result.append(dist_collect)
dist_result=np.array(dist_result)
print(dist_result.shape)
np.savetxt('droplet_polyD10_distr.xvg',dist_result,delimiter='	')


    
    



