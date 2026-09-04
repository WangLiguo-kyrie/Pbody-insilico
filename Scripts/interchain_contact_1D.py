#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os
    
def contacts_within_cutoff(group_a, group_b,box_dim, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions,box=box_dim)
    #print(dist.shape)
    # determine which distances <= radius
    n_contacts = contacts.contact_matrix(dist, radius).sum(axis=1)
    #print(n_contacts.shape)
    return n_contacts    

#trajectory pbc and protein extraction
trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate-replica3/'

u = mda.Universe(trajectory_dir+'production9-protein.tpr',trajectory_dir+'production9-Verlet-pbc.xtc')
condensate=u.select_atoms('name BB and index 0-68742') #excluding client polyD here 
print('condensate size {}'.format(len(condensate.atoms)))
# virtual stie in Martini couldn't be recognized as part of protein chain,so must use BB to avoid the fragments of virtual site bead. 
ag = u.select_atoms('name BB and index 41587-52524') #Pat1
print(len(ag.atoms))
fragments = ag.atoms.fragments
print('copy number {}'.format(len(fragments)))

#discard initial 2us 
start = 4000
end = -1
#every 1ns
step = 10

#dilutes=np.loadtxt(trajectory_dir+'phase_Exchange/Phase_Exchange.txt')
#Dilute_monomer_list=np.unique(dilutes[:,1])
#print(Dilute_monomer_list)
Dilute_monomer_list=[]

#loop over each fragment
for i, frag in enumerate(fragments): 
    print('***** Protein Monomer {} *****'.format(i))       
    results_contact=[]
    others=condensate-frag
    others_BB=others.select_atoms('name BB')
    frag_BB=frag.select_atoms('name BB')
    print(frag_BB.atoms)
    print(others_BB.atoms)
    timeseries = []
      
    if i in Dilute_monomer_list:
        dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
        for index, ts in enumerate(u.trajectory[start:end:step]): 
            if ts.time not in dilute_interval[:,0]:
                ca=contacts_within_cutoff(frag_BB, others_BB, ts.dimensions)
                timeseries.append(ca)
                if ts.time%100000==0:
                    print('Time {}ps, Overall Contacts {}'.format(ts.time,ca.sum()))

    else:   
        for index, ts in enumerate(u.trajectory[start:end:step]):    
            ca=contacts_within_cutoff(frag_BB, others_BB, ts.dimensions)
            timeseries.append(ca)
            if ts.time%100000==0:
                print('Time {}ps, Overall Contacts {}'.format(ts.time,ca.sum())) 
    
    timeseries=np.array(timeseries)
    print(timeseries.shape)
    frames=timeseries.shape[0]
    print('Overall {} frames within condensate'.format(frames))
    contacts_strength=timeseries.mean(axis=0)
    np.savetxt('replica3-interchain_contact_monomer{}.xvg'.format(i),contacts_strength.T,delimiter='	')
  

    


