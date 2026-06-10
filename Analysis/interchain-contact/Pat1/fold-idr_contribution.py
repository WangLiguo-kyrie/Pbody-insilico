#! /usr/bin/python
#collect the interchain contact contribution from folded and idr
import sys
import os
import numpy as np
import MDAnalysis as mda
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

def fold_filter(fold_array):
    #filter the buried residues in folded domains, using max(replicas)<1 
    print('contributing folded domain residues: {}'.format(fold_array.shape[0]))
    fold_array_filter=fold_array[fold_array.max(axis=1)>1]
    print('contributing folded domain residues after filter: {}'.format(fold_array_filter.shape[0]))
    return fold_array_filter

#Pat1: Idr1 1-255, idr2 375-475, MD 256-374, HEAT 476-796
contact_collect=[]
for i in range(6):    
    contact = np.loadtxt('interchain_contact_monomer{}.xvg'.format(i))   
    if i==0:contact_collect=contact
    else:contact_collect=np.vstack((contact_collect,contact))
print(contact_collect.shape)  
contact_collect_trans=contact_collect.T
print(contact_collect_trans.shape) 

idr1=contact_collect_trans[:255,:]
idr2=contact_collect_trans[374:475,:]
fold1=contact_collect_trans[255:374,:]
fold2=contact_collect_trans[475:,:]
fold1_filter=fold_filter(fold1)
fold2_filter=fold_filter(fold2)

#idr_contribution=np.vstack((idr1,idr2))
idr_contribution=idr1
fold_contribution=np.vstack((fold1_filter,fold2_filter))
print(idr_contribution.shape)
print(fold_contribution.shape)
np.savetxt('idr_contribution.xvg',idr_contribution,delimiter='	')
np.savetxt('fold_contribution.xvg',fold_contribution,delimiter='	')


#solvation
#Pat1: Idr1 1-255, idr2 375-475, MD 256-374, HEAT 476-796
contact_collect=[]
for i in range(6):    
    contact = np.loadtxt('solvation_contact_monomer{}.xvg'.format(i))   
    if i==0:contact_collect=contact
    else:contact_collect=np.vstack((contact_collect,contact))
print(contact_collect.shape)  
contact_collect_trans=contact_collect.T
print(contact_collect_trans.shape) 

idr1=contact_collect_trans[:255,:]
idr2=contact_collect_trans[374:475,:]
fold1=contact_collect_trans[255:374,:]
fold2=contact_collect_trans[475:,:]
fold1_filter=fold_filter(fold1)
fold2_filter=fold_filter(fold2)

#idr_contribution=np.vstack((idr1,idr2))
idr_contribution=idr1
fold_contribution=np.vstack((fold1_filter,fold2_filter))
print(idr_contribution.shape)
print(fold_contribution.shape)
np.savetxt('idr_solvation.xvg',idr_contribution,delimiter='	')
np.savetxt('fold_solvation.xvg',fold_contribution,delimiter='	')
