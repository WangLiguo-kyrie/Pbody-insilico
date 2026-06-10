#! /usr/bin/python
#collect the interchain contact contribution from folded and idr
import sys
import os
import numpy as np
import MDAnalysis as mda
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager


contact_collect=[]
for i in range(36):    
    contact = np.loadtxt('interchain_contact_monomer{}.xvg'.format(i))   
    if i==0:contact_collect=contact
    else:contact_collect=np.vstack((contact_collect,contact))
print(contact_collect.shape)  
contact_collect_trans=contact_collect.T
print(contact_collect_trans.shape) 

idr_contribution=contact_collect_trans
print(idr_contribution.shape)
np.savetxt('idr_contribution.xvg',idr_contribution,delimiter='	')

#solvation
contact_collect=[]
for i in range(36):    
    contact = np.loadtxt('solvation_contact_monomer{}.xvg'.format(i))   
    if i==0:contact_collect=contact
    else:contact_collect=np.vstack((contact_collect,contact))
print(contact_collect.shape)  
contact_collect_trans=contact_collect.T
print(contact_collect_trans.shape) 

idr_contribution=contact_collect_trans
print(idr_contribution.shape)
np.savetxt('idr_solvation.xvg',idr_contribution,delimiter='	')
