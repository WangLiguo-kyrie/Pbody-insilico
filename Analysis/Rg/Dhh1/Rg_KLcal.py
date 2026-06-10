#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate KL divergence of individual monomers to the overall Rg distribution
import sys
import os
import numpy as np
import argparse
from statistics import mean
from scipy.stats import entropy

Rg_all=[]
KLs=[]

for i in range(5):    
    data = np.loadtxt('Condensate_phase_monomer{}.xvg'.format(i))
    Rg = data[:,3]/10
    Rg_all=np.hstack((Rg_all,Rg))
    print(mean(Rg))   
print(mean(Rg_all))  
hist_all,edge_all=np.histogram(Rg_all,bins=np.arange(1.5,6.5,0.05))   #define the range 1.5~4.5 nm
prob_sll=hist_all/len(Rg_all)

for i in range(5):    
    data = np.loadtxt('Condensate_phase_monomer{}.xvg'.format(i))
    Rg_i = data[:,3]/10
    hist_i,edge_i=np.histogram(Rg_i,bins=np.arange(1.5,6.5,0.05))   #define the range 1.5~4.5 nm
    prob_i=hist_i/len(Rg_i)
    KL_i=entropy(prob_i,prob_sll)
    print(KL_i)
    KLs.append([i,KL_i])

KLs=np.array(KLs)
print(KLs)
np.savetxt('Rg_KLs.xvg',KLs,delimiter='	')
print(KLs.mean(axis=0))
print('Mean KL:{}'.format(KLs.mean(axis=0)[1]))




