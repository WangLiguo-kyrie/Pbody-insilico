#! /usr/bin/python
# -*- coding: utf-8 -*-
# test the KL and shannon entropy sensitivity to bin width
import sys
import os
import numpy as np
import argparse
from statistics import mean
from scipy.stats import entropy

results=[]
binwidths=[0.025,0.05,0.075,0.1]
for binwidth in binwidths:
    #combine all replicas
    Rg_all=[]
    for replica_id in range(1,4):
        for i in range(6):    
            data = np.loadtxt('replica{}-Condensate_phase_monomer{}.xvg'.format(replica_id, i))
            Rg = data[:,3]/10
            Rg_all=np.hstack((Rg_all,Rg))
            #print(mean(Rg))   
    print(mean(Rg_all))  
    print('***')
    print(Rg_all.min())
    print(Rg_all.max())
    hist_all,edge_all=np.histogram(Rg_all,bins=np.arange(Rg_all.min(),Rg_all.max()+binwidth, binwidth))   
    prob_all=hist_all/len(Rg_all)

    KLs=[]
    entropys=[]
    for replica_id in range(1,4):
        for i in range(6):    
            data = np.loadtxt('replica{}-Condensate_phase_monomer{}.xvg'.format(replica_id, i))
            Rg_i = data[:,3]/10
            hist_i,edge_i=np.histogram(Rg_i,bins=np.arange(Rg_all.min(),Rg_all.max()+binwidth, binwidth))   
            prob_i=hist_i/len(Rg_i)
            KL_i=entropy(prob_i,prob_all)
            Rg_entropy=entropy(prob_i)
            print(KL_i)
            print(Rg_entropy)
            KLs.append(KL_i)
            entropys.append(Rg_entropy)
    print(KLs)
    print(entropys)
    results.append([binwidth, mean(KLs), mean(entropys)])
results=np.array(results)
np.savetxt('binwidth_sensitivity.xvg', results,delimiter='	')


