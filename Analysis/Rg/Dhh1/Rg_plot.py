#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate gyrate distribution in different trajectory interval to analyze convergence and conformation ensemble
import sys
import os
import numpy as np
import argparse
import pandas as pd
import seaborn as sns
#sns.set(color_codes=True)
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from statistics import mean
from scipy.stats import entropy

Rgs=[]
Rg_all=[]
Rg_entropy_collect=[]
KLs=np.loadtxt('Rg_KLs.xvg')

fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=16) 
cm = plt.cm.get_cmap('Set2')
for i in range(5):    
    data = np.loadtxt('Condensate_phase_monomer{}.xvg'.format(i))
    Rg = data[:,3]/10
    #use consistent 0.05nm bin for entropy
    hist_i,edge_i=np.histogram(Rg,bins=np.arange(Rg.min(),Rg.max()+0.05,0.05)) 
    prob_i=hist_i/len(Rg)
    Rg_entropy=entropy(prob_i)
    Rg_all=np.hstack((Rg_all,Rg))
    print(mean(Rg))   
    Rgs.append([i,mean(Rg)])    
    Rg_entropy_collect.append(Rg_entropy)
    
    color_index=i%8
    sns.kdeplot(Rg, color=cm.colors[color_index],alpha=0.9)  
sns.kdeplot(Rg_all, color='black', linestyle='dashed',linewidth=3)  
print(Rg_all.shape) 
ax.text( 1.5,0.90,'KL divergence: {:.2f}'.format(KLs.mean(axis=0)[1]),horizontalalignment='left',fontproperties=legend_prop)
ax.text( 1.5,0.95,'Entropy: {:.2f}'.format(mean(Rg_entropy_collect)),horizontalalignment='left',fontproperties=legend_prop)
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Density',fontproperties=font_prop)
ax.set_xlabel('Dhh1 C-ter Rg (nm)',fontproperties=font_prop)
plt.ylim(0,1)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('Condensate_Rg_distribution.png',dpi=600,bbox_inches='tight')
plt.show()  
np.savetxt('Condensate_Rg_monomers.xvg',Rgs,delimiter='	')

