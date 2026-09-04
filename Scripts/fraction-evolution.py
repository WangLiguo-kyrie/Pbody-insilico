#! /usr/bin/python
# -*- coding: utf-8 -*-
#plot the fraction of  polyD10 monomers in different phase
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



dilute=np.loadtxt('Phase_Exchange.txt')
print(dilute.shape)
radial_dis=np.loadtxt('droplet_polyD10_distr.xvg')
#25ns step 
result=[]
for time in range(1000000, 30000000, 25000):
    dilute_monomers=dilute[dilute[:,0]==time][:,1].tolist()
    print(dilute_monomers)
    print(len(dilute_monomers))
    
    others=[i for i in range(30) if i not in dilute_monomers]
    all_dist=radial_dis[radial_dis[:,0]==time].flatten()
    #use raidal dist>9nm as interface criterion 
    indices = np.where(all_dist > 90)[0] 
    index_update=indices[1:]-1
    interface_monomers=list(set(index_update) - set(dilute_monomers))
    print(interface_monomers)
    print(len(interface_monomers))
    
    core_number=30-len(dilute_monomers)-len(interface_monomers)
    result.append([time, len(dilute_monomers), len(interface_monomers), core_number])
result=np.array(result)
print(result.shape)
np.savetxt('state_fraction_record.xvg',result,delimiter='	')

fig, ax = plt.subplots(figsize=(8,6))
font_path =  '/grain/liguo/Script/python-plot/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=15) 
cm = plt.cm.get_cmap('tab20')

ax.plot(result[:,0]/1000, result[:,1]/30*100,  color=cm.colors[0], label='Dilute')  
ax.plot(result[:,0]/1000, result[:,2]/30*100,  color=cm.colors[1], label='Interface')  
ax.plot(result[:,0]/1000, result[:,3]/30*100,  color=cm.colors[2], label='Condensate Interior')  
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('Time (ns)',fontproperties=font_prop)
ax.set_ylabel('Fraction (%)',fontproperties=font_prop)
plt.ylim(0,100)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
plt.legend(prop=legend_prop, frameon=False)
plt.savefig('fraction_evolution.png',dpi=600,bbox_inches='tight')
plt.show()  



