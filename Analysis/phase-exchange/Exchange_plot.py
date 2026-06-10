#! /usr/bin/python
# -*- coding: utf-8 -*-
#plot the condensate-dilute exchagne of different polyD10 monomers
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

fig, ax = plt.subplots(figsize=(12,8))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=22)
tick_prop = font_manager.FontProperties(fname=font_path, size=18)
legend_prop = font_manager.FontProperties(fname=font_path, size=18) 
cm = plt.cm.get_cmap('tab20')
for i in range(30):    
    print('Monomer{}'.format(i))
    if os.path.exists('Dilute_monoer{}.txt'.format(i)):
        data = np.loadtxt('Dilute_monoer{}.txt'.format(i))
        state=[]
        for time in range(1000000, 30000000, 10000):
            if time in data[:,0]:
                state.append([time,i])
        state=np.array(state)
        print(state.shape)

        color_index=i%20
        ax.plot(state[:,0]/1000, state[:,1], '.', color=cm.colors[color_index])  
ax.spines[['right', 'top']].set_visible(False)
ax.set_xlabel('Time (ns)',fontproperties=font_prop)
ax.set_ylabel('Monomer Index',fontproperties=font_prop)
#plt.xlim(0.2,0.45)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('Exchange_plot.png',dpi=600,bbox_inches='tight')
plt.show()  



