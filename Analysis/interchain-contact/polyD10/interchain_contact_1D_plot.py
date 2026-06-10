#! /usr/bin/python
# -*- coding: utf-8 -*-
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

exchange_record=np.unique(np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate_small/phase_Exchange/Phase_Exchange.txt')[:,1])
exchange=[int(i) for i in exchange_record]
non_exchange=[i for i in range(30) if i not in exchange]
print(non_exchange)
print(exchange)


fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=16) 
cm = plt.cm.get_cmap('Accent')
non_exchange_overall=[]
for idx,i in enumerate(non_exchange):    
    data = np.loadtxt('interchain_contact_monomer{}.xvg'.format(i))
    if idx==0:
        non_exchange_overall=data
    else:non_exchange_overall=np.vstack((non_exchange_overall,data))
    ax.plot(data, color=cm.colors[4],alpha=0.6,label='Residing') 
ax.plot(non_exchange_overall.mean(axis=0), color=cm.colors[4],linewidth=3.5,linestyle='dashdot')     
exchange_overall=[]
for idx,i in enumerate(exchange): 
    data = np.loadtxt('interchain_contact_monomer{}.xvg'.format(i))
    if idx==0:
        exchange_overall=data
    else:exchange_overall=np.vstack((exchange_overall,data))
    ax.plot(data, color=cm.colors[5],alpha=0.6,label='Exchanging')
ax.plot(exchange_overall.mean(axis=0), color=cm.colors[5],linewidth=3.5,linestyle='dashdot')  
ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Interchain Contact Number',fontproperties=font_prop)
ax.set_xlabel('Residue',fontproperties=font_prop)
plt.ylim(1,18)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#remove duplicates in legend
handles, labels = plt.gca().get_legend_handles_labels()
newLabels, newHandles = [], []
for handle, label in zip(handles, labels):
  if label not in newLabels:
    newLabels.append(label)
    newHandles.append(handle)
plt.legend(newHandles, newLabels,prop=legend_prop,frameon=False)
plt.savefig('interchain_contact_1D_polyD10.png',dpi=600,bbox_inches='tight')
plt.show()  


#polyD10  polyD60
fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=16) 
cm = plt.cm.get_cmap('Accent')
non_exchange_overall=[]
for idx,i in enumerate(non_exchange):    
    data = np.loadtxt('interchain_contact_monomer{}.xvg'.format(i))
    if idx==0:
        non_exchange_overall=data
    else:non_exchange_overall=np.hstack((non_exchange_overall,data))
print(non_exchange_overall)
polyD60_overall=[]
for i in range(36):
    data = np.loadtxt('../polyD60/interchain_contact_monomer{}.xvg'.format(i))
    if i==0:
        polyD60_overall=data
    else:polyD60_overall=np.hstack((polyD60_overall,data))
print(polyD60_overall)

sns.violinplot([non_exchange_overall,polyD60_overall],palette="Accent",inner_kws=dict(box_width=10, whis_width=1.5))
ax.spines[['right', 'top']].set_visible(False)       
ax.set_ylabel('Interchain Contact Number',fontproperties=font_prop)
ax.set_xlabel('Client',fontproperties=font_prop)
plt.xticks(ticks=[0,1],labels=[r'Residing (polyD)$_{10}$',r'(polyD)$_{60}$'],fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
#plt.savefig('interchain_contact_1D_diff_client.png',dpi=600,bbox_inches='tight')
#plt.show()
