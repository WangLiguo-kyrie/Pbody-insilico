#! /usr/bin/python
#cluster the pore within condensates
#each grid voxel not containing protein beads, viewed as pore; and then cluster the pore grids
import sys
import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from sklearn.cluster import AgglomerativeClustering

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate_small/'

u = mda.Universe(trajectory_dir+'production8-Verlet-pbc-protein.gro',trajectory_dir+'production8-Verlet-pbc.xtc')
condensate=u.select_atoms('name BB or name SC*')  

#based on condensate COG, define a core cubix box with side=4.8*2nm, avoiding sampling too much region outside condensate, because condensate is not a perfect spherical droplet
box_side=49.5
#define the grid with 0.8nm size
grid_size=9
voxels=int(box_side*2/grid_size)

#discard initial 2us 
start = 4000
end = -1
#every 250ns
step = 500

pore_size_collect=[]
for index, ts in enumerate(u.trajectory[start:end:step]): 
    condensate_cog=np.around(condensate.center_of_geometry())
    print(condensate_cog)
    #minim = condensate_cog - box_side
    #maxim = condensate_cog + box_side
    min_X=condensate_cog[0]-box_side
    max_X=condensate_cog[0]+box_side
    min_Y=condensate_cog[1]-box_side
    max_Y=condensate_cog[1]+box_side
    min_Z=condensate_cog[2]-box_side
    max_Z=condensate_cog[2]+box_side
    core_protein=u.select_atoms('prop {}<x and prop x<{} and prop  {}<y and prop y<{} and prop {}<z and prop z<{}'.format(min_X,max_X,min_Y,max_Y,min_Z,max_Z))
    print(len(core_protein))

    #initialize number_grid every frame
    number_grid = np.zeros((voxels, voxels, voxels))
    protein_coords = core_protein.positions 
    # Assign each protein atom to a voxel
    for atom in protein_coords:
        #IDX = ((atom - minim) / grid_size).astype(int)
        x_idx = int((atom[0] - min_X) / grid_size)
        y_idx = int((atom[1] - min_Y) / grid_size)
        z_idx = int((atom[2] - min_Z) / grid_size)
        number_grid[x_idx, y_idx, z_idx] += 1
    print(number_grid.shape)
    pore_pos=[]
    for ix in range(voxels):
        for iy in range(voxels):
            for iz in range(voxels):
                #not contain protein
                if number_grid[ix, iy, iz] ==0:
                    pore_pos.append([ix, iy, iz])
    pore_pos=np.array(pore_pos)
    print(pore_pos.shape)
    #np.savetxt('Pore_grid.xvg',pore_pos,delimiter='	')

    #cluster for pore determination        
    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=1.1, linkage='single').fit(pore_pos)  
    #clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=1.1, linkage='ward').fit(pore_pos)  
    print(clustering.labels_)
    unique, counts=np.unique(clustering.labels_, return_counts=True)
    pore_size_collect=pore_size_collect+counts.tolist()

pore_size_collect=np.array(pore_size_collect)
print(pore_size_collect.shape)
frame=len(u.trajectory[start:end:step])
voxel_size=0.729

fig, ax = plt.subplots(figsize=(6,4.7))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)   
binwidth=1
bin_min=3
bin_max=int(max(pore_size_collect))+2
#bin_max=50
bins_edge=[x-0.5 for x in range(bin_min, bin_max+binwidth, binwidth)]
counts, bins, bars=ax.hist(pore_size_collect,bins=bins_edge,edgecolor='white',density=False,alpha=0.8,color='#3d5a9c')
ax.spines[['right', 'top']].set_visible(False)  
ax.set_ylabel('Frequency',fontproperties=font_prop)
ax.set_xlabel('Water cavity volume $(nm^3)$',fontproperties=font_prop)
#plt.xticks([7,50],labels=[5*voxel_size,10*voxel_size,20*voxel_size,30*voxel_size,40*voxel_size,50*voxel_size],fontproperties=tick_prop)
#plt.xticks([7,50,100,150],labels=[2.401,17.15,34.3],fontproperties=tick_prop)
plt.xticks([3,10,20],labels=[2.187,7.29,14.58],fontproperties=tick_prop)
#plt.xticks([3,5,10,15,20],labels=[3*voxel_size,5*voxel_size,10*voxel_size,15*voxel_size,20*voxel_size],fontproperties=tick_prop)
plt.yticks([0.2*frame,0.4*frame,0.6*frame,0.8*frame,1.0*frame],labels=[0.2,0.4,0.6,0.8,1.0],fontproperties=tick_prop) #rescale the frame cumulation
plt.savefig('pore_cluster-size_Voxel9.png',dpi=600,bbox_inches='tight')
plt.show()

#calculate the overall cavity vol
bins_reshape=np.array([(bins[i]+bins[i+1])/2 for i in range(len(bins)-1)])
print(bins_reshape)
print(bins_reshape.shape)
print(counts)
print(counts.shape)
print(bars)
vol=0
for x in range(len(bins_reshape)):
    vol+=bins_reshape[x]*counts[x]
print('ensemble average cavity volume:{}nm3'.format(vol*voxel_size/frame))
with open ('cavity_vol_Voxel9.txt','w') as g:
    print('ensemble average cavity volume:{}nm3'.format(vol*voxel_size/frame),file=g)

