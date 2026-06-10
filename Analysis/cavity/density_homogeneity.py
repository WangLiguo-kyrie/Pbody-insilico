#! /usr/bin/python

import sys
import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate_small/'

u = mda.Universe(trajectory_dir+'production8-Verlet-30us.gro',trajectory_dir+'production8-Verlet-pbc-skip50.xtc')
condensate=u.select_atoms('name BB or name SC*')  

#based on condensate COG, define a core cubix box with side=5.0*2nm, avoiding sampling too much region outside condensate, because condensate is not a perfect spherical droplet
box_side=49
#define the voxels with 1nm size
grid_size=7
voxels=int(box_side*2/grid_size)

#discard initial 2us 
start = 80
end = -1
#every 250ns
step = 10

number_grid_collect=[]
density_grid_collect=[]
for index, ts in enumerate(u.trajectory[start:end:step]): 
    condensate_cog=np.around(condensate.center_of_geometry())
    print(condensate_cog)
    #based on condensate COG, define a core cubix box with side=5.5*2nm, avoiding sampling too much region outside condensate, because condensate is not a perfect spherical droplet
    min_X=condensate_cog[0]-box_side
    max_X=condensate_cog[0]+box_side
    min_Y=condensate_cog[1]-box_side
    max_Y=condensate_cog[1]+box_side
    min_Z=condensate_cog[2]-box_side
    max_Z=condensate_cog[2]+box_side
    core_water=u.select_atoms('name W and prop {}<x and prop x<{} and prop  {}<y and prop y<{} and prop {}<z and prop z<{}'.format(min_X,max_X,min_Y,max_Y,min_Z,max_Z))
    print(len(core_water))
    print(core_water)

    #initialize number_grid every frame
    number_grid = np.zeros((voxels, voxels, voxels))
    coords = core_water.positions 
    # Assign each protein atom to a voxel
    for atom in coords:
        x_idx = int((atom[0] - min_X) / grid_size)
        y_idx = int((atom[1] - min_Y) / grid_size)
        z_idx = int((atom[2] - min_Z) / grid_size)
        number_grid[x_idx, y_idx, z_idx] += 1
    voxel_volume = (grid_size/10)**3  
    # Calculate density (beads per cubic Angstrom) in each voxel
    density_grid = number_grid/voxel_volume 
    
    density_flat = density_grid.flatten()
    number_flat = number_grid.flatten()
    density_grid_collect.append(density_flat)
    number_grid_collect.append(number_flat)

density_grid_collect=np.array(density_grid_collect)
number_grid_collect=np.array(number_grid_collect)
print(number_grid_collect.shape)
print(density_grid_collect.shape)
density_flat_all = density_grid_collect.flatten()
print(density_flat_all)
number_flat_all = number_grid_collect.flatten()
density_mean = np.mean(density_flat_all)
density_std = np.std(density_flat_all)
# Calculate the coefficient of variation (CV) to quantify homogeneity
cv = density_std / density_mean  # A smaller CV indicates more homogeneous density
print(f"Density Mean: {density_mean} beads/nm³")
print(f"Density Standard Deviation: {density_std} beads/nm³")
print(f"Density Coefficient of Variation: {cv}")
#with open('density_homogeneity.txt','w') as g:
#    print(f"Density Mean: {density_mean} beads/nm³",file=g)
#    print(f"Density Standard Deviation: {density_std} beads/nm³",file=g)
#    print(f"Density Coefficient of Variation: {cv}",file=g)


fig, ax = plt.subplots(figsize=(6,4.7))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)             
binwidth=1
bin_min=int(min(number_flat_all))
bin_max=int(max(number_flat_all))+1
print('max water number per voxel:{}'.format(int(max(number_flat_all))))
bins_edge=[x-0.5 for x in range(bin_min, bin_max+binwidth, binwidth)]
ax.hist(number_grid.flatten(),bins=bins_edge,edgecolor='white',density=True,alpha=0.8,color='#3d5a9c')
ax.spines[['right', 'top']].set_visible(False)  
ax.set_ylabel('Probability',fontproperties=font_prop)
ax.set_xlabel('Water beads per voxel',fontproperties=font_prop)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop) 
plt.savefig('water_hist_Voxel7.png',dpi=600,bbox_inches='tight')
plt.show()
