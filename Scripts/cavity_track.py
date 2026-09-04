import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
import pickle
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from statistics import mean

from collections import defaultdict
from typing import List, Set, Tuple, Dict, Any

def match_cavities(prev_cavities: List[Set[int]],
                   curr_cavities: List[Set[int]],
                   overlap_frac: float = 0.3) :
    """
    Match cavities between two consecutive frames using voxel overlap.

    Parameters
    ----------
    prev_cavities : list of sets
        Each set contains voxel indices of a cavity in the previous frame.
    curr_cavities : list of sets
        Each set contains voxel indices of a cavity in the current frame.
    overlap_frac : float
        Minimum fractional overlap relative to the smaller cavity size to consider a match.
    """
    n_prev = len(prev_cavities)
    n_curr = len(curr_cavities)

    # Compute all overlaps
    overlaps = {}
    for i, prev in enumerate(prev_cavities):
        for j, curr in enumerate(curr_cavities):
            overlap = len(prev & curr)
            if overlap > 0:
                # Fraction relative to the smaller cavity size
                min_size = min(len(prev), len(curr))
                frac = overlap / min_size if min_size > 0 else 0.0
                if frac >= overlap_frac:
                    overlaps[(i, j)] = overlap

    # Greedy one-to-one matching: sort by overlap descending, assign best available
    matches = []
    used_prev = set()
    used_curr = set()
    for (i, j), overlap in sorted(overlaps.items(), key=lambda x: -x[1]):
        if i not in used_prev and j not in used_curr:
            matches.append((i, j))
            used_prev.add(i)
            used_curr.add(j)

    unmatched_prev = [i for i in range(n_prev) if i not in used_prev]
    unmatched_curr = [j for j in range(n_curr) if j not in used_curr]

    return matches, unmatched_prev, unmatched_curr


def track_cavities(all_frames_cavities: List[List[Set[int]]],
                   overlap_frac: float = 0.3) :
    """
    Track cavities across all frames.

    Parameters
    ----------
    all_frames_cavities : list of list of sets
        For each frame, a list of cavities (each cavity is a set of voxel indices).
    overlap_frac : float
        Overlap fraction for matching.
    """
    n_frames = len(all_frames_cavities)
    if n_frames == 0:
        return {'tracks': [], 'formation_events': [], 'disappearance_events': [], 'n_frames': 0}

    # Map from persistent track ID to list of (frame_idx, local cavity index)
    track_dict = defaultdict(list)
    # Map from (frame_idx, local_cavity_idx) -> track ID
    frame_cavity_to_track = {}
    next_track_id = 0

    # Initialize with first frame: all cavities are formations
    first_frame_cavities = all_frames_cavities[0]
    for j, cavity in enumerate(first_frame_cavities):
        track_id = next_track_id
        next_track_id += 1
        track_dict[track_id].append((0, j))
        frame_cavity_to_track[(0, j)] = track_id

    formation_events = [(0, j) for j in range(len(first_frame_cavities))]  # all appear at frame 0
    disappearance_events = []

    # Process consecutive frames
    for frame_idx in range(1, n_frames):
        prev_cavities = all_frames_cavities[frame_idx - 1]
        curr_cavities = all_frames_cavities[frame_idx]

        matches, unmatched_prev, unmatched_curr = match_cavities(prev_cavities, curr_cavities, overlap_frac)

        # Handle matches: continue the track of the previous cavity
        for prev_local_idx, curr_local_idx in matches:
            prev_global_key = (frame_idx - 1, prev_local_idx)
            if prev_global_key in frame_cavity_to_track:
                track_id = frame_cavity_to_track[prev_global_key]
                track_dict[track_id].append((frame_idx, curr_local_idx))
                frame_cavity_to_track[(frame_idx, curr_local_idx)] = track_id
            else:
                # Should not happen if tracking is consistent, but fallback: create new track
                track_id = next_track_id
                next_track_id += 1
                track_dict[track_id].append((frame_idx, curr_local_idx))
                frame_cavity_to_track[(frame_idx, curr_local_idx)] = track_id
                formation_events.append((frame_idx, curr_local_idx))

        # Handle unmatched previous cavities: they disappeared
        for prev_local_idx in unmatched_prev:
            prev_global_key = (frame_idx - 1, prev_local_idx)
            if prev_global_key in frame_cavity_to_track:
                disappearance_events.append((frame_idx - 1, prev_local_idx))
                # Track will be marked with death_frame later

        # Handle unmatched current cavities: they are new formations
        for curr_local_idx in unmatched_curr:
            track_id = next_track_id
            next_track_id += 1
            track_dict[track_id].append((frame_idx, curr_local_idx))
            frame_cavity_to_track[(frame_idx, curr_local_idx)] = track_id
            formation_events.append((frame_idx, curr_local_idx))

    # Build track summaries
    tracks = []
    for track_id, occurrences in track_dict.items():
        # Sort by frame index (should already be sorted)
        occurrences.sort(key=lambda x: x[0])
        birth_frame = occurrences[0][0]
        death_frame = occurrences[-1][0]
        lifetime_frames = death_frame - birth_frame + 1
        frames_present = [occ[0] for occ in occurrences]
        tracks.append({
            'id': track_id,
            'birth_frame': birth_frame,
            'death_frame': death_frame,
            'lifetime_frames': lifetime_frames,
            'frames_present': frames_present,
            'occurrences': occurrences  # list of (frame_idx, local_idx)
        })

    # Sort tracks by birth_frame for convenience
    tracks.sort(key=lambda t: t['birth_frame'])

    return {
        'tracks': tracks,
        'formation_events': formation_events,
        'disappearance_events': disappearance_events,
        'n_frames': n_frames
    }


def compute_lifetimes(tracks: List[Dict[str, Any]], dt: float = 1.0):
    """
    Convert lifetime_frames to physical time using frame interval dt.

    Parameters
    ----------
    tracks : list of track dicts from track_cavities()
    dt : frame time step (unit: ns)
    """
    return [t['lifetime_frames'] * dt for t in tracks]


def compute_spatial_persistence(all_frames_cavities: List[List[Set[int]]],
                                grid_shape: Tuple[int, int, int]):
    """
    Compute voxel occupancy probability (spatial persistence map).

    Parameters
    ----------
    all_frames_cavities : list of list of sets, each set contains voxel linear indices.
    grid_shape : tuple (nx, ny, nz) of the voxel grid dimensions.
    """
    nx, ny, nz = grid_shape
    total_voxels = nx * ny * nz
    occupancy_count = np.zeros(grid_shape, dtype=np.float32)
    frame_count = 0

    for frame_cavities in all_frames_cavities:
        frame_count += 1
        # Mark all voxels that are part of any cavity in this frame
        for cavity in frame_cavities:
            for voxel_idx in cavity:
                occupancy_count[voxel_idx] += 1.0

    occupancy_map = occupancy_count.reshape(grid_shape) / frame_count
    return occupancy_map



trajectory_dir='/grain/liguo/MYC/IDP-conf-change/Project4_Condensate_reconstitution/Pbody_yeast/condensate-replica2/'

u = mda.Universe(trajectory_dir+'production9-Verlet-pbc-protein.gro',trajectory_dir+'production9-Verlet-pbc-fit.xtc')
condensate=u.select_atoms('name BB or name SC*')  

#based on condensate COG, define a core cubix box with side=4.8*2nm, avoiding sampling too much region outside condensate, because condensate is not a perfect spherical droplet
box_side=48
#define the grid with 0.8nm size
grid_size=8
voxels=int(box_side*2/grid_size)

#discard initial 2us 
start = 80
end = -1
#every 25ns
step = 1

all_frames=[]
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

    clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=1.1, linkage='single').fit(pore_pos)  
    print(clustering.labels_)
    unique, counts=np.unique(clustering.labels_, return_counts=True)
    clusters=np.asarray((unique, counts)).T
    #filter the small clusters with threshold size 5
    cluster_filter=clusters[clusters[:,1]>4]
    print(cluster_filter)

    cavity_set=[]
    for cluster_index in cluster_filter[:,0]:
        cluster_pos=[]
        for i, x in enumerate(clustering.labels_):
            if x == cluster_index:
                grid_pos=pore_pos[i]
                cluster_pos.append(tuple(grid_pos))
        cavity_set.append(set(cluster_pos))
    print(cavity_set)
    all_frames.append(cavity_set)
print(all_frames)
print(len(all_frames))
pickle.dump(all_frames, open( "replica2/cavity_frames.p", "wb" ))


# Track with overlap fraction 0.3
results = track_cavities(all_frames, overlap_frac=0.3)
pickle.dump(results, open( "replica2/cavity_track.p", "wb" ))
print("Number of tracks:", len(results['tracks']))
for track in results['tracks']:
    print(f"Track {track['id']}: birth frame {track['birth_frame']}, "
              f"death frame {track['death_frame']}, lifetime {track['lifetime_frames']} frames")

print("\nFormation events (frame, local_idx):", results['formation_events'])
print('Formation events number:{}'.format(len(results['formation_events'])))
print("Disappearance events (frame, local_idx):", results['disappearance_events'])
print('Disappearance events number:{}'.format(len(results['disappearance_events'])))

# Lifetimes, step size 25ns
lifetimes = compute_lifetimes(results['tracks'], dt=25.0)
print("Lifetimes:", lifetimes)
np.savetxt('replica2/lifetimes.xvg', np.array(lifetimes),delimiter='	')

occupancy = compute_spatial_persistence(all_frames, (12,12,12))
pickle.dump(occupancy, open( "replica2/spatial_occupancy.p", "wb" ))
print(occupancy.shape)

coords = np.argwhere(occupancy != 0)  # Exclude value = 0
values = occupancy[occupancy != 0]
print(values)
# Extract X, Y, Z coordinates
x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(projection='3d')
font_path = '/grain/liguo/Script/python-plot/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=16)
tick_prop = font_manager.FontProperties(fname=font_path, size=14)
legend_prop = font_manager.FontProperties(fname=font_path, size=14)  
sc = ax.scatter(x, y, z, c=values, s=30,cmap='viridis')
cbar=fig.colorbar(sc, ax=ax)
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_title('Fraction',fontproperties=tick_prop)
ax.set_xlabel('X Axis',fontproperties=font_prop)
ax.set_ylabel('Y Axis',fontproperties=font_prop)
ax.set_zlabel('Z Axis',fontproperties=font_prop)
plt.yticks(fontproperties=tick_prop)    
plt.xticks(fontproperties=tick_prop)   
ax.set_zticks(np.arange(0,12,2),np.arange(0,12,2),fontproperties=tick_prop)
plt.savefig('replica2/spatial-persistance.png',dpi=600,bbox_inches='tight')
plt.show()
