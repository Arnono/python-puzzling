# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 13:53:30 2020

@author: arnow
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import os
import gdal
from gdalconst import GA_ReadOnly
import tables
from h5py import File
import scipy.ndimage
from scipy import stats
from PIL import Image

     
base_dir = 'C:/school_data/spat_numerical/project/'


# =============================================================================
# Generate sample data and Visualize it
# =============================================================================

dx, dy = 0.05, 0.05
y, x = np.mgrid[slice(1, 10 + dy, dy),
                slice(1, 10 + dx, dx)]
z = np.sin(x)**2 + np.cos(10 + y*x/4) * np.cos(x)
plt.imshow(z, cmap = 'terrain')
plt.imshow(z)
template = z[1:181,1:181]
new_template =(template*3)**2
new_template_rot = np.rot90(new_template)
new_template_rot2 = np.rot90(new_template_rot)
new_template_merge = (new_template_rot + new_template_rot2/2 + new_template)*2
template =new_template_merge

#define tile sizes
split_size = int(template.shape[0]/4)
temp_size = (template.shape[0],template.shape[0])
puzzle_size = (template.shape[0]/split_size,template.shape[0]/split_size)
index = np.arange(0,int(puzzle_size[0]**2))

#View
plt.imshow(new_template_merge, cmap = 'terrain')
plt.colorbar()


# =============================================================================
# Stack tile files for easy access
# =============================================================================
  
FILTERS = tables.Filters(complib='zlib', complevel=5)
h5file = tables.open_file(base_dir +'puzzle_stack8.hdf', mode='w', 
                          title='puzzle', filters=FILTERS)

#Load hdf file with data  
SIZE = (split_size, split_size, index.shape[0])
ds = h5file.create_carray(h5file.root, 'tiles', tables.Float32Atom(),(SIZE))
tile_index = h5file.create_carray(h5file.root, 'index', tables.Int16Atom(),(SIZE[2],))
tile_index[:] = index 

#load simulated image
row = 0
stack_index = 0
for i in range(int(puzzle_size[0])): 
    col = 0
    for j in range(int(puzzle_size[1])):
        tile = template[row:row+45,col:col+45]
        ds[:,:,stack_index] = tile
        print(row,col,stack_index)
        col+=45
        stack_index += 1
    row+=45
      
h5file.close()    

#read stored stack file
puzzle_in = File(base_dir + 'puzzle_stack4.hdf' , 'r')
index_tle = puzzle_in['index'][:]

#Check data
samp_tle = puzzle_in['tiles'][:,:,1]
plt.imshow(samp_tle , cmap = 'terrain')

mini = np.min(samp_tle)
maxi = np.max(samp_tle)


# =============================================================================
# Create tile boarders array
# =============================================================================

sides_arr = np.ones((samp_tle.shape[0],(index_tle.shape[0]*4)))
n = 0
for i in range ((16)):
    print(i)
    samp_tle = puzzle_in['tiles'][:,:,i]
    s1,s2,s3,s4 = (samp_tle[0:1,:].T,samp_tle[:,-1:],samp_tle[-1:,:].T,
                   samp_tle[:,0:1]) #from top clockwise
    s1_mean,s2_mean,s3_mean,s4_mean = np.mean(s1,axis=1),np.mean(s2,axis=1),np.mean(s3,axis=1),np.mean(s4,axis=1)
    sides_arr[:,n] = s1_mean
    sides_arr[:,n+1] = s2_mean
    sides_arr[:,n+2] = s3_mean
    sides_arr[:,n+3] = s4_mean
    n+=4
#visualise boarders
plt.imshow(sides_arr , cmap = 'terrain')
plt.tight_layout

# =============================================================================
# Check possible allignments
# =============================================================================

link_side = -1
n = 0
ref_line_ind = 0

#for i in range ((sides_arr.shape[1])):  
for i in range ((20)): 
    side = i%4
    if side == 0:
        n +=1       
    highest_corr_coeff = 0
    for j in range ((64)):
        side2 = j%4
        if (side ==0 and side2 ==2) or (side ==2 and side2 ==0) or (side ==1 and side2 ==3) or (side ==3 and side2 ==1):
            if (j != ref_line_ind):
                ref_line = sides_arr[:,ref_line_ind] 
                check_line = sides_arr[:,j]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.94 and corr_coeff > highest_corr_coeff :      #if corr_coeff > 0.94
                    print('')
                    print(side,side2)        
                    print('tile',n,'side',side,'matches Tile',int(j/4)+1,
                          'side',(j%4))
                    print(ref_line_ind, corr_coeff,'match sighted')                  
                    highest_corr_coeff = corr_coeff 
    ref_line_ind +=1


# =============================================================================
# Initialize reconstruction
# =============================================================================

tle_zero = puzzle_in['tiles'][:,:,0]
recon_array = np.zeros(((sides_arr.shape[0]*4),(sides_arr.shape[0]*4)))
recon_array[0:split_size,0:split_size] = tle_zero #insert reference tile

#Check image
plt.imshow(recon_array, cmap = 'terrain')

n = 0
ref_line_ind = 0 #from boarders array
recon_tile_number = 0
recon_row = 0
recon_col = 0
direction = 1
dir_set = -1
 
ref_line = tle_zero[:,-1] 

#for i in range ((sides_arr.shape[1])):  
for i in range ((index.shape[0]-1)):                       
    highest_corr_coeff = 0
    print(int(i/index.shape[0]-1*100),'% ' ,end='')
    j=0
    k=0
    fig,axs = plt.subplots(figsize=(7.68, 7.68) ) 
    axs.imshow(recon_array, cmap = 'terrain', vmin = mini, vmax = maxi )
    axs.set_title('Tiling', fontsize=16, fontweight='bold')
    axs.axes.get_xaxis().set_visible(False)
    axs.axes.get_yaxis().set_visible(False)
    axs.axes.get_xaxis().set_visible(False)
    axs.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig('%04d.png' % (i+1))
    plt.close('all')
    if dir_set ==1:
        direction = 1
    if dir_set ==0:
        direction = 0       
    if direction == 1:       
        for j in range ((64)):              #search through the boarders array
            recon_row = int(j/16)+1
            recon_col = int((j%16)/4)+1
            y_mod = j%4
            if (y_mod ==3 and direction ==1):
                check_line = sides_arr[:,j]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.94 and corr_coeff > highest_corr_coeff :  #0.9
                    highest_corr_coeff = corr_coeff
                    fit_tile_indx = int(j/4)
                    new_ref_tle = puzzle_in['tiles'][:,:,fit_tile_indx]
                    ref_line = new_ref_tle[:,split_size-1]
                    recon_array[(recon_row-1)*split_size :recon_row*split_size , 
                                split_size*(recon_col-1) : split_size*recon_col] = new_ref_tle
                    if recon_col == 4:
                        dir_set =0
                        ref_line = recon_array[(split_size*recon_row)-1,0:split_size]                        
    if direction == 0:       
        for k in range ((64)):             #search through the boarders array
            recon_row = int(k/16)+1
            recon_col = int((k%16)/4)+1
            y_mod = k%4
            if (y_mod ==0 and direction ==0 and recon_col ==1):        
                check_line = sides_arr[:,k]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.94 and corr_coeff > highest_corr_coeff :   #0.9
                    highest_corr_coeff = corr_coeff
                    fit_tile_indx = int(k/4)
                    new_ref_tle = puzzle_in['tiles'][:,:,fit_tile_indx]
                    recon_array[(recon_row-1)*split_size :recon_row*split_size , 
                                split_size*(recon_col-1) : split_size*recon_col] = new_ref_tle
            if recon_row == 4 and k == 63:
                dir_set =1
                ref_line = new_ref_tle[:,split_size-1]






# =============================================================================
#
# Real data 
#
# =============================================================================

base_dir = 'C:/school_data/spat_numerical/project/'

#define tile sizes
split_size = 2000
index = np.arange(0,16)

# =============================================================================
# Stack tile files for easy access
# =============================================================================
  
#Initiate stack file
FILTERS = tables.Filters(complib='zlib', complevel=5)
h5file = tables.open_file(base_dir +'Berlin_puzzle_stack7.hdf', mode='w', 
                          title='puzzle', filters=FILTERS)

#Load hdf file with data  
SIZE = (split_size, split_size, index.shape[0])
ds = h5file.create_carray(h5file.root, 'tiles', tables.Float32Atom(),(SIZE))
tile_index = h5file.create_carray(h5file.root, 'index', tables.Int16Atom(),(SIZE[2],))
tile_index[:] = index 


#load real image
image_tiles_dir= 'C:/school_data/spat_numerical/outs_ims/s2a_tiles_berlin_arranged/'
for i, f in enumerate(os.listdir(image_tiles_dir)):
    path = os.path.join(image_tiles_dir,f)
    im_tile = Image.open(path).convert('L')          
    np_image = np.array(im_tile)
    ds[:,:,i] = np_image
    print(path,' ',i)
        
h5file.close()  
   
#read stored stack file
puzzle_in = File(base_dir + 'Berlin_puzzle_stack7.hdf' , 'r')
index_tle = puzzle_in['index'][:]

#Check data
samp_tle = puzzle_in['tiles'][:,:,1]
plt.imshow(samp_tle , cmap = 'terrain')

mini = np.min(samp_tle)
maxi = np.max(samp_tle)

# =============================================================================
# Create tile boarders array
# =============================================================================

sides_arr = np.ones((samp_tle.shape[0],(index_tle.shape[0]*4)))
n = 0
for i in range ((16)):
    print(i)
    samp_tle = puzzle_in['tiles'][:,:,i]
    s1,s2,s3,s4 = (samp_tle[0:1,:].T,samp_tle[:,-1:],samp_tle[-1:,:].T,
                   samp_tle[:,0:1]) #from top clockwise
    s1_mean,s2_mean,s3_mean,s4_mean = np.mean(s1,axis=1),np.mean(s2,axis=1),np.mean(s3,axis=1),np.mean(s4,axis=1)
    sides_arr[:,n] = s1_mean
    sides_arr[:,n+1] = s2_mean
    sides_arr[:,n+2] = s3_mean
    sides_arr[:,n+3] = s4_mean
    n+=4
#visualise boarders
plt.imshow(sides_arr , cmap = 'terrain')
plt.tight_layout

# =============================================================================
# Check possible allignments
# =============================================================================

link_side = -1
n = 0
ref_line_ind = 0

#for i in range ((sides_arr.shape[1])):  
for i in range ((20)): 
    side = i%4
    if side == 0:
        n +=1       
    highest_corr_coeff = 0
    for j in range ((64)):
        side2 = j%4
        if (side ==0 and side2 ==2) or (side ==2 and side2 ==0) or (side ==1 and side2 ==3) or (side ==3 and side2 ==1):
            if (j != ref_line_ind):
                ref_line = sides_arr[:,ref_line_ind] 
                check_line = sides_arr[:,j]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.94 and corr_coeff > highest_corr_coeff :      
                    print('')
                    print(side,side2)        
                    print('tile',n,'side',side,'matches Tile',int(j/4)+1,
                          'side',(j%4))
                    print(ref_line_ind, corr_coeff,'match sighted')                  
                    highest_corr_coeff = corr_coeff 
    ref_line_ind +=1


# =============================================================================
# Initialie reconstruction
# =============================================================================

tle_zero = puzzle_in['tiles'][:,:,0]
recon_array = np.zeros(((sides_arr.shape[0]*4),(sides_arr.shape[0]*4)))
recon_array[0:split_size,0:split_size] = tle_zero #insert reference tile

#Check image
plt.imshow(recon_array, cmap = 'terrain')

n = 0
ref_line_ind = 0 #from boarders array
recon_tile_number = 0
recon_row = 0
recon_col = 0
direction = 1
dir_set = -1
 
ref_line = tle_zero[:,-1] 

#for i in range ((sides_arr.shape[1])):  
for i in range ((index.shape[0]-1)):                       
    highest_corr_coeff = 0
    print(int(i/index.shape[0]-1*100),'% ' ,end='')
    j=0
    k=0
    fig,axs = plt.subplots(figsize=(7.68, 7.68) ) 
    axs.imshow(recon_array, cmap = 'terrain', vmin = mini, vmax = maxi )
    axs.set_title('Tiling', fontsize=16, fontweight='bold')
    axs.axes.get_xaxis().set_visible(False)
    axs.axes.get_yaxis().set_visible(False)
    axs.axes.get_xaxis().set_visible(False)
    axs.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    plt.savefig('%04d.png' % (i+1))
    plt.close('all')
    if dir_set ==1:
        direction = 1
    if dir_set ==0:
        direction = 0       
    if direction == 1:       
        for j in range ((64)):              #search through the boarders array
            recon_row = int(j/16)+1
            recon_col = int((j%16)/4)+1
            y_mod = j%4
            if (y_mod ==3 and direction ==1):
                check_line = sides_arr[:,j]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.9 and corr_coeff > highest_corr_coeff :  #0.9
                    highest_corr_coeff = corr_coeff
                    fit_tile_indx = int(j/4)
                    new_ref_tle = puzzle_in['tiles'][:,:,fit_tile_indx]
                    ref_line = new_ref_tle[:,split_size-1]
                    recon_array[(recon_row-1)*split_size :recon_row*split_size , 
                                split_size*(recon_col-1) : split_size*recon_col] = new_ref_tle
                    if recon_col == 4:
                        dir_set =0
                        ref_line = recon_array[(split_size*recon_row)-1,0:split_size]                        
    if direction == 0:       
        for k in range ((64)):             #search through the boarders array
            recon_row = int(k/16)+1
            recon_col = int((k%16)/4)+1
            y_mod = k%4
            if (y_mod ==0 and direction ==0 and recon_col ==1):        
                check_line = sides_arr[:,k]
                corr_coeff = stats.pearsonr(ref_line,check_line)[0]
                if corr_coeff > 0.9 and corr_coeff > highest_corr_coeff :   #0.9
                    highest_corr_coeff = corr_coeff
                    fit_tile_indx = int(k/4)
                    new_ref_tle = puzzle_in['tiles'][:,:,fit_tile_indx]
                    recon_array[(recon_row-1)*split_size :recon_row*split_size , 
                                split_size*(recon_col-1) : split_size*recon_col] = new_ref_tle
            if recon_row == 4 and k == 63:
                dir_set =1
                ref_line = new_ref_tle[:,split_size-1]





















