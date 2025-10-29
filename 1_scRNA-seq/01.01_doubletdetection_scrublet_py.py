
# %% [markdown]
# # 01.01_Doublet detection and exclusion
# Author: GM

# %% [markdown]
# ## Notebook setup

# %%
# %matplotlib widget
# %load_ext watermark
 
import warnings
warnings.filterwarnings('ignore')
 
import os, sys, json, operator, getpass
from pathlib import Path
from datetime import datetime
 
import numpy as np
import pandas as pd
import scanpy as sc
 
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
 
import ipywidgets as widgets

import scrublet as scr
from collections import Counter
import anndata as ann
 
sc.settings.verbosity = 3             # show some output
sc.settings.file_format_figs = 'svg'  # set this to 'svg' (notebook) or 'pdf' (files) if you want vector graphics
sc.settings.savefigs = False
 
# plt.rcParams['font.family'] = 'sans-serif'
# plt.rcParams['font.sans-serif'] = 'Arial'
# plt.rc('font', size=14)
 
home = str(Path.home())
user = getpass.getuser()
 
# %watermark --iversions

# %%
basedir = os.path.join(laurenti, user, project)
 
sc.settings.writedir = os.path.join(basedir, 'output/objects/')

# %%
datadir = os.path.join(basedir, 'processed')

# %%
basename = '_CHIP27_'

# %% [markdown]
# ### Concatenate raw data
#
# ---
#

# %%
now = datetime.now()
prefix = now.strftime('%Y%m%d')

with open('/.singularity.d/labels.json') as fh:
    singularity = json.load(fh)
    
singularity['Version']
os.environ['SLURM_JOB_PARTITION'
os.environ['SLURM_NPROCS']

# %%
samples = ['SITTA2', 'SITTB2', 'SITTC2', 'SITTD2', 'SITTE2', 'SITTF2']

# %%
# %%capture

libraries = {}

for sample in samples:
    libraries[sample] = sc.read_10x_h5(os.path.join(datadir, sample, 'outs', 
                                                    'filtered_feature_bc_matrix.h5'), 'hg19')
    libraries[sample].var_names_make_unique()
    libraries[sample].obs['library'] = sample

# %% [markdown]
# explore data from one library:

# %%
data = libraries['SITTA2'].concatenate( libraries['SITTB2'], libraries['SITTC2'], libraries['SITTD2'], libraries['SITTE2'], libraries['SITTF2'])

# %%
data.obs.index = data.obs.apply(lambda x : x.library + '_' + x.name[0:-4] , axis = 1 )

# %%
print("Raw data has %d genes in %d cells" %(data.X.shape[1], data.X.shape[0]))

sample_count = Counter(data.obs['library'])

for sample in samples:
    print(sample, 'has %d cells' %sample_count[sample])

# %% [markdown]
# #### export concat anndata

# %%
# %%time

filename = prefix + basename + 'concatenated_filtered_gene_bc_expression'
print(filename)

sc.write(filename, data)

# %% [markdown]
# ### automatic Scrublet detection

# %% [markdown]
# #### import concat anndata

# %%
data = sc.read(os.path.join(basedir, 'output','objects', '20240320_CHIP27_concatenated_filtered_gene_bc_expression'))

# %%
samples = ['SITTA2', 'SITTB2', 'SITTC2', 'SITTD2', 'SITTE2', 'SITTF2']

# %%
# %%time

scrub = {}
doublet_scores = {}
predicted_doublets = {}

# emr = {'SITTA2': 0.058, #emr removed because it performs worse for SITTA2 (way less doublets removed), as the unbiased approach estimate higher number of duplets for this cell n
#        'SITTB2': 0.058,
#        'SITTC2': 0.027,
#        'SITTD2': 0.058,
#        'SITTE2': 0.043,
#        'SITTF2': 0.048,
#       }

for sample in samples:
    print(sample)

    scrub[sample] = scr.Scrublet( data[np.array(data.obs['library'] == sample), :].X#, 
                                  # expected_doublet_rate=emr[sample] 
                                )

        
    doublet_scores[sample], predicted_doublets[sample] = scrub[sample].scrub_doublets()
    print('\n\n')

# %% [markdown]
# ### Manual threshold adjustment

# %% [markdown]
# #### SITTA2

# %%
sample = 'SITTA2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.43

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# #### SITTB2

# %%
sample = 'SITTB2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.38

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %%
sample = 'SITTB2'

# %%
scrub[sample].threshold_

# %%
#adj thresh
thres = 0.38

# %%
x = scrub[sample].plot_histogram()

# %%
thres = scrub[sample].threshold_

# %% jupyter={"outputs_hidden": true}
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')
axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )


# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %% jupyter={"outputs_hidden": true}
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# #### SITTC2

# %%
sample = 'SITTC2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
#with different binning
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.43

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# #### SITTD2

# %%
sample = 'SITTD2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.37

# %%

fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# #### SITTE2

# %%
sample = 'SITTE2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.3

# %%

fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# #### SITTF2

# %%
sample = 'SITTF2'

# %% [markdown]
# **automatic threshold**

# %%
scrub[sample].threshold_ 

# %%
x = scrub[sample].plot_histogram()

# %% [markdown]
# **manually adjusted**

# %%
#visualise automatic threshold first
thres = scrub[sample].threshold_

# %%
fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
#adj thresh
thres = 0.35

# %%

fig_size = (8,3)
fig, axs = plt.subplots(1, 2, figsize = fig_size)
# axes[0].plot(x1, y1)
# axes[1].plot(x2, y2)
# fig.tight_layout()
x = scrub[sample].doublet_scores_obs_
y = scrub[sample].doublet_scores_sim_
num_bins = 200

axs[0].hist(x, num_bins, facecolor='blue', alpha=0.5)
axs[0].axvline(thres, color='red')

axs[1].hist(y, num_bins, facecolor='blue', alpha=0.5)
axs[1].axvline(thres, color='red')

fig.tight_layout()

# %%
sum( scrub[sample].call_doublets(threshold=thres) )

# %%
x = scrub[sample].plot_histogram()

# %%
scrub[sample].set_embedding('UMAP', scr.get_umap(scrub[sample].manifold_obs_, 10, min_dist=0.3))

# %%
scrub[sample].plot_embedding('UMAP', order_points=True)

# %%
predicted_doublets[sample] = scrub[sample].call_doublets(threshold=thres)

# %% [markdown]
# ---

# %% [markdown]
# ### Filter anndata

# %%
data_doublets = os.path.join( basedir , 'output', 'doublets')

if not os.path.exists( data_doublets ):
    os.makedirs( data_doublets )

# %%
for key in doublet_scores:
    np.savetxt(os.path.join(data_doublets, prefix+'_'+key+'_doublet_scores.txt'), doublet_scores[key])

# %%
doublet_scores_list = []
for key in doublet_scores:
    #print(key)
    doublet_scores_list += list(doublet_scores[key])
    
data.obs['doublet_score'] = doublet_scores_list

# %%
len(doublet_scores_list)

# %%
data.shape

# %%
predicted_doublets_mask = []

for key in predicted_doublets:
    #print(key)
    predicted_doublets_mask += list(predicted_doublets[key])

# %%
predicted_singletons_mask = [not i for i in predicted_doublets_mask]

# %%
data = data[np.array(predicted_singletons_mask), :].copy()

doublet_count = len(predicted_singletons_mask) - sum(predicted_singletons_mask)

print( 'Removing %d cells due to doublet scoring' %(doublet_count) )

# %%
for sample in samples:
    print(sample, ':', sum(predicted_doublets[sample]))

# %% [markdown]
# ### export anndata

# %%
print("Raw data has %d genes in %d cells" %(data.X.shape[1], data.X.shape[0]))

sample_count = Counter(data.obs['library'])

for sample in samples:
    print(sample, 'has %d cells' %sample_count[sample])

# %%
# %%time

filename = prefix + basename + 'filtered_gene_bc_expression_minus_'+str(doublet_count)+'_putative_doublets'

print(filename)

sc.write(filename, data)
