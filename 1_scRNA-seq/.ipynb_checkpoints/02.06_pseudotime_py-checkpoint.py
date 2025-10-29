# %% [markdown]
# # 02.06_Pseudotime
#
# Author:GM/HB <br>
# useful links: https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.dpt.html

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
 
 
sc.settings.verbosity = 3             
sc.settings.file_format_figs = 'svg'  
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

basename = '_CHIP27_'

refdir = os.path.join(basedir, 'references')

os.makedirs(os.path.join(basedir, 'output/figures', '09.01'), exist_ok=True)
figdir = os.path.join(basedir, 'output/figures', '09.01')

os.makedirs(os.path.join(basedir, 'output/tables', '09.01'), exist_ok=True)
tabdir = os.path.join(basedir, 'output/tables', '09.01')


# %% [markdown]
# ### def functions

# %%
def plot_root_cell(dat, emb, idx, ax, color='cyan'):
    x = dat.obsm[emb][idx][0]
    y = dat.obsm[emb][idx][1]
    
    # root cell name
    # dat.obs.iloc[idx].name
    
    ax.scatter(x=x, y=y, s=20, c=color, marker='*')
    #ax.annotate('root', xy=(x, y), xytext=(x, y+1), color='red', size=50)
    plt.draw()


# %% [markdown]
# ### def bin parameter

# %%
n_bins = 18

# %% [markdown]
# ### Load data

# %%
# %%time
data= sc.read('20241120_CHIP27_after_lognorm_BBKNN_Kwok_Zeng_cellcycle_clustering1115_1120')

# %%
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))

ax = ax.ravel()
sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data, color='cluster_1120', ax=ax[1], show=False, legend_loc="on data" )

# %% [markdown]
# ## Define root cell

# %%
np.argmin( data[data.obs['cluster_1120'] == '2_5'].obsm['X_umap'][:,1] )

# %%
name2 = data[data.obs['cluster_1120'] == '2_5'].obs.iloc[3354].name

ind2 = data.obs.index.get_loc(name2)

# %%
basis = 'X_umap' 
clustering = 'cluster_1120'
root_cluster = '2_5'

ax = sc.pl.embedding(data, basis=basis, show=False)

sc.pl.embedding(data[data.obs[clustering] == root_cluster],
                color=clustering,
                basis=basis, 
                ax=ax, show=False)
plt.draw()

plot_root_cell(data, basis, ind2, ax, color='cyan')

# %%
data.obs.iloc[ind2]

# %%
data.uns['iroot'] = ind2

# %% [markdown]
# ## Compute pseudotime

# %%
sc.tl.dpt(data)

# %%
data.obs['dpt_pseudotime'].min(), data.obs['dpt_pseudotime'].max()

# %%
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 3))

ax = ax.ravel()

sc.pl.umap(data, color='cluster_1120', ax=ax[0], show=False, legend_loc="on data" )
sc.pl.umap(data, color='dpt_pseudotime', ax=ax[1], show=False, cmap='YlOrRd', vmax=0.6)#max is 1

# %%
sc.pl.violin(
    data,
    keys=["dpt_pseudotime"],
    groupby="cluster_1120",
    rotation=45
)

# %% [markdown]
# ### export adata with pseudotime values
# To do pretty plots in R

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok_Zeng_cellcycle_clustering1115_1120_pseudotime'
print(filename)
print('\n')

sc.write(filename, data)

# %%
data.obs.to_csv(os.path.join(tabdir, prefix+basename+"dataobs_with_pseudotime.csv"))
