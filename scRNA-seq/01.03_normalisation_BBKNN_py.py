# %% [markdown]
# # B_01.03_CHIP27_Normalisation
# Author: GM <br>

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
project = "CHIP27/"
 
basedir = os.path.join(laurenti, user, project)
 
sc.settings.writedir = os.path.join(basedir, 'output/objects/')


# %%
basename = '_CHIP27_'

# %% [markdown]
# ---

# %% [markdown]
# ### Environment info

# %%
now = datetime.now()
prefix = now.strftime('%Y%m%d')
print(prefix)

# %%
with open('/.singularity.d/labels.json') as fh:
    singularity = json.load(fh)
    
singularity['Version']
os.environ['SLURM_JOB_PARTITION']
os.environ['SLURM_NPROCS']

# %% [markdown]
# ---

# %% [markdown]
# ### Import data

# %%
# %%time
data=sc.read('20241003_CHIP27_postQC_beforenorm')

# %%
data.layers['counts'] = data.X.copy()

# %% [markdown]
# ---

# %% [markdown]
# #### 1. normalise (to target 1e4: 10,000)

# %%
p1 = sns.histplot(data.obs["total_counts"], bins=100, kde=False)

# %% [markdown]
# Visualise top 25 expressed genes:

# %%
fig, ax = plt.subplots(1,1, figsize=(6,4),dpi=150 )
sc.pl.highest_expr_genes(data, n_top=25, ax=ax, show=False)

# %%
# %%time
sc.pp.normalize_total(data, target_sum=1e4, exclude_highly_expressed=True)

# %% [markdown]
# #### 2. log transform (base2)

# %%
sc.pp.log1p(data, base=2)

# %%
data.layers["lognorm"] = data.X.copy()

# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
p1 = sns.histplot(data.layers["counts"].sum(1), bins=100, kde=False, ax=axes[0])
axes[0].set_title("Total counts")
p2 = sns.histplot(data.layers["lognorm"].sum(1), bins=100, kde=False, ax=axes[1])
axes[1].set_title("Total counts after shifted log norm")
plt.show()

# %% [markdown]
# ### Feature selection - highly variable genes
#

# %%
# %%time
#parameter changed:  and go lower in the dispersion to select fewer genes. 
sc.pp.highly_variable_genes(data, 
                            min_mean= 0.05, 
                            max_mean=7, min_disp=0.25)

# %%
sc.pl.highly_variable_genes(data)

# %%
data.var['highly_variable'].sum()

# %% [markdown]
# ### regressing out unwanted sources of variation = total counts, % mito counts

# %%
# %%time
sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'], n_jobs=5)

# %%
data.layers["regressed"] = data.X.copy()

# %% [markdown]
# ### PCA

# %% [markdown]
# #### scaling to unit variance and zero mean

# %%
sc.pp.scale(data)

# %% [markdown]
# #### compute PCA

# %%
# %%time
sc.tl.pca(data, svd_solver='arpack', n_comps = 50)

# %% [markdown]
# #### plot PCA

# %%
fig, ax = plt.subplots(1,1, figsize=(3,2))
fig.tight_layout()
sc.pl.pca_scatter(data, ax=ax, color=['library'])

fig.subplots_adjust(right=0.8)

# %%
samples =  data.obs.library.unique().to_list()

# %%
fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(6, 10))

ax = ax.ravel()

i=0

for sample in samples:
    sc.pl.pca_scatter(data, ax=ax[i], show=False)
    sc.pl.pca_scatter(data[data.obs['library'] == sample], color=['library'], title=sample, legend_loc=None, ncols=1, ax=ax[i], show=False)
    i=i+1


# %% [markdown]
# #### elbow plot + first PCs

# %%
sc.pl.pca_variance_ratio(data, log=True, n_pcs=50)

# %% jupyter={"outputs_hidden": true}
sc.pl.pca_variance_ratio(data, log=True)

# %%
sc.pl.pca_loadings(data)

# %% [markdown]
# ### UMAP

# %%
# %%time
sc.pp.neighbors(data)

# %%
sc.tl.umap(data, n_components=2)

# %%
sc.pl.umap(data, color=['library'], ncols=1)

# %%
fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(6, 8.5))

ax = ax.ravel()

i=0

for sample in samples:
    sc.pl.umap(data, ax=ax[i], show=False)
    sc.pl.umap(data[data.obs['library'] == sample], color=['library'], title=sample, legend_loc=None, ncols=1, ax=ax[i], show=False)
    i=i+1
    


# %% [markdown]
# **COMMENT** this does not consider batch effect between libraries.
# In iteration_A we integrated with BKNN. It makes sense BUT I don't think the segregation here is crazy:
# * non-trated conditions are similar
#   * SITTB2 includes also one little "horn" that could be due to the fact that CD14-CD15-CD66b- could be neut progenitors in SITTA2 (MUT) and mono progenitors in SITTB2 (WT)
#   * SITTD2 and SITTE2 are enriched and have an extra "blob" of cells each, that is different between mut and wt. Could this be a real difference in transcriptional classes of monocytes that are present?
# * IFN-y treated conditions are similar: there is a clear IFNy switch but they are quite matched between the two 

# %% [markdown]
# ---

# %% [markdown]
# ### Batch effect removal with  Batch Balanced KNN
# https://github.com/Teichlab/bbknn

# %%
import bbknn
#bknn version 1.6.0

# %% [markdown]
# ##### clustering overview before bbkkn

# %%
sc.tl.leiden(data, resolution=0.5, key_added='leiden.0.5')

# %%
fig, ax = plt.subplots(1,2, figsize=(11,4))

sc.pl.umap(data, color='library', ax=ax[0], show=False )
sc.pl.umap(data, color='leiden.0.5', ax=ax[1], show=False )
plt.subplots_adjust( wspace=0.45, right=0.8 )


# %%
# %%time

num_pcs = 18 #see elbow plot

#create a new object
adata_bbknn  = bbknn.bbknn(data, 
                                 neighbors_within_batch=3, #default
                                 n_pcs=num_pcs,
                                 approx=False, 
                                 copy=True, 
                                 batch_key='library')



sc.tl.umap(adata_bbknn)
sc.tl.leiden(adata_bbknn, resolution=0.5, key_added='bb.leiden.0.5')

# %%
fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(6, 8.5))

ax = ax.ravel()

i=0

for sample in samples:
    sc.pl.umap(adata_bbknn, ax=ax[i], show=False)
    sc.pl.umap(adata_bbknn[adata_bbknn.obs['library'] == sample], color=['library'], title=sample, legend_loc=None, ncols=1, ax=ax[i], show=False)
    i=i+1
    


# %%
fig, ax = plt.subplots(1,2, figsize=(11, 4))

sc.pl.umap(adata_bbknn, color='library', ax=ax[0], show=False )
sc.pl.umap(adata_bbknn, color='bb.leiden.0.5', ax=ax[1], show=False)

plt.subplots_adjust( wspace=0.45, right=0.8 )

# %% [markdown]
# #### export adata

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN'
print(filename)
print('\n')

sc.write(filename, adata_bbknn)

# %% [markdown]
# ---
