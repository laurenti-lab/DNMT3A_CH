# %% [markdown]
# # 02.03_Label transfer visualisation
# Visualise label tx from:
# 1. dataset from Kwok 2023 paper Nat Immunology
# 2. Azimuth PBMC reference map
# 3. (on top) Andy Zeng BM database
#
# Author:GM <br>

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
os.makedirs(os.path.join(basedir, 'output/figures', '02.03'), exist_ok=True)
figdir = os.path.join(basedir, 'output/figures', '02.03')
refdir=os.path.join(basedir, 'output/tables/02.02_labeltx/')
refdir_andy=os.path.join(basedir, 'output/tables/02.04/')

# %% [markdown]
# ---

# %% [markdown]
# ### Load data

# %%
# %%time
data= sc.read('20241004_CHIP27_after_lognorm_BBKNN')

# %% [markdown]
# ### Kwok Nat Immunology 2023

# %%
meta = pd.read_csv(os.path.join(refdir, '20241008_CHIP27_vs_kwokWB_azimuth_metadata.txt'), sep='\t')

# %%
meta = meta[['predicted.fine_annot.score',
       'predicted.fine_annot', 'predicted.broad_annot.score',
       'predicted.broad_annot', 'mapping.score']]

# %%
newmeta= pd.merge(data.obs, meta, left_index=True, right_index=True)

# %%
newmeta.loc[data.obs.index]

# %%
data.obs = newmeta.loc[data.obs.index].copy()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 7))


sc.pl.umap(data, color='predicted.broad_annot', ax=axes[0], show=False)
sc.pl.umap(data, color='predicted.broad_annot.score', ax=axes[1], show=False)
sc.pl.umap(data, color='mapping.score', ax=axes[2], show=False)


plt.tight_layout()
plt.show()

# %% jupyter={"outputs_hidden": true}
fig, ax = plt.subplots(nrows=1, ncols=5, figsize=(14, 3))

ax = ax.ravel()

sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted.broad_annot'] == 'Cycling_TNK'], color='predicted.broad_annot', ax=ax[0], show=False, size=20)

sc.pl.umap(data, ax=ax[0+1], show=False)
sc.pl.umap(data[data.obs['predicted.broad_annot'] == 'Neutrophil_progenitors'], color='predicted.broad_annot', ax=ax[0+1], show=False)

sc.pl.umap(data, ax=ax[1+1], show=False)
sc.pl.umap(data[data.obs['predicted.broad_annot'] == 'Immature_neutrophils'], color='predicted.broad_annot', ax=ax[1+1], show=False,size=20)

sc.pl.umap(data, ax=ax[2+1], show=False)
sc.pl.umap(data[data.obs['predicted.broad_annot'] == 'Mature_neutrophils'], color='predicted.broad_annot', ax=ax[2+1], show=False,size=20)

sc.pl.umap(data, ax=ax[3+1], show=False)
sc.pl.umap(data[data.obs['predicted.broad_annot'] == 'Classical_monocytes'], color='predicted.broad_annot', ax=ax[3+1], show=False)




for axis in ax:
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)

plt.tight_layout()
plt.show()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 7))

vector=['Cycling_TNK','HSPCs','Mast_cells/eosiniophils', 'Non-classical_monocytes', 'cDCs', 'Degranulating_neutrophils']
dotsize=11

sc.pl.umap(data[data.obs['predicted.broad_annot'].isin(vector)], color='predicted.broad_annot', ax=axes[0], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.broad_annot'].isin(vector)], color='predicted.broad_annot.score', ax=axes[1], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.broad_annot'].isin(vector)], color='mapping.score', ax=axes[2], show=False, size=dotsize)


plt.tight_layout()
plt.show()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 9))


sc.pl.umap(data, color='predicted.fine_annot', ax=axes[0], show=False)
sc.pl.umap(data, color='predicted.fine_annot.score', ax=axes[1], show=False)
sc.pl.umap(data, color='mapping.score', ax=axes[2], show=False)


plt.tight_layout()
plt.show()

# %%
data.obs.groupby('predicted.fine_annot').size()

# %% jupyter={"outputs_hidden": true}
fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(13, 6))

ax = ax.ravel()

sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Cycling_TNK'], color='predicted.fine_annot', ax=ax[0], show=False)

sc.pl.umap(data, ax=ax[1], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Cycling_neutrophil_progenitors'], color='predicted.fine_annot', ax=ax[1], show=False,size=50)

sc.pl.umap(data, ax=ax[2], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'MPO+_immature_neutrophils_or_progenitors'], color='predicted.fine_annot', ax=ax[2], show=False,size=50)

sc.pl.umap(data, ax=ax[3], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'IL1R2+_immature_neutrophils'], color='predicted.fine_annot', ax=ax[3], show=False,size=50)

sc.pl.umap(data, ax=ax[4], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'PADI4+_immature_neutrophils'], color='predicted.fine_annot', ax=ax[4], show=False,size=50)

sc.pl.umap(data, ax=ax[5], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'S100A8-9_hi_neutrophils'], color='predicted.fine_annot', ax=ax[5], show=False,size=50)

sc.pl.umap(data, ax=ax[6], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Mature_neutrophils'], color='predicted.fine_annot', ax=ax[6], show=False,size=50)

sc.pl.umap(data, ax=ax[7], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Degranulating_neutrophils'], color='predicted.fine_annot', ax=ax[7], show=False,size=50)


for axis in ax:
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)

plt.tight_layout()
plt.show()

# %%
figname = os.path.join(figdir, prefix + '_kwok_neuts.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')


# %% jupyter={"outputs_hidden": true}
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(8.5, 3))

ax = ax.ravel()

sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Classical_monocytes'], color='predicted.fine_annot', ax=ax[0], show=False)

sc.pl.umap(data, ax=ax[1], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'Non-classical_monocytes'], color='predicted.fine_annot', ax=ax[1], show=False,size=50)

sc.pl.umap(data, ax=ax[2], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'] == 'cDCs'], color='predicted.fine_annot', ax=ax[2], show=False,size=50)


for axis in ax:
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)

plt.tight_layout()
plt.show()

# %%
figname = os.path.join(figdir, prefix + '_kwok_mono.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')


# %% [markdown]
# #### export anndata

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok'
print(filename)
print('\n')

sc.write(filename, data)

# %% [markdown]
# ### Andy Zeng BM dataset

# %%
meta = pd.read_csv(os.path.join(refdir_andy, '20241010_CHIP27_vs_andyBM.csv'), sep=',')

# %%
meta = meta.set_index('Cell')

# %%
meta = meta[['mapping_error_score', 'predicted_CellType', 'predicted_CellType_prob']]

# %%
newmeta= pd.merge(data.obs, meta, left_index=True, right_index=True)

# %%
newmeta.loc[data.obs.index]

# %%
data.obs = newmeta.loc[data.obs.index].copy()

# %%
data.obs.groupby('predicted_CellType').size()

# %%
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 7))

vector=['CD14 Mono','CD16 Mono','Early GMP', 'Early ProMono', 'GMP-Cycle', 'GMP-Mono', 'GMP-Neut', 'Late ProMono', 'Pre-cDC', 'cDC2']
dotsize=2

sc.pl.umap(data[data.obs['predicted_CellType'].isin(vector)], color='predicted_CellType', ax=axes[0], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted_CellType'].isin(vector)], color='predicted_CellType_prob', ax=axes[1], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted_CellType'].isin(vector)], color='mapping_error_score', ax=axes[2], show=False, size=dotsize)


plt.tight_layout()
plt.show()

# %%
figname = os.path.join(figdir, prefix + '_andy_scores.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')

# %%
fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(12, 5))

ax = ax.ravel()

# custom_colors = {
#     'Pre-cDC': 'red'}

sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'GMP-Cycle'], color='predicted_CellType', ax=ax[0], show=False, size=20)

sc.pl.umap(data, ax=ax[1], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'Early GMP'], color='predicted_CellType', ax=ax[1], show=False, size=20)

sc.pl.umap(data, ax=ax[1+1], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'GMP-Neut'], color='predicted_CellType', ax=ax[1+1], show=False,size=20)

sc.pl.umap(data, ax=ax[2+1], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'GMP-Mono'], color='predicted_CellType', ax=ax[2+1], show=False,size=20)

sc.pl.umap(data, ax=ax[3+1], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'Early ProMono'], color='predicted_CellType', ax=ax[3+1], show=False)

sc.pl.umap(data, ax=ax[5], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'Late ProMono'], color='predicted_CellType', ax=ax[5], show=False)

sc.pl.umap(data, ax=ax[6], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'CD14 Mono'], color='predicted_CellType', ax=ax[6], show=False)

sc.pl.umap(data, ax=ax[7], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'CD16 Mono'], color='predicted_CellType', ax=ax[7], show=False, size=20)

sc.pl.umap(data, ax=ax[8], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'Pre-cDC'], color='predicted_CellType', palette=custom_colors, ax=ax[8], show=False, size=20)

sc.pl.umap(data, ax=ax[9], show=False)
sc.pl.umap(data[data.obs['predicted_CellType'] == 'cDC2'], color='predicted_CellType', ax=ax[9], show=False, size=20)


for axis in ax:
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)

plt.tight_layout()
plt.show()

# %%
cell_types = [
    ('GMP-Cycle', 'GMP-Cycle'),
    ('Early GMP', 'Early GMP'),
    ('GMP-Neut', 'GMP-Neut'),
    ('GMP-Mono', 'GMP-Mono'),
    ('Early ProMono', 'Early ProMono'),
    ('Late ProMono', 'Late ProMono'),
    ('CD14 Mono', 'CD14 Mono'),
    ('CD16 Mono', 'CD16 Mono'),
    ('Pre-cDC', 'Pre-cDC'),
    ('cDC2', 'cDC2')
]

fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(11, 4))
ax = ax.ravel()

for i, (cell_type, title) in enumerate(cell_types):
    sc.pl.umap(data, ax=ax[i], show=False)
    sc.pl.umap(data[data.obs['predicted_CellType'] == cell_type], color='predicted_CellType_prob', ax=ax[i], show=False, size=20)
    ax[i].set_title(title)

fig.suptitle('predicted_CellType_prob', fontsize=10, ha='left')
plt.tight_layout()
plt.show()


# %%
figname = os.path.join(figdir, prefix + '_andy_probability.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')

# %% [markdown]
# #### export anndata - appending Andy to Kwok

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok_Zeng'
print(filename)
print('\n')

sc.write(filename, data)

# %% [markdown]
# ### Azimuth human pbmcs

# %%
meta = pd.read_csv(os.path.join(refdir, '20241008_CHIP27_vs_pbmc_azimuth_metadata.txt'), sep='\t')

# %%
meta = meta[['predicted.celltype.l1.score',
       'predicted.celltype.l1', 'predicted.celltype.l2.score',
       'predicted.celltype.l2', 'predicted.celltype.l3.score',
       'predicted.celltype.l3', 'mapping.score']]

# %%
data.obs=data.obs[['library', 'batch', 'DNMT3A_genotype', 'treatment', 'FACS',
       'doublet_score', 'n_genes_by_counts', 'log1p_n_genes_by_counts',
       'total_counts', 'log1p_total_counts', 'total_counts_mt',
       'log1p_total_counts_mt', 'pct_counts_mt', 'total_counts_ribo',
       'log1p_total_counts_ribo', 'pct_counts_ribo', 'total_counts_hb',
       'log1p_total_counts_hb', 'pct_counts_hb', 'n_genes', 'bb.leiden.0.5']]

# %%
data.obs.head()

# %%
newmeta= pd.merge(data.obs, meta, left_index=True, right_index=True) 

# %%
newmeta.loc[data.obs.index]

# %%
data.obs = newmeta.loc[data.obs.index].copy()

# %%
data.obs.groupby('predicted.celltype.l1').size()

# %%
data.obs.groupby('predicted.celltype.l2').size()

# %%
data.obs.groupby('predicted.celltype.l3').size()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(5, 7))


sc.pl.umap(data, color='predicted.celltype.l2', ax=axes[0], show=False)
sc.pl.umap(data, color='predicted.celltype.l2.score', ax=axes[1], show=False)
sc.pl.umap(data, color='mapping.score', ax=axes[2], show=False)


plt.tight_layout()
plt.show()

# %%
figname = os.path.join(figdir, prefix + '_azimuth_pbmc_l2.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')
#figname

# %%
data.obs.groupby('predicted.celltype.l2').size()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(4,8))

vector=['ILC','cDC2', 'Eryth']
dotsize=11

sc.pl.umap(data, ax=axes[0], show=False)
sc.pl.umap(data, ax=axes[1], show=False)
sc.pl.umap(data, ax=axes[2], show=False)

sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='predicted.celltype.l2', ax=axes[0], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='predicted.celltype.l2.score', ax=axes[1], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='mapping.score', ax=axes[2], show=False, size=dotsize)


plt.tight_layout()
plt.show()

# %%
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(3, 8))

ax = ax.ravel()

sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted.celltype.l2'] == 'ILC'], color='predicted.celltype.l2', ax=ax[0], show=False, size=50)

sc.pl.umap(data, ax=ax[1], show=False)
sc.pl.umap(data[data.obs['predicted.celltype.l2'] == 'cDC2'], color='predicted.celltype.l2', ax=ax[1], show=False,size=50)

sc.pl.umap(data, ax=ax[2], show=False)
sc.pl.umap(data[data.obs['predicted.celltype.l2'] == 'Eryth'], color='predicted.celltype.l2', ax=ax[2], show=False,size=50)

for axis in ax:
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, ncol=2)
    
plt.tight_layout()
plt.show()

# %% [markdown]
# ### Azimuth human bmmc

# %%
# %%time
data= sc.read('20241004_CHIP27_after_lognorm_BBKNN')

# %%
meta = pd.read_csv(os.path.join(refdir, '20241009_CHIP27_vs_bmmc_azimuth_metadata.txt'), sep='\t')

# %%
meta = meta[['predicted.celltype.l2.score',
       'predicted.celltype.l2', 'predicted.celltype.l1.score',
       'predicted.celltype.l1', 'mapping.score']]

# %%
data.obs.columns

# %%
newmeta= pd.merge(data.obs, meta, left_index=True, right_index=True)

# %%
newmeta.loc[data.obs.index]

# %%
data.obs = newmeta.loc[data.obs.index].copy()

# %%
data.obs.groupby('predicted.celltype.l1').size()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(4.5, 7))


sc.pl.umap(data, color='predicted.celltype.l1', ax=axes[0], show=False)
sc.pl.umap(data, color='predicted.celltype.l1.score', ax=axes[1], show=False)
sc.pl.umap(data, color='mapping.score', ax=axes[2], show=False)


plt.tight_layout()
plt.show()

# %%
data.obs.groupby('predicted.celltype.l2').size()

# %% jupyter={"outputs_hidden": true}
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(4.5, 7))


sc.pl.umap(data, color='predicted.celltype.l2', ax=axes[0], show=False)
sc.pl.umap(data, color='predicted.celltype.l2.score', ax=axes[1], show=False)
sc.pl.umap(data, color='mapping.score', ax=axes[2], show=False)


plt.tight_layout()
plt.show()

# %%
figname = os.path.join(figdir, prefix + '_azimuth_bmmc_l2.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')
#figname

# %%
data.obs.groupby('predicted.celltype.l2').size()

# %%
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(4,8))

vector=['Macrophage','HSC', 'cDC2', 'pre-mDC']
dotsize=11

sc.pl.umap(data, ax=axes[0], show=False)
sc.pl.umap(data, ax=axes[1], show=False)
sc.pl.umap(data, ax=axes[2], show=False)

sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='predicted.celltype.l2', ax=axes[0], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='predicted.celltype.l2.score', ax=axes[1], show=False, size=dotsize)
sc.pl.umap(data[data.obs['predicted.celltype.l2'].isin(vector)], color='mapping.score', ax=axes[2], show=False, size=dotsize)


plt.tight_layout()
plt.show()
