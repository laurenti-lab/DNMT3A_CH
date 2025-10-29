# %% [markdown]
# # 02.05_Clustering
#
# Author:GM / HB <br>

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
os.makedirs(os.path.join(basedir, 'output/figures', '03.01'), exist_ok=True)
figdir = os.path.join(basedir, 'output/figures', '03.01')
os.makedirs(os.path.join(basedir, 'output/tables', '03.01'), exist_ok=True)
tabdir = os.path.join(basedir, 'output/tables', '03.01')


# %% [markdown]
# #### define functions and dictionaries

# %%
def save_dex_results(adata, key, grouping, basename, outpath):
    
    statcols =['gene_name', 'scores', 'logfoldchanges','pvals', 'pvals_adj' ]
    
    groups = adata.obs[grouping].unique()   
    n = adata.shape[1]
    
    marker_genes_full = {}   
    for g in groups:
        marker_genes = np.reshape( np.array(adata.uns[key]['names'][g]), (n,1) )
       
        marker_genes_full.setdefault(g, None )
        marker_genes_full[g] = np.reshape( np.array(adata.uns[key]['names'][g]), (n,1) )
        marker_genes_full[g] = np.concatenate( [marker_genes_full[g], 
                                                   np.reshape(adata.uns[key]['scores'][g], (n,1) ) ], axis=1)
        marker_genes_full[g] = np.concatenate( [marker_genes_full[g], 
                                                   np.reshape(adata.uns[key]['logfoldchanges'][g], (n,1) ) ], axis=1)
        marker_genes_full[g] = np.concatenate( [marker_genes_full[g], 
                                                   np.reshape(adata.uns[key]['pvals'][g], (n,1) ) ], axis=1)
        marker_genes_full[g] = np.concatenate( [marker_genes_full[g], 
                                                   np.reshape(adata.uns[key]['pvals_adj'][g], (n,1) ) ], axis=1)   

    for c in marker_genes_full.keys():
        temp = pd.DataFrame(marker_genes_full[c], columns=statcols)   
        temp.to_csv( os.path.join(outpath,  
                    prefix+basename+key+'_group_'+c+'.txt'), 
                    index=None, sep='\t')


# %%
def export_top_n_diff_genes(adata, key, groups, n, basename, outpath):
  
    marker_genes = np.reshape( np.array(adata.uns[key]['names'][groups[0]][:n]  ), (n,1) )

    for c in groups[1:]:
        marker_genes = np.concatenate( [marker_genes, np.reshape(adata.uns[key]['names'][c][:n], (n,1) ) ], axis=1)

    markers_frame = pd.DataFrame( marker_genes, columns=groups )

    markers_frame.to_csv(os.path.join(outpath,
                         prefix+basename+key+'_top'+str(n)+'_gene_markers.txt'), 
                         index=None, sep='\t')


# %%
def bundled_DEx_shenanigans(adata, group, key, bname, outpath, n_top_markers=100):
    
    # beware of hard-coded variables (backtracked also)
    
    sc.tl.rank_genes_groups(adata, 
                        groupby=group,
                        key_added=key,
                        method='wilcoxon', 
                        tie_correct=False, #True,     
                        use_raw=False)
    
    save_dex_results(adata, key, group, bname, outpath)
    
    targets = adata.obs[group].unique()
    export_top_n_diff_genes(adata, key, targets, n_top_markers, bname, outpath)


# %%
marker_genes = {
    "fig": ['DEFA3', 'ELANE', 'MPO', 'SMC4', 'LTF', 'C5AR1', 'S100A12', 'CYTIP', 'HLA-DRA', 'CD14', 'CD36', 'IRF8']

}

# %%
#run this after loading data
marker_genes_in_data = dict()
for ct, markers in marker_genes.items():
    markers_found = list()
    for marker in markers:
        if marker in data.var.index:
            markers_found.append(marker)
    marker_genes_in_data[ct] = markers_found


# %%
def markerplots(data, celltype, size=30, sort_order=False):
    print(celltype)
    sc.pl.umap(
        data,
        color=marker_genes_in_data[celltype],
        vmin=0,
        vmax="p99",  
        sort_order=sort_order,  
        frameon=False,
        cmap="Reds",
        size=size
    )
    


# %% [markdown]
# ### Load data

# %%
# %%time
data= sc.read('20241009_CHIP27_after_lognorm_BBKNN_Kwok')

# %% [markdown]
# ## Leiden clustering

# %%
sc.tl.leiden(data, resolution=0.25, key_added='bb.leiden.0.25')

# %%
sc.tl.leiden(data, resolution=1, key_added='bb.leiden.1')

# %% [markdown]
# ### clusters per label tx

# %%
data.obs.groupby('predicted.fine_annot').size()

# %% jupyter={"outputs_hidden": true}
fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(7, 12))

ax=ax.ravel()

vector=['Classical_monocytes','Cycling_TNK','Cycling_neutrophil_progenitors', 'IL1R2+_immature_neutrophils','MPO+_immature_neutrophils_or_progenitors', 'Degranulating_neutrophils', 'Mature_neutrophils','S100A8-9_hi_neutrophils']


sc.pl.umap(data, ax=ax[0], show=False)
sc.pl.umap(data[data.obs['predicted.fine_annot'].isin(vector)], color='predicted.fine_annot', ax=ax[0], show=False )
sc.pl.umap(data, color='bb.leiden.0.5', ax=ax[1], show=False )
sc.pl.umap(data, color='bb.leiden.0.25', ax=ax[2], show=False )
sc.pl.umap(data, color='bb.leiden.1', ax=ax[3], show=False )

# for axis in ax:
#     handles, labels = axis.get_legend_handles_labels()
#     axis.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, -0.05), fancybox=True, shadow=True, frameon=False, ncol=2)

plt.tight_layout() 

# %%
plt.gcf()

# %%
samples =  data.obs.library.unique().to_list()

fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(5, 7))
ax=ax.ravel()
i=0
for sample in samples:
    sc.pl.umap(data, ax=ax[i], show=False)
    sc.pl.umap(data[data.obs['library'] == sample], color=['library'], title=sample, legend_loc=None, ncols=1, ax=ax[i], show=False)
    i=i+1
plt.tight_layout() 

# %%
plt.gcf()

# %% [markdown]
# ### export anndata

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok_clustering'
print(filename)
print('\n')

sc.write(filename, data)

# %% [markdown]
# ## Marker gene viz

# %%
data.X = data.layers["lognorm"] 

# %%
counter=0

# %%
celltype= list(marker_genes.keys())[counter]
markerplots(data, celltype=celltype)

# %%
counter += 1
plt.gcf().suptitle(celltype, fontsize=16)
plt.gcf()

# %%
plt.gcf().savefig(f"{figdir}/{prefix}_markergenes.pdf", format="pdf", dpi=300, bbox_inches="tight")

# %%
#RASTERISE UMAP
fig = plt.gcf() 
ax = fig.axes[0]
ax.get_children()[0].set_zorder(-10)
ax.set_rasterization_zorder(0)
fig.savefig(f"{figdir}/{prefix}_markergenes_raster.pdf", format="pdf", dpi=300, bbox_inches="tight")

# %% [markdown]
# ## Final clustering (id: 1120)
# Informed by leiden clustering and label transfer

# %%
mapping = {
    '0': '0',
    '1': '1_7',
    '2': '2_5',
    '3': '3_4',
    '4':'3_4',
    '5': '2_5',
    '6': '6',
    '7': '1_7',
    '8': '8',
    '9': '9',
    '10' :'10'
}


# %%
data.obs['cluster_1120'] = data.obs['bb.leiden.0.5'].map(mapping)


# %%
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 3.5))

ax=ax.ravel()

sc.pl.umap(data, color='bb.leiden.0.5', ax=ax[0], show=False, legend_loc="on data", size=5)
sc.pl.umap(data, color='cluster_1115', ax=ax[1], show=False, legend_loc="on data", size=5)
sc.pl.umap(data, color='cluster_1120', ax=ax[2], show=False, legend_loc="on data", size=5)

plt.tight_layout() 

# %%
plt.gcf()

# %%
figname = os.path.join(figdir, prefix + '_cluster_1115_1120.pdf')
plt.gcf().savefig(figname, dpi=200, bbox_inches='tight')


# %% [markdown]
# #### export anndata

# %%
# %%time

filename = prefix + basename + 'after_lognorm_BBKNN_Kwok_Zeng_cellcycle_clustering1115_1120'
print(filename)
print('\n')

sc.write(filename, data)

# %% [markdown]
# ## export top DEGs per cluster
# important: use lognorm counts

# %%
plt.rcParams['figure.figsize'] = (4, 3)
sc.pl.umap(data, color='cluster_1120', show=False, legend_loc="on data" )

# %%
plt.gcf()

# %%
data.X= data.layers['lognorm'].copy()

# %%
# %%time
bundled_DEx_shenanigans(adata=data, 
                        group='cluster_1120', 
                        key='cluster_1120',
                        bname=basename, 
                        outpath=output_path,
                        n_top_markers=50)
