# %% [markdown]
# # 03.01 Create pseudobulks for DGE
# Author:GM <br>
# useful links: https://www.sc-best-practices.org/conditions/differential_gene_expression.html

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
import anndata as ann

# %%
import random
import sc_toolbox

# #%load_ext rpy2.ipython


# %%

basedir = os.path.join(laurenti, user, project)
 
sc.settings.writedir = os.path.join(basedir, 'output/objects/')

basename = '_CHIP27_'
figdir = os.path.join(basedir, 'output/figures', '05.01')

if not os.path.exists(figdir):
    os.makedirs(figdir)
    print(figdir, "created successfully")
else:
    print(figdir, "already exists")

# %% [markdown]
# ### define functions and dictionaries

# %%
#function edited from vignette in link
NUM_OF_CELL_PER_DONOR = 25

def aggregate_and_filter(
    adata, #data_DEG
    cell_identity, #this is defined just before running function
    donor_key="library", #library
    condition_key="DNMT3A_MUT", 
    cell_identity_key="cell_type",
    obs_to_keep=[], 
    replicates_per_patient=1,
):
    # subset adata to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR - this is not going to apply to our case because we have many more cells but I left it in the function
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing library {i+1} out of {len(donors)}...", end="\r") 
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # create replicates for each donor
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first" #overrides the default aggregation function and instructs the aggregation process to retain the first encountered value for that observation instead of performing any aggregation operation
                
                df_replicate = pd.DataFrame(adata_replicate.X.A)
                df_replicate.index = adata_replicate.obs_names
                df_replicate.columns = adata_replicate.var_names
                df_replicate = df_replicate.join(adata_replicate.obs[obs_to_keep]) #NB obs_to_keep SHOULD NOT CONTAIN ANY GENE NAME that might cause joining issues due to duplicate values
                
                df_replicate = df_replicate.groupby(donor_key).agg(agg_dict) 
                df_replicate[donor_key] = donor
                df.loc[f"replicate_{donor}_{i}"] = df_replicate.loc[donor]
    print("\n")
    
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop



# %%
genotype_dict={'SITTA2': 'MUT', 'SITTD2': 'MUT', 'SITTC2': 'MUT', 'SITTB2': 'WT', 'SITTE2': 'WT', 'SITTF2': 'WT'}

# %%
treatment_dict={'SITTA2':'CTR', 'SITTB2':'CTR', 'SITTD2':'CTR', 'SITTE2':'CTR', 'SITTC2':'IFNy', 'SITTF2':'IFNy'}

# %%
sort_dict={'SITTA2':'notenriched', 'SITTB2':'notenriched', 'SITTD2':'enriched', 'SITTE2':'enriched', 'SITTC2':'enriched', 'SITTF2':'enriched'}

# %% [markdown]
# ### Import data and cluster numerosity overview

# %%
# %%time
data = sc.read('20241120_CHIP27_after_lognorm_BBKNN_Kwok_Zeng_cellcycle_clustering1115_1120')

# %%
samples =  data.obs.library.unique().to_list()

# %% jupyter={"outputs_hidden": true}
fig, ax = plt.subplots(nrows=2, ncols=5, figsize=(16, 6))

ax = ax.ravel()

i=0

for sample in samples:
    sc.pl.umap(data, ax=ax[i], show=False)
    sc.pl.umap(data[data.obs['library'] == sample], color=['library'], title=sample, legend_loc=None, ncols=1, ax=ax[i], show=False)
    i=i+1

sc.pl.umap(data, color='bb.leiden.0.5', ax=ax[6], show=False, legend_loc="on data" )
sc.pl.umap(data, color='cluster_1115', ax=ax[7], show=False, legend_loc="on data" )
sc.pl.umap(data, color='cluster_1120', ax=ax[8], show=False, legend_loc="on data" )

# %%
plt.gcf()

# %% [markdown]
# #### overview of all number of cells per cluster

# %%
obs_table = (
    data.obs.groupby(['library', 'cluster_1120'])
    .size()
    .reset_index(name='count')
    .pivot(index='library', columns='cluster_1120', values='count')
    .fillna(0)
)

obs_table

# %%
data_DEG = ann.AnnData(data.layers["counts"], 
                       obs=data.obs[['library', 'cluster_1120']],
                       var=data.var[['gene_ids']])

# %%
del data

# %% jupyter={"outputs_hidden": true}
genotype = {
    'SITTA2': 'MUT',
    'SITTB2': 'WT',
    'SITTC2': 'MUT',
    'SITTD2': 'MUT',
    'SITTE2': 'WT',
    'SITTF2': 'WT'
}

stimulation = {
    'SITTA2': 'CTR',
    'SITTB2': 'CTR',
    'SITTC2': 'IFNy',
    'SITTD2': 'CTR',
    'SITTE2': 'CTR',
    'SITTF2': 'IFNy'
}

sort = {
    'SITTA2': 'broad',
    'SITTB2': 'broad',
    'SITTC2': 'enriched',
    'SITTD2': 'enriched',
    'SITTE2': 'enriched',
    'SITTF2': 'enriched'
}


data_DEG.obs['genotype'] = data_DEG.obs['library'].map(genotype)
data_DEG.obs['treatment'] = data_DEG.obs['library'].map(stimulation)
data_DEG.obs['FACS'] = data_DEG.obs['library'].map(sort)


data_DEG.obs.tail()

# %%
columns=['library','cluster_1120','genotype','treatment','FACS']

for col in columns:
    data_DEG.obs[col] = data_DEG.obs[col].astype('category')

# %% [markdown]
# ### create pseudo-bulks
#

# %% [markdown]
# #### run function

# %%
data_DEG.obs_keys()

# %%
obs_to_keep=['library', 'cluster_1120', 'genotype', 'treatment', 'FACS']

# %%
data_DEG.obs['cluster_1120'].cat.categories[0:5]

# %%
data_DEG.obs['cluster_1120'].cat.categories[5]

# %%
# %%time

cell_type = data_DEG.obs["cluster_1120"].cat.categories[5]
print(
    f'Processing cluster {cell_type} (1 out of {len(data_DEG.obs["cluster_1120"].cat.categories)})...'
)
adata_pb = aggregate_and_filter(data_DEG, cell_type, obs_to_keep=obs_to_keep, 
                                cell_identity_key="cluster_1120", replicates_per_patient=1) # too few cells in cluster 8, then removed from the analysis

for i, cell_type in enumerate(data_DEG.obs["cluster_1120"].cat.categories[0:5]): 
    print(
        f'Processing cluster {cell_type} ({i+2} out of {len(data_DEG.obs["cluster_1120"].cat.categories)})...'
    )
    adata_cell_type = aggregate_and_filter(data_DEG, cell_type, obs_to_keep=obs_to_keep, 
                                           cell_identity_key="cluster_1120",replicates_per_patient=2)
    adata_pb = adata_pb.concatenate(adata_cell_type) #this concatenation will add batch=1 to the last anndata added, ignore it

# %%
adata_pb

# %%
adata_pb.obs.drop(columns=['batch'], inplace=True)

# %%
adata_pb.obs['pseudoreplicate'] = adata_pb.obs.groupby(['library', 'cluster_1120']).cumcount()

# %%
adata_pb.obs['pseudor_pair'] = (
    adata_pb.obs['library'].astype(str) + '_' +
    adata_pb.obs['cluster_1120'].astype(str)
)

# %%
adata_pb.obs['pseudor_pair'] = adata_pb.obs['pseudor_pair'].astype('category')

# %% [markdown]
# #### basic EDA to check for outliers
# We perform very basic EDA on the created pseudo-replicates to check if some patients/pseudobulks are outliers that we need to exclude so as not to bias the DE results:

# %%
adata_pb.layers['counts'] = adata_pb.X.copy()

# %%
sc.pp.normalize_total(adata_pb, target_sum=1e5)
sc.pp.log1p(adata_pb)
sc.pp.pca(adata_pb)

# %%
adata_pb.obs["lib_size"] = np.sum(adata_pb.layers["counts"], axis=1)

# %%
adata_pb

# %%
sc.pl.pca(adata_pb, color=adata_pb.obs, ncols=1, size=100)

# %%
plt.gcf()

# %% [markdown]
# 

# %% jupyter={"outputs_hidden": true}
sc.pl.pca(adata_pb[adata_pb.obs['pseudor_pair'].isin(['SITTB2_3_4','SITTA2_0'])], color=adata_pb.obs, ncols=1, size=100)


# %%
plt.gcf()

# %% [markdown]
# #### clean and export pseudobulk anndata
#

# %%
adata_pb.obs.drop(columns=['pseudor_pair'], inplace= True)

# %%
adata_pb.obs['replicate'] = adata_pb.obs['library'].astype(str) + '_' + adata_pb.obs['pseudoreplicate'].astype(str)

# %%
adata_pb.obs.head()

# %%
adata_pb.obs['lib_size'] = adata_pb.obs['lib_size'].astype(float)

# %%
type(adata_pb.layers['counts'])

# %%
np.matrix(adata_pb.layers['counts'])

# %%
adata_pb.layers['counts'].astype(float)

# %%
import scipy.sparse as sp

adata_pb.layers['counts'] = sp.csr_matrix(adata_pb.layers['counts'].astype(float))


# %%
adata_pb.X=adata_pb.layers['counts'].copy()

# %%
# %%time

filename = prefix + basename + 'pseudobulks_cluster_1120'
print(filename)
print('\n')

sc.write(filename, adata_pb)
